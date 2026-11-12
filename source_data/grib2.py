"""GRIB2消息与场级采集:严格边界、位图归属、简单打包解码与扫描反演。

仅按Section 0声明的64位总长度切分顺序相接的edition=2消息,消息体内的
GRIB字节不构成新边界;错误长度、截断、非法节序与越界节一律抛出GribError,
不静默跳过。数值解码限于网格模板3.0、产品模板4.0与简单打包模板5.0;未知
模板保留完整原始节与undecoded状态,不沿用前一场的参数。消息边界重置网格、
局部使用节与位图继承。
"""
import hashlib
import struct
from fractions import Fraction

END_MARK = b'7777'


class GribError(ValueError):
    """GRIB2结构或资料约定校验失败。"""


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _u16(data, offset):
    return int.from_bytes(data[offset:offset + 2], 'big')


def _u32(data, offset):
    return int.from_bytes(data[offset:offset + 4], 'big')


def _s32(data, offset):
    return int.from_bytes(data[offset:offset + 4], 'big', signed=True)


def _scale16(data, offset):
    """E/D尺度:GRIB规范的符号位绝对值表示,不是二进制补码。"""
    raw = _u16(data, offset)
    return -(raw & 0x7FFF) if raw & 0x8000 else raw


def _bit(data, index):
    return (data[index >> 3] >> (7 - (index & 7))) & 1


def _bits(data, start, length):
    value = 0
    for k in range(start, start + length):
        value = (value << 1) | _bit(data, k)
    return value


def _public(record):
    return {k: v for k, v in record.items() if not k.startswith('_')}


def split_messages(raw, max_message_bytes=None):
    """按Section 0声明的总长度切出消息,原封保留字节,不重新同步。"""
    if not raw:
        raise GribError('输入为空,没有GRIB2消息')
    messages = []
    pos = 0
    while pos < len(raw):
        rest = len(raw) - pos
        if rest < 16:
            raise GribError(f'偏移{pos}处剩余{rest}字节,不足Section 0')
        if raw[pos:pos + 4] != b'GRIB':
            raise GribError(f'偏移{pos}处不是GRIB标识,声明的长度链已断裂')
        if raw[pos + 7] != 2:
            raise GribError(f'偏移{pos}处edition={raw[pos + 7]},仅支持2')
        total = int.from_bytes(raw[pos + 8:pos + 16], 'big')
        if total < 20:
            raise GribError(f'偏移{pos}处声明总长度{total}非法')
        if max_message_bytes is not None and total > max_message_bytes:
            raise GribError(f'偏移{pos}处消息{total}字节超出上限{max_message_bytes}')
        if total > rest:
            raise GribError(f'偏移{pos}处消息声明{total}字节,输入截断')
        message = raw[pos:pos + total]
        if message[-4:] != END_MARK:
            raise GribError(f'偏移{pos}处消息未以7777结束')
        messages.append({'index': len(messages), 'offset': pos, 'raw': message})
        pos += total
    return messages


def _walk_sections(message, message_offset):
    """顺序读出消息内全部节,节必须相接且不得越过7777。"""
    sections = []
    pos = 16
    end = len(message) - 4
    while pos < end:
        if end - pos < 5:
            raise GribError(f'消息偏移{message_offset + pos}处节头不足5字节')
        length = _u32(message, pos)
        if length < 5:
            raise GribError(f'消息偏移{message_offset + pos}处节长度{length}非法')
        if pos + length > end:
            raise GribError(f'消息偏移{message_offset + pos}处节越过7777结束标记')
        sections.append((message[pos + 4], pos, length))
        pos += length
    if not sections:
        raise GribError(f'消息偏移{message_offset}处没有任何节')
    return sections


def _parse_section1(sec):
    if len(sec) != 21:
        raise GribError('Section 1长度必须为21')
    return {
        'center': _u16(sec, 5),
        'subcenter': _u16(sec, 7),
        'master_version': sec[9],
        'local_version': sec[10],
        'significance': sec[11],
        'reference_time': '%04d-%02d-%02dT%02d:%02d:%02dZ' % (
            _u16(sec, 12), sec[14], sec[15], sec[16], sec[17], sec[18]),
        'production_status': sec[19],
        'data_type': sec[20],
    }


def _parse_section3(sec, offset):
    """网格定义;只支持无缩减列表的模板3.0,其余保留点数并标记未支持。"""
    if len(sec) < 14:
        raise GribError('Section 3长度不足14字节')
    grid = {'template': _u16(sec, 12), 'points': _u32(sec, 6), 'offset': offset,
            'sha256': _sha(sec), 'supported': False}
    if grid['template'] != 0:
        grid['unsupported_reason'] = f"网格模板3.{grid['template']}超出解码范围"
        return grid
    if sec[5] != 0:
        grid['unsupported_reason'] = '网格来源不是模板给定'
        return grid
    if sec[10]:
        grid['unsupported_reason'] = '不接受带缩减列表的网格'
        return grid
    if len(sec) != 72:
        raise GribError('模板3.0的Section 3长度必须为72')
    ni = _u32(sec, 30)
    nj = _u32(sec, 34)
    if ni == 0 or nj == 0:
        raise GribError('Ni/Nj必须为正')
    if ni * nj != grid['points']:
        raise GribError(f'格点数{grid["points"]}与Ni*Nj={ni * nj}不符')
    basic = _u32(sec, 38)
    subdivisions = _u32(sec, 42)
    if basic == 0 or subdivisions in (0, 0xFFFFFFFF):
        factor = Fraction(1, 1000000)
    else:
        factor = Fraction(basic, subdivisions)
    la1 = _s32(sec, 46)
    lo1 = _s32(sec, 50)
    la2 = _s32(sec, 55)
    lo2 = _s32(sec, 59)
    di = _u32(sec, 63)
    dj = _u32(sec, 67)
    if abs(la2 - la1) != (nj - 1) * dj:
        raise GribError('La1/La2与Nj、Dj不一致')
    turn = Fraction(360) / factor
    if Fraction(lo2 - lo1) % turn != Fraction((ni - 1) * di) % turn:
        raise GribError('Lo1/Lo2与Ni、Di不一致')
    grid.update({'supported': True, 'ni': ni, 'nj': nj, 'scan': sec[71],
                 'la1': la1, 'lo1': lo1, 'la2': la2, 'lo2': lo2,
                 'di': di, 'dj': dj, 'factor': factor, 'resolution': sec[54]})
    return grid


def _surface_value(surface_type, scale_byte, value):
    if surface_type == 0xFF or (scale_byte == 0xFF and value == 0xFFFFFFFF):
        return None
    scale = int.from_bytes(bytes([scale_byte]), 'big', signed=True)
    return float(Fraction(value) * Fraction(10) ** (-scale))


def _parse_section4(sec, offset):
    """产品定义;模板4.0之外只取参数类别/编号,层面与预报时长不猜测。"""
    if len(sec) < 9:
        raise GribError('Section 4长度不足9字节')
    product = {'template': _u16(sec, 7), 'offset': offset,
               'sha256': _sha(sec), 'supported': False}
    if len(sec) >= 11:
        product['parameter_category'] = sec[9]
        product['parameter_number'] = sec[10]
    if product['template'] != 0:
        return product
    coordinates = _u16(sec, 5)
    if len(sec) != 34 + 4 * coordinates:
        raise GribError('模板4.0的Section 4长度与坐标值数不符')
    product.update({
        'supported': True,
        'generating_process_type': sec[11],
        'background_process': sec[12],
        'generating_process_id': sec[13],
        'cutoff_hours': _u16(sec, 14),
        'cutoff_minutes': sec[16],
        'forecast_unit': sec[17],
        'forecast_time': _s32(sec, 18),
        'surface1_type': sec[22],
        'surface1_value': _surface_value(sec[22], sec[23], _u32(sec, 24)),
        'surface2_type': sec[28],
        'surface2_value': _surface_value(sec[28], sec[29], _u32(sec, 30)),
    })
    return product


def _parse_section5(sec, offset):
    """数据表示;模板5.0之外只取打包点数,不猜测打包参数。"""
    if len(sec) < 11:
        raise GribError('Section 5长度不足11字节')
    packing = {'template': _u16(sec, 9), 'packed_count': _u32(sec, 5),
               'offset': offset, 'sha256': _sha(sec), 'supported': False}
    if packing['template'] != 0:
        return packing
    if len(sec) != 21:
        raise GribError('模板5.0的Section 5长度必须为21')
    reference = struct.unpack('>f', sec[11:15])[0]
    packing.update({
        'supported': True,
        '_reference': Fraction(reference),
        'reference_value': reference,
        'binary_scale': _scale16(sec, 15),
        'decimal_scale': _scale16(sec, 17),
        'bits_per_value': sec[19],
        'original_type': sec[20],
    })
    return packing


def _parse_bitmap(sec, offset, points, packed, state, field_index, message_index):
    """位图节;0本节位图,255无位图,254复用本消息内先前明示位图。"""
    if len(sec) < 6:
        raise GribError(f'消息{message_index}的Section 6长度不足6字节')
    indicator = sec[5]
    if indicator == 0:
        need = (points + 7) // 8
        if len(sec) != 6 + need:
            raise GribError(f'消息{message_index}的位图字节数与格点数{points}不符')
        bits = sec[6:]
        pad = need * 8 - points
        if pad and _bits(bits, points, pad):
            raise GribError(f'消息{message_index}的位图尾部补齐位非零')
        ones = sum(_bit(bits, k) for k in range(points))
        if ones != packed:
            raise GribError(f'消息{message_index}的位图1数量{ones}与打包值数{packed}不符')
        info = {'indicator': 0, 'source': 'explicit', 'section_offset': offset,
                'points': points, 'ones': ones, '_bits': bits, '_field': field_index}
        return info, info
    if indicator == 255:
        if len(sec) != 6:
            raise GribError(f'消息{message_index}的无位图Section 6长度必须为6')
        if packed != points:
            raise GribError(f'消息{message_index}无位图时打包值数{packed}必须等于格点数{points}')
        return {'indicator': 255, 'source': 'none'}, state
    if indicator == 254:
        if len(sec) != 6:
            raise GribError(f'消息{message_index}的复用Section 6长度必须为6')
        if state is None:
            raise GribError(f'消息{message_index}的254复用前,当前消息内没有明示位图')
        if state['points'] != points:
            raise GribError(
                f'消息{message_index}的254位图点数{state["points"]}与当前格点数{points}不兼容')
        if state['ones'] != packed:
            raise GribError(
                f'消息{message_index}的位图1数量{state["ones"]}与打包值数{packed}不符')
        return {'indicator': 254, 'source': 'reused',
                'source_section_offset': state['section_offset'],
                'source_field': state['_field'],
                'points': state['points'], 'ones': state['ones'],
                '_bits': state['_bits']}, state
    raise GribError(f'消息{message_index}的预定义位图编号{indicator}不支持')


def _scan_map(ni, nj, mode):
    """按扫描标志生成扫描顺序到网格坐标(gi,gj)的映射。

    gi沿网格定义的i方向(Lo1向Lo2),gj沿La1向La2;标志位依次处理i负向、
    j正向、邻接j与相邻行反向四类。
    """
    i_negative = bool(mode & 0x80)
    j_positive = bool(mode & 0x40)
    j_consecutive = bool(mode & 0x20)
    alternate = bool(mode & 0x10)
    order = []
    outer, inner = (ni, nj) if j_consecutive else (nj, ni)
    for a in range(outer):
        for b in range(inner):
            if j_consecutive:
                gi_scan = a
                gj_scan = nj - 1 - b if alternate and a % 2 else b
            else:
                gj_scan = a
                gi_scan = ni - 1 - b if alternate and a % 2 else b
            gi = ni - 1 - gi_scan if i_negative else gi_scan
            gj = gj_scan if j_positive else nj - 1 - gj_scan
            order.append((gi, gj))
    return order


def _decode_values(grid, packing, bitmap, sec7, sec7_offset):
    """按Y=(R+X*2^E)*10^(-D)解码,位流高位在先,缺测点不参与解码。"""
    bits = packing['bits_per_value']
    count = packing['packed_count']
    data = sec7[5:]
    expected = (count * bits + 7) // 8
    if len(data) != expected:
        raise GribError(f'Section 7数据{len(data)}字节与打包参数期望{expected}不符')
    if count and bits:
        pad = expected * 8 - count * bits
        if pad and _bits(data, count * bits, pad):
            raise GribError('Section 7尾部补齐位非零')
    xs = [_bits(data, k * bits, bits) for k in range(count)] if bits else [0] * count
    binary = Fraction(2) ** packing['binary_scale']
    decimal = Fraction(10) ** (-packing['decimal_scale'])
    ni = grid['ni']
    nj = grid['nj']
    flip = grid['la2'] < grid['la1']
    factor = grid['factor']
    turn = Fraction(360) / factor
    la_south = min(grid['la1'], grid['la2'])
    bitmap_bits = bitmap.get('_bits')
    entries = {}
    packed_index = 0
    for k, (gi, gj) in enumerate(_scan_map(ni, nj, grid['scan'])):
        j = nj - 1 - gj if flip else gj
        lat = float((la_south + j * grid['dj']) * factor)
        lon = float(((grid['lo1'] + gi * grid['di']) % turn) * factor)
        if bitmap_bits is not None and not _bit(bitmap_bits, k):
            entries[j, gi] = {'j': j, 'i': gi, 'lat': lat, 'lon': lon,
                              'missing': True, 'value': None}
            continue
        x = xs[packed_index]
        y = (packing['_reference'] + x * binary) * decimal
        bit_range = None
        if bits:
            bit_range = {'section7_bit_offset': packed_index * bits,
                         'bit_length': bits,
                         'file_bit_offset': (sec7_offset + 5) * 8 + packed_index * bits}
        entries[j, gi] = {'j': j, 'i': gi, 'lat': lat, 'lon': lon,
                          'missing': False, 'x': x, 'value': float(y),
                          'bit_range': bit_range}
        packed_index += 1
    if packed_index != count:
        raise GribError('位图1数量与打包值数不符')
    return [entries[j, i] for j in range(nj) for i in range(ni)], packed_index


def _grid_report(grid, source):
    report = {'template': f"3.{grid['template']}", 'supported': grid['supported'],
              'points': grid['points'], 'section_offset': grid['offset'],
              'sha256': grid['sha256'], 'source': source}
    if grid['supported']:
        factor = grid['factor']
        report.update({
            'ni': grid['ni'], 'nj': grid['nj'],
            'scanning_mode': grid['scan'],
            'scanning_flags': {
                'i_negative': bool(grid['scan'] & 0x80),
                'j_positive': bool(grid['scan'] & 0x40),
                'j_consecutive': bool(grid['scan'] & 0x20),
                'alternate_rows': bool(grid['scan'] & 0x10)},
            'la1': float(grid['la1'] * factor), 'lo1': float(grid['lo1'] * factor),
            'la2': float(grid['la2'] * factor), 'lo2': float(grid['lo2'] * factor),
            'di': float(grid['di'] * factor), 'dj': float(grid['dj'] * factor),
            'resolution_flags': grid['resolution']})
    elif 'unsupported_reason' in grid:
        report['unsupported_reason'] = grid['unsupported_reason']
    return report


def _product_report(product):
    report = {'template': f"4.{product['template']}", 'supported': product['supported'],
              'section_offset': product['offset'], 'sha256': product['sha256']}
    for key in ('parameter_category', 'parameter_number', 'generating_process_type',
                'background_process', 'generating_process_id', 'cutoff_hours',
                'cutoff_minutes', 'forecast_unit', 'forecast_time',
                'surface1_type', 'surface1_value', 'surface2_type', 'surface2_value'):
        if key in product:
            report[key] = product[key]
    return report


def _packing_report(packing):
    report = {'template': f"5.{packing['template']}", 'supported': packing['supported'],
              'section_offset': packing['offset'], 'sha256': packing['sha256'],
              'packed_count': packing['packed_count']}
    for key in ('reference_value', 'binary_scale', 'decimal_scale',
                'bits_per_value', 'original_type'):
        if key in packing:
            report[key] = packing[key]
    return report


def _finish_field(message_index, message_offset, field_index, identification,
                  pending, sec7, sec7_offset, accepted):
    grid = pending['grid']
    product = pending['product']
    packing = pending['packing']
    reasons = []
    if not grid['supported']:
        reasons.append(grid.get('unsupported_reason', '网格模板超出解码范围'))
    elif grid['template'] not in accepted['grid']:
        reasons.append('网格模板未被资料配置接受')
    if not product['supported']:
        reasons.append(f"产品模板4.{product['template']}超出解码范围")
    elif product['template'] not in accepted['product']:
        reasons.append('产品模板未被资料配置接受')
    if not packing['supported']:
        reasons.append(f"打包模板5.{packing['template']}超出解码范围")
    elif packing['template'] not in accepted['packing']:
        reasons.append('打包模板未被资料配置接受')
    field = {
        'index': field_index,
        'identity': {
            'message_index': message_index,
            'message_offset': message_offset,
            'field_index': field_index,
            'parameter_category': product.get('parameter_category'),
            'parameter_number': product.get('parameter_number'),
            'surface_type': product.get('surface1_type'),
            'level_value': product.get('surface1_value'),
            'reference_time': identification['reference_time'],
            'forecast_time': product.get('forecast_time'),
            'forecast_unit': product.get('forecast_unit'),
        },
        'identification': identification,
        'sections': {
            'grid': grid['offset'],
            'product': product['offset'],
            'packing': packing['offset'],
            'bitmap': pending['bitmap_offset'],
            'data': sec7_offset,
            'local_use': pending['local_use']['offset'] if pending['local_use'] else None,
        },
        'local_use': pending['local_use'],
        'grid': _grid_report(grid, 'defined' if pending['grid_fresh'] else 'inherited'),
        'product': _product_report(product),
        'packing': _packing_report(packing),
        'bitmap': _public(pending['bitmap']),
    }
    if reasons:
        field['status'] = 'undecoded'
        field['undecoded_reason'] = ';'.join(reasons)
        field['values'] = None
        field['raw_sections_hex'] = {
            'product': pending['product_sec'].hex(),
            'packing': pending['packing_sec'].hex(),
            'bitmap': pending['bitmap_sec'].hex(),
            'data': sec7.hex(),
        }
    else:
        values, present = _decode_values(grid, packing, pending['bitmap'], sec7, sec7_offset)
        field['status'] = 'decoded'
        field['constant'] = packing['bits_per_value'] == 0
        field['value_count'] = {'packed': packing['packed_count'], 'present': present,
                                'missing': grid['points'] - present}
        field['values'] = values
    return field


def _parse_message(message, index, offset, accepted):
    sections = _walk_sections(message, offset)
    if sections[0][0] != 1:
        raise GribError(f'消息{index}的首节必须是Section 1')
    identification = _parse_section1(message[sections[0][1]:sections[0][1] + sections[0][2]])
    grid = None
    grid_fresh = False
    local_use = None
    bitmap_state = None
    pending = None
    fields = []
    for number, at, length in sections[1:]:
        sec = message[at:at + length]
        absolute = offset + at
        if number == 1:
            raise GribError(f'消息{index}的Section 1重复出现')
        if number == 2:
            if pending is not None:
                raise GribError(f'消息{index}的Section 2位于场组内部')
            local_use = {'offset': absolute, 'length': length, 'sha256': _sha(sec)}
            continue
        if number == 3:
            if pending is not None:
                raise GribError(f'消息{index}的Section 3位于场组内部')
            grid = _parse_section3(sec, absolute)
            grid_fresh = True
            continue
        if number == 4:
            if pending is not None:
                raise GribError(f'消息{index}的Section 4位于场组内部')
            if grid is None:
                raise GribError(f'消息{index}的首个场组缺少Section 3')
            pending = {'grid': grid, 'grid_fresh': grid_fresh, 'local_use': local_use,
                       'product': _parse_section4(sec, absolute), 'product_sec': sec}
            grid_fresh = False
            continue
        if number == 5:
            if pending is None or 'packing' in pending:
                raise GribError(f'消息{index}的Section 5节序非法')
            pending['packing'] = _parse_section5(sec, absolute)
            pending['packing_sec'] = sec
            continue
        if number == 6:
            if pending is None or 'packing' not in pending or 'bitmap' in pending:
                raise GribError(f'消息{index}的Section 6节序非法')
            info, bitmap_state = _parse_bitmap(
                sec, absolute, pending['grid']['points'],
                pending['packing']['packed_count'], bitmap_state, len(fields), index)
            pending['bitmap'] = info
            pending['bitmap_sec'] = sec
            pending['bitmap_offset'] = absolute
            continue
        if number == 7:
            if pending is None or 'bitmap' not in pending:
                raise GribError(f'消息{index}的Section 7节序非法')
            fields.append(_finish_field(index, offset, len(fields), identification,
                                        pending, sec, absolute, accepted))
            pending = None
            continue
        raise GribError(f'消息{index}出现非法节号{number}')
    if pending is not None:
        raise GribError(f'消息{index}的场组在7777前不完整')
    if not fields:
        raise GribError(f'消息{index}不含任何场')
    return {
        'index': index,
        'offset': offset,
        'total_bytes': len(message),
        'sha256': _sha(message),
        'discipline': message[6],
        'identification': identification,
        'sections': [{'number': n, 'offset': offset + at, 'length': length,
                      'sha256': _sha(message[at:at + length])}
                     for n, at, length in sections],
        'fields': fields,
    }


def scan(raw, max_message_bytes=None, accepted_grid=(0,), accepted_product=(0,),
         accepted_packing=(0,)):
    """采集raw中顺序相接的全部GRIB2消息与场,返回可JSON序列化的报告。

    报告给出文件与每条消息的SHA-256、全部节的文件偏移、每场的位图来源,
    以及每个解码值在Section 7中的打包位范围;消息字节按声明长度原封切片,
    可由offset与total_bytes逐字节取回。
    """
    accepted = {'grid': set(accepted_grid), 'product': set(accepted_product),
                'packing': set(accepted_packing)}
    messages = split_messages(raw, max_message_bytes)
    return {
        'edition': 2,
        'file_bytes': len(raw),
        'file_sha256': _sha(raw),
        'message_count': len(messages),
        'messages': [_parse_message(item['raw'], item['index'], item['offset'], accepted)
                     for item in messages],
    }
