"""人工构造GRIB2字节的拼装助手,仅供测试使用;不依赖任何打包库。"""
import struct


def section(number, body):
    return (5 + len(body)).to_bytes(4, 'big') + bytes([number]) + body


def sec1(year=2026, month=1, day=1, hour=0, minute=0, second=0):
    body = (0).to_bytes(2, 'big') + (0).to_bytes(2, 'big') + bytes([2, 1, 1])
    body += year.to_bytes(2, 'big') + bytes([month, day, hour, minute, second])
    body += bytes([0, 1])
    return section(1, body)


def sec2(payload=b'local-use'):
    return section(2, payload)


def sec3(ni, nj, la1, lo1, la2, lo2, di, dj, scan, points=None, template=0):
    """经纬度以微度为单位;template非0时填充未解释模板体。"""
    if points is None:
        points = ni * nj
    body = bytes([0]) + points.to_bytes(4, 'big') + bytes([0, 0]) + template.to_bytes(2, 'big')
    if template == 0:
        body += bytes([6, 0xFF]) + b'\x00' * 4 + bytes([0xFF]) + b'\x00' * 4
        body += bytes([0xFF]) + b'\x00' * 4
        body += ni.to_bytes(4, 'big') + nj.to_bytes(4, 'big')
        body += b'\x00' * 4 + b'\xFF' * 4
        body += la1.to_bytes(4, 'big', signed=True)
        body += lo1.to_bytes(4, 'big', signed=True)
        body += bytes([0x30])
        body += la2.to_bytes(4, 'big', signed=True)
        body += lo2.to_bytes(4, 'big', signed=True)
        body += di.to_bytes(4, 'big') + dj.to_bytes(4, 'big') + bytes([scan])
    else:
        body += b'\x00' * 10
    return section(3, body)


def sec4(category, number, unit=1, forecast=0, surface=1, level=0, template=0):
    body = (0).to_bytes(2, 'big') + template.to_bytes(2, 'big')
    if template == 0:
        body += bytes([category, number, 2, 0, 0])
        body += (0).to_bytes(2, 'big') + bytes([0])
        body += bytes([unit]) + forecast.to_bytes(4, 'big', signed=True)
        body += bytes([surface, 0]) + level.to_bytes(4, 'big')
        body += bytes([0xFF, 0xFF]) + b'\xFF' * 4
    else:
        body += bytes([category, number]) + b'\x00' * 8
    return section(4, body)


def _scale(value):
    return ((0x8000 | abs(value)) if value < 0 else value).to_bytes(2, 'big')


def sec5(reference, binary_scale, decimal_scale, bits, count, template=0):
    body = count.to_bytes(4, 'big') + template.to_bytes(2, 'big')
    if template == 0:
        body += struct.pack('>f', reference)
        body += _scale(binary_scale) + _scale(decimal_scale)
        body += bytes([bits, 0])
    else:
        body += b'\x00' * 6
    return section(5, body)


def sec6_explicit(flags):
    """flags为扫描顺序上每个格点的0/1。"""
    return section(6, bytes([0]) + pack_bits(flags, 1))


def sec6_none():
    return section(6, bytes([255]))


def sec6_reuse():
    return section(6, bytes([254]))


def sec7(data):
    return section(7, data)


def message(*sections, discipline=0, edition=2):
    body = b''.join(sections) + b'7777'
    head = b'GRIB' + b'\x00\x00' + bytes([discipline, edition])
    return head + (16 + len(body)).to_bytes(8, 'big') + body


def pack_bits(values, bits):
    """高位在先打包;最后一个字节剩余低位补零。"""
    out = bytearray()
    acc = 0
    n = 0
    for value in values:
        for k in range(bits - 1, -1, -1):
            acc = (acc << 1) | ((value >> k) & 1)
            n += 1
            if n == 8:
                out.append(acc)
                acc = 0
                n = 0
    if n:
        out.append(acc << (8 - n))
    return bytes(out)


def grid2x2(scan=0x40):
    return sec3(2, 2, 0, 0, 1000000, 1000000, 1000000, 1000000, scan)


def grid3x2(scan=0x00):
    return sec3(3, 2, 0, 0, 1000000, 2000000, 1000000, 1000000, scan)
