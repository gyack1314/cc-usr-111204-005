"""GRIB2消息与场级采集的约定检查;全部样本为人工构造的GRIB2字节。"""
import hashlib
import io
import json
import unittest
from contextlib import redirect_stdout
from pathlib import Path

from source_data import grib2
from source_data.cli import main as cli_main
from source_data.grib2 import GribError
from tests.grib2_samples import (grid2x2, grid3x2, message, pack_bits, sec1,
                                 sec2, sec3, sec4, sec5, sec6_explicit,
                                 sec6_none, sec6_reuse, sec7, section)


def good_message():
    return message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 4),
                   sec6_none(), sec7(pack_bits([0, 1, 2, 3], 2)))


def only_field(report, message_index=0, field_index=0):
    return report['messages'][message_index]['fields'][field_index]


def flat_values(field):
    return [v['value'] for v in field['values']]


class SampleFileTests(unittest.TestCase):
    def test_fields_file(self):
        raw = Path('data/fields.grib2').read_bytes()
        report = grib2.scan(raw)
        json.dumps(report)  # 报告必须可JSON序列化
        self.assertEqual(report['file_bytes'], 248)
        self.assertEqual(report['file_sha256'], hashlib.sha256(raw).hexdigest())
        self.assertEqual(report['message_count'], 1)
        msg = report['messages'][0]
        self.assertEqual(msg['offset'], 0)
        self.assertEqual(msg['total_bytes'], 248)
        self.assertEqual(msg['sha256'], hashlib.sha256(raw[:248]).hexdigest())
        self.assertEqual(msg['identification']['reference_time'],
                         '2030-01-01T00:00:00Z')
        self.assertEqual(
            [(s['number'], s['offset'], s['length']) for s in msg['sections']],
            [(1, 16, 21), (3, 37, 72), (4, 109, 34), (5, 143, 21),
             (6, 164, 7), (7, 171, 6), (4, 177, 34), (5, 211, 21),
             (6, 232, 6), (7, 238, 6)])
        first, second = msg['fields']
        self.assertEqual(first['status'], 'decoded')
        self.assertEqual(first['identity']['message_offset'], 0)
        self.assertEqual(first['identity']['field_index'], 0)
        self.assertEqual(first['identity']['parameter_category'], 0)
        self.assertEqual(first['identity']['parameter_number'], 0)
        self.assertEqual(first['identity']['forecast_time'], 0)
        self.assertEqual(first['identity']['forecast_unit'], 1)
        self.assertEqual(first['identity']['surface_type'], 1)
        self.assertEqual(first['identity']['level_value'], 0.0)
        self.assertEqual(first['grid']['source'], 'defined')
        self.assertTrue(first['grid']['scanning_flags']['j_positive'])
        self.assertEqual(first['bitmap']['source'], 'explicit')
        self.assertEqual(first['bitmap']['section_offset'], 164)
        self.assertEqual(first['bitmap']['ones'], 3)
        self.assertEqual(first['packing']['bits_per_value'], 2)
        self.assertEqual(flat_values(first), [1.0, 2.0, None, 3.0])
        self.assertEqual(first['value_count'],
                         {'packed': 3, 'present': 3, 'missing': 1})
        ranges = [v.get('bit_range') for v in first['values']]
        self.assertEqual([r['section7_bit_offset'] if r else None for r in ranges],
                         [0, 2, None, 4])
        self.assertEqual([r['file_bit_offset'] if r else None for r in ranges],
                         [1408, 1410, None, 1412])
        # 两个场的产品定义逐字节相同,靠场序号区分身份,不能只按参数名覆盖
        self.assertEqual(second['identity']['parameter_number'], 0)
        self.assertEqual(second['identity']['forecast_time'], 0)
        self.assertEqual(second['identity']['field_index'], 1)
        self.assertNotEqual(tuple(sorted(first['identity'].items())),
                            tuple(sorted(second['identity'].items())))
        self.assertEqual(second['grid']['source'], 'inherited')
        self.assertEqual(second['grid']['section_offset'], 37)
        self.assertEqual(second['bitmap']['source'], 'reused')
        self.assertEqual(second['bitmap']['source_field'], 0)
        self.assertEqual(second['bitmap']['source_section_offset'], 164)
        self.assertEqual(flat_values(second), [4.0, 3.0, None, 2.0])

    def test_two_messages_file(self):
        raw = Path('data/two-messages.grib2').read_bytes()
        report = grib2.scan(raw)
        self.assertEqual(report['message_count'], 2)
        self.assertEqual([m['offset'] for m in report['messages']], [0, 248])
        for index, msg in enumerate(report['messages']):
            self.assertEqual(msg['total_bytes'], 248)
            self.assertEqual(msg['sha256'],
                             hashlib.sha256(raw[msg['offset']:msg['offset'] + 248]).hexdigest())
            self.assertEqual(len(msg['fields']), 2)
            for field in msg['fields']:
                self.assertEqual(field['identity']['message_index'], index)
                self.assertEqual(field['identity']['message_offset'], msg['offset'])
        second_message = report['messages'][1]
        explicit = second_message['fields'][0]['bitmap']
        reused = second_message['fields'][1]['bitmap']
        self.assertEqual(explicit['source'], 'explicit')
        self.assertEqual(explicit['section_offset'], 248 + 164)
        self.assertEqual(reused['source'], 'reused')
        self.assertEqual(reused['source_section_offset'], 248 + 164)


class DecodeTests(unittest.TestCase):
    def test_constant_field(self):
        raw = message(sec1(), grid2x2(), sec4(0, 0), sec5(2.5, 0, 1, 0, 4),
                      sec6_none(), sec7(b''))
        field = only_field(grib2.scan(raw))
        self.assertEqual(field['status'], 'decoded')
        self.assertTrue(field['constant'])
        self.assertEqual(flat_values(field), [0.25, 0.25, 0.25, 0.25])
        self.assertEqual([v['x'] for v in field['values']], [0, 0, 0, 0])
        self.assertEqual([v['bit_range'] for v in field['values']],
                         [None, None, None, None])

    def test_negative_scales_are_sign_magnitude(self):
        raw = message(sec1(), grid2x2(), sec4(0, 0), sec5(0.0, -1, -2, 4, 4),
                      sec6_none(), sec7(pack_bits([1, 2, 3, 4], 4)))
        field = only_field(grib2.scan(raw))
        self.assertEqual(field['packing']['binary_scale'], -1)
        self.assertEqual(field['packing']['decimal_scale'], -2)
        # Y=(0+X*2^-1)*10^2=50X;若误按二进制补码解释尺度,结果将完全不同
        self.assertEqual(flat_values(field), [50.0, 100.0, 150.0, 200.0])

    def test_scanning_modes(self):
        expected = {
            0x00: [4, 5, 6, 1, 2, 3],  # i正向、j负向、i邻接
            0x40: [1, 2, 3, 4, 5, 6],  # j正向
            0x80: [6, 5, 4, 3, 2, 1],  # i负向
            0x20: [2, 4, 6, 1, 3, 5],  # 邻接j
            0x10: [6, 5, 4, 1, 2, 3],  # 相邻行反向
            0x30: [2, 3, 6, 1, 4, 5],  # 邻接j且相邻行反向
        }
        for mode, want in expected.items():
            with self.subTest(mode=hex(mode)):
                raw = message(sec1(), grid3x2(mode), sec4(0, 0),
                              sec5(0.0, 0, 0, 4, 6), sec6_none(),
                              sec7(pack_bits([1, 2, 3, 4, 5, 6], 4)))
                field = only_field(grib2.scan(raw))
                self.assertEqual(flat_values(field),
                                 [float(v) for v in want])
                self.assertEqual([v['j'] for v in field['values']],
                                 [0, 0, 0, 1, 1, 1])
                self.assertEqual([v['i'] for v in field['values']],
                                 [0, 1, 2, 0, 1, 2])
        raw = message(sec1(), grid3x2(0x40), sec4(0, 0), sec5(0.0, 0, 0, 4, 6),
                      sec6_none(), sec7(pack_bits([1, 2, 3, 4, 5, 6], 4)))
        values = only_field(grib2.scan(raw))['values']
        self.assertEqual((values[0]['lat'], values[0]['lon']), (0.0, 0.0))
        self.assertEqual((values[-1]['lat'], values[-1]['lon']), (1.0, 2.0))

    def test_north_first_grid_definition(self):
        # La1在北、La2在南的网格定义同样映射到j北向递增
        raw = message(sec1(), sec3(2, 2, 1000000, 0, 0, 1000000, 1000000,
                                   1000000, 0x40),
                      sec4(0, 0), sec5(0.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([0, 1, 2, 3], 2)))
        field = only_field(grib2.scan(raw))
        self.assertEqual(flat_values(field), [2.0, 3.0, 0.0, 1.0])
        self.assertEqual([v['lat'] for v in field['values']],
                         [0.0, 0.0, 1.0, 1.0])

    def test_longitude_wrap_keeps_grid_definition(self):
        raw = message(sec1(), sec3(2, 1, 0, 350000000, 0, 10000000, 20000000,
                                   0, 0x40),
                      sec4(0, 0), sec5(0.0, 0, 0, 4, 2), sec6_none(),
                      sec7(pack_bits([5, 7], 4)))
        field = only_field(grib2.scan(raw))
        self.assertEqual(flat_values(field), [5.0, 7.0])
        self.assertEqual([v['lon'] for v in field['values']], [350.0, 10.0])

    def test_bit_ranges_reported(self):
        raw = message(sec1(), grid2x2(), sec4(0, 0), sec5(0.0, 0, 0, 4, 4),
                      sec6_none(), sec7(pack_bits([1, 2, 3, 4], 4)))
        field = only_field(grib2.scan(raw))
        data_offset = field['sections']['data']
        for k, entry in enumerate(field['values']):
            self.assertEqual(entry['bit_range']['section7_bit_offset'], k * 4)
            self.assertEqual(entry['bit_range']['bit_length'], 4)
            self.assertEqual(entry['bit_range']['file_bit_offset'],
                             (data_offset + 5) * 8 + k * 4)


class BitmapTests(unittest.TestCase):
    def test_reuse_within_message(self):
        raw = message(sec1(), grid2x2(),
                      sec4(0, 0), sec5(1.0, 0, 0, 2, 3),
                      sec6_explicit([1, 1, 0, 1]), sec7(pack_bits([0, 1, 2], 2)),
                      sec4(0, 2, forecast=1), sec5(1.0, 0, 0, 2, 3),
                      sec6_reuse(), sec7(pack_bits([3, 2, 1], 2)))
        report = grib2.scan(raw)
        first, second = report['messages'][0]['fields']
        self.assertEqual(first['bitmap']['source'], 'explicit')
        self.assertEqual(flat_values(first), [1.0, 2.0, None, 3.0])
        self.assertEqual(second['bitmap']['source'], 'reused')
        self.assertEqual(second['bitmap']['source_field'], 0)
        self.assertEqual(flat_values(second), [4.0, 3.0, None, 2.0])

    def test_none_indicator_does_not_clear_explicit_bitmap(self):
        raw = message(sec1(), grid2x2(),
                      sec4(0, 0), sec5(1.0, 0, 0, 2, 1),
                      sec6_explicit([1, 0, 0, 0]), sec7(pack_bits([0], 2)),
                      sec4(0, 1), sec5(1.0, 0, 0, 2, 4),
                      sec6_none(), sec7(pack_bits([0, 1, 2, 3], 2)),
                      sec4(0, 2), sec5(1.0, 0, 0, 2, 1),
                      sec6_reuse(), sec7(pack_bits([1], 2)))
        fields = grib2.scan(raw)['messages'][0]['fields']
        self.assertEqual(fields[1]['bitmap']['source'], 'none')
        self.assertEqual(flat_values(fields[1]), [1.0, 2.0, 3.0, 4.0])
        self.assertEqual(fields[2]['bitmap']['source'], 'reused')
        self.assertEqual(fields[2]['bitmap']['source_field'], 0)
        self.assertEqual(flat_values(fields[2]), [2.0, None, None, None])

    def test_reuse_without_prior_explicit_bitmap(self):
        raw = message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 3),
                      sec6_reuse(), sec7(pack_bits([0, 1, 2], 2)))
        with self.assertRaises(GribError):
            grib2.scan(raw)

    def test_reuse_across_messages_is_illegal(self):
        first = message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 3),
                        sec6_explicit([1, 1, 0, 1]), sec7(pack_bits([0, 1, 2], 2)))
        second = message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 3),
                         sec6_reuse(), sec7(pack_bits([0, 1, 2], 2)))
        with self.assertRaises(GribError):
            grib2.scan(first + second)

    def test_reuse_with_incompatible_points(self):
        raw = message(sec1(), grid2x2(),
                      sec4(0, 0), sec5(1.0, 0, 0, 2, 3),
                      sec6_explicit([1, 1, 0, 1]), sec7(pack_bits([0, 1, 2], 2)),
                      grid3x2(), sec4(0, 1), sec5(1.0, 0, 0, 2, 3),
                      sec6_reuse(), sec7(pack_bits([0, 1, 2], 2)))
        with self.assertRaises(GribError):
            grib2.scan(raw)

    def test_bitmap_ones_must_match_packed_count(self):
        raw = message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 4),
                      sec6_explicit([1, 1, 0, 1]),
                      sec7(pack_bits([0, 1, 2, 3], 2)))
        with self.assertRaises(GribError):
            grib2.scan(raw)

    def test_no_bitmap_requires_full_count(self):
        raw = message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 3),
                      sec6_none(), sec7(pack_bits([0, 1, 2], 2)))
        with self.assertRaises(GribError):
            grib2.scan(raw)

    def test_predefined_bitmap_indicator_unsupported(self):
        raw = message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 4),
                      section(6, bytes([1])), sec7(pack_bits([0, 1, 2, 3], 2)))
        with self.assertRaises(GribError):
            grib2.scan(raw)

    def test_bitmap_padding_bits_must_be_zero(self):
        grid = sec3(5, 1, 0, 0, 0, 4000000, 1000000, 0, 0x40)
        body = bytes([0]) + bytes([0b11011001])  # 5个有效位后补齐位非零
        raw = message(sec1(), grid, sec4(0, 0), sec5(1.0, 0, 0, 2, 4),
                      section(6, body), sec7(pack_bits([0, 1, 2, 3], 2)))
        with self.assertRaises(GribError):
            grib2.scan(raw)
        clean = message(sec1(), grid, sec4(0, 0), sec5(1.0, 0, 0, 2, 4),
                        section(6, bytes([0, 0b11011000])),
                        sec7(pack_bits([0, 1, 2, 3], 2)))
        field = only_field(grib2.scan(clean))
        self.assertEqual(flat_values(field), [1.0, 2.0, None, 3.0, 4.0])

    def test_data_padding_bits_must_be_zero(self):
        raw = message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 3),
                      sec6_explicit([1, 1, 0, 1]), sec7(bytes([0b00011001])))
        with self.assertRaises(GribError):
            grib2.scan(raw)


class StructureTests(unittest.TestCase):
    def test_grib_bytes_inside_message_are_not_boundary(self):
        payload = b'local GRIB and GRIB\x00\x00 again'
        raw = message(sec1(), sec2(payload), grid2x2(), sec4(0, 0),
                      sec5(1.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([0, 1, 2, 3], 2)))
        report = grib2.scan(raw)
        self.assertEqual(report['message_count'], 1)
        field = only_field(report)
        self.assertEqual(field['local_use']['sha256'],
                         hashlib.sha256(sec2(payload)).hexdigest())
        self.assertEqual(flat_values(field), [1.0, 2.0, 3.0, 4.0])

    def test_local_use_updates_at_field_boundary(self):
        raw = message(sec1(), sec2(b'alpha'), grid2x2(),
                      sec4(0, 0), sec5(1.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([0, 1, 2, 3], 2)),
                      sec2(b'beta'),
                      sec4(0, 1), sec5(1.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([3, 2, 1, 0], 2)))
        first, second = grib2.scan(raw)['messages'][0]['fields']
        self.assertEqual(first['local_use']['sha256'],
                         hashlib.sha256(sec2(b'alpha')).hexdigest())
        self.assertEqual(second['local_use']['sha256'],
                         hashlib.sha256(sec2(b'beta')).hexdigest())

    def test_broken_message_boundaries(self):
        good = good_message()
        cases = {
            'empty': b'',
            'not_grib': b'not-a-grib-message-padded',
            'trailing_garbage': good + b'\x00',
            'truncated': good[:-1],
            'no_end_mark': good[:-4] + b'XXXX',
            'bad_edition': message(sec1(), grid2x2(), sec4(0, 0),
                                   sec5(1.0, 0, 0, 2, 4), sec6_none(),
                                   sec7(pack_bits([0, 1, 2, 3], 2)), edition=1),
        }
        tiny = bytearray(good)
        tiny[15] = 10
        cases['declared_length_too_small'] = bytes(tiny)
        overrun = bytearray(good)
        overrun[16] = 0x7F
        cases['section_length_overrun'] = bytes(overrun)
        for name, raw in cases.items():
            with self.subTest(case=name):
                with self.assertRaises(GribError):
                    grib2.scan(raw)

    def test_illegal_section_order(self):
        cases = [
            message(sec1(), sec5(1.0, 0, 0, 2, 4)),
            message(sec1(), sec4(0, 0), sec5(1.0, 0, 0, 2, 4), sec6_none(),
                    sec7(pack_bits([0, 1, 2, 3], 2))),
            message(sec1(), grid2x2(), sec4(0, 0), sec6_none()),
            message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 4),
                    sec6_none()),
            message(sec1(), sec1(), grid2x2(), sec4(0, 0),
                    sec5(1.0, 0, 0, 2, 4), sec6_none(),
                    sec7(pack_bits([0, 1, 2, 3], 2))),
            message(grid2x2(), sec1(), sec4(0, 0), sec5(1.0, 0, 0, 2, 4),
                    sec6_none(), sec7(pack_bits([0, 1, 2, 3], 2))),
            message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 4),
                    sec5(1.0, 0, 0, 2, 4)),
            message(sec1(), section(9, b'xx')),
            message(sec1(), grid2x2(), sec4(0, 0), sec5(1.0, 0, 0, 2, 4),
                    sec6_none(), sec7(pack_bits([0, 1, 2, 3], 2)),
                    sec7(pack_bits([0, 1, 2, 3], 2))),
        ]
        for index, raw in enumerate(cases):
            with self.subTest(case=index):
                with self.assertRaises(GribError):
                    grib2.scan(raw)

    def test_max_message_bytes(self):
        raw = good_message()
        self.assertEqual(grib2.scan(raw, max_message_bytes=len(raw))['message_count'], 1)
        with self.assertRaises(GribError):
            grib2.scan(raw, max_message_bytes=len(raw) - 1)


class TemplateTests(unittest.TestCase):
    def test_unknown_packing_template_kept_raw(self):
        unknown5 = sec5(0, 0, 0, 0, 4, template=2)
        raw = message(sec1(), grid2x2(),
                      sec4(0, 0), sec5(1.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([0, 1, 2, 3], 2)),
                      sec4(0, 1), unknown5, sec6_none(), sec7(b'\x00' * 4))
        first, second = grib2.scan(raw)['messages'][0]['fields']
        self.assertEqual(first['status'], 'decoded')
        self.assertEqual(second['status'], 'undecoded')
        self.assertIn('5.2', second['undecoded_reason'])
        self.assertIsNone(second['values'])
        self.assertNotIn('reference_value', second['packing'])
        self.assertEqual(second['packing']['packed_count'], 4)
        self.assertEqual(second['raw_sections_hex']['packing'], unknown5.hex())

    def test_unknown_product_template_kept_raw(self):
        raw = message(sec1(), grid2x2(), sec4(3, 7, template=8),
                      sec5(1.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([0, 1, 2, 3], 2)))
        field = only_field(grib2.scan(raw))
        self.assertEqual(field['status'], 'undecoded')
        self.assertIn('4.8', field['undecoded_reason'])
        self.assertEqual(field['product']['template'], '4.8')
        self.assertEqual(field['identity']['parameter_category'], 3)
        self.assertEqual(field['identity']['parameter_number'], 7)
        self.assertIsNone(field['identity']['forecast_time'])
        self.assertIsNone(field['identity']['level_value'])

    def test_unknown_grid_template_kept_raw(self):
        raw = message(sec1(),
                      sec3(2, 2, 0, 0, 0, 0, 0, 0, 0, template=1),
                      sec4(0, 0), sec5(1.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([0, 1, 2, 3], 2)))
        field = only_field(grib2.scan(raw))
        self.assertEqual(field['status'], 'undecoded')
        self.assertIn('3.1', field['undecoded_reason'])
        self.assertFalse(field['grid']['supported'])
        self.assertEqual(field['grid']['points'], 4)

    def test_accepted_template_lists(self):
        raw = good_message()
        self.assertEqual(only_field(grib2.scan(raw))['status'], 'decoded')
        for kwargs in ({'accepted_grid': ()}, {'accepted_product': ()},
                       {'accepted_packing': ()}):
            with self.subTest(kwargs=kwargs):
                field = only_field(grib2.scan(raw, **kwargs))
                self.assertEqual(field['status'], 'undecoded')
                self.assertIn('未被资料配置接受', field['undecoded_reason'])


class IdentityTests(unittest.TestCase):
    def test_identity_not_keyed_by_parameter_alone(self):
        raw = message(sec1(), grid2x2(),
                      sec4(0, 0, forecast=0), sec5(1.0, 0, 0, 2, 4),
                      sec6_none(), sec7(pack_bits([0, 1, 2, 3], 2)),
                      sec4(0, 0, forecast=6), sec5(1.0, 0, 0, 2, 4),
                      sec6_none(), sec7(pack_bits([0, 1, 2, 3], 2)),
                      sec4(0, 0, forecast=0, surface=100, level=5000),
                      sec5(1.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([0, 1, 2, 3], 2)))
        fields = grib2.scan(raw)['messages'][0]['fields']
        self.assertEqual(len(fields), 3)
        keys = {tuple(sorted(f['identity'].items())) for f in fields}
        self.assertEqual(len(keys), 3)
        self.assertEqual(fields[1]['identity']['forecast_time'], 6)
        self.assertEqual(fields[2]['identity']['surface_type'], 100)
        self.assertEqual(fields[2]['identity']['level_value'], 5000.0)

    def test_grid_redefinition_and_inheritance(self):
        raw = message(sec1(), grid2x2(),
                      sec4(0, 0), sec5(1.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([0, 1, 2, 3], 2)),
                      sec4(0, 1), sec5(1.0, 0, 0, 2, 4), sec6_none(),
                      sec7(pack_bits([0, 1, 2, 3], 2)),
                      grid3x2(),
                      sec4(0, 2), sec5(0.0, 0, 0, 4, 6), sec6_none(),
                      sec7(pack_bits([1, 2, 3, 4, 5, 6], 4)))
        first, second, third = grib2.scan(raw)['messages'][0]['fields']
        self.assertEqual(first['grid']['source'], 'defined')
        self.assertEqual(second['grid']['source'], 'inherited')
        self.assertEqual(second['grid']['ni'], 2)
        self.assertEqual(third['grid']['source'], 'defined')
        self.assertEqual(third['grid']['ni'], 3)
        self.assertNotEqual(first['grid']['section_offset'],
                            third['grid']['section_offset'])


class CliTests(unittest.TestCase):
    def run_cli(self, *argv):
        buffer = io.StringIO()
        with redirect_stdout(buffer):
            cli_main(list(argv))
        return json.loads(buffer.getvalue())

    def test_scan_command(self):
        report = self.run_cli('scan', 'data/normal.json')
        self.assertEqual(report['message_count'], 1)
        self.assertEqual(len(report['messages'][0]['fields']), 2)
        raw = Path('data/fields.grib2').read_bytes()
        self.assertEqual(report['source_sha256'], hashlib.sha256(raw).hexdigest())
        report = self.run_cli('scan', 'data/two-messages.json')
        self.assertEqual(report['message_count'], 2)

    def test_inspect_command_still_works(self):
        report = self.run_cli('inspect', 'data/normal.json')
        self.assertEqual(report['summary']['first_declared_message_bytes'], 248)

    def test_scan_rejects_invalid_signature(self):
        with self.assertRaises(ValueError):
            self.run_cli('scan', 'data/invalid-domain.json')


if __name__ == '__main__':
    unittest.main()
