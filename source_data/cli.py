"""显示资料库存与原始字段摘要,并采集GRIB2消息与场。"""
import argparse
import json
from pathlib import Path
from .io import load, asset, integer
from .domain import summarize
from . import grib2


def _templates(value, key):
    items = value.get(key, [])
    if not isinstance(items, list) or any(type(x) is not int or x < 0 or x > 65535 for x in items):
        raise ValueError('模板配置非法')
    return tuple(items)


def scan_config(value, root):
    raw, digest = asset(root, value['file'])
    maximum = value.get('max_message_bytes')
    if maximum is not None:
        integer(maximum, 1)
    report = grib2.scan(raw, max_message_bytes=maximum,
                        accepted_grid=_templates(value, 'accepted_grid_templates'),
                        accepted_product=_templates(value, 'accepted_product_templates'),
                        accepted_packing=_templates(value, 'accepted_packing_templates'))
    report['source_sha256'] = digest
    return report


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument('command', choices=['inspect', 'scan'])
    p.add_argument('file')
    a = p.parse_args(argv)
    value, digest = load(a.file)
    if a.command == 'inspect':
        print(json.dumps({'source_sha256': digest, 'summary': summarize(value, Path(a.file).parent)}, ensure_ascii=False, indent=2))
    else:
        print(json.dumps(scan_config(value, Path(a.file).parent), ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
