#!/usr/bin/env python3
"""Validate params/**/*.ini against NA6PRecoParam.h without building the project.

This checks declarations and INI values, not reconstruction behavior or physics.
"""

import argparse
import configparser
import math
from pathlib import Path
import re
import sys


REPOSITORY = Path(__file__).resolve().parents[2]


def read_schema(header):
    source = re.sub(r'/\*.*?\*/|//[^\n]*', '', header.read_text(), flags=re.S)
    constants = dict(re.findall(r'static\s+constexpr\s+int\s+(\w+)\s*=\s*(\d+)\s*;', source))
    schema = {}
    for kind, name, size in re.findall(
            r'^\s*(float|int|bool|std::string)\s+(\w+)\s*(?:\[\s*(\w+)\s*\])?\s*=',
            source, flags=re.M):
        length = int(constants.get(size, size)) if size else None
        schema[name] = (kind, length)
    if not schema:
        raise ValueError(f'No reconstruction parameter declarations found in {header}')
    return schema


def validate_value(kind, value):
    if kind == 'bool':
        if value not in ('true', 'false', '0', '1'):
            raise ValueError('expected true, false, 0, or 1')
    elif kind == 'int':
        if not re.fullmatch(r'[+-]?[0-9]+', value):
            raise ValueError('expected an integer')
        if not -(2**31) <= int(value) < 2**31:
            raise ValueError('integer exceeds the 32-bit range')
    elif kind == 'float':
        if not re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', value):
            raise ValueError('expected a finite decimal number')
        number = float(value)
        if not math.isfinite(number) or abs(number) > 3.4028234663852886e38:
            raise ValueError('number exceeds the finite float range')
    elif not value:
        raise ValueError('expected a nonempty string')


def validate_file(path, schema):
    parser = configparser.ConfigParser(interpolation=None, delimiters=('=',),
                                       inline_comment_prefixes=('#', ';'),
                                       empty_lines_in_values=False)
    parser.optionxform = str  # Parameter names are case sensitive.
    try:
        with path.open() as stream:
            parser.read_file(stream)
    except (configparser.Error, UnicodeError) as error:
        return [str(error)]
    errors = []
    if parser.defaults():
        errors.append('[DEFAULT] values are not supported')
    if parser.sections() != ['reco']:
        errors.append('expected exactly one [reco] section')
    if not parser.has_section('reco'):
        return errors
    if not parser.items('reco'):
        errors.append('[reco] has no parameter assignments')
    for key, value in parser.items('reco'):
        match = re.fullmatch(r'(\w+)(?:\[([0-9]+)\])?', key)
        if not match or match[1] not in schema:
            errors.append(f'{key}: unknown parameter')
            continue
        kind, length = schema[match[1]]
        index = match[2]
        if length is None and index is not None:
            errors.append(f'{key}: scalar parameter cannot have an index')
        elif length is not None and (index is None or int(index) >= length):
            errors.append(f'{key}: expected an array index from 0 to {length - 1}')
        else:
            try:
                validate_value(kind, value)
            except ValueError as error:
                errors.append(f'{key}: {error} (got {value!r})')
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--params', type=Path, default=REPOSITORY / 'params')
    parser.add_argument('--header', type=Path,
                        default=REPOSITORY / 'rec/include/NA6PRecoParam.h')
    args = parser.parse_args()
    schema = read_schema(args.header)
    configs = sorted(args.params.rglob('*.ini'))
    if not configs:
        parser.error(f'No .ini files found in {args.params}')
    failures = 0
    for config in configs:
        errors = validate_file(config, schema)
        if errors:
            failures += 1
            for error in errors:
                print(f'FAIL {config}: {error}', file=sys.stderr)
        else:
            print(f'PASS {config}')
    print(f'{len(configs) - failures}/{len(configs)} parameter files passed.')
    return bool(failures)


if __name__ == '__main__':
    try:
        sys.exit(main())
    except (OSError, ValueError) as error:
        sys.exit(str(error))
