#!/usr/bin/env python3
"""Validate params/**/*.ini against parameter headers without a build.

Checks declarations, array indices, and basic scalar types, not physics choices.
"""

import argparse
import configparser
import math
from pathlib import Path
import re
import sys


REPOSITORY = Path(__file__).resolve().parents[2]
HEADERS = {
    'reco': 'rec/include/NA6PRecoParam.h',
    'layout': 'base/include/NA6PLayoutParam.h',
    'beam': 'base/include/NA6PBeamParam.h',
    'keyval': 'utils/include/KeyValParam.h',
    'gencut': 'sim/include/NA6PGenCutParam.h',
}


def read_schema(header):
    source = re.sub(r'/\*.*?\*/|//[^\n]*', '', header.read_text(), flags=re.S)
    constants = dict(re.findall(r'static\s+constexpr\s+int\s+(\w+)\s*=\s*(\d+)\s*;', source))
    schema = {}
    for kind, name, size in re.findall(
            r'^\s*(float|int|bool|size_t|std::string)\s+(\w+)\s*(?:\[\s*(\w+)\s*\])?\s*=',
            source, flags=re.M):
        length = int(constants.get(size, size)) if size else None
        schema[name] = (kind, length)
    if not schema:
        raise ValueError(f'No parameter declarations found in {header}')
    return schema


def validate_value(kind, value):
    if kind == 'bool':
        if value not in ('true', 'false', '0', '1'):
            raise ValueError('expected true, false, 0, or 1')
    elif kind in ('int', 'size_t'):
        if not re.fullmatch(r'[+-]?[0-9]+', value):
            raise ValueError('expected an integer')
        if kind == 'size_t' and not 0 <= int(value) < 2**64:
            raise ValueError('expected an unsigned 64-bit integer')
        if kind == 'int' and not -(2**31) <= int(value) < 2**31:
            raise ValueError('integer exceeds the 32-bit range')
    elif kind == 'float':
        if not re.fullmatch(r'[+-]?(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?', value):
            raise ValueError('expected a finite decimal number')
        number = float(value)
        if not math.isfinite(number) or abs(number) > 3.4028234663852886e38:
            raise ValueError('number exceeds the finite float range')
    # Strings may be empty, e.g. unused medAbsorberPlug array entries.


def validate_file(path, schemas):
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
    if not parser.sections():
        errors.append('expected a parameter section')
    for section in parser.sections():
        if section not in schemas:
            errors.append(f'[{section}]: unknown parameter section')
            continue
        if not parser.items(section):
            errors.append(f'[{section}] has no parameter assignments')
        for key, value in parser.items(section):
            name = f'{section}.{key}'
            match = re.fullmatch(r'(\w+)(?:\[([^\]]*)\])?', key)
            if not match or match[1] not in schemas[section]:
                errors.append(f'{name}: unknown parameter')
                continue
            kind, length = schemas[section][match[1]]
            index = match[2]
            if length is None and index is not None:
                errors.append(f'{name}: scalar parameter cannot have an index')
            elif length is not None and (index is None or
                                         not re.fullmatch(r'[0-9]+', index) or
                                         int(index) >= length):
                errors.append(f'{name}: expected an array index from 0 to {length - 1}')
            else:
                try:
                    validate_value(kind, value)
                    if section == 'reco' and kind == 'std::string' and not value:
                        raise ValueError('expected a nonempty string')
                except ValueError as error:
                    errors.append(f'{name}: {error} (got {value!r})')
    return errors


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--params', type=Path, default=REPOSITORY / 'params',
                        help='INI file or directory to inspect recursively (default: params/)')
    parser.add_argument('--header-root', type=Path, default=REPOSITORY,
                        help='repository containing the parameter headers')
    parser.add_argument('--header', type=Path,
                        help='override the reconstruction parameter header')
    args = parser.parse_args()
    schemas = {section: read_schema(args.header if section == 'reco' and args.header
                                   else args.header_root / header)
               for section, header in HEADERS.items()}
    configs = [args.params] if args.params.is_file() else sorted(args.params.rglob('*.ini'))
    if not configs:
        parser.error(f'No .ini files found in {args.params}')
    failures = 0
    for config in configs:
        errors = validate_file(config, schemas)
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
