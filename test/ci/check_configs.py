#!/usr/bin/env python3
"""Validate config/**/*.ini against parameter headers without a build.

Checks that every section present explicitly assigns all its parameters, including
all array elements, and validates the names and values of those assignments.
Discovers parameter classes in the source headers and reads their INI section names
from NA6PParamDef.
"""

import argparse
import configparser
import math
import os
from pathlib import Path
import re
import sys


REPOSITORY = Path(__file__).resolve().parents[2]
CONFIG_DIR = REPOSITORY / 'config'
HEADER_SUFFIXES = {'.h', '.hh', '.hpp', '.hxx'}
CPP_LITERAL = r'"(?:\\.|[^"\\])*"' r"|'(?:\\.|[^'\\])*'"
CPP_COMMENTS_AND_LITERALS = re.compile(r'/\*.*?\*/|//[^\n]*|' + CPP_LITERAL, re.S)
CPP_BRACES_AND_LITERALS = re.compile(CPP_LITERAL + r'|[{}]', re.S)
CLASS_DECLARATION = re.compile(
    r'\b(?:class|struct)\s+(?P<name>\w+)\s*(?:final\s*)?:\s*(?P<bases>[^;{}]*)\{')
PARAMETER_BASE = re.compile(
    r'(?:^|,)\s*(?:(?:public|protected|private|virtual)\s+)*'
    r'(?:::\s*)?na6p\s*::\s*conf\s*::\s*ConfigurableParamHelper\s*'
    r'<\s*(?:::\s*)?(?:\w+\s*::\s*)*(?P<name>\w+)\s*>\s*(?=,|$)')
SECTION_DECLARATION = re.compile(
    r'\bNA6PParamDef\s*\(\s*(?P<name>\w+)\s*,\s*"(?P<section>[^"\n]+)"\s*\)')


# Read a header and remove C++ comments while preserving quoted values.
def read_header(header):
    return CPP_COMMENTS_AND_LITERALS.sub(
        lambda match: ' ' if match[0].startswith('/') else match[0], header.read_text())


# Find the closing brace of a class, accounting for nested braces and literals.
def class_body(source, start, header):
    depth = 1
    for token in CPP_BRACES_AND_LITERALS.finditer(source, start):
        if token[0] == '{':
            depth += 1
        elif token[0] == '}':
            depth -= 1
            if depth == 0:
                return source[start:token.start()]
    raise ValueError(f'Unclosed configurable parameter class in {header}')


# Read configurable member types and array sizes from the C++ declarations.
def read_parameter_specs(source, header):
    constants = dict(re.findall(r'static\s+constexpr\s+int\s+(\w+)\s*=\s*(\d+)\s*;', source))
    parameter_specs = {}
    for kind, name, size in re.findall(
            r'^\s*(float|int|bool|size_t|std::string)\s+(\w+)\s*(?:\[\s*(\w+)\s*\])?\s*=',
            source, flags=re.M):
        length = int(constants.get(size, size)) if size else None
        parameter_specs[name] = (kind, length)
    if not parameter_specs:
        raise ValueError(f'No parameter declarations found in {header}')
    return parameter_specs


# Find configurable classes and their parameter specs by INI section.
def discover_parameter_specs(root, reco_header=None):
    parameter_specs_by_section = {}
    origins = {}
    for directory, subdirectories, filenames in os.walk(root):
        if Path(directory) != root and 'CMakeCache.txt' in filenames:
            subdirectories.clear()
            continue
        subdirectories[:] = sorted(name for name in subdirectories
                                   if not name.startswith('.') and name not in {'build', 'install'})
        for filename in sorted(filenames):
            header = Path(directory) / filename
            if header.suffix not in HEADER_SUFFIXES:
                continue
            source = read_header(header)
            for declaration in CLASS_DECLARATION.finditer(source):
                base = PARAMETER_BASE.search(declaration['bases'])
                name = declaration['name']
                if base is None or base['name'] != name:
                    continue
                body = class_body(source, declaration.end(), header)
                registration = next((match for match in SECTION_DECLARATION.finditer(body)
                                     if match['name'] == name), None)
                if registration is None:
                    raise ValueError(f'{header}: {name} is missing its NA6PParamDef declaration')
                section = registration['section']
                if section in origins:
                    raise ValueError(f'{header}: duplicate section [{section}], also in {origins[section]}')
                origins[section] = header
                if section != 'reco' or reco_header is None:
                    parameter_specs_by_section[section] = read_parameter_specs(body, header)
    if not origins:
        raise ValueError(f'No ConfigurableParamHelper classes found under {root}')
    if reco_header is not None:
        parameter_specs_by_section['reco'] = read_parameter_specs(read_header(reco_header), reco_header)
    return parameter_specs_by_section


# Check that an INI value fits the type declared in C++.
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


# Check parameter assignments for each section present in one INI file.
def validate_file(path, parameter_specs_by_section):
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
        if section not in parameter_specs_by_section:
            errors.append(f'[{section}]: unknown parameter section')
            continue
        parameter_specs = parameter_specs_by_section[section]
        options = set(parser.options(section))
        for key, (_, length) in parameter_specs.items():
            expected = (key,) if length is None else (f'{key}[{i}]' for i in range(length))
            for name in expected:
                if name not in options:
                    errors.append(f'{section}.{name}: missing explicit assignment')
        for key, value in parser.items(section):
            name = f'{section}.{key}'
            match = re.fullmatch(r'(\w+)(?:\[([^\]]*)\])?', key)
            if not match or match[1] not in parameter_specs:
                errors.append(f'{name}: unknown parameter')
                continue
            kind, length = parameter_specs[match[1]]
            index = match[2]
            if length is None and index is not None:
                errors.append(f'{name}: scalar parameter cannot have an index')
            elif length is not None and (index is None or
                                         not re.fullmatch(r'0|[1-9][0-9]*', index) or
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


# Discover parameter specs and validate the selected INI files.
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--params', type=Path, default=CONFIG_DIR,
                        help='INI file or directory to inspect recursively (default: config/)')
    parser.add_argument('--header-root', type=Path, default=REPOSITORY,
                        help='repository containing the parameter headers')
    parser.add_argument('--header', type=Path,
                        help='override the reconstruction parameter header')
    args = parser.parse_args()
    parameter_specs_by_section = discover_parameter_specs(args.header_root, args.header)
    configs = [args.params] if args.params.is_file() else sorted(args.params.rglob('*.ini'))
    if not configs:
        parser.error(f'No .ini files found in {args.params}')
    failures = 0
    for config in configs:
        errors = validate_file(config, parameter_specs_by_section)
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
