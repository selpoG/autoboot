#!/usr/bin/env python3
"""Check the release version in metadata and, when supplied, a Git tag.

Only the top-level CFF version scalar is read here; this is not a CFF schema
validator. Keeping this check in the standard library avoids CI dependencies.
"""
import argparse
import json
import os
from pathlib import Path
import re


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--tag', help='release tag, such as v1.0.0')
    args = parser.parse_args()
    root = Path(__file__).resolve().parents[1]
    version = (root / 'VERSION').read_text().strip()
    number = r'(?:0|[1-9][0-9]*)'
    if not re.fullmatch(rf'{number}\.{number}\.{number}', version):
        parser.error('VERSION must contain a stable MAJOR.MINOR.PATCH version')
    cff = re.findall(r'^version:\s*["\x27]?([0-9]+\.[0-9]+\.[0-9]+)["\x27]?\s*$',
                     (root / 'CITATION.cff').read_text(), re.MULTILINE)
    if cff != [version]:
        parser.error('CITATION.cff version must match VERSION exactly once')
    if json.loads((root / '.zenodo.json').read_text()).get('version') != version:
        parser.error('.zenodo.json version must match VERSION')
    tag = args.tag
    if tag is None and os.environ.get('GITHUB_REF_TYPE') == 'tag':
        tag = os.environ.get('GITHUB_REF_NAME')
    if tag is not None and tag != 'v' + version:
        parser.error(f'release tag must be v{version}, got {tag!r}')
    print(f'RELEASE_VERSION={version} PASS')


if __name__ == '__main__':
    main()
