#!/usr/bin/env python3
"""Fail if the release tree contains private paths, secrets, or large binaries."""
from __future__ import annotations
from pathlib import Path
import re
import sys

ROOT = Path(__file__).resolve().parents[1]
EXCLUDE_DIRS = {'.git', '.venv', 'venv', '__pycache__', '.pytest_cache'}
TEXT = {'.md', '.py', '.json', '.txt', '.yml', '.yaml', '.cff', '.toml', '.csv'}
BAD_TEXT = [
    re.compile(r'/Users/[A-Za-z0-9._-]+'),
    re.compile(r'/home/[A-Za-z0-9._-]+'),
    re.compile(r'C:\\\\Users\\\\'),
    re.compile(r'gho_[A-Za-z0-9_]+'),
    re.compile(r'github_pat_[A-Za-z0-9_]+'),
    re.compile(r'sk-[A-Za-z0-9]{16,}'),
    re.compile(r'AKIA[0-9A-Z]{16}'),
    re.compile(r'BEGIN (RSA |EC |OPENSSH )?PRIVATE KEY'),
]
BAD_NAMES = {'.DS_Store'}
BAD_SUFFIXES = {'.pdf', '.zip', '.tar', '.gz', '.7z', '.dmg'}
MAX_BYTES = 8_000_000

errors = []
for p in ROOT.rglob('*'):
    if any(part in EXCLUDE_DIRS for part in p.parts):
        continue
    if p.is_dir():
        continue
    if p.name in BAD_NAMES:
        errors.append(f'bad file name: {p.relative_to(ROOT)}')
    if p.suffix.lower() in BAD_SUFFIXES:
        errors.append(f'bad suffix: {p.relative_to(ROOT)}')
    if p.stat().st_size > MAX_BYTES:
        errors.append(f'large file ({p.stat().st_size} bytes): {p.relative_to(ROOT)}')
    if p.suffix.lower() in TEXT:
        try:
            text = p.read_text(encoding='utf-8', errors='ignore')
        except OSError:
            continue
        for rx in BAD_TEXT:
            if rx.search(text):
                errors.append(f'private/secret pattern {rx.pattern!r}: {p.relative_to(ROOT)}')
                break

print(f'files scanned: {sum(1 for p in ROOT.rglob("*") if p.is_file())}')
print(f'errors: {len(errors)}')
for e in errors:
    print(' -', e)
raise SystemExit(1 if errors else 0)
