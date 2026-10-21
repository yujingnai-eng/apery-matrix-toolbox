#!/usr/bin/env python3
"""Run the public toolbox verifiers. Allowlist only; no shell."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent
MANIFEST = HERE / 'manifest.json'


def dig(data, dotted):
    cur = data
    for part in dotted.split('.'):
        if not isinstance(cur, dict) or part not in cur:
            return None
        cur = cur[part]
    return cur


def run_one(tool, python):
    cwd = HERE / tool['cwd']
    cmd = [python if x == '{python}' else x for x in tool['command']]
    try:
        proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True,
                              timeout=tool.get('timeout', 300))
        out = cwd / tool['output']
        data = json.loads(out.read_text()) if out.exists() else None
        problems = []
        if proc.returncode != 0:
            problems.append(f'exit {proc.returncode}')
        if data is None:
            problems.append(f'missing {tool["output"]}')
        elif data.get('status') != tool['status']:
            problems.append(f'status={data.get("status")!r} != {tool["status"]!r}')
        if data is not None:
            for path, want in tool.get('counts', {}).items():
                got = dig(data, path)
                if got != want:
                    problems.append(f'{path}={got!r} != {want!r}')
        return dict(id=tool['id'], status='pass' if not problems else 'fail',
                    exit=proc.returncode, problems=problems,
                    output=str(out.relative_to(HERE)) if out.exists() else None,
                    stdout=proc.stdout[-4000:], stderr=proc.stderr[-4000:])
    except subprocess.TimeoutExpired:
        return dict(id=tool['id'], status='fail', exit=None,
                    problems=[f'timeout after {tool.get("timeout", 300)}s'])
    except OSError as exc:
        return dict(id=tool['id'], status='fail', exit=None,
                    problems=[f'OSError: {exc}'])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--list', action='store_true')
    ap.add_argument('--run')
    ap.add_argument('--all', action='store_true')
    ap.add_argument('--python', default=sys.executable)
    args = ap.parse_args()
    manifest = json.loads(MANIFEST.read_text())
    tools = manifest['tools']
    if args.list:
        for t in tools:
            print(f"{t['id']:14} {t['description']}")
        return 0
    if args.all:
        selected = tools
    elif args.run:
        selected = [t for t in tools if t['id'] == args.run]
        if not selected:
            print(f'unknown tool: {args.run}', file=sys.stderr)
            return 2
    else:
        ap.error('use --list, --run ID, or --all')
    results = [run_one(t, args.python) for t in selected]
    print(f"{'id':14} {'status':6} {'exit':>4}  detail")
    print('-' * 60)
    for r in results:
        detail = r['problems'][0] if r['problems'] else (r.get('output') or '')
        print(f"{r['id']:14} {r['status']:6} {str(r['exit']):>4}  {detail}")
    overall = 'pass' if all(r['status'] == 'pass' for r in results) else 'fail'
    print(json.dumps({'overall': overall, 'results': results}, ensure_ascii=False))
    return 0 if overall == 'pass' else 1


if __name__ == '__main__':
    raise SystemExit(main())
