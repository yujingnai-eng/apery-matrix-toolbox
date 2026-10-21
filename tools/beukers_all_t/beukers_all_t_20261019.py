#!/usr/bin/env python3
"""All Beukers Appendix C families and all listed t.

For every family/point and every good prime in [5, 200], compute the first
truncated block H_p(t)=sum_{n<p} A_n t^n and its p-adic valuation, and compare
with chi_p = (Delta/p).

This is a finite consistency scan.  It is not a proof of the Hasse/Deuring
bridge, and the higher-order inert cases are recorded as exceptions rather
than smoothed away.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as Fr
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[4]
REG = ROOT / '00_当前主线/外部审稿_20261017/勘误版注册表'
sys.path.insert(0, str(REG))
import appendix_c_registry as R  # noqa: E402


def is_prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d*d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


def vp(x, p):
    if x == 0:
        return None
    k = 0
    while x % p == 0:
        x //= p
        k += 1
    return k


def vp_frac(x, p):
    x = Fr(x)
    return vp(x.numerator, p) - vp(x.denominator, p)


def kron(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p-1)//2, p) == 1 else -1


def family_sequence(fam, p):
    """A_0..A_{p-1} for the family recurrence, or None if nonintegral."""
    a = [1]
    for n in range(1, p):
        num = -(R.poly(fam["p1"], n-1)*a[n-1] +
                (R.poly(fam["p2"], n-2)*a[n-2] if n >= 2 else 0))
        q, rem = divmod(num, n**3)
        if rem:
            return None
        a.append(q)
    return a


def main():
    primes = [p for p in range(5, 201) if is_prime(p)]
    cases = []
    skipped = Counter()
    for fam in R.FAMILIES:
        if not fam["points"]:
            skipped["empty_family"] += 1
            continue
        for p in primes:
            if (fam["level"]*6) % p == 0:
                skipped["bad_reduction"] += len(fam["points"])
                continue
            seq = family_sequence(fam, p)
            if seq is None:
                skipped["nonintegral"] += len(fam["points"])
                continue
            for idx, pt in enumerate(fam["points"], 1):
                t = Fr(pt["t"])
                if t.denominator % p == 0 or t.numerator % p == 0:
                    skipped["non_p_integral_t"] += 1
                    continue
                H = sum(seq[n]*(t**n) for n in range(p))
                v = vp_frac(H, p)
                chi = kron(pt["delta"], p)
                cases.append(dict(
                    section=fam["section"], family=fam.get("name", ""),
                    point=idx, t=str(t), d=pt["d"], delta=pt["delta"],
                    p=p, chi=chi, v=v,
                ))
    by_chi = defaultdict(Counter)
    by_family = defaultdict(lambda: defaultdict(Counter))
    exceptions = []
    for c in cases:
        by_chi[c["chi"]][c["v"]] += 1
        by_family[c["section"]][c["chi"]][c["v"]] += 1
        bad = False
        if c["chi"] == 1 and c["v"] != 0:
            bad = True
        if c["chi"] == 0 and c["v"] != 1:
            bad = True
        if c["chi"] == -1 and (c["v"] is None or c["v"] < 2):
            bad = True
        if bad:
            exceptions.append(c)
    out = dict(
        schema_version=1,
        date='2026-10-19',
        status='finite consistency scan; not a proof',
        scope=dict(
            families=len([f for f in R.FAMILIES if f["points"]]),
            points=sum(len(f["points"]) for f in R.FAMILIES),
            primes=[primes[0], primes[-1]],
            prime_count=len(primes),
            cases=len(cases),
            skipped=dict(skipped),
        ),
        by_chi={str(k): {str(v): n for v, n in sorted(c.items(), key=lambda z: -1 if z[0] is None else z[0])}
                for k, c in sorted(by_chi.items(), key=lambda z: z[0])},
        by_family={sec: {str(ch): {str(v): n for v, n in sorted(c.items(), key=lambda z: -1 if z[0] is None else z[0])}
                         for ch, c in sorted(dd.items(), key=lambda z: z[0])}
                   for sec, dd in sorted(by_family.items())},
        exceptions=exceptions,
        exception_count=len(exceptions),
        higher_inert=[c for c in cases if c["chi"] == -1 and c["v"] is not None and c["v"] >= 3],
        limits=[
            'Finite primes 5..200 and the 118 listed Appendix C points only.',
            'H_p is the first truncated block, not the whole depth tower.',
            'The Hasse/Deuring identification is a proposed bridge, not proved here.',
            'Higher-order inert valuations are preserved as exceptions.',
            'No originality or literature check.',
        ],
        producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    (HERE / 'beukers_all_t_20261019.json').write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + '\n')
    print('cases', len(cases), 'exceptions', len(exceptions))
    print('by chi', dict(out['by_chi']))
    print('skipped', dict(skipped))
    for c in exceptions[:10]:
        print('EXC', c)


if __name__ == '__main__':
    main()
