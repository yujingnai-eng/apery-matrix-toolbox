#!/usr/bin/env python3
"""Independent check of the all-t first-block scan.

Recomputes the aggregates with a separate loop and Horner evaluation.
Does not import the producer.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from fractions import Fraction as Fr
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REG = HERE
sys.path.insert(0, str(REG))
import appendix_c_registry as R  # noqa: E402

DATA = HERE / 'beukers_all_t_20261019.json'


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


def H_direct(fam, t, p):
    a = [1]
    for n in range(1, p):
        num = -(R.poly(fam["p1"], n-1)*a[n-1] +
                (R.poly(fam["p2"], n-2)*a[n-2] if n >= 2 else 0))
        q, rem = divmod(num, n**3)
        if rem:
            return None
        a.append(q)
    tt = Fr(t)
    s = 0
    for n in range(p-1, -1, -1):
        s = s*tt + a[n]
    return s


def main():
    data = json.loads(DATA.read_text())
    by_chi = defaultdict(Counter)
    by_family = defaultdict(lambda: defaultdict(Counter))
    cases = 0
    exceptions = []
    primes = [p for p in range(5, 201) if is_prime(p)]
    for fam in R.FAMILIES:
        if not fam["points"]:
            continue
        for p in primes:
            if (fam["level"]*6) % p == 0:
                continue
            for idx, pt in enumerate(fam["points"], 1):
                t = Fr(pt["t"])
                if t.denominator % p == 0 or t.numerator % p == 0:
                    continue
                h = H_direct(fam, t, p)
                if h is None:
                    continue
                v = vp_frac(h, p)
                chi = kron(pt["delta"], p)
                cases += 1
                by_chi[chi][v] += 1
                by_family[fam["section"]][chi][v] += 1
                bad = ((chi == 1 and v != 0) or (chi == 0 and v != 1) or
                       (chi == -1 and (v is None or v < 2)))
                if bad:
                    exceptions.append(dict(section=fam["section"], point=idx,
                                           t=str(t), p=p, chi=chi, v=v))
    recomputed_chi = {str(k): {str(v): n for v, n in c.items()}
                      for k, c in by_chi.items()}
    recomputed_family = {sec: {str(ch): {str(v): n for v, n in c.items()}
                               for ch, c in dd.items()}
                         for sec, dd in by_family.items()}
    assert cases == data["scope"]["cases"], (cases, data["scope"]["cases"])
    assert recomputed_chi == data["by_chi"], (recomputed_chi, data["by_chi"])
    assert recomputed_family == data["by_family"]
    assert exceptions == data["exceptions"]
    assert data["exception_count"] == 0
    out = dict(
        status='pass',
        cases=cases,
        by_chi=recomputed_chi,
        exceptions=len(exceptions),
        hashes={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (HERE/'beukers_all_t_20261019.py', DATA, Path(__file__))},
        limits=['Same finite primes and registry points; no all-t/all-p proof.'],
    )
    (HERE / 'verify_beukers_all_t_20261019.json').write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
