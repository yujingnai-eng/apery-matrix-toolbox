#!/usr/bin/env python3
"""Independent checks for the deep probes. No producer import.

Covers the finite claims P01, P02, P03 (A1=5), P06 (r=3,4), P07 and P08.
P04, P05 and P09 are exploratory and are not independently certified here.
"""
from fractions import Fraction as Q
from math import comb
import hashlib
import itertools
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
DATA = HERE / 'deep_probes_20261019.json'
CHECKS = {}


def require(condition, group):
    if not condition:
        raise AssertionError(group)
    CHECKS[group] = CHECKS.get(group, 0) + 1


def v2(x):
    x = Q(x)
    if x == 0:
        return None
    n, d = x.numerator, x.denominator
    c = 0
    while n % 2 == 0:
        n //= 2
        c += 1
    while d % 2 == 0:
        d //= 2
        c -= 1
    return c


def source(family, n):
    if family == 'Apery':
        return sum((comb(n, k)*comb(n+k, k))**2 for k in range(n+1))
    p = 3 if family == 'Franel3' else 4
    return sum(comb(n, k)**p for k in range(n+1))


def block(A, a, s, signed=True):
    lo, hi = a*(1 << s), (a+1)*(1 << s)
    if signed:
        return sum(((-1)**n)*A[n] for n in range(lo, hi))
    return sum(A[n] for n in range(lo, hi))


def verify01(e):
    A = [source('Apery', n) for n in range(129)]
    for s in range(1, 7):
        bs = 1 << s
        for a in range(64//bs):
            require(v2(block(A, a, s)) == 3*s-1, 'P01_signed_block')
            require(v2(block(A, a, s, False)) == s, 'P01_unsigned_block')
    for r in range(1, 7):
        f = block(A, 0, r)
        b = block(A, 1, r)
        require(v2(f) == 3*r-1, 'P01_origin_f')
        require(v2(b) == 3*r-1, 'P01_origin_b')
        require(v2(b-7*f) == 4*r+2, 'P01_origin_congruence')
    for family in ('Franel3', 'Franel4'):
        B = [source(family, n) for n in range(65)]
        for s in range(1, 5):
            bs = 1 << s
            vv = {v2(block(B, a, s)) for a in range(64//bs)}
            require(len(vv) > 1, 'P01_franel_not_uniform')


def verify02(e):
    for n in range(33):
        hs = [comb(n, k)**2*comb(n+k, k)**2 for k in range(n+1)]
        hm = [Q(hs[k]-hs[n-k], 2) for k in range(n+1)]
        require(sum(hm) == 0, 'P02_total_hm_zero')
    A = [source('Apery', n) for n in range(65)]
    for s in range(1, 5):
        bs = 1 << s
        for a in range(64//bs):
            shm = sum(((-1)**n)*sum(Q(comb(n, k)**2*comb(n+k, k)**2 -
                                      comb(n, n-k)**2*comb(n+n-k, n-k)**2, 2)
                                    for k in range(n+1))
                      for n in range(a*bs, (a+1)*bs))
            require(shm == 0, 'P02_block_hm_zero')
            shp = sum(((-1)**n)*sum(Q(comb(n, k)**2*comb(n+k, k)**2 +
                                      comb(n, n-k)**2*comb(n+n-k, n-k)**2, 2)
                                    for k in range(n+1))
                      for n in range(a*bs, (a+1)*bs))
            require(shp == block(A, a, s), 'P02_block_hp_matches')


def verify03(e):
    A = [source('Apery', n) for n in range(65)]
    vals = []
    for r in range(2, 7):
        n = 1 << r
        P = sum(A[0:n:2])
        O = sum(A[1:n:2])
        vals.append((v2(P), v2(O), v2(P-O)))
    require([x[0] for x in vals] == [1, 2, 3, 4, 5], 'P03_special_vP')
    require([x[2] for x in vals] == [5, 8, 11, 14, 17], 'P03_special_vPO')
    require(all(x[0] == x[1] for x in vals), 'P03_special_pairing')


def verify06(e):
    A = [source('Apery', n) for n in range(64)]
    Fm, Fp = [0], [0]
    for n in range(1, 65):
        Fm.append(Fm[-1] + ((-1)**(n-1))*A[n-1])
        Fp.append(Fp[-1] + A[n-1])
    for r in (3, 4):
        n = 1 << r
        rows = [[(-1)**i*A[i], ((-1)**(i-1)*A[i-1] if i else 0), Fm[i], Fp[i]]
                for i in range(n)]
        g = 0
        for idx in itertools.combinations(range(n), 4):
            d = int(sp.Matrix([rows[i] for i in idx]).det())
            g = sp.gcd(g, abs(d))
        require(int(g) == 10, 'P06_lattice_index_10')


def verify07(e):
    A = [source('Apery', n) for n in range(129)]
    F = [0]
    for n in range(1, 129):
        F.append(F[-1] + ((-1)**(n-1))*A[n-1])
    V = {r: sp.Matrix([sp.Rational(F[1 << (r+1)], 1 << (3*r+2)),
                       sp.Rational(F[1 << r], 1 << (3*r-1))]) for r in range(1, 7)}

    def mod128(q):
        q = sp.Rational(q)
        n, d = q.as_numer_denom()
        return int(n % 128)*pow(int(d % 128), -1, 128) % 128

    expected = {1: [[83, 110], [1, 0]], 2: [[35, 94], [1, 0]],
                3: [[3, 126], [1, 0]], 4: [[3, 126], [1, 0]]}
    for r in range(1, 5):
        M = sp.Matrix.hstack(V[r], V[r+1])
        Nm = sp.Matrix.hstack(V[r+1], V[r+2])
        D = Nm*M.inv()
        got = [[mod128(D[i, j]) for j in range(2)] for i in range(2)]
        require(got == expected[r], 'P07_transfer_shadow')


def verify08(e):
    for n in range(1, 33):
        Hp = lambda m: sum(k*comb(m, k)**2*comb(m+k, k)**2 for k in range(m+1))
        P = 34*n**3 + 51*n**2 + 27*n + 5
        R = (n+1)**3*Hp(n+1) - P*Hp(n) + n**3*Hp(n-1)
        require(v2(R) == 2, 'P08_jet_v2_flat')


def main():
    data = json.loads(DATA.read_text())
    source_file = HERE / 'deep_probes_20261019.py'
    require(hashlib.sha256(source_file.read_bytes()).hexdigest() ==
            data['producer_sha256'], 'producer_hash')
    verify01(data['probes']['P01_dyadic'])
    verify02(data['probes']['P02_kernel_split'])
    verify03(data['probes']['P03_initial_scan'])
    verify06(data['probes']['P06_lattice'])
    verify07(data['probes']['P07_transfer'])
    verify08(data['probes']['P08_jet'])
    out = dict(
        status='pass within declared scope',
        covered=['P01(s<=6,63 blocks)', 'P02(n<=32,s<=4)', 'P03(A1=5)',
                 'P06(r=3,4)', 'P07(r<=4)', 'P08(n<=32)'],
        not_independently_certified=['P04 Mahler', 'P05 2-kernel', 'P09 spectrum'],
        check_counts=CHECKS,
        limits=['Finite ranges only; no all-r proof of G1/G2.',
                'P04/P05/P09 are exploratory and not certified here.',
                'Producer ranges are larger than verifier ranges; see the report boundary section.',
                'No originality or literature check.'],
        hashes={p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in (source_file, DATA, Path(__file__))},
    )
    (HERE / 'verify_deep_probes_20261019.json').write_text(
        json.dumps(out, ensure_ascii=False, indent=2) + '\n')
    print(json.dumps(out, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    main()
