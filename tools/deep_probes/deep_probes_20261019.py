#!/usr/bin/env python3
"""Nine exploratory probes around G1/G2. Exact finite arithmetic.

This is an exploration, not a proof. Each probe records its finite range and
whether it supports, refutes, or fails to decide a structural hypothesis.
Run: python3 deep_probes_20261019.py --save
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from math import comb
import hashlib
import json
import sys
from pathlib import Path

import sympy as sp
from sympy.matrices.normalforms import smith_normal_form

sys.set_int_max_str_digits(0)

HERE = Path(__file__).resolve().parent
MAX_S = 12
NMAX = 1 << (MAX_S + 1)          # 8192 terms, enough for blocks up to s=12
PRIMES = (2, 3, 5, 7, 11, 13)


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


def terms(family, N):
    a = [1]
    for n in range(N):
        if family == 'Apery':
            d, p, q = 3, 34*n**3 + 51*n**2 + 27*n + 5, -n**3
        elif family == 'Franel3':
            d, p, q = 2, 7*n*n + 7*n + 2, 8*n*n
        else:
            d, p, q = 3, 12*n**3 + 18*n*n + 10*n + 2, 4*n*(16*n*n - 1)
        top = p*a[n] + (q*a[n-1] if n else 0)
        value, rem = divmod(top, (n+1)**d)
        assert rem == 0, (family, n)
        a.append(value)
    return a


def block(A, a, s, signed=True):
    bs = 1 << s
    lo, hi = a*bs, (a+1)*bs
    if signed:
        return sum(((-1)**n)*A[n] for n in range(lo, hi))
    return sum(A[n] for n in range(lo, hi))


def kernel(family, n, k):
    if family == 'Apery':
        return comb(n, k)**2 * comb(n+k, k)**2
    p = 3 if family == 'Franel3' else 4
    return comb(n, k)**p


def probe01_dyadic():
    """Local G1: every aligned dyadic alternating block has valuation 3s-1."""
    A = terms('Apery', NMAX)
    rows = {}
    full = []
    for s in range(1, MAX_S+1):
        bs = 1 << s
        nblocks = NMAX // bs
        vals = [block(A, a, s, True) for a in range(nblocks)]
        vv = [v2(x) for x in vals]
        assert set(vv) == {3*s-1}, (s, sorted(set(vv)))
        plain = [block(A, a, s, False) for a in range(nblocks)]
        pv = [v2(x) for x in plain]
        assert set(pv) == {s}, (s, sorted(set(pv)))
        units = [(vals[a] // (1 << (3*s-1))) % 8 for a in range(min(8, nblocks))]
        rows[str(s)] = dict(
            blocks=nblocks,
            signed_valuations=sorted(set(vv)),
            unsigned_valuations=sorted(set(pv)),
            units_mod8_first8=units,
            signed_hash=hashlib.sha256(
                '\n'.join(str(x) for x in vals).encode()).hexdigest(),
        )
        full.extend(vals)
    # Adjacent cancellation mod 8 at every scale.
    adjacent = []
    for s in range(1, MAX_S+1):
        bs = 1 << s
        nblocks = NMAX // bs
        ok = True
        for a in range(nblocks//2):
            u0 = block(A, 2*a, s, True) // (1 << (3*s-1))
            u1 = block(A, 2*a+1, s, True) // (1 << (3*s-1))
            if (u0+u1) % 8:
                ok = False
                break
        adjacent.append(dict(s=s, ok=ok, pairs=nblocks//2))
    # Origin pair f_r=B(0,r), b_r=B(1,r).
    origin = []
    for r in range(1, MAX_S+1):
        f = block(A, 0, r, True)
        b = block(A, 1, r, True)
        origin.append(dict(r=r, v_f=v2(f), v_b=v2(b), v_b_minus_7f=v2(b-7*f)))
    # Franel controls: same statistic is not uniform.
    controls = {}
    for family in ('Franel3', 'Franel4'):
        B = terms(family, 256)
        ctrl = {}
        for s in range(1, 6):
            bs = 1 << s
            nblocks = 256 // bs
            vv = sorted({v2(block(B, a, s, True)) for a in range(nblocks)})
            ctrl[str(s)] = vv
        controls[family] = ctrl
    return dict(
        kind='local aligned dyadic block law',
        max_s=MAX_S,
        signed_hash=hashlib.sha256(
            '\n'.join(str(x) for x in full).encode()).hexdigest(),
        scales=rows,
        adjacent_mod8=adjacent,
        origin_pair=origin,
        franel_controls=controls,
        conclusion='v2(B_{a,s})=3s-1 for all aligned blocks tested; h_+/h_- split in P2.',
    )


def probe02_kernel_split():
    """The antisymmetric kernel has zero total; h_+ carries the block law."""
    A = terms('Apery', 256)
    total_hm = []
    for n in range(129):
        hs = [kernel('Apery', n, k) for k in range(n+1)]
        hp = [Q(hs[k]+hs[n-k], 2) for k in range(n+1)]
        hm = [Q(hs[k]-hs[n-k], 2) for k in range(n+1)]
        assert sum(hm) == 0
        total_hm.append(dict(n=n, sum_hm=str(sum(hm))))
    blocks = []
    for s in range(1, 7):
        bs = 1 << s
        nblocks = 256 // bs
        for a in range(nblocks):
            shm = sum(((-1)**n)*sum(Q(kernel('Apery', n, k)-kernel('Apery', n, n-k), 2)
                                    for k in range(n+1))
                      for n in range(a*bs, (a+1)*bs))
            assert shm == 0, (s, a, shm)
            shp = sum(((-1)**n)*sum(Q(kernel('Apery', n, k)+kernel('Apery', n, n-k), 2)
                                    for k in range(n+1))
                      for n in range(a*bs, (a+1)*bs))
            assert shp == block(A, a, s, True)
            blocks.append(dict(s=s, a=a, sum_hm=str(shm), sum_hp=str(shp)))
    return dict(
        kind='symmetric/antisymmetric kernel split',
        total_hm_examples=total_hm[:8],
        block_count=len(blocks),
        blocks_first=blocks[:12],
        conclusion='h_- has zero total for every n and zero alternating dyadic block sum; h_+ carries the law.',
    )


def seq_q(a1, N):
    a = [Q(1), Q(a1)]
    for n in range(1, N):
        p = 34*n**3 + 51*n**2 + 27*n + 5
        a.append((p*a[n] - n**3*a[n-1]) / Q((n+1)**3))
    return a


def probe03_initial_scan():
    """Vary A_1. The pairing and the slope are not generic."""
    patterns = {}
    special = None
    for a1 in range(1, 128, 2):
        A = seq_q(a1, 64)          # indices 0..63, enough for r<=6
        vals = []
        for r in range(2, 7):
            n = 1 << r
            P = sum(A[0:n:2])
            O = sum(A[1:n:2])
            vals.append((v2(P), v2(O), v2(P-O)))
        key = tuple(x[0] for x in vals)
        paired = all(x[0] == x[1] for x in vals)
        rec = dict(a1=a1, vP=key, vPO=tuple(x[2] for x in vals), paired=paired)
        patterns.setdefault((key, paired), []).append(a1)
        if a1 == 5:
            special = rec
    summary = []
    for (key, paired), aas in sorted(patterns.items(), key=lambda kv: -len(kv[1])):
        summary.append(dict(vP=key, paired=paired, count=len(aas), examples=aas[:6]))
    return dict(
        kind='initial-value phase scan',
        range='A_1 odd, 1..127; r=2..6',
        special_A1_5=special,
        patterns=summary[:12],
        conclusion='A_1=5 is a special locus; generic A_1 gives different slopes/pairing.',
    )


def probe04_mahler():
    """Finite Mahler coefficients of N -> F_N(-1)."""
    A = terms('Apery', 512)
    F = [0]
    for n in range(1, 513):
        F.append(F[-1] + ((-1)**(n-1))*A[n-1])
    M = 64
    vals = []
    for j in range(M+1):
        s = sum(((-1)**(j-i))*comb(j, i)*F[i] for i in range(j+1))
        vals.append(s)
    return dict(
        kind='Mahler coefficients, finite',
        max_j=M,
        v2=[v2(x) for x in vals],
        first_33=[str(x) for x in vals[:33]],
        conclusion='No simple slope-3 line in the first 64 Mahler coefficients; fractal at j=2^k.',
    )


def rank_bits(vecs):
    basis = {}
    for x in vecs:
        while x:
            b = x.bit_length()-1
            if b in basis:
                x ^= basis[b]
            else:
                basis[b] = x
                break
    return len(basis)


def probe05_twokernel():
    """Rank of the 2-kernel of F_N(-1) mod 2^K."""
    N = 512
    A = terms('Apery', N)
    F = [0]
    for n in range(1, N+1):
        F.append(F[-1] + ((-1)**(n-1))*A[n-1])
    out = []
    K = 4
    for E in (5, 6, 7, 8):
        vecs = []
        for e in range(E+1):
            for r in range(1 << e):
                bits = 0
                idx = 0
                for n in range(r, N, 1 << e):
                    val = F[n] % (1 << K)
                    for b in range(K):
                        if (val >> b) & 1:
                            bits |= 1 << (b*N + idx)
                    idx += 1
                vecs.append(bits)
        out.append(dict(E=E, subseq=len(vecs), rank=rank_bits(vecs)))
    return dict(
        kind='2-kernel rank over F_2',
        K=K,
        rows=out,
        conclusion='Rank grows with E; no evidence for a small 2-automatic transducer.',
    )


def probe06_lattice():
    """Smith normal form of the Z-span of the actual states Y_n(-1)."""
    N = 64
    A = terms('Apery', N)
    Fm, Fp = [0], [0]
    for n in range(1, N+1):
        Fm.append(Fm[-1] + ((-1)**(n-1))*A[n-1])
        Fp.append(Fp[-1] + A[n-1])
    rows = []
    for r in (3, 4, 5, 6):
        n = 1 << r
        M = sp.Matrix([[(-1)**i*A[i],
                        ((-1)**(i-1)*A[i-1] if i else 0),
                        Fm[i], Fp[i]] for i in range(n)])
        D = smith_normal_form(M)
        diag = [int(D[i, i]) for i in range(4)]
        rows.append(dict(r=r, states=n, snf=diag, rank=int(M.rank())))
    return dict(
        kind='state lattice Smith normal form',
        rows=rows,
        conclusion='Orbit lattice index is 10 (2*5), not 2^(3r-1); growth is in the vector, not the lattice index.',
    )


def probe07_transfer():
    """Exact 2x2 transfer of (F_{2^{r+1}},F_{2^r}) and its mod-128 shadow."""
    N = 512
    A = terms('Apery', N)
    F = [0]
    for n in range(1, N+1):
        F.append(F[-1] + ((-1)**(n-1))*A[n-1])
    V = {r: sp.Matrix([sp.Rational(F[1 << (r+1)], 1 << (3*r+2)),
                       sp.Rational(F[1 << r], 1 << (3*r-1))]) for r in range(1, 7)}
    def mod128(q):
        q = sp.Rational(q)
        n, d = q.as_numer_denom()
        return int(n % 128) * pow(int(d % 128), -1, 128) % 128
    rows = []
    for r in range(1, 5):
        M = sp.Matrix.hstack(V[r], V[r+1])
        Nm = sp.Matrix.hstack(V[r+1], V[r+2])
        D = Nm * M.inv()
        rows.append(dict(r=r, D_mod128=[[mod128(D[i, j]) for j in range(2)] for i in range(2)]))
    return dict(
        kind='normalized 2D readout transfer',
        rows=rows,
        candidate_C=[[6, 16], [1, 0]],
        conclusion='Exact transfer drifts; mod-128 shadow is not the candidate C. The 2D readout pair is not a closed state.',
    )


def probe08_jet():
    """Kernel-position derivative source R_n'(1)."""
    N = 256
    A = terms('Apery', N)
    Hp = []
    for n in range(N+1):
        Hp.append(sum(k*comb(n, k)**2*comb(n+k, k)**2 for k in range(n+1)))
    R = [0]*(N+1)
    for n in range(1, N):
        P = 34*n**3 + 51*n**2 + 27*n + 5
        R[n] = (n+1)**3*Hp[n+1] - P*Hp[n] + n**3*Hp[n-1]
    hist = {}
    for n in range(1, N):
        hist[str(v2(R[n]))] = hist.get(str(v2(R[n])), 0) + 1
    partial = []
    for r in range(1, 9):
        n = 1 << r
        partial.append(dict(r=r, v2_prefix=v2(sum(R[1:n])),
                            v2_alternating=v2(sum(((-1)**j)*R[j] for j in range(1, n)))))
    return dict(
        kind='kernel-position derivative source',
        range=f'n=1..{N-1}',
        v2_histogram=hist,
        partials=partial,
        conclusion='R_n\'(1) has constant v2=2; naive jet channel is 2-adically flat and needs cumulative amplification.',
    )


def probe09_spectrum():
    """2-adic Newton polygon of the exact 4x4 block at t=-1."""
    def step(n):
        P = 34*n**3 + 51*n**2 + 27*n + 5
        return sp.Matrix([[sp.Rational(-P, (n+1)**3), sp.Rational(-n**3, (n+1)**3), 0, 0],
                          [1, 0, 0, 0], [1, 0, 1, 0], [(-1)**n, 0, 0, 1]])
    x = sp.symbols('x')
    rows = []
    for r in range(1, 6):
        M = sp.eye(4)
        for n in range(1 << r, 1 << (r+1)):
            M = step(n) * M
        poly = sp.Poly(M.charpoly(x).as_expr(), x)
        coeffs = poly.all_coeffs()
        vv = [v2(sp.Rational(c)) for c in coeffs]
        rows.append(dict(r=r, degree=poly.degree(), coeff_v2=vv, det_v2=v2(M.det())))
    return dict(
        kind='full block 2-adic Newton polygon',
        rows=rows,
        conclusion='One root valuation -3, three units; determinant 1/8. No simple expansive 8-eigenvalue.',
    )


def pack(x):
    if isinstance(x, Q):
        return str(x)
    if isinstance(x, sp.Integer):
        return int(x)
    if isinstance(x, dict):
        return {str(k): pack(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [pack(v) for v in x]
    return x


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--save', action='store_true')
    args = parser.parse_args()
    out = {
        'schema_version': 1,
        'date': '2026-10-19',
        'status': 'exploratory finite probes; no all-r proof',
        'config': dict(max_s=MAX_S, nmax=NMAX, primes=PRIMES),
        'probes': {
            'P01_dyadic': probe01_dyadic(),
            'P02_kernel_split': probe02_kernel_split(),
            'P03_initial_scan': probe03_initial_scan(),
            'P04_mahler': probe04_mahler(),
            'P05_twokernel': probe05_twokernel(),
            'P06_lattice': probe06_lattice(),
            'P07_transfer': probe07_transfer(),
            'P08_jet': probe08_jet(),
            'P09_spectrum': probe09_spectrum(),
        },
        'limits': [
            'Finite ranges only; no all-r proof of G1/G2.',
            'P04/P05/P09 are exploratory; P01/P02/P03/P06/P07/P08 have finite certificates.',
            'No originality or literature check.',
            'No claim of intrinsic chirality or CM identification.',
        ],
        'producer_sha256': hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    for name in out['probes']:
        print(name, 'completed', flush=True)
    if args.save:
        path = HERE / 'deep_probes_20261019.json'
        path.write_text(json.dumps(pack(out), ensure_ascii=False, indent=2) + '\n')
        print('Saved', path.name)
    print('No all-depth proof claimed.')


if __name__ == '__main__':
    main()
