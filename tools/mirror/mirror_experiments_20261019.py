#!/usr/bin/env python3
"""Ten source-backed mirror probes. Exact arithmetic; no fitted matrices.

Run: python3 mirror_experiments_20261019.py --save
Source recurrences generate terms; the separate verifier uses binomial sums.
An expected negative control is a successful falsification, not a theorem.
No experiment proves G1/G2, geometric chirality, or a CM identification.
"""
from __future__ import annotations

import argparse
from fractions import Fraction as Q
from functools import lru_cache
import hashlib
import json
from math import comb, isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
FAMILIES = ('Apery', 'Franel3', 'Franel4')
PRIMES = (2, 3, 5, 7, 11, 13)
MAX_N = 512


def eye(n):
    return [[Q(i == j) for j in range(n)] for i in range(n)]


def mm(a, b):
    return [[sum((x*y for x, y in zip(row, col)), Q(0))
             for col in zip(*b)] for row in a]


def mv(a, x):
    return [sum((v*w for v, w in zip(row, x)), Q(0)) for row in a]


def add(a, b, sign=1):
    return [[x+sign*y for x, y in zip(r, s)] for r, s in zip(a, b)]


def scale(a, c):
    return [[c*x for x in row] for row in a]


def zero(a):
    return all(not x for row in a for x in row)


def inv(a):
    n = len(a)
    z = [list(map(Q, row)) + e for row, e in zip(a, eye(n))]
    for j in range(n):
        k = next((k for k in range(j, n) if z[k][j]), None)
        if k is None:
            raise ValueError('singular; no modular inverse has been assumed')
        z[j], z[k] = z[k], z[j]
        d = z[j][j]
        z[j] = [x/d for x in z[j]]
        for k in range(n):
            if k != j:
                d = z[k][j]
                z[k] = [x-d*y for x, y in zip(z[k], z[j])]
    return [r[n:] for r in z]


def rank(a):
    a = [list(map(Q, r)) for r in a]
    at = 0
    for j in range(len(a[0])):
        k = next((k for k in range(at, len(a)) if a[k][j]), None)
        if k is None:
            continue
        a[at], a[k] = a[k], a[at]
        d = a[at][j]
        a[at] = [x/d for x in a[at]]
        for k in range(at+1, len(a)):
            d = a[k][j]
            if d:
                a[k] = [x-d*y for x, y in zip(a[k], a[at])]
        at += 1
        if at == len(a):
            break
    return at


def det2(a):
    return a[0][0]*a[1][1]-a[0][1]*a[1][0]


def vp(x, p):
    x = Q(x)
    if not x:
        return None  # Exact zero, NOT a censored lower bound.
    def v(n):
        n, k = abs(n), 0
        while n % p == 0:
            n //= p
            k += 1
        return k
    return v(x.numerator)-v(x.denominator)


def matrix_vp(a, p):
    z = [vp(x, p) for row in a for x in row if x]
    return min(z) if z else None


def unit(x, p, digits):
    x = Q(x)
    if not x:
        return None
    v = vp(x, p)
    u = x / Q(p)**v
    mod = p**digits
    assert u.denominator % p
    return u.numerator*pow(u.denominator, -1, mod) % mod


def smith2(a, p):
    if not det2(a):
        raise ValueError('Smith exponents here require an invertible matrix')
    lo = matrix_vp(a, p)
    return [lo, vp(det2(a), p)-lo]


def coefficients(family, n):
    if family == 'Apery':
        return 3, 34*n**3+51*n*n+27*n+5, -n**3
    if family == 'Franel3':
        return 2, 7*n*n+7*n+2, 8*n*n
    return 3, 12*n**3+18*n*n+10*n+2, 4*n*(16*n*n-1)


@lru_cache(None)
def terms(family):
    a = [1]
    for n in range(MAX_N):
        d, p, q = coefficients(family, n)
        top = p*a[n]+(q*a[n-1] if n else 0)
        value, rem = divmod(top, (n+1)**d)
        assert rem == 0
        a.append(value)
    return a


def kernel(family, n, k):
    if family == 'Apery':
        return comb(n, k)**2*comb(n+k, k)**2
    return comb(n, k)**(3 if family == 'Franel3' else 4)


def endpoint(family, n, t):
    d, p, q = coefficients(family, n)
    return [[Q(t)*Q(p, (n+1)**d), Q(t)**2*Q(q, (n+1)**d)],
            [Q(1), Q(0)]]


def step(family, n, t):
    b = endpoint(family, n, t)
    return [b[0]+[Q(0), Q(0)], b[1]+[Q(0), Q(0)],
            [Q(1), Q(0), Q(1), Q(0)],
            [Q((-1)**n), Q(0), Q(0), Q(1)]]


def step_prime(family, n, t):
    d, p, q = coefficients(family, n)
    z = scale(eye(4), 0)
    z[0][0] = Q(p, (n+1)**d)
    z[0][1] = 2*Q(t)*Q(q, (n+1)**d)
    return z


@lru_cache(None)
def block(family, a, b, t):
    out = eye(4)
    for n in range(a, b):
        out = mm(step(family, n, t), out)
    return out


def state(family, n, t, derivative=False):
    a = terms(family)
    t = Q(t)
    if not derivative:
        return [t**n*a[n], t**(n-1)*a[n-1] if n else Q(0),
                sum((t**j*a[j] for j in range(n)), Q(0)),
                sum(((-t)**j*a[j] for j in range(n)), Q(0))]
    return [n*t**(n-1)*a[n] if n else Q(0),
            (n-1)*t**(n-2)*a[n-1] if n > 1 else Q(0),
            sum((j*t**(j-1)*a[j] for j in range(1, n)), Q(0)),
            sum(((-1)**j*j*t**(j-1)*a[j] for j in range(1, n)), Q(0))]


def mirror(n):
    sign=1 if n % 2 == 0 else -1
    return [[Q(sign), 0, 0, 0], [0, Q(-sign), 0, 0],
            [0, 0, 0, 1], [0, 0, 1, 0]]


def diagonal_blocks(a, b):
    n = len(a)
    z = [[Q(0)]*n for _ in range(n)]
    return [r+s for r, s in zip(a, z)] + [r+s for r, s in zip(z, b)]


def lengths():
    out = []
    for p in PRIMES:
        r, n = 0, 1
        while n <= MAX_N:
            out.append((p, r, n))
            n *= p
            r += 1
    return out


def probe01():
    rows = []
    fixed = mirror(0)
    for f in FAMILIES:
        for t in (Q(1), Q(-1)):
            for a in (1, 2, 3, 4, 5, 7, 8, 11, 13, 16, 32):
                b = 2*a
                k, km = block(f, a, b, t), block(f, a, b, -t)
                good = add(mm(km, mirror(a)), mm(mirror(b), k), -1)
                assert zero(good)
                assert mv(k, state(f, a, t)) == state(f, b, t)
                assert mv(mirror(a), state(f, a, t)) == state(f, a, -t)
                bad = add(mm(km, fixed), mm(fixed, k), -1)
                rows.append(dict(family=f, t=t, a=a, b=b, K=k,
                                 moving_mirror_residual=good,
                                 constant_mirror_min_v2=matrix_vp(bad, 2),
                                 constant_mirror_fails=not zero(bad)))
    assert any(x['constant_mirror_fails'] for x in rows)
    return dict(kind='source-backed exact transport',
                statement='K(-t) S_a = S_b K(t); S_n has endpoint parity phases',
                state='(t^n g_n,t^(n-1)g_(n-1),F_n(t),F_n(-t))',
                ring='Q; no reduction of nonunit denominators', cases=rows)


def probe02():
    j, s, h = [[0,-1],[1,0]], [[1,0],[0,-1]], [[1,1],[1,-1]]
    assert mm(j,j) == scale(eye(2),-1)
    assert mm(s,s) == eye(2)
    assert mm(mm(s,j),s) == scale(j,-1)
    rows=[]
    for f in FAMILIES:
        for n in (1,2,3,4,8,16,32,64,128,256,512):
            x=state(f,n,1)[2:]
            half=mv(h,x)
            assert all(z.denominator == 1 and z % 2 == 0 for z in half)
            even=sum(terms(f)[0:n:2]); odd=sum(terms(f)[1:n:2])
            assert [z/2 for z in half] == [even,odd]
            rows.append(dict(family=f,N=n,paired=x,even=even,odd=odd))
    return dict(kind='representation control plus actual prefix lattice',
                J=j,S=s,H=h,det_H=-2,
                smith_H={str(p):smith2(h,p) for p in PRIMES},
                cases=rows,negative_control=dict(ambient_pair=[1,0],
                half_H_pair=[Q(1,2),Q(1,2)]),
                conclusion='No unimodular split on Z_2^2; actual prefixes lie in x=y mod2 and DO split integrally there.')


def probe03():
    data={}; moment_arrays={}
    for f in FAMILIES:
        rows=[]; cc=[]
        for n in range(129):
            h=[kernel(f,n,k) for k in range(n+1)]
            assert sum(h) == terms(f)[n]
            c=sum((2*k-n)*z for k,z in enumerate(h))
            assert sum((2*k-n)*z for k,z in enumerate(reversed(h))) == -c
            if f != 'Apery':
                assert c == 0
            cc.append(c)
            rows.append(dict(n=n,g=terms(f)[n],mirror_moment=c,
                             valuations={str(p):vp(c,p) for p in PRIMES},
                             marked_coefficients=h if n<=8 else None))
        data[f]=rows; moment_arrays[f]=cc
    tower=[]
    for r in range(1,8):
        n=2**r
        total=sum((-1)**j*moment_arrays['Apery'][j] for j in range(n))
        tower.append(dict(r=r,N=n,moment=moment_arrays['Apery'][n],
                          signed_moment_sum=total,v2=vp(total,2),
                          original_F=state('Apery',n,-1)[2],
                          original_v2=vp(state('Apery',n,-1)[2],2)))
    # The marked polynomial cannot simply reuse the unmarked four-state transport.
    defects=[]
    for f in FAMILIES:
        for z in (Q(2),Q(1,2)):
            marked=[sum((kernel(f,n,k)*z**k for k in range(n+1)),Q(0))
                    for n in range(10)]
            for n in range(1,8):
                d,p,q=coefficients(f,n)
                rem=(n+1)**d*marked[n+1]-p*marked[n]-q*marked[n-1]
                defects.append(dict(family=f,z=z,n=n,recurrence_residual=rem))
    assert any(x['recurrence_residual'] for x in defects)
    return dict(kind='new marked observable, not a geometric chirality claim',
                definition='c_n=sum_k(2k-n)h(n,k)',families=data,tower=tower,
                marked_state_obstruction=defects,
                conclusion='Franel odd moments vanish by kernel symmetry; marked kernels require additional source channels.')


def probe04():
    rows=[]
    for f in FAMILIES:
        a=terms(f)
        for p,r,n in lengths():
            differences=[a[j]-a[n-1-j] for j in range(n)]
            z=[vp(x,p) for x in differences if x]
            value=state(f,n,-1)[2]
            reflected=sum(((-1)**j*a[n-1-j] for j in range(n)),Q(0))
            assert reflected == (-1)**(n-1)*value
            assert list(reversed(list(reversed(a[:n])))) == a[:n]
            rows.append(dict(family=f,p=p,r=r,N=n,
                             coefficient_difference_min_vp=min(z) if z else None,
                             original_at_minus1=value,mirror_at_minus1=reflected))
    controls=[]
    for name,a in [('constant',[1]*16),('linear',list(range(16))),
                   ('irregular',[j*j+3*j+1 for j in range(16)])]:
        for n in (2,4,8,16):
            x=sum((-1)**j*a[j] for j in range(n))
            y=sum((-1)**j*a[n-1-j] for j in range(n))
            assert y == -x
            controls.append(dict(sequence=name,N=n,F=x,mirror=y))
    return dict(kind='exact finite-polynomial reversal',cases=rows,
                controls=controls,conclusion='Mirror sign at -1 for even N is universal, not evidence for 3r-1.')


def deck_jet(a,n,j):
    # [x^j] (1-t(x))*F_N(t(x)), t(x)=-(1+x)/(1-x).
    return 2*sum((-1)**k*a[k]*sum(comb(k,i)*comb(k+j-i,j-i)
                 for i in range(min(k,j)+1)) for k in range(n))


def probe05():
    a,b=terms('Apery'),[comb(2*k,k)*terms('Franel3')[k] for k in range(33)]
    for n in range(33):
        assert a[n] == sum(b[k]*comb(n+k,2*k) for k in range(n+1))
    coordinates=[]
    for t in (Q(-1),Q(2),Q(3),Q(-2),Q(1,2)):
        y=t/(1-t)**2; x=(t+1)/(t-1)
        assert (1/t)/(1-1/t)**2 == y
        assert ((1/t+1)/(1/t-1)) == -x
        assert y == (x*x-1)/4
        coordinates.append(dict(t=t,x=x,y=y))
    finite=[]
    for n in (2,4,8):
        for t in (Q(2),Q(1,2),Q(-1)):
            y=t/(1-t)**2
            left=(1-t)*sum((a[k]*t**k for k in range(n)),Q(0))
            right=sum((b[k]*y**k for k in range(n)),Q(0))
            finite.append(dict(N=n,t=t,finite_A=left,finite_B=right,boundary=left-right))
    assert any(x['boundary'] for x in finite)
    rows=[]
    for r in range(1,9):
        n=2**r
        co=[deck_jet(a,n,j) for j in range(6)]
        assert co[0] == 2*state('Apery',n,-1)[2]
        rows.append(dict(r=r,N=n,deck_jet_coefficients=co,
                         valuations_v2=[vp(x,2) for x in co]))
    return dict(kind='exact rational parameter cover plus finite-boundary controls',
                coordinate_identity='y=t/(1-t)^2=(x^2-1)/4; t -> 1/t sends x -> -x',
                formal_Strehl_check_n=[0,32],coordinates=coordinates,
                finite_boundary_counterexamples=finite,deck_jets=rows,
                conclusion='Finite B_N(y(x)) is even in x, but the finite A-side has nonzero odd boundary jets.')


def probe06():
    rows=[]
    for f in FAMILIES:
        for t in (Q(1),Q(-1),Q(2,3)):
            jet=eye(8)
            for n in range(16):
                z=[[Q(0)]*4 for _ in range(4)]
                m,mp=step(f,n,t),step_prime(f,n,t)
                jj=[r+s for r,s in zip(m,z)]+[r+s for r,s in zip(mp,m)]
                jet=mm(jj,jet)
                source=state(f,0,t)+state(f,0,t,True)
                assert mv(jet,source) == state(f,n+1,t)+state(f,n+1,t,True)
                mir=diagonal_blocks(mirror(n),scale(mirror(n),-1))
                mir_next=diagonal_blocks(mirror(n+1),scale(mirror(n+1),-1))
                m2,mp2=step(f,n,-t),step_prime(f,n,-t)
                other=[r+s for r,s in zip(m2,z)]+[r+s for r,s in zip(mp2,m2)]
                assert mm(other,mir) == mm(mir_next,jj)
            rows.append(dict(family=f,t=t,N=16,jet_matrix=jet,
                             readout=mv(jet,source),all_step_checks=16))
    freeze=[]
    for f in FAMILIES:
        a=terms(f)
        for n in (2,3,4,8):
            good=sum(a[2*j]-a[2*j+1] for j in range(n))
            bad=sum((-1)**j*(a[2*j]-a[2*j+1]) for j in range(n))
            assert good == state(f,2*n,-1)[2]
            freeze.append(dict(family=f,N=n,correct_t_squared_equals1=good,
                               incorrectly_frozen_at_minus1=bad,residual=bad-good))
    assert any(x['residual'] for x in freeze)
    return dict(kind='exact differentiated COMPLETE source transport',
                dimensions=dict(endpoint_jet=4,complete_paired_jet=8),
                formula='[[T,0],[dT/dt,T]]; derivative reflection diag(S_n,-S_n)',
                cases=rows,wrong_binary_freeze=freeze,
                conclusion='Complete differentiated state is eight-dimensional here; no four-dimensional closure is assumed.')


def probe07():
    swap=[[0,1],[1,0]]; identities=[]; protocols=[]
    for f in ('Apery','Franel3'):
        for n in range(1,13):
            for t in (Q(1),Q(-1),Q(2,3)):
                target=1/t if f == 'Apery' else -1/(8*t)
                r=endpoint(f,n,t)
                reflected=endpoint(f,-n-1,target)
                assert reflected == mm(mm(swap,inv(r)),swap)
                identities.append(dict(family=f,n=n,t=t,reflected_t=target))
    for f in FAMILIES:
        for a in (1,2,3,4,5,7,8):
            for t in (Q(1),Q(-1)):
                k1,k2=block(f,a,2*a,t),block(f,2*a,3*a,t)
                forward,reverse=mm(k2,k1),mm(k1,k2)
                assert forward == block(f,a,3*a,t)
                bk=[r[:2] for r in forward[:2]]
                hh=[r[:2] for r in forward[2:]]
                bi=inv(bk)
                formula=[r+[Q(0),Q(0)] for r in bi]
                formula += [r+e for r,e in zip(scale(mm(hh,bi),-1),eye(2))]
                assert formula == inv(forward)
                x=state(f,a,t)
                assert mv(inv(forward),mv(forward,x)) == x
                diff=[x-y for x,y in zip(mv(forward,x),mv(reverse,x))]
                protocols.append(dict(family=f,a=a,t=t,forward=forward,reverse=reverse,
                                      source_readout_difference=diff[2:],
                                      commutator_min_v2=matrix_vp(add(forward,reverse,-1),2)))
    # Franel4 does not inherit Apéry's reciprocal parameter reflection.
    f4=add(endpoint('Franel4',-2,Q(1)),
           mm(mm(swap,inv(endpoint('Franel4',1,Q(1)))),swap),-1)
    assert not zero(f4)
    return dict(kind='source recurrence reversibility and ordered block controls',
                endpoint_reflection=identities,protocols=protocols,
                Franel4_wrong_reciprocal_residual=f4,
                conclusion='Correct reverse transport includes -H B^(-1); reversing a schedule is not the original forward recurrence.')


def character(p,d):
    x=pow((-d)%p,(p-1)//2,p)
    return 0 if not x else (1 if x == 1 else -1)


def probe08():
    rows=[]; i2=eye(2)
    for p in (5,7,11,13,17,19):
        for t,d in ((1,2),(-1,3)):
            ch=character(p,d)
            if ch == 1:
                candidates=[a for a in range(isqrt(p)+1) for b in range(isqrt(p//d)+1)
                            if a*a+d*b*b == p]
                assert candidates
                tr=4*candidates[0]**2-2*p
            else:
                tr=0
            c=[[Q(tr),Q(-ch*p*p)],[Q(1),Q(0)]]
            base=[[0,-d],[1,0]]
            row=dict(p=p,t=t,d=d,chi_CM=ch,T_CM=tr,C_CM=c,
                     naive_iota_anticommutator=add(mm(c,base),mm(base,c)))
            if ch == -1:
                # Engineered partner exists for ANY d, not a CM identification.
                fake=[[Q(1-d,2),Q(p*(d+1),2)],
                      [Q(-(d+1),2*p),Q(d-1,2)]]
                assert mm(fake,fake) == scale(i2,-d)
                assert zero(add(mm(c,fake),mm(fake,c)))
                g=[[Q(p),0],[0,Q(1)]]
                adapted=mm(mm(inv(g),fake),g)
                assert matrix_vp(fake,p) == -1
                assert matrix_vp(adapted,p) >= 0
                assert all((a*a+d)%p for a in range(p))
                row.update(engineered_iota=fake,engineered_iota_on_G_lattice=adapted,
                           lattice_G=g,standard_lattice_obstruction='a^2=-d mod p has no solution',
                           conclusion='Engineered generator works on G lattice; it is noncanonical and not geometric evidence.')
            else:
                assert tr
                row['invertible_anticommutator_obstruction']='trace C != 0 in characteristic zero'
            rows.append(row)
    c=[[6,16],[1,0]]
    assert (2*6)%128 != 0
    return dict(kind='candidate arithmetic matrices; not identified geometric Frobenius',
                cases=rows,fixed_mod128_obstruction=dict(C=c,
                necessary_trace_test_2tr_mod128=12,
                conclusion='No invertible anticommuting partner exists modulo128.'))


def qmul(x,y):
    a,b,c,d=x; e,f,g,h=y
    return [a*e-b*f-c*g-d*h,a*f+b*e+c*h-d*g,
            a*g-b*h+c*e+d*f,a*h+b*g-c*f+d*e]


def qconj(q):
    return [q[0]]+[-x for x in q[1:]]


def qmatrix(q,right=False):
    basis=eye(4)
    cols=[qmul(x,q) if right else qmul(q,x) for x in basis]
    return [list(z) for z in zip(*cols)]


def commutant_dimension(matrices):
    n=len(matrices[0]); equations=[]
    for a in matrices:
        for i in range(n):
            for j in range(n):
                row=[Q(0)]*(n*n)
                for k in range(n):
                    row[k*n+j]+=a[i][k]
                    row[i*n+k]-=a[k][j]
                equations.append(row)
    return n*n-rank(equations)


def probe09():
    conjugation=[[Q(i==j)*(1 if i==0 else -1) for j in range(4)] for i in range(4)]
    qi,qj=[0,1,0,0],[0,0,1,0]
    assert qconj(qmul(qi,qj)) == qmul(qconj(qj),qconj(qi))
    lhs=mm(mm(conjugation,qmatrix(qi)),conjugation)
    assert lhs == qmatrix(qconj(qi),True)
    wrong=add(lhs,qmatrix(qconj(qi)),-1)
    assert not zero(wrong)
    model_dimension=commutant_dimension([qmatrix(qi),qmatrix(qj)])
    assert model_dimension == 4
    rows=[]; common=[]
    for f in FAMILIES:
        for t in (Q(1),Q(-1)):
            for a in (1,2,3,5,7):
                k=block(f,a,2*a,t)
                b=[r[:2] for r in k[:2]]
                tr=b[0][0]+b[1][1]; dt=det2(b)
                discriminant=tr*tr-4*dt
                assert discriminant
                rows.append(dict(family=f,t=t,a=a,b=2*a,tr_B=tr,det_B=dt,
                                 discriminant_B=discriminant,
                                 obstruction='charpoly(K)=(X-1)^2 charpoly(B) is not a polynomial square'))
            dim=commutant_dimension([step(f,n,t) for n in range(1,5)])
            common.append(dict(family=f,t=t,common_commutant_dimension=dim,
                               cannot_be_one_common_quaternion_left_algebra=dim<4))
    return dict(kind='quaternion algebra control and actual complete-block obstruction',
                quaternion_mirror='C L(q) C = R(conjugate(q)), not L(conjugate(q))',
                wrong_left_only_residual=wrong,model_commutant_dimension=model_dimension,
                actual_block_obstructions=rows,actual_common_commutants=common,
                conclusion='Mirror switches left/right actions; actual complete blocks are not single quaternion left multiplications.')


def wedge(x,y):
    return x[0]*y[1]-x[1]*y[0]


def probe10():
    s=[[1,0],[0,-1]]; rows=[]
    for f in FAMILIES:
        for p in PRIMES:
            a,b=p,2*p
            k=block(f,a,b,Q(1)); bb=[r[:2] for r in k[:2]]
            x=state(f,a,Q(1))[:2]; y=[Q(1),Q(0)]
            initial=wedge(x,y); out=wedge(mv(bb,x),mv(bb,y))
            assert out == det2(bb)*initial
            mirrored=wedge(mv(s,mv(bb,x)),mv(s,mv(bb,y)))
            assert mirrored == -out
            mb=mm(mm(s,bb),s)
            assert smith2(bb,p) == smith2(mb,p)
            assert vp(out,p) == vp(mirrored,p)
            assert unit(out,p,2) != unit(mirrored,p,2)
            rows.append(dict(family=f,p=p,a=a,b=b,endpoint_B=bb,
                             oriented_readout=out,mirror_oriented_readout=mirrored,
                             determinant=det2(bb),smith=smith2(bb,p),
                             v=vp(out,p),unit_mod_p2=unit(out,p,2),
                             mirror_unit_mod_p2=unit(mirrored,p,2)))
    return dict(kind='exterior-square readout with unit-sensitive valuations',cases=rows,
                conclusion='Valuation/Smith exponents miss orientation; wedge transport itself is already determined by det(B).')


def pack(x):
    if isinstance(x,Q):
        return str(x)
    if isinstance(x,dict):
        return {str(k):pack(v) for k,v in x.items()}
    if isinstance(x,(list,tuple)):
        return [pack(v) for v in x]
    return x


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--save',action='store_true')
    args=parser.parse_args()
    results={}
    for k,fn in enumerate((probe01,probe02,probe03,probe04,probe05,
                          probe06,probe07,probe08,probe09,probe10),1):
        results[f'E{k:02d}']=fn()
        print(f'E{k:02d}: completed, including negative controls',flush=True)
    output=dict(schema_version=1,date='2026-10-19',status='pass within bounded scope',
                configuration=dict(families=FAMILIES,primes=PRIMES,max_term_index=MAX_N,
                                   moment_n_range=[0,128],zero_valuation='null means exact zero'),
                source_terms={f:terms(f) for f in FAMILIES},experiments=results,
                limits=['No all-r proof of G1/G2','No intrinsic chirality or geometric CM identification',
                        'Q-valued complete blocks are not assumed Z_p-integral',
                        'Expected failed hypotheses are preserved as negative controls',
                        'No regression or reconstruction of allowed lattices from desired output'],
                producer_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest())
    if args.save:
        path=Path(__file__).with_suffix('.json')
        path.write_text(json.dumps(pack(output),ensure_ascii=False,indent=2)+'\n')
        print(f'Saved {path.name}',flush=True)
    print('No all-depth proof claimed.',flush=True)


if __name__ == '__main__':
    main()
