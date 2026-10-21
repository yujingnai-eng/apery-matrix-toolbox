#!/usr/bin/env python3
"""Independent verification: no producer imports and no shared matrix helpers.

Requires sympy. Terms are recomputed from defining binomial sums. Complete
blocks are reconstructed from two scalar fundamental solutions, not products
of producer step matrices. Exact symbolic checks supplement finite records.
"""
from functools import lru_cache
import hashlib
import json
from math import comb
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent
DATA = HERE / 'mirror_experiments_20261019.json'
CHECKS = {}


def require(condition, group):
    if not condition:
        raise AssertionError(group)
    CHECKS[group] = CHECKS.get(group, 0)+1


def mat(a):
    return sp.Matrix([[sp.Rational(x) for x in row] for row in a])


def vec(a):
    return sp.Matrix([sp.Rational(x) for x in a])


def order(x,p):
    x=sp.Rational(x)
    if x == 0:
        return None
    num,den=x.as_numer_denom()
    z=0
    while num % p == 0:
        num//=p; z+=1
    while den % p == 0:
        den//=p; z-=1
    return z


def unit2(x,p):
    x=sp.Rational(x)/sp.Rational(p)**order(x,p)
    num,den=x.as_numer_denom()
    return int(num)*pow(int(den),-1,p*p)%(p*p)


def source_term(f,n):
    if f == 'Apery':
        return sum((comb(n,k)*comb(n+k,k))**2 for k in range(n+1))
    power={'Franel3':3,'Franel4':4}[f]
    return sum(comb(n,k)**power for k in range(n+1))


def recurrence(f,n):
    if f == 'Apery':
        return (n+1)**3, (2*n+1)*(17*n*n+17*n+5), -n**3
    if f == 'Franel3':
        return (n+1)**2, 7*n*(n+1)+2, 8*n*n
    return (n+1)**3, 2*(2*n+1)*(3*n*n+3*n+1), 4*n*(4*n-1)*(4*n+1)


@lru_cache(None)
def other_solution(f):
    b=[sp.Rational(0),sp.Rational(1)]
    for n in range(1,96):
        denominator,middle,back=recurrence(f,n)
        b.append((middle*b[-1]+back*b[-2])/denominator)
    return b


def prefix(a,n,t):
    return sum((a[j]*t**j for j in range(n)),sp.Rational(0))


def source_state(a,n,t):
    return sp.Matrix([a[n]*t**n,a[n-1]*t**(n-1),prefix(a,n,t),prefix(a,n,-t)])


@lru_cache(None)
def fundamental_block(f,a,b,t):
    g=VALUES[f]; h=other_solution(f)
    xa=sp.Matrix([[t**a*g[a],t**a*h[a]],
                  [t**(a-1)*g[a-1],t**(a-1)*h[a-1]]])
    xb=sp.Matrix([[t**b*g[b],t**b*h[b]],
                  [t**(b-1)*g[b-1],t**(b-1)*h[b-1]]])
    sums=sp.Matrix([[sum(z[j]*tt**j for j in range(a,b)) for z in (g,h)]
                    for tt in (t,-t)])
    return (xb*xa.inv()).row_join(sp.zeros(2)).col_join(
           (sums*xa.inv()).row_join(sp.eye(2)))


def step(f,n,t,parity=None):
    d,p,q=recurrence(f,n)
    e=(-1)**int(n) if parity is None else parity
    return sp.Matrix([[t*p/d,t*t*q/d,0,0],[1,0,0,0],
                      [1,0,1,0],[e,0,0,1]])


def mirror(n):
    sign=1 if n % 2 == 0 else -1
    return sp.diag(sign,-sign,sp.Matrix([[0,1],[1,0]]))


def verify_symbolic():
    n,t=sp.symbols('n t', nonzero=True)
    for f in VALUES:
        for parity in (1,-1):
            sn=sp.diag(parity,-parity,sp.Matrix([[0,1],[1,0]]))
            sn1=sp.diag(-parity,parity,sp.Matrix([[0,1],[1,0]]))
            residual=step(f,n,-t,parity)*sn-sn1*step(f,n,t,parity)
            require(residual.applyfunc(sp.cancel)==sp.zeros(4),'symbolic_source_mirror')
            j=step(f,n,t,parity).diff(t).row_join(step(f,n,t,parity))
            jet=step(f,n,t,parity).row_join(sp.zeros(4)).col_join(j)
            jminus=jet.subs(t,-t)
            jm=sp.diag(sn,-sn); jm1=sp.diag(sn1,-sn1)
            require((jminus*jm-jm1*jet).applyfunc(sp.cancel)==sp.zeros(8),
                    'symbolic_differentiated_mirror')
    e=sp.Matrix([[0,1],[1,0]])
    for f,tt in [('Apery',1/t),('Franel3',-1/(8*t))]:
        a=step(f,n,t,1)[:2,:2]
        b=step(f,-n-1,tt,1)[:2,:2]
        require((b-e*a.inv()*e).applyfunc(sp.cancel)==sp.zeros(2),'symbolic_recurrence_reflection')
    x=(t+1)/(t-1); y=t/(1-t)**2
    require(sp.cancel(y-(x*x-1)/4)==0,'symbolic_parameter_cover')
    require(sp.cancel(y-y.subs(t,1/t))==0,'symbolic_parameter_cover')


def verify01(e):
    for row in e['cases']:
        f,a,b,t=row['family'],row['a'],row['b'],sp.Rational(row['t'])
        k=fundamental_block(f,a,b,t)
        require(k==mat(row['K']),'E01_independent_fundamental_blocks')
        require(fundamental_block(f,a,b,-t)*mirror(a)==mirror(b)*k,'E01_moving_mirror')
        bad=fundamental_block(f,a,b,-t)*mirror(0)-mirror(0)*k
        require(bool(bad!=sp.zeros(4))==row['constant_mirror_fails'],'E01_negative_control')


def verify02(e):
    j,s,h=mat(e['J']),mat(e['S']),mat(e['H'])
    require(j*j==-sp.eye(2) and s*s==sp.eye(2) and s*j*s==-j,'E02_dihedral_algebra')
    for row in e['cases']:
        f,n=row['family'],row['N']; a=VALUES[f]
        x=sp.Matrix([sum(a[:n]),sum((-1)**i*a[i] for i in range(n))])
        require(x==vec(row['paired']),'E02_actual_pairs')
        require(h*x/2==sp.Matrix([sum(a[:n:2]),sum(a[1:n:2])]),'E02_actual_integral_split')
    require(h.det()==-2 and h*sp.Matrix([1,0])/2==vec(e['negative_control']['half_H_pair']),
            'E02_ambient_lattice_counterexample')


def verify03(e):
    moments={}
    for f,rows in e['families'].items():
        cc=[]
        for row in rows:
            n=row['n']
            if f=='Apery':
                kernel=[comb(n,k)**2*comb(n+k,k)**2 for k in range(n+1)]
            else:
                kernel=[comb(n,k)**({'Franel3':3,'Franel4':4}[f]) for k in range(n+1)]
            val=2*sum(k*v for k,v in enumerate(kernel))-n*sum(kernel)
            require(val==row['mirror_moment'],'E03_direct_marked_moment')
            require(sum(kernel)==row['g'],'E03_source_sum')
            if f!='Apery':
                require(kernel==kernel[::-1] and val==0,'E03_symmetric_control')
            for p in (2,3,5,7,11,13):
                require(order(val,p)==row['valuations'][str(p)],'E03_valuation')
            cc.append(val)
        moments[f]=cc
    for row in e['tower']:
        n=row['N']; total=sum((-1)**j*moments['Apery'][j] for j in range(n))
        require(total==row['signed_moment_sum'] and order(total,2)==row['v2'],'E03_tower_counterexample')
        require(order(prefix(VALUES['Apery'],n,sp.Rational(-1)),2)==row['original_v2'],'E03_original_readout')
    for row in e['marked_state_obstruction']:
        f,n,z=row['family'],row['n'],sp.Rational(row['z'])
        def h(m):
            return sum(((comb(m,k)*comb(m+k,k))**2 if f=='Apery'
                       else comb(m,k)**({'Franel3':3,'Franel4':4}[f]))*z**k for k in range(m+1))
        den,mid,back=recurrence(f,n)
        require(den*h(n+1)-mid*h(n)-back*h(n-1)==sp.Rational(row['recurrence_residual']),
                'E03_marked_closure_obstruction')


def verify04(e):
    for row in e['cases']:
        a,n,p=VALUES[row['family']],row['N'],row['p']
        vals=[order(a[i]-a[n-1-i],p) for i in range(n) if a[i]!=a[n-1-i]]
        require((min(vals) if vals else None)==row['coefficient_difference_min_vp'],'E04_coefficient_reversal')
        v=sum((-1)**j*a[j] for j in range(n))
        require(v==sp.Rational(row['original_at_minus1']),'E04_direct_prefix')
        require(sp.Rational(row['mirror_at_minus1'])==(-1)**(n-1)*v,'E04_universal_sign')
    for row in e['controls']:
        require(row['mirror']==-row['F'],'E04_nonspecial_controls')


def independent_deck_jet(a,N,d):
    # A common denominator (1-x)^N, not the producer's termwise denominator.
    numerator=[]
    for j in range(d+1):
        numerator.append(2*sum((-1)**n*a[n]*sum(
            comb(n,i)*(-1)**(j-i)*comb(N-1-n,j-i)
            for i in range(j+1) if i<=n and j-i<=N-1-n) for n in range(N)))
    return [sum(numerator[k]*comb(N+j-k-1,j-k) for k in range(j+1)) for j in range(d+1)]


def verify05(e):
    a=VALUES['Apery']; b=[comb(2*k,k)*VALUES['Franel3'][k] for k in range(33)]
    for n in range(33):
        require(a[n]==sum(b[k]*comb(n+k,2*k) for k in range(n+1)),'E05_formal_bridge')
    for row in e['finite_boundary_counterexamples']:
        n,t=row['N'],sp.Rational(row['t']); y=t/(1-t)**2
        left=(1-t)*prefix(a,n,t); right=prefix(b,n,y)
        require(left==sp.Rational(row['finite_A']) and right==sp.Rational(row['finite_B']),
                'E05_finite_boundary_values')
        require(left-right==sp.Rational(row['boundary']),'E05_finite_boundary_obstruction')
    for row in e['deck_jets']:
        cc=independent_deck_jet(a,row['N'],5)
        require(cc==row['deck_jet_coefficients'],'E05_independent_common_denominator_jets')
        require([order(z,2) for z in cc]==row['valuations_v2'],'E05_jet_valuations')


def verify06(e):
    t=sp.Symbol('t')
    for row in e['cases']:
        a,n=VALUES[row['family']],row['N']
        # K_(0,N) from the defining scalar solution and cumulative polynomials.
        k=sp.Matrix([[a[n]*t**n,0,0,0],[a[n-1]*t**(n-1),0,0,0],
                     [prefix(a,n,t),0,1,0],[prefix(a,n,-t),0,0,1]])
        jet=k.row_join(sp.zeros(4)).col_join(k.diff(t).row_join(k))
        require(jet.subs(t,sp.Rational(row['t']))==mat(row['jet_matrix']),
                'E06_independent_polynomial_derivative')
    for row in e['wrong_binary_freeze']:
        a,n=VALUES[row['family']],row['N']
        good=prefix(a,2*n,sp.Rational(-1))
        bad=sum((-1)**j*(a[2*j]-a[2*j+1]) for j in range(n))
        require(good==row['correct_t_squared_equals1'] and bad==row['incorrectly_frozen_at_minus1'],
                'E06_binary_parameter_counterexample')


def verify07(e):
    for row in e['protocols']:
        f,a,t=row['family'],row['a'],sp.Rational(row['t'])
        k1=fundamental_block(f,a,2*a,t); k2=fundamental_block(f,2*a,3*a,t)
        require(k2*k1==mat(row['forward']) and k1*k2==mat(row['reverse']),
                'E07_independent_ordered_blocks')
        initial=source_state(VALUES[f],a,t)
        difference=(k2*k1-k1*k2)*initial
        require(difference[2:]==list(vec(row['source_readout_difference'])),
                'E07_ordered_readout')
    require(mat(e['Franel4_wrong_reciprocal_residual'])!=sp.zeros(2),'E07_wrong_family_reflection')


def verify08(e):
    for row in e['cases']:
        p,d=row['p'],row['d']; c=mat(row['C_CM'])
        leg=1 if pow(-d%p,(p-1)//2,p)==1 else -1
        require(leg==row['chi_CM'],'E08_character')
        if leg==-1:
            ii=mat(row['engineered_iota']); g=mat(row['lattice_G'])
            require(ii*ii==-d*sp.eye(2) and c*ii+ii*c==sp.zeros(2),'E08_engineered_algebra_only')
            require(min(order(x,p) for x in ii if x)==-1,'E08_nonintegral_generator')
            require(g.inv()*ii*g==mat(row['engineered_iota_on_G_lattice']),'E08_adapted_lattice')
            require(all((a*a+d)%p for a in range(p)),'E08_standard_lattice_obstruction')
        else:
            require(c.trace()!=0,'E08_trace_obstruction')
    require(2*mat(e['fixed_mod128_obstruction']['C']).trace()%128==12,'E08_fixed128_obstruction')


def verify09(e):
    c=sp.diag(1,-1,-1,-1)
    li=sp.Matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,-1],[0,0,1,0]])
    ri=sp.Matrix([[0,-1,0,0],[1,0,0,0],[0,0,0,1],[0,0,-1,0]])
    require(c*li*c==-ri and c*li*c!=-li,'E09_left_right_switch')
    x=sp.Symbol('x')
    a,b,cc,d,alpha,beta=sp.symbols('a b c d alpha beta')
    generic=sp.Matrix([[a,alpha*b,beta*cc,-alpha*beta*d],
                       [b,a,beta*d,-beta*cc],[cc,-alpha*d,a,alpha*b],[d,-cc,b,a]])
    norm=a*a-alpha*b*b-beta*cc*cc+alpha*beta*d*d
    require(sp.expand(generic.charpoly(x).as_expr()-(x*x-2*a*x+norm)**2)==0,
            'symbolic_general_quaternion_left_characteristic_polynomial')
    for row in e['actual_block_obstructions']:
        f,a,b=row['family'],row['a'],row['b']; k=fundamental_block(f,a,b,sp.Rational(row['t']))
        bb=k[:2,:2]
        disc=bb.trace()**2-4*bb.det()
        require(disc==sp.Rational(row['discriminant_B']) and disc!=0,'E09_quaternion_charpoly_obstruction')
        require(sp.expand(k.charpoly(x).as_expr()-(x-1)**2*bb.charpoly(x).as_expr())==0,
                'E09_complete_characteristic_polynomial')
    for row in e['actual_common_commutants']:
        matrices=[step(row['family'],n,sp.Rational(row['t'])) for n in range(1,5)]
        equations=sp.Matrix.vstack(*[sp.kronecker_product(sp.eye(4),a)-
                                    sp.kronecker_product(a.T,sp.eye(4)) for a in matrices])
        dimension=16-equations.rank()
        require(dimension==row['common_commutant_dimension']==1,'E09_independent_kronecker_commutant')


def verify10(e):
    s=sp.diag(1,-1)
    for row in e['cases']:
        f,p,a,b=row['family'],row['p'],row['a'],row['b']
        bb=fundamental_block(f,a,b,sp.Rational(1))[:2,:2]
        v=sp.Matrix([VALUES[f][a],VALUES[f][a-1]]); w=sp.Matrix([1,0])
        val=sp.Matrix.hstack(bb*v,bb*w).det()
        require(val==sp.Rational(row['oriented_readout']),'E10_independent_wedge')
        require(val==-sp.Rational(row['mirror_oriented_readout']),'E10_orientation')
        require(order(val,p)==row['v'] and unit2(val,p)==row['unit_mod_p2'], 'E10_unit_sensitive_measurement')
        aa=min(order(z,p) for z in bb if z)
        require([aa,order(bb.det(),p)-aa]==row['smith'],'E10_smith')
        mb=s*bb*s
        require(min(order(z,p) for z in mb if z)==aa,'E10_smith_mirror_blindness')


def main():
    global VALUES
    data=json.loads(DATA.read_text())
    source=HERE/'mirror_experiments_20261019.py'
    require(hashlib.sha256(source.read_bytes()).hexdigest()==data['producer_sha256'],'producer_hash')
    VALUES={}
    for f,expected in data['source_terms'].items():
        actual=[source_term(f,n) for n in range(len(expected))]
        require(actual==expected,'independent_binomial_terms_513_per_family')
        VALUES[f]=actual
    verify_symbolic()
    for number,fn in enumerate((verify01,verify02,verify03,verify04,verify05,
                               verify06,verify07,verify08,verify09,verify10),1):
        fn(data['experiments'][f'E{number:02d}'])
        print(f'E{number:02d}: independent checks passed',flush=True)
    out=dict(status='pass within declared scope',source_method='defining binomial sums, n=0..512 for three families',
             block_method='two scalar fundamental solutions and direct cumulative sums',
             symbolic_checks=['source mirror, both parities','complete derivative mirror',
                              'Apéry and Franel3 recurrence reflections','rational parameter cover',
                              'general quaternion left-regular characteristic polynomial is a square'],
             check_counts=CHECKS,all_ten_experiments_covered=True,
             limits=['No all-r proof of G1/G2','No proof of intrinsic chirality or CM identification',
                     'Finite observed valuation patterns remain finite observations'],
             hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (source,DATA,Path(__file__))})
    path=HERE/'verify_mirror_experiments_20261019.json'
    path.write_text(json.dumps(out,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(out,ensure_ascii=False,indent=2))


if __name__=='__main__':
    main()
