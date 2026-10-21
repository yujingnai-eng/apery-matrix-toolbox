"""Wilson--Gessel mod-p^2 four-state scheme, no fitting or external packages.

Input theorem: Gessel's Apéry congruence, using the derivative defined in
Rowland--Yassawi--Krattenthaler (2021), equations (3)--(4).
Scope is the unit-factorial D tower, not the unnormalized E or original F tower.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent
PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


def apery(limit):
    a = [1, 5]
    for n in range(1, limit):
        q, r = divmod((34*n**3+51*n*n+27*n+5)*a[-1]-n**3*a[-2], (n+1)**3)
        assert r == 0
        a.append(q)
    return a


def residue(q, modulus):
    q = Q(q)
    return q.numerator*pow(q.denominator, -1, modulus) % modulus


def coefficients(p, a):
    h = [Q(0)]
    for n in range(1, 2*p):
        h.append(h[-1]+Q(1, n))
    aprime = [sum((2*comb(ell, j)**2*comb(ell+j, j)**2*(h[ell+j]-h[ell-j])
                   for j in range(ell+1)), Q(0)) for ell in range(p)]
    w = (factorial(p-1)+1)//p
    e = [Q(a[ell], factorial(ell)) for ell in range(p)]
    b0 = [(aprime[ell]-h[ell]*a[ell])/factorial(ell) for ell in range(p)]
    beta = [b0[ell]+w*e[ell] for ell in range(p)]
    return aprime, e, beta, b0, w


def unit_coefficients(a, p):
    modulus, unit, out = p*p, 1, []
    for n, value in enumerate(a):
        if n:
            x = n
            while x % p == 0:
                x //= p
            unit = unit*x % modulus
        out.append(value*pow(unit, -1, modulus) % modulus)
    return out


def eval_poly(poly, t, modulus):
    value = 0
    for c in reversed(poly):
        value = (value*t+c) % modulus
    return value


def inner_matrix(e, beta, t, p):
    m = p*p
    e, b = [residue(x, m) for x in e], [residue(x, m) for x in beta]
    E, B = eval_poly(e, t, m), eval_poly(b, t, m)
    thetaE = eval_poly([j*x for j, x in enumerate(e)], t, m)
    thetaB = eval_poly([j*x for j, x in enumerate(b)], t, m)
    return [[E, p*B % m], [thetaE, p*(E+thetaB) % m]]


def outer_matrix(plus, minus):
    return [[0, 0]+plus[0], [0, 0]+plus[1],
            minus[0]+[0, 0], minus[1]+[0, 0]]


def mv(matrix, state, modulus):
    return [sum(x*y for x, y in zip(row, state)) % modulus for row in matrix]


def state(d, n, t, modulus):
    answer = []
    for arg in [t, -t]:
        answer.extend([eval_poly(d[:n], arg, modulus),
                       eval_poly([j*x for j, x in enumerate(d[:n])], arg, modulus)])
    return answer


def teichmueller(a, p, K):
    x, modulus = a, p**K
    for _ in range(K-1):
        x = pow(x, p, modulus)
    assert pow(x, p, modulus) == x
    return x


def wilson_root(p, K):
    x = 1
    for j in range(1, K):
        modulus = p**(j+1)
        candidates = [x+c*p**j for c in range(p)
                      if (pow(x+c*p**j, p-1, modulus)+factorial(p-1)) % modulus == 0]
        assert len(candidates) == 1
        x = candidates[0]
    return x


def build():
    a = apery(2047)
    records, fixed, maps = [], [], []
    naive_counter = None
    for p in PRIMES:
        m = p*p
        aprime, e, beta, b0, w = coefficients(p, a)
        d = unit_coefficients(a, p)
        digit_matrices = []
        for ell in range(p):
            E, B = residue(e[ell], m), residue(beta[ell], m)
            inner = [[E, p*B % m], [ell*E % m, p*(E+ell*B) % m]]
            negative = [[(-1)**ell*x % m for x in row] for row in inner]
            digit_matrices.append(outer_matrix(inner, negative))
        for n in range(1025):
            k, ell = divmod(n, p)
            source = [d[k], k*d[k] % m, (-1)**k*d[k] % m, (-1)**k*k*d[k] % m]
            target = [d[n], n*d[n] % m, (-1)**n*d[n] % m, (-1)**n*n*d[n] % m]
            assert mv(digit_matrices[ell], source, m) == target
        # The Wilson-free inner scheme moves the parameter by sigma_W(t).
        sigma = lambda t: pow(t, p, m)*pow(factorial(p-1), -1, m) % m
        semilinear = []
        for t in range(p):
            r = 0
            while p**(r+1) <= 2048:
                n = p**r
                lhs = state(d, p*n, t, m)[:2]
                source = state(d, n, sigma(t), m)[:2]
                matrix = inner_matrix(e, b0, t, p)
                assert lhs == mv(matrix, source, m)
                semilinear.append(dict(t=t, r=r, sigma=sigma(t), source=source, target=lhs))
                r += 1
        for base in range(p):
            tau = teichmueller(base, p, 2)
            for kind, t, b in [('teichmueller', tau, beta),
                               ('wilson_adapted', wilson_root(p, 2)*tau % m, b0)]:
                matrix = outer_matrix(inner_matrix(e, b, t, p), inner_matrix(e, b, -t, p))
                if kind == 'teichmueller':
                    assert pow(t, p, m) == t
                else:
                    assert sigma(t) == (-t) % m
                states, r = [], 0
                while p**r <= 2048:
                    value = state(d, p**r, t, m)
                    if states:
                        assert mv(matrix, states[-1], m) == value
                    states.append(value)
                    r += 1
                fixed.append(dict(p=p, base=base, kind=kind, t=t, matrix=matrix, states=states))
        # A general integer parameter is not necessarily fixed by t -> t^p mod p^2.
        if naive_counter is None:
            for t in range(2, p):
                matrix = outer_matrix(inner_matrix(e, beta, t, p), inner_matrix(e, beta, -t, p))
                prediction = mv(matrix, state(d, p, t, m), m)
                actual = state(d, p*p, t, m)
                if prediction[0] != actual[0]:
                    naive_counter = dict(p=p, t=t, r=1, incorrectly_frozen=prediction,
                                         actual=actual, t_to_p=pow(t, p, m))
                    break
        records.append(dict(p=p, wilson_quotient_mod_p=w % p,
                            derivative=[str(x) for x in aprime],
                            e_mod_p2=[residue(x, m) for x in e],
                            beta_mod_p=[residue(x, p) for x in beta],
                            beta0_mod_p=[residue(x, p) for x in b0],
                            digit_matrices=digit_matrices,
                            inner_plus=inner_matrix(e, beta, 1, p),
                            inner_minus=inner_matrix(e, beta, -1, p),
                            digit_indices=1025, semilinear=semilinear))
        for K in range(1, 7):
            modulus = p**K
            eta = wilson_root(p, K)
            assert pow(eta, p-1, modulus) == -factorial(p-1) % modulus
            lifts = []
            for base in range(p):
                t = eta*teichmueller(base, p, K) % modulus
                assert pow(t, p, modulus)*pow(factorial(p-1), -1, modulus) % modulus == (-t) % modulus
                lifts.append(t)
            maps.append(dict(p=p, K=K, eta=eta, adapted_lifts=lifts))
    assert naive_counter is not None
    return dict(schema_version=1, status='fixed_precision_consequence_of_gessel_and_wilson',
                scope=dict(primes=PRIMES, digit_max_index=1024, tower_max_length=2048,
                           precision=2, wilson_parameter_precision_max=6,
                           note='All-index mod-p^2 scheme uses Gessel theorem. No mod-p^3 or G1/G2 upgrade.'),
                records=records, fixed_towers=fixed, wilson_parameter_lifts=maps,
                negative_controls=dict(non_teichmueller_frozen_parameter=naive_counter,
                                       old_scheme_mod49=dict(p=7, n=7, old=44, actual=16)),
                counts=dict(digit_indices=len(PRIMES)*1025, fixed_parameter_cases=len(fixed),
                            fixed_transitions=sum(len(x['states'])-1 for x in fixed),
                            semilinear_transitions=sum(len(x['semilinear']) for x in records),
                            wilson_root_cases=len(maps),
                            wilson_sign_flip_cases=sum(len(x['adapted_lifts']) for x in maps)),
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                theorem_source='https://ericrowland.github.io/papers/Lucas_congruences_for_the_Apery_numbers_modulo_p%5E2.pdf')


if __name__ == '__main__':
    c = build()
    output = HERE/'wilson_lift2_certificate.json'
    output.write_text(json.dumps(c, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(c['counts']))
    print('counterexample', c['negative_controls']['non_teichmueller_frozen_parameter'])
    print('certificate_sha256', sha256(output.read_bytes()).hexdigest())
