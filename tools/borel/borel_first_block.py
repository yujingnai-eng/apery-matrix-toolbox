"""Apéry Borel first blocks and the Wilson-signed unit-factorial tower.

Standard library only. Algebraic proofs are in borel_first_block.md.
Finite certificates do not establish CM identifications or G1/G2.
"""
from fractions import Fraction as Q
from hashlib import sha256
from math import factorial
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent


def primes(limit):
    return [n for n in range(3, limit+1)
            if all(n % d for d in range(2, int(n**0.5)+1))]


def apery(limit):
    a = [1, 5]
    for n in range(1, limit):
        num = (2*n+1)*(17*n*n+17*n+5)*a[n]-n**3*a[n-1]
        z, remainder = divmod(num, (n+1)**3)
        assert remainder == 0
        a.append(z)
    return a[:limit+1]


def trim(a):
    while a and not a[-1]:
        a.pop()
    return a


def remainder(a, b, p):
    a = a[:]
    inv = pow(b[-1], -1, p)
    while len(a) >= len(b):
        c, j = a[-1]*inv % p, len(a)-len(b)
        for k, x in enumerate(b):
            a[k+j] = (a[k+j]-c*x) % p
        trim(a)
    return a


def gcd(a, b, p):
    while b:
        a, b = b, remainder(a, b, p)
    return [x*pow(a[-1], -1, p) % p for x in a]


def evaluate(c, t, p):
    z = 0
    for x in reversed(c):
        z = (z*t+x) % p
    return z


def matrix_multiply(a, b):
    return [[sum(x*y for x, y in zip(row, col))
             for col in zip(*b)] for row in a]


def full_step(n, t):
    p = (2*n+1)*(17*n*n+17*n+5)
    return [[Q(t*p, (n+1)**4), Q(-t*t*n*n, (n+1)**4), 0, 0],
            [1, 0, 0, 0], [1, 0, 1, 0], [(-1)**n, 0, 0, 1]]


def encode_matrix(m):
    return [[str(x) for x in row] for row in m]


def polynomial_add(*polys):
    out = [0]*max(map(len, polys))
    for poly in polys:
        for k, value in enumerate(poly):
            out[k] += value
    return trim(out) or [0]


def polynomial_multiply(a, b):
    out = [0]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i+j] += x*y
    return trim(out) or [0]


def polynomial_mm(a, b):
    return [[polynomial_add(*(polynomial_multiply(a[i][k], b[k][j])
                              for k in range(4))) for j in range(4)] for i in range(4)]


def quantum_gauge():
    # Low-to-high polynomial coefficients. z is the Borel variable, not depth r.
    c = [[[0], [1], [0], [0]], [[0], [0], [1], [0]],
         [[0], [0], [0], [1]],
         [[0, 5, -1], [0, 27, -2], [0, 51, -1], [0, 34]]]
    g = [[[0, -5, -12, -180], [0, -22, -48], [0, -29], [1]],
         [[0, -5, -36], [0, -17], [1], [0]],
         [[0, -5], [1], [0], [0]], [[1], [0], [0], [0]]]
    m = [[[0, 5], [0, 0, 96], [0, 0, 0, 1692], [0, 0, 0, 0, 12816]],
         [[1], [0, 12], [0, 0, 216], [0, 0, 0, 1692]],
         [[0], [1], [0, 12], [0, 0, 96]], [[0], [0], [1], [0, 5]]]
    gc, mg = polynomial_mm(g, c), polynomial_mm(m, g)
    for i in range(4):
        for j in range(4):
            theta_g = [k*x for k, x in enumerate(g[i][j])]
            assert polynomial_add(theta_g, gc[i][j], [-x for x in mg[i][j]]) == [0]
    return dict(companion=c, change_of_basis=g, quantum=m,
                change_of_basis_determinant=1,
                truncated_source_direction=[1, 0, 0, 0],
                reference='Golyshev-Zagier, Gamma Conjecture, section 2, first displayed 4x4 matrix',
                reference_url='https://people.mpim-bonn.mpg.de/zagier/files/tex/GammaConjecture/GammaConj.pdf')


def ode_boundary(a, n):
    c = [Q(a[k], factorial(k)) for k in range(n)]
    out = [Q(0)]*(n+2)
    for k, x in enumerate(c):
        out[k] += k**4*x
        out[k+1] -= (2*k+1)*(17*k*k+17*k+5)*x
        out[k+2] += (k+1)**2*x
    expected = [Q(0)]*(n+2)
    expected[n] = Q(-n**3*a[n], factorial(n-1))
    expected[n+1] = Q(n*n*a[n-1], factorial(n-1))
    assert out == expected
    return [str(x) for x in out]


def build():
    a = apery(2048)
    first = []
    for p in primes(199):
        unit, c = 1, []
        for n in range(p):
            if n:
                unit = unit*n % p
            c.append(a[n]*pow(unit, -1, p) % p)
        h = [x % p for x in a[:p]]
        g = gcd(c, [i*c[i] % p for i in range(1, p)], p)
        partner = [a[n]*factorial(n) % p for n in range(p)]
        assert c[::-1] == [(-(-1)**n*partner[n]) % p for n in range(p)]
        roots = [t for t in range(p) if evaluate(c, t, p) == 0]
        first.append(dict(p=p, coefficients=c, ordinary_coefficients=h,
                          roots=roots, repeated_gcd=g,
                          h_pm=[evaluate(h, 1, p), evaluate(h, -1, p)],
                          e_pm=[evaluate(c, 1, p), evaluate(c, -1, p)],
                          factorial_partner=partner))
    digits = []
    unit_corrections = []
    for p in [3, 5, 7, 11, 13, 17, 19]:
        unit, d = 1, []
        for n in range(1025):
            if n:
                x = n
                while x % p == 0:
                    x //= p
                unit = unit*x % p
            d.append(a[n]*pow(unit, -1, p) % p)
        c = next(row['coefficients'] for row in first if row['p'] == p)
        for n, value in enumerate(d):
            k, ell = divmod(n, p)
            assert value == ((-1)**k*c[ell]*d[k]) % p
        digits.append(dict(p=p, residues=d))
        wp = (factorial(p-1)+1)//p
        harmonics = [sum(pow(j, -1, p) for j in range(1, ell+1)) % p for ell in range(p)]
        units = []
        value = 1
        for n in range(513):
            if n % p:
                value = value*n % (p*p)
            k, ell = divmod(n, p)
            predicted = (-1)**k*factorial(ell)*(1+p*k*(harmonics[ell]-wp)) % (p*p)
            assert value == predicted
            units.append(value)
        unit_corrections.append(dict(p=p, wilson_quotient_mod_p=wp % p,
                                     harmonic_mod_p=harmonics, p_free_product_mod_p2=units))
    towers = []
    for p in primes(31):
        unit, d = 1, []
        for n in range(2048):
            if n:
                x = n
                while x % p == 0:
                    x //= p
                unit = unit*x % p
            d.append(a[n]*pow(unit, -1, p) % p)
        c = next(row['coefficients'] for row in first if row['p'] == p)
        for t in range(p):
            eplus, eminus = evaluate(c, t, p), evaluate(c, -t, p)
            previous, states, r = [1, 1], [[1, 1]], 1
            while p**r <= 2048:
                n = p**r
                direct = [sum(d[j]*pow(t, j, p) for j in range(n)) % p,
                          sum(d[j]*pow(-t, j, p) for j in range(n)) % p]
                assert direct == [eplus*previous[1] % p, eminus*previous[0] % p]
                states.append(direct)
                previous, r = direct, r+1
            towers.append(dict(p=p, t=t, matrix=[[0, eplus], [eminus, 0]],
                               states=states))
    blocks = []
    for t in [1, -1, 2]:
        for start, end in [(1, 2), (1, 4), (2, 8), (1, 12)]:
            m = [[int(i == j) for j in range(4)] for i in range(4)]
            for n in range(start, end):
                m = matrix_multiply(full_step(n, t), m)
            state = [Q(t**start*a[start], factorial(start)),
                     Q(t**(start-1)*a[start-1], factorial(start-1)),
                     sum(Q(t**k*a[k], factorial(k)) for k in range(start)),
                     sum(Q((-t)**k*a[k], factorial(k)) for k in range(start))]
            out = [sum(x*y for x, y in zip(row, state)) for row in m]
            target = [Q(t**end*a[end], factorial(end)),
                      Q(t**(end-1)*a[end-1], factorial(end-1)),
                      sum(Q(t**k*a[k], factorial(k)) for k in range(end)),
                      sum(Q((-t)**k*a[k], factorial(k)) for k in range(end))]
            assert out == target
            blocks.append(dict(t=t, start=start, end=end, matrix=encode_matrix(m),
                               input=[str(x) for x in state],
                               output=[str(x) for x in out]))
    boundaries = [dict(p=p, coefficients=ode_boundary(a, p)) for p in primes(31)]
    e8 = sum(Q(a[k], factorial(k)) for k in range(8))
    assert e8 == Q(6477413, 42)
    return dict(schema_version=1, status='finite_observation_with_algebraic_derivations',
                scope=dict(first_block_prime_max=199, first_block_prime_min=3,
                           digit_max_index=1024, tower_max_length=2048,
                           excludes='CM identification, growing precision, G1/G2 proof'),
                first_blocks=first, digital_relations=digits, towers=towers,
                unit_factorial_corrections=unit_corrections,
                full_blocks=blocks, ode_boundaries=boundaries,
                differential_gauge=quantum_gauge(),
                counterexamples=dict(nonintegral_E8_at_7=str(e8),
                                     E4_minus=str(sum(Q((-1)**k*a[k], factorial(k)) for k in range(4))),
                                     naive_mod49_digit_p7_n7=dict(actual=16, predicted_without_correction=44),
                                     split_p17_t1_matrix=[[0, 15], [0, 0]],
                                     p5_repeated_factor=[3, 0, 1]),
                counts=dict(first_blocks=len(first),
                            digital_indices=sum(len(x['residues']) for x in digits),
                            tower_parameter_pairs=len(towers),
                            tower_transitions=sum(len(x['states'])-1 for x in towers),
                            full_blocks=len(blocks), ode_boundaries=len(boundaries)),
                source_sha256=sha256(Path(__file__).read_bytes()).hexdigest())


if __name__ == '__main__':
    cert = build()
    path = HERE/'borel_first_block_certificate.json'
    path.write_text(json.dumps(cert, ensure_ascii=False, indent=2)+'\n')
    print(json.dumps(cert['counts']))
    print('certificate_sha256', sha256(path.read_bytes()).hexdigest())
