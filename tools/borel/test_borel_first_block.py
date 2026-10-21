"""Independent binomial-kernel verification; imports no project calculation code."""
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial
from pathlib import Path
import json
import unittest

HERE = Path(__file__).resolve().parent


def direct_apery(n):
    # u_k = binom(n,k) binom(n+k,k); exact k-axis updates, not n recurrence.
    u, total = 1, 1
    for k in range(n):
        num = u*(n-k)*(n+k+1)
        u, rest = divmod(num, (k+1)**2)
        assert rest == 0
        total += u*u
    return total


def ps(limit):
    return [n for n in range(3, limit+1)
            if not any(n % d == 0 for d in range(2, n))]


def legendre(n, p):
    result = 0
    while n:
        n //= p
        result += n
    return result


def poly_gcd(a, b, p):
    def clean(x):
        while x and x[-1] % p == 0:
            x.pop()
        return x
    a, b = clean(a[:]), clean(b[:])
    while b:
        while len(a) >= len(b):
            i, multiplier = len(a)-len(b), a[-1]*pow(b[-1], -1, p) % p
            for j in range(len(b)):
                a[i+j] = (a[i+j]-multiplier*b[j]) % p
            clean(a)
        a, b = b, a
    return [int(x*pow(a[-1], -1, p) % p) for x in a]


def evaluate(c, t, p):
    return sum(x*pow(t, n, p) for n, x in enumerate(c)) % p


def check(cert, a, facts):
    assert cert['schema_version'] == 1
    assert cert['source_sha256'] == sha256((HERE/'borel_first_block.py').read_bytes()).hexdigest()
    assert [row['p'] for row in cert['first_blocks']] == ps(199)
    for row in cert['first_blocks']:
        p = row['p']
        c = [a[n]*pow(facts[n] % p, -1, p) % p for n in range(p)]
        h = [a[n] % p for n in range(p)]
        assert row['coefficients'] == c
        assert row['ordinary_coefficients'] == h
        assert row['h_pm'] == [evaluate(h, 1, p), evaluate(h, -1, p)]
        assert row['e_pm'] == [evaluate(c, 1, p), evaluate(c, -1, p)]
        assert row['roots'] == [t for t in range(p) if evaluate(c, t, p) == 0]
        assert row['repeated_gcd'] == poly_gcd(c, [j*c[j] % p for j in range(1, p)], p)
        assert row['factorial_partner'] == [a[n]*facts[n] % p for n in range(p)]
        assert c[::-1] == [(-(-1)**n*a[n]*facts[n]) % p for n in range(p)]
    assert [row['p'] for row in cert['digital_relations']] == [3, 5, 7, 11, 13, 17, 19]
    for row in cert['digital_relations']:
        p = row['p']
        d = [a[n]*pow((facts[n]//p**legendre(n, p)) % p, -1, p) % p
             for n in range(1025)]
        assert row['residues'] == d
        for n in range(1025):
            k, ell = divmod(n, p)
            assert d[n] == ((-1)**k*a[ell]*pow(facts[ell] % p, -1, p)*d[k]) % p
    assert [row['p'] for row in cert['unit_factorial_corrections']] == [3, 5, 7, 11, 13, 17, 19]
    for row in cert['unit_factorial_corrections']:
        p = row['p']
        assert row['wilson_quotient_mod_p'] == (facts[p-1]+1)//p % p
        assert row['harmonic_mod_p'] == [sum(pow(j, -1, p) for j in range(1, ell+1)) % p for ell in range(p)]
        # Exact quotient of integer factorials, not the producer's running product.
        expected = [facts[n]//(p**(n//p)*facts[n//p]) % (p*p) for n in range(513)]
        assert row['p_free_product_mod_p2'] == expected
        for n, value in enumerate(expected):
            k, ell = divmod(n, p)
            h, w = row['harmonic_mod_p'][ell], row['wilson_quotient_mod_p']
            assert value == (-1)**k*facts[ell]*(1+p*k*(h-w)) % (p*p)
    assert [(row['p'], row['t']) for row in cert['towers']] == [
        (p, t) for p in ps(31) for t in range(p)]
    ds = {p: [a[n]*pow((facts[n]//p**legendre(n, p)) % p, -1, p) % p
              for n in range(2048)] for p in ps(31)}
    for row in cert['towers']:
        p, t = row['p'], row['t']
        c = [a[n]*pow(facts[n] % p, -1, p) % p for n in range(p)]
        v, w = evaluate(c, t, p), evaluate(c, -t, p)
        assert row['matrix'] == [[0, v], [w, 0]]
        expected, r = [], 0
        while p**r <= 2048:
            n = p**r
            expected.append([sum(ds[p][k]*pow(sign*t, k, p) for k in range(n)) % p
                             for sign in (1, -1)])
            r += 1
        assert row['states'] == expected
    assert [(x['t'], x['start'], x['end']) for x in cert['full_blocks']] == [
        (t, i, j) for t in [1, -1, 2] for i, j in [(1, 2), (1, 4), (2, 8), (1, 12)]]
    for row in cert['full_blocks']:
        t, start, end = row['t'], row['start'], row['end']
        matrix = [[Fraction(x) for x in r] for r in row['matrix']]
        for column in range(4):
            now, before, plus, minus = [Fraction(i == column) for i in range(4)]
            for n in range(start, end):
                plus, minus = plus+now, minus+(-1)**n*now
                now, before = (t*(34*n**3+51*n*n+27*n+5)*now-t*t*n*n*before)/(n+1)**4, now
            assert [now, before, plus, minus] == [matrix[i][column] for i in range(4)]
        for n, key in [(start, 'input'), (end, 'output')]:
            expected = [Fraction(t**n*a[n], facts[n]), Fraction(t**(n-1)*a[n-1], facts[n-1])]
            expected += [sum(Fraction((s*t)**k*a[k], facts[k]) for k in range(n)) for s in [1, -1]]
            assert [Fraction(x) for x in row[key]] == expected
    assert [r['p'] for r in cert['ode_boundaries']] == ps(31)
    for row in cert['ode_boundaries']:
        p = row['p']
        def c(j):
            return Fraction(a[j], facts[j]) if 0 <= j < p else Fraction(0)
        coefficients = []
        for n in range(p+2):
            coeff = n**4*c(n)-(34*(n-1)**3+51*(n-1)**2+27*(n-1)+5)*c(n-1)+(n-1)**2*c(n-2)
            coefficients.append(coeff)
        assert [str(x) for x in coefficients] == row['coefficients']
        assert coefficients[:p] == [0]*p
        assert coefficients[p:] == [Fraction(-p**3*a[p], facts[p-1]), Fraction(p*p*a[p-1], facts[p-1])]
    gauge = cert['differential_gauge']
    assert gauge['change_of_basis_determinant'] == 1
    assert gauge['truncated_source_direction'] == [1, 0, 0, 0]
    # The degree caps below imply degree <=8. Nine exact evaluations
    # certify its polynomial vanishing, independently of producer convolution.
    g, c, m = [gauge[k] for k in ('change_of_basis', 'companion', 'quantum')]
    assert max(len(poly) for mat in (g, c, m) for row in mat for poly in row) <= 5
    for z in range(9):
        def ev(poly):
            return sum(value*z**k for k, value in enumerate(poly))
        G, C, M = [[[ev(poly) for poly in row] for row in mat] for mat in (g, c, m)]
        assert [row[3] for row in G] == [1, 0, 0, 0]
        for i in range(4):
            for j in range(4):
                derivative = sum(k*value*z**k for k, value in enumerate(g[i][j]))
                assert derivative+sum(G[i][k]*C[k][j]-M[i][k]*G[k][j] for k in range(4)) == 0
        # Determinant via permutation expansion, no symbolic algebra dependency.
        from itertools import permutations
        determinant = 0
        for perm in permutations(range(4)):
            inversions = sum(perm[i] > perm[j] for i in range(4) for j in range(i+1, 4))
            term = (-1)**inversions
            for i in range(4):
                term *= G[i][perm[i]]
            determinant += term
        assert determinant == 1
    assert cert['counterexamples']['nonintegral_E8_at_7'] == str(sum(Fraction(a[k], facts[k]) for k in range(8)))
    assert cert['counterexamples']['E4_minus'] == '-625/3'
    assert a[7]*pow(facts[7]//7, -1, 49) % 49 == 16
    assert cert['counterexamples']['naive_mod49_digit_p7_n7'] == dict(actual=16, predicted_without_correction=44)
    assert cert['counterexamples']['split_p17_t1_matrix'] == [[0, 15], [0, 0]]
    assert cert['counterexamples']['p5_repeated_factor'] == [3, 0, 1]
    assert cert['counts'] == dict(first_blocks=len(ps(199)), digital_indices=7*1025,
                                 tower_parameter_pairs=sum(ps(31)),
                                 tower_transitions=sum(len(x['states'])-1 for x in cert['towers']),
                                 full_blocks=12, ode_boundaries=len(ps(31)))


class Verification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a = [direct_apery(n) for n in range(2048)]
        cls.facts = [factorial(n) for n in range(2048)]
        cls.cert = json.loads((HERE/'borel_first_block_certificate.json').read_text())
        for n in range(60):
            assert cls.a[n] == sum(comb(n, k)**2*comb(n+k, k)**2 for k in range(n+1))

    def test_certificate(self):
        check(self.cert, self.a, self.facts)

    def test_negative_controls(self):
        mutations = [
            lambda c: c['first_blocks'].pop(),
            lambda c: c['first_blocks'][1]['coefficients'].__setitem__(1, 1),
            lambda c: c['first_blocks'][1].__setitem__('repeated_gcd', [1]),
            lambda c: c['digital_relations'][0]['residues'].__setitem__(3, 0),
            lambda c: c['towers'][0]['matrix'][0].__setitem__(0, 1),
            lambda c: c['ode_boundaries'][0]['coefficients'].__setitem__(-1, '0'),
            lambda c: c['full_blocks'][0]['matrix'][0].__setitem__(0, '0'),
            lambda c: c['differential_gauge']['change_of_basis'][0][0].__setitem__(3, -181),
        ]
        for mutation in mutations:
            bad = deepcopy(self.cert)
            mutation(bad)
            with self.assertRaises(AssertionError):
                check(bad, self.a, self.facts)


if __name__ == '__main__':
    suite = unittest.defaultTestLoader.loadTestsFromTestCase(Verification)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    output = dict(schema_version=1, status='pass' if result.wasSuccessful() else 'fail',
                  source='independent k-axis binomial kernel, no generator imports',
                  binomial_indices=2048, literal_comb_crosschecks=60,
                  certificate_sha256=sha256((HERE/'borel_first_block_certificate.json').read_bytes()).hexdigest(),
                  generator_sha256=sha256((HERE/'borel_first_block.py').read_bytes()).hexdigest(),
                  verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                  tests_run=result.testsRun, mutation_cases=8,
                  limits='Finite checks; universal statements require the separate algebraic derivations. No G1/G2 claim.')
    (HERE/'borel_first_block_verification.json').write_text(json.dumps(output, ensure_ascii=False, indent=2)+'\n')
    raise SystemExit(not result.wasSuccessful())
