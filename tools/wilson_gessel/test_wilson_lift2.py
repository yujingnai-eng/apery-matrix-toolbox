"""Independent source, derivative, factorial, and matrix checks for Wilson lift.

No imports from the producer or the project library. Universal conclusions use
the Wilson derivation and Gessel input theorem, not these finite checks alone.
"""
from copy import deepcopy
from fractions import Fraction as Q
from hashlib import sha256
from math import comb, factorial
from pathlib import Path
import json
import unittest

HERE = Path(__file__).resolve().parent
PRIMES = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]


def direct_a(n):
    v, result = 1, 1
    for k in range(n):
        v, rem = divmod(v*(n-k)*(n+k+1), (k+1)**2)
        assert rem == 0
        result += v*v
    return result


def derivative(n):
    # Differentiate finite polynomial products, without harmonic-number input.
    total = Q(0)
    for k in range(1, n+1):
        value, tangent = 1, 0
        for j in range(-k+1, k+1):
            value, tangent = value*(n+j), tangent*(n+j)+value
        total += Q(2*value*tangent, factorial(k)**4)
    return total


def md(q, modulus):
    q = Q(q)
    return q.numerator*pow(q.denominator, -1, modulus) % modulus


def vp_factorial(n, p):
    out = 0
    while n:
        n //= p
        out += n
    return out


def read_state(d, length, t, modulus):
    result = []
    for z in (t, -t):
        powers = [pow(z, n, modulus) for n in range(length)]
        result.extend([sum(x*y for x, y in zip(d, powers)) % modulus,
                       sum(n*d[n]*powers[n] for n in range(length)) % modulus])
    return result


def action(m, x, modulus):
    return [sum(m[i][j]*x[j] for j in range(len(x))) % modulus for i in range(len(m))]


def check(cert, a, facts, derivatives):
    assert cert['schema_version'] == 1
    assert cert['source_sha256'] == sha256((HERE/'wilson_lift2.py').read_bytes()).hexdigest()
    assert [x['p'] for x in cert['records']] == PRIMES
    all_d, matrices = {}, {}
    for row in cert['records']:
        p, m = row['p'], row['p']**2
        w = (facts[p-1]+1)//p
        e = [md(Q(a[l], facts[l]), m) for l in range(p)]
        H = [sum((Q(1, j) for j in range(1, l+1)), Q(0)) for l in range(p)]
        b0 = [md((derivatives[l]-H[l]*a[l])/facts[l], p) for l in range(p)]
        b = [(b0[l]+w*e[l]) % p for l in range(p)]
        assert row['derivative'] == [str(x) for x in derivatives[:p]]
        assert row['wilson_quotient_mod_p'] == w % p
        assert row['e_mod_p2'] == e and row['beta_mod_p'] == b and row['beta0_mod_p'] == b0
        expected = []
        for l in range(p):
            f, g, s = e[l], p*b[l] % m, (-1)**l
            expected.append([[0, 0, f, g], [0, 0, l*f % m, p*(f+l*b[l]) % m],
                             [s*f % m, s*g % m, 0, 0],
                             [s*l*f % m, s*p*(f+l*b[l]) % m, 0, 0]])
        assert row['digit_matrices'] == expected
        matrices[p] = expected
        d = [a[n]*pow((facts[n]//p**vp_factorial(n, p)) % m, -1, m) % m for n in range(2048)]
        all_d[p] = d
        assert row['digit_indices'] == 1025
        for n in range(1025):
            k, ell = divmod(n, p)
            src = [d[k], k*d[k] % m, (-1)**k*d[k] % m, (-1)**k*k*d[k] % m]
            dst = [d[n], n*d[n] % m, (-1)**n*d[n] % m, (-1)**n*n*d[n] % m]
            assert action(expected[ell], src, m) == dst
        for sign, key in [(1, 'inner_plus'), (-1, 'inner_minus')]:
            full = [[sum(sign**l*expected[l][i][j] for l in range(p)) % m
                     for j in range(4)] for i in range(4)]
            assert row[key] == [r[2:] for r in full[:2]]
        valid_pairs = [(t, r) for t in range(p) for r in range(8) if p**(r+1) <= 2048]
        assert [(x['t'], x['r']) for x in row['semilinear']] == valid_pairs
        for x in row['semilinear']:
            t, r = x['t'], x['r']
            sigma = pow(t, p, m)*pow(facts[p-1], -1, m) % m
            source = read_state(d, p**r, sigma, m)[:2]
            target = read_state(d, p**(r+1), t, m)[:2]
            B = [[sum(e[l]*pow(t, l, m) for l in range(p)) % m,
                  p*sum(b0[l]*pow(t, l, m) for l in range(p)) % m],
                 [sum(l*e[l]*pow(t, l, m) for l in range(p)) % m,
                  p*sum((e[l]+l*b0[l])*pow(t, l, m) for l in range(p)) % m]]
            assert x['source'] == source and x['target'] == target and x['sigma'] == sigma
            assert action(B, source, m) == target
    assert [(x['p'], x['base'], x['kind']) for x in cert['fixed_towers']] == [
        (p, base, kind) for p in PRIMES for base in range(p) for kind in ['teichmueller', 'wilson_adapted']]
    for row in cert['fixed_towers']:
        p, t, m = row['p'], row['t'], row['p']**2
        assert 0 <= t < m and t % p == row['base']
        digits = matrices[p]
        if row['kind'] == 'teichmueller':
            assert pow(t, p, m) == t
        else:
            assert pow(t, p, m)*pow(facts[p-1], -1, m) % m == (-t) % m
            # Remove just the Wilson quotient term from both inner blocks.
            w = (facts[p-1]+1)//p
            source_row = next(x for x in cert['records'] if x['p'] == p)
            digits = deepcopy(digits)
            for ell, f in enumerate(source_row['e_mod_p2']):
                correction = p*w*f
                for i, j, factor in [(0, 3, 1), (1, 3, ell), (2, 1, (-1)**ell), (3, 1, (-1)**ell*ell)]:
                    digits[ell][i][j] = (digits[ell][i][j]-factor*correction) % m
        M = [[sum(pow(t, ell, m)*digits[ell][i][j] for ell in range(p)) % m
              for j in range(4)] for i in range(4)]
        assert row['matrix'] == M
        expected = [read_state(all_d[p], p**r, t, m) for r in range(8) if p**r <= 2048]
        assert row['states'] == expected
        for left, right in zip(expected, expected[1:]):
            assert action(M, left, m) == right
    assert [(x['p'], x['K']) for x in cert['wilson_parameter_lifts']] == [
        (p, K) for p in PRIMES for K in range(1, 7)]
    previous = {}
    for row in cert['wilson_parameter_lifts']:
        p, K, eta = row['p'], row['K'], row['eta']
        modulus = p**K
        assert eta % p == 1 and pow(eta, p-1, modulus) == -facts[p-1] % modulus
        if K > 1:
            assert eta % (p**(K-1)) == previous[p]
        previous[p] = eta
        assert len(row['adapted_lifts']) == p
        for base, t in enumerate(row['adapted_lifts']):
            assert t % p == base
            assert pow(t, p, modulus)*pow(facts[p-1], -1, modulus) % modulus == (-t) % modulus
    counter = cert['negative_controls']['non_teichmueller_frozen_parameter']
    p, t, m = counter['p'], counter['t'], counter['p']**2
    M = [[sum(pow(t, ell, m)*matrices[p][ell][i][j] for ell in range(p)) % m
          for j in range(4)] for i in range(4)]
    actual = read_state(all_d[p], p*p, t, m)
    wrong = action(M, read_state(all_d[p], p, t, m), m)
    assert counter['incorrectly_frozen'] == wrong and counter['actual'] == actual
    assert wrong[0] != actual[0] and counter['t_to_p'] == pow(t, p, m)
    assert cert['negative_controls']['old_scheme_mod49'] == dict(p=7, n=7, old=44, actual=16)
    assert all_d[7][7] == 16
    assert cert['counts'] == dict(digit_indices=10250, fixed_parameter_cases=316,
                                 fixed_transitions=712, semilinear_transitions=356,
                                 wilson_root_cases=60, wilson_sign_flip_cases=948)


class Verification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.a = [direct_a(n) for n in range(2048)]
        cls.facts = [factorial(n) for n in range(2048)]
        cls.derivatives = [derivative(n) for n in range(31)]
        cls.cert = json.loads((HERE/'wilson_lift2_certificate.json').read_text())
        for n in range(60):
            assert cls.a[n] == sum(comb(n, k)**2*comb(n+k, k)**2 for k in range(n+1))

    def test_certificate(self):
        check(self.cert, self.a, self.facts, self.derivatives)

    def test_negative_controls(self):
        mutations = [lambda c: c['records'].pop(),
                     lambda c: c['records'][2]['beta_mod_p'].__setitem__(0, 0),
                     lambda c: c['records'][0]['digit_matrices'][0][0].__setitem__(2, 0),
                     lambda c: c['fixed_towers'].pop(),
                     lambda c: c['fixed_towers'][4]['states'][-1].__setitem__(0, 0),
                     lambda c: c['wilson_parameter_lifts'][1].__setitem__('eta', 1),
                     lambda c: c['records'][0]['semilinear'][0].__setitem__('sigma', 1)]
        for mutation in mutations:
            bad = deepcopy(self.cert)
            mutation(bad)
            with self.assertRaises(AssertionError):
                check(bad, self.a, self.facts, self.derivatives)


if __name__ == '__main__':
    result = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Verification))
    receipt = dict(schema_version=1, status='pass' if result.wasSuccessful() else 'fail',
                   independent_methods=['binomial k-axis kernel', 'polynomial product differentiation',
                                        'integer factorial quotient', 'digit-matrix weighted sums'],
                   certificate_sha256=sha256((HERE/'wilson_lift2_certificate.json').read_bytes()).hexdigest(),
                   generator_sha256=sha256((HERE/'wilson_lift2.py').read_bytes()).hexdigest(),
                   verifier_sha256=sha256(Path(__file__).read_bytes()).hexdigest(),
                   source_binomial_terms=2048, literal_comb_checks=60, tests_run=result.testsRun,
                   mutation_cases=7, scope='fixed mod-p^2 D tower; no growing precision or G1/G2 certification')
    (HERE/'wilson_lift2_verification.json').write_text(json.dumps(receipt, ensure_ascii=False, indent=2)+'\n')
    raise SystemExit(not result.wasSuccessful())
