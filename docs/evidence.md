# Evidence and boundaries

## 1. Status vocabulary

| Status | Meaning |
|---|---|
| exact / proved | exact identity or proved statement |
| fixed precision | fixed \(p^K\), declared range |
| finite observation | finite data, no all-\(r\) theorem |
| conjecture | open target |
| refuted | explicit counterexample |
| superseded | replaced by a later audit |

## 2. Beukers Appendix C

- 15 families; 9.4 and 9.6 have no points; 13 non-empty sections.
- 118 CM points; 44 primes \(5\le p\le199\).
- 5,192 point–prime grid; 5,130 legacy eligible; 5,092 strict non-ramified.

First-block classification on 5,130 cases:

| type | cases | \(v_p(H_p)\) |
|---|---:|---|
| split | 2,433 | 0 |
| ramified | 38 | 1 |
| inert | 2,659 | 2: 2,592; 3: 60; 4: 6; 5: 1 |

The 38 legacy ramified cases also vanish. Therefore the naive statement

> \(H_p(t_0)=0\) if and only if \(p\) is inert

is false without non-ramified/good-reduction hypotheses.

## 3. Geometric vs rational roots

Let \(H_D^{\rm CM}(X)\) be the actual order class polynomial and

\[
\mathrm{SS}_p(X)=\prod_{j\text{ supersingular}}(X-j)
\]

the supersingular polynomial over \(\mathbb F_p\). The geometric condition is

\[
\gcd(H_D^{\rm CM}\bmod p,\mathrm{SS}_p)\neq1.
\]

The rational condition is

\[
\gcd\bigl(\gcd(H_D^{\rm CM},\mathrm{SS}_p),X^p-X\bigr)\neq1.
\]

They are not equivalent. In the strict non-ramified grid:

- 5,092 cases: first-block vanishing agrees with the geometric condition, 0 mismatches;
- **561** cases have a geometric common root but no \(\mathbb F_p\) root.

The 561 cases are the main reason “supersingular” cannot be read as “\(j\in\mathbb F_p\)”.

## 4. Order class numbers

Field class numbers: 66/40/12 for \(h=1/2/4\).

Actual order class numbers:

\[
\boxed{43/54/16/5}
\qquad(1/2/4/8).
\]

Among 118 points:

- 46 non-maximal orders;
- 29 class-number changes;
- 23 of the old 66 “class number one” points have actual \(j(\tau)\) degree \(>1\);
- 29 old table \(j\) values are correct; 14 are wrong integers.

The old “2,884 pointwise \(j\) test” is superseded: it used a field representative \(j\), not the actual \(j(\tau)\).

## 5. Local dyadic blocks

- 8,190 aligned blocks with \(s\le12\);
- \(v_2(B_{a,s})=3s-1\), 0 exceptions;
- \(v_2(b_r-7f_r)=4r+2\) for \(r\le12\);
- \(C=[[6,16],[1,0]]\bmod128\) for \(r\ge1\); \(r=0\) fails with residual 40.

A fixed matrix modulo 128 does not determine the integer valuation. The repository keeps the distinction between congruence and valuation.

## 6. G2

Candidate:

\[
C=\begin{pmatrix}T&-D\\1&0\end{pmatrix},
\qquad
D=\chi p^2.
\]

For inert \(p\), the candidate has \(T=0\); for split \(p=a^2+db^2\), \(T=4a^2-2p\).

The target is

\[
\Pi_{r+1}\mathcal B_r x
\equiv
C\Pi_r x
\pmod{p^{3r+5}}.
\]

What is available:

- candidate form;
- fixed-precision digital states;
- \(p=3\) K18 certificate for two families, all lengths, \(j=0,\ldots,5\), modulo \(3^{18}\).

What is missing:

- derivation of \(T\);
- the 2-dimensional symmetric-square quotient and readout;
- the growing modulus \(p^{3r+5}\);
- an all-\(r\) allowed lattice.

K18 is not G2.

## 7. What is not claimed

- G1 and G2 are not proved.
- The first block is not proved to equal the Hasse invariant for every CM point.
- The 561 geometric-only, 38 ramified, and 67 higher-inert cases are not classified.
- The CM audit is finite (118 points), not a classification theorem.
- The Wilson–Gessel and Borel constructions are different towers from the original \(F_N\).

## 8. Why keep negative results

The repository keeps:

- the \(r=0\) failure of the mod-128 matrix;
- the 38 ramified cases;
- the 561 geometric/rational mismatch;
- the old \(j\) test as superseded;
- the distinction between fixed precision and growing precision.

These are not noise. They are the boundary of the current theorem.
