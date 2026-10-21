# Overview

## 1. The problem

The Apéry numbers

\[
A_n=\sum_{k=0}^{n}\binom nk^2\binom{n+k}{k}^2
\]

satisfy a three-term recurrence and their generating function is D-finite. Beukers placed them in a supercongruence framework: for a truncated sum

\[
F_N(t)=\sum_{n<N}A_n t^n,
\]

one asks how \(F_{p^r}(t)\) behaves \(p\)-adically as the depth \(r\) grows.

The project studies a stronger structural question:

> Can the exact \(n\)-axis transport be compressed, on an allowed lattice and after a readout, to a small depth matrix whose modulus grows with \(r\)?

The central comparison is

\[
\Pi_{r+1}\mathcal B_r x
\equiv
C_r\Pi_r x
\pmod{p^{K(r)}}
\qquad(x\in\Lambda_r),
\]

where \(\mathcal B_r\) is the complete block from length \(p^r\) to \(p^{r+1}\), \(\Lambda_r\) is an allowed lattice, \(\Pi_r\) is a readout, and \(C_r\) is the compressed matrix.

## 2. Two open targets

**G1 (dyadic).** For classical Apéry, \(p=2,t=-1\),

\[
v_2\bigl(F_{2^r}(-1)\bigr)=3r-1
\qquad(r\ge1).
\]

**G2 (odd CM).** For classical Apéry, \(p\ge5,t=\pm1\), the candidate

\[
C=\begin{pmatrix}T&-D\\1&0\end{pmatrix},
\qquad D=\chi p^2,
\]

is expected to satisfy the comparison modulo \(p^{3r+5}\).

Both are open. The repository provides exact and finite evidence, not proofs.

## 3. Exact objects

The \(n\)-axis state is

\[
X_n=(u_n,u_{n-1},F_n)^T,
\qquad
X_{n+1}=\mathcal R_n(t)X_n,
\]

with

\[
\mathcal R_n(t)=
\begin{pmatrix}
\dfrac{tP(n)}{(n+1)^3}&-\dfrac{t^2n^3}{(n+1)^3}&0\\
1&0&0\\
1&0&1
\end{pmatrix},
\qquad
P(n)=34n^3+51n^2+27n+5.
\]

Adding the \(F_n(-t)\) channel gives a four-state exact transport. Complete blocks

\[
\mathcal K_{a,b}=\mathcal R_{b-1}\cdots\mathcal R_a
=
\begin{pmatrix}
\mathsf B_{a,b}&0\\
\rho_{a,b}&1
\end{pmatrix}
\]

compose by

\[
(\mathsf B_2,\rho_2)\circ(\mathsf B_1,\rho_1)
=
(\mathsf B_2\mathsf B_1,\ \rho_1+\rho_2\mathsf B_1).
\]

The mirror covariance

\[
\mathcal K_{a,b}(-t)S_a
=
S_b\mathcal K_{a,b}(t)
\]

is exact over \(\mathbb Q(t)\).

## 4. The arithmetic side

For a CM point \(t=t_0\), the first block is

\[
H_p(t)=\sum_{n<p}A_n t^n.
\]

The observed trichotomy is:

| \(\chi=(-d/p)\) | type | \(v_p(H_p)\) |
|---|---|---|
| \(1\) | split | \(0\) |
| \(0\) | ramified | \(1\) |
| \(-1\) | inert | \(\ge2\) |

The geometric explanation is the classical Deuring/Hasse picture: \(p\) inert in the CM field means the CM elliptic curve is supersingular at \(p\), so the Hasse invariant vanishes.

The hard part is the comparison theorem: why the first block reads the actual order class polynomial, and how the depth matrix sees it. This is not proved in general.

## 5. What is exact, fixed, finite, open

| Layer | Examples | Status |
|---|---|---|
| Exact | \(\mathcal R_n,\mathcal K_{a,b}\), mirror covariance, Borel four-state transport, \(\eta_p\) | proved/exact |
| Fixed precision | \(C=[[6,16],[1,0]]\bmod128\), \(p=3\) K18, Wilson–Gessel \(K=2\) | fixed range |
| Finite observation | local dyadic blocks, first-block trichotomy, CM order audit | finite, not theorem |
| Open | G1, G2, \(T\), \(t\leftrightarrow j\) comparison | conjecture |

## 6. What this repository contains

- a small allowlisted runner (`run_tool.py`);
- five independent verifiers with frozen inputs;
- exact certificates for mirror, deep probes, Beukers all-\(t\), Wilson–Gessel, and Borel;
- a small CM audit summary;
- documentation of the boundaries.

It is deliberately not the full research archive. The full archive contains large audit bundles, historical reports, and third-party material that are not needed to reproduce the included checks.

## 7. Reading order

1. This overview.
2. [`toolbox.md`](toolbox.md) for the matrix tools.
3. [`evidence.md`](evidence.md) for the numbers and status.
4. `run_tool.py --all` to reproduce.
