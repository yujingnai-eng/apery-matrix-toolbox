# Matrix toolbox

## 1. N matrices

The project uses a recursive \(2\times2\) block representation. For a coefficient ring \(R\),

\[
\mathscr N_0(R)=R,
\qquad
\mathscr N_{h+1}(R)=M_2(\mathscr N_h(R)).
\]

A \(2\times2\) matrix of \(2\times2\) matrices is \(\mathscr N_2(R)\), expanded as \(4\times4\). The nesting depth \(h\) is not the truncation depth \(r\).

All layers use the standard ordered block multiplication

\[
\begin{pmatrix}A&B\\C&D\end{pmatrix}
\begin{pmatrix}E&F\\G&H\end{pmatrix}
=
\begin{pmatrix}AE+BG&AF+BH\\CE+DG&CF+DH\end{pmatrix}.
\]

A scalar \(c\) embeds as \(cI\). Not every block has a scalar representative.

## 2. Mirror covariance

For the four-state \(Y_n=(u_n,u_{n-1},F_n(t),F_n(-t))^T\), the exact mirror is

\[
S_n=\operatorname{diag}\bigl((-1)^n,(-1)^{n-1},J\bigr),
\qquad
J=\begin{pmatrix}0&1\\1&0\end{pmatrix},
\]

with

\[
\mathcal K_{a,b}(-t)S_a
=
S_b\mathcal K_{a,b}(t).
\]

The mirror tool checks this exact covariance and ten finite experiments. It does not prove G1/G2, and it does not identify an intrinsic chirality or a geometric Frobenius.

## 3. Local dyadic blocks

For aligned dyadic blocks,

\[
B_{a,s}=\sum_{n=a2^s}^{(a+1)2^s-1}(-1)^nA_n.
\]

The observed law is

\[
v_2(B_{a,s})=3s-1
\qquad(s\le12),
\]

over 8,190 blocks. For the origin pair

\[
f_r=B_{0,r},
\qquad
b_r=B_{1,r},
\]

one observes

\[
v_2(f_r)=v_2(b_r)=3r-1,
\qquad
v_2(b_r-7f_r)=4r+2
\qquad(r\le12).
\]

The recursion

\[
f_{r+1}=8f_r+(b_r-7f_r)
\]

explains why a proof of the local block law would imply G1. The missing step is the all-\(r\) induction, not more data.

## 4. Wilson–Gessel unit-factorial tower

Define

\[
u_p(n!)=\frac{n!}{p^{v_p(n!)}},
\qquad
d_n=\frac{A_n}{u_p(n!)},
\]

and

\[
D_N(t)=\sum_{n<N}d_nt^n,
\qquad
J_N(t)=\sum_{n<N}n d_nt^n.
\]

The unit-factorial family is not the original \(A_n\) tower and not the Borel sum.

Wilson’s unit block and Gessel’s congruence give, for \(n=pk+\ell\),

\[
d_{pk+\ell}\equiv(-1)^k(e_\ell+pkb_\ell)d_k\pmod{p^2},
\]

with

\[
e_\ell=\frac{A_\ell}{\ell!},
\qquad
b_\ell=\frac{\dot A_\ell-H_\ell A_\ell}{\ell!}+W_pe_\ell,
\]

\[
W_p=\frac{(p-1)!+1}{p},
\qquad
\dot A_\ell
=
2\sum_{j=0}^{\ell}\binom{\ell}{j}^2\binom{\ell+j}{j}^2
\bigl(H_{\ell+j}-H_{\ell-j}\bigr).
\]

The \(p k\) term reads \(k d_k\), so two states are insufficient. The four-state matrix is

\[
M_\ell=
\begin{pmatrix}
0_2&\mathsf B_\ell\\
(-1)^\ell\mathsf B_\ell&0_2
\end{pmatrix},
\qquad
\mathsf B_\ell=
\begin{pmatrix}
e_\ell&pb_\ell\\
\ell e_\ell&p(e_\ell+\ell b_\ell)
\end{pmatrix}.
\]

This is a fixed \(K=2\) construction. It is a template, not G2.

## 5. Borel first block

The exponential generating function is

\[
E_N(t)=\sum_{n<N}\frac{A_n}{n!}t^n.
\]

It has an exact four-state transport and satisfies

\[
L_E\mathcal E=0,
\qquad
L_E=\theta^4-zP(\theta)+z^2(\theta+1)^2.
\]

The truncated boundary is

\[
L_E E_N(z)
=
-\frac{N^3A_N}{(N-1)!}z^N
+\frac{N^2A_{N-1}}{(N-1)!}z^{N+1}.
\]

This shows that dividing by \(n!\) does **not** simply turn the differential equation into a difference equation. The Borel object is still D-finite, but it is a different tower from the original \(F_N\).

## 6. Running the tools

```sh
python3 run_tool.py --list
python3 run_tool.py --run mirror
python3 run_tool.py --run deep_probes
python3 run_tool.py --run beukers_all_t
python3 run_tool.py --run wilson_gessel
python3 run_tool.py --run borel
python3 run_tool.py --all
```

Each tool writes a JSON certificate in its own directory. The runner checks exit code, status, and key counts. It does not check the mathematics beyond the verifier’s declared scope.
