# Apéry Matrix Toolbox

Exact matrix-transport and independent-verification tools for Apéry supercongruences.

**Status: research toolbox, not a proof of G1/G2.**  
The two target statements remain open. This repository contains reproducible finite certificates, exact structural checks, and tools that make the boundary between *proved*, *fixed precision*, and *finite observation* explicit.

[Overview](docs/overview.md) · [Toolbox](docs/toolbox.md) · [Evidence](docs/evidence.md) · [Publishing](PUBLISHING.md)

---

## 中文简介

这是一个 Apéry 超同余的**矩阵工具箱**：把原数列的精确推进、完整块、数字矩阵、镜像、局部二进块、Wilson–Gessel 单位阶乘族和 Borel 首块做成可独立复算的脚本与证书。

它**不证明** G1（\(v_2(F_{2^r}(-1))=3r-1\)）或 G2（\(p^{3r+5}\) 的二维 CM 候选）。  
它做的是：把已经能精确验证的部分固定下来，把有限观察与定理分开，把反例保留下来。

---

## Quick start

```sh
git clone https://github.com/yujingnai-eng/apery-matrix-toolbox.git
cd apery-matrix-toolbox
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

python3 run_tool.py --list
python3 run_tool.py --all
```

Expected: five tools pass, about 30 seconds on a laptop.

A single tool:

```sh
python3 run_tool.py --run mirror
python3 run_tool.py --run beukers_all_t
```

All tools are allowlisted in [`manifest.json`](manifest.json). `run_tool.py` never uses a shell and never accepts arbitrary commands.

---

## What is inside

| Tool | What it checks | Scope |
|---|---|---|
| `mirror` | Ten mirror experiments; exact \(K(-t)S_a=S_bK(t)\) covariance plus finite controls | finite + symbolic |
| `deep_probes` | Local dyadic blocks, kernel splitting, initial phase diagram | finite ranges |
| `beukers_all_t` | 118 Beukers Appendix C points, \(p\le199\), 5,130 eligible cases | finite |
| `wilson_gessel` | Unit-factorial Wilson–Gessel fixed \(p^2\) four-state matrix | fixed \(K=2\) |
| `borel` | Borel first block: exact four-state differential transport | exact + finite checks |

The repository also contains a small CM audit summary:

- `data/cm_audit/summary.json`
- `data/cm_audit/points_table.md`
- `data/cm_audit/class_polynomials.json`

The full 38 MB audit bundle is not included; the summary records its main numbers and boundaries.

---

## Key numbers

Beukers Appendix C:

- 15 families; 9.4 and 9.6 have no points; **13 non-empty sections**.
- **118** CM points; **44** primes \(5\le p\le199\).
- **5,192** point–prime grid; **5,130** legacy eligible; **5,092** strict non-ramified.
- First-block trichotomy on 5,130 cases: **2,433 split**, **38 ramified**, **2,659 inert**.
- Inert valuations: 2,592 at \(v=2\); 60 at 3; 6 at 4; 1 at 5.
- **561** cases have a geometric supersingular common root but no \(\mathbb F_p\) root.
- **38** legacy ramified cases also vanish; the naive “iff inert” needs non-ramified/good-reduction hypotheses.

CM audit:

- Field class numbers 66/40/12; **actual order class numbers 43/54/16/5**.
- 46 non-maximal orders; 29 class-number changes.
- 23 of the old 66 “class number one” points have actual \(j(\tau)\) of degree \(>1\).
- 54 class polynomials; 2,376 gcd certificates.

Local dyadic blocks:

- \(v_2(B_{a,s})=3s-1\) for \(s\le12\) (8,190 blocks, 0 exceptions).
- \(v_2(b_r-7f_r)=4r+2\) for \(r\le12\).
- G1 itself is **not proved**.

---

## Repository layout

```text
.
├─ run_tool.py              # allowlisted runner
├─ manifest.json            # tool list, commands, expected outputs
├─ tools/
│  ├─ mirror/               # verifier + producer + frozen input
│  ├─ deep_probes/
│  ├─ beukers_all_t/
│  ├─ wilson_gessel/
│  └─ borel/
├─ data/cm_audit/           # small audit summary (not the 38 MB bundle)
├─ docs/
│  ├─ overview.md
│  ├─ toolbox.md
│  └─ evidence.md
├─ scripts/check_public_tree.py
├─ CITATION.cff
├─ LICENSE                  # MIT for code
└─ LICENSE-docs             # CC BY 4.0 for docs
```

---

## What this is not

- Not a proof of G1 or G2.
- Not a proof that the first block equals the Hasse invariant for every CM point.
- Not a classification of the 561 geometric-only cases, the 38 ramified cases, or the 67 higher inert cases.
- Not a replacement for the full research archive; it is a curated, reproducible subset.
- Not a claim of originality for Deuring/Hasse, class polynomials, or Beukers’ supercongruences.

The intended use is: **reproduce the finite certificates, inspect the exact structural identities, and build on them without confusing finite consistency with a theorem.**

---

## Reproducibility

- Exact integer/rational arithmetic; no floating-point decisions.
- Each tool writes a JSON certificate; `run_tool.py` checks status and key counts.
- Frozen inputs are included next to their verifiers.
- `scripts/check_public_tree.py` scans for absolute paths, secrets, PDFs, large files, and private archive names before publishing.

```sh
python3 scripts/check_public_tree.py
```

---

## Data and third-party material

The repository contains **project-generated code and data only**. It does not redistribute third-party PDFs or referee reports. References to Beukers, Deuring, Hasse, and Golyshev–Zagier are citations, not bundled copies.

If you reuse the data, cite this repository and the relevant papers.

---

## Citation

See [`CITATION.cff`](CITATION.cff). Replace the placeholder author metadata with your preferred name/ORCID before a public release.

---

## License

- **Code**: MIT, see [`LICENSE`](LICENSE).
- **Documentation and data**: CC BY 4.0, see [`LICENSE-docs`](LICENSE-docs).

---

## Contributing

Issues and pull requests are welcome, especially:

- independent re-runs of a verifier;
- counterexamples to a claimed finite pattern;
- corrections to the CM/order data;
- clearer separation between exact, fixed-precision, and finite claims.

Please do not add a “proof” of G1/G2 without a complete, checkable argument. Finite agreement is not a proof.

---

## Related literature

- F. Beukers, *Supercongruences using modular forms*, arXiv:2403.03301.
- F. Beukers, *p-linear schemes*, arXiv:2211.15240.
- G. Golyshev and D. Zagier, *Proof of the Gamma Conjecture for Fano 3-folds of Picard rank one*.
- Classical Deuring/Hasse theory of CM elliptic curves and supersingular reduction.
