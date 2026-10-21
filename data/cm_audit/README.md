# CM audit summary

This directory contains a small summary of the 118-point order / class-polynomial audit:

- `summary.json`: aggregate counts and distributions;
- `points_table.md`: the 118 corrected points, discriminants, class numbers, and actual \(j\);
- `class_polynomials.json`: the class polynomials used in the audit.

Key numbers:

- field class numbers 66/40/12;
- actual order class numbers **43/54/16/5**;
- 46 non-maximal orders; 29 class-number changes;
- 23 of the old 66 “class number one” points have actual \(j(\tau)\) degree \(>1\);
- 54 class polynomials; 2,376 gcd certificates.

The full 38 MB audit bundle (including all independent certificates and first-block cases)
is not included in this repository. The summary is enough to check the main numbers and
to see why the old “field representative \(j\)” test was superseded.
