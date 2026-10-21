#!/usr/bin/env python3
"""Machine-readable registry for Beukers Appendix C (v3, sections 9.1--9.15).

配套修订版（2026-10-17）。Beukers arXiv:2403.03301v3 Appendix C 有20条
参数/坐标/CM域配套疑点；这不等于20个参数都不是CM值。
  * 9.7、9.8 两节：系数为 (-1)^n 扭曲序列，但所列 t 值属于未扭曲变量；
    对扭曲序列应取 -t。佐证：9.7 原表的 1/4、1/16 是未扭曲 Domb 方程的奇点，
    而 9.7 自身算子 (1+20t+64t^2) 的奇点是 -1/4、-1/16。
  * 6 个孤立点的 t 值有误（见 ERRATA）。
有限一致性检验比较H_p(t)是否为零与所配域的惰性字符。
将它识别为Deuring判据，仍需逐族的Hasse识别和良好约化前提；不是非CM判定器。
至少4条原参数已识别另一CM配对；其余16条在报告的有限域搜索内尚未找到匹配，非CM性未证。
修正值先由算术反求得到，后用原文自己的 mirror map（Appendix A 的 h_N）在原文 α 处高精度求值直接确认：
6 个孤立点 t(α) 即修正值；9.7、9.8 的 14 点满足 t(α+1/2) = -t_表（未扭曲→扭曲即 τ→τ+1/2）。
alpha始终保留原文字符串；alpha_published为同义溯源字段。
alpha_corrected与tau_evaluated明确给出当前t的求值坐标，14条记录的tau_shift为1/2，其余为0。
高精度求值只是数值核验，不冒充特殊值恒等式证明；修正值尚未经原作者确认。
"""

from __future__ import annotations

import math
from fractions import Fraction
from typing import Iterable


SOURCE_URL = "https://arxiv.org/html/2403.03301v3"


def point(alpha: str, t: str, d: int, fixed: str | None = None, note: str = "", *, tau_shift: str = "0") -> dict:
    """Keep legacy alpha and published t; expose a separate evaluation coordinate."""
    if any(d % (q * q) == 0 for q in range(2, math.isqrt(d) + 1)):
        raise ValueError(f"CM radicand must be square-free: {d}")
    delta = -d if d % 4 == 3 else -4 * d
    published = Fraction(t)
    value = Fraction(fixed) if fixed is not None else published
    shift = Fraction(tau_shift)
    if shift not in (Fraction(0), Fraction(1, 2)):
        raise ValueError("only the explicitly checked zero/half coordinate shifts are supported")
    evaluated = f"({alpha})+1/2" if shift else alpha
    return {"alpha": alpha, "t": value, "d": d, "delta": delta,
            "alpha_published": alpha, "alpha_corrected": evaluated,
            "tau_evaluated": evaluated, "tau_shift": str(shift),
            "t_published": published, "t_corrected": value,
            "correction_status": "numerically_checked_pending_author_confirmation" if fixed is not None else "published_unmodified",
            "published_cm_match_status": "no_match_found_in_reported_finite_search" if fixed is not None else "not_assessed",
            "published_cm_d_alternative": None,
            "erratum": note if fixed is not None else ""}


SIGN = "所列14条记录在扭曲mirror map下数值满足t(alpha_published+1/2)=-t_published；不是任意族的通用换号规则，非CM性未由此推出"
TYPO = "原文t与所列坐标不匹配；原文mirror map在原坐标处的数值支持修正值；原参数在有限搜索中尚未找到CM匹配，非CM性未证"


FAMILIES = [
    {
        "section": "9.1", "id": "hyp_1_2", "name": "3F2(1/2,1/2,1/2)",
        "oeis": "A002897", "level": 4, "p1": [-8, -48, -96, -64], "p2": [0, 0, 0, 0],
        "formula": "C(2n,n)^3",
        "points": [
            point("i/2", "1/64", 1), point("(1+sqrt(-2))/2", "-1/64", 2),
            point("sqrt(-3)/2", "1/256", 3), point("(1+sqrt(-3))/4", "1/16", 3),
            point("(1+i)/2", "-1/2", 1, fixed="-1/8", note=TYPO), point("1/2+i", "-1/512", 1),
            point("(3+sqrt(-7))/8", "1", 7), point("sqrt(-7)/2", "1/4096", 7),
        ],
    },
    {
        "section": "9.2", "id": "hyp_1_4", "name": "3F2(1/4,1/2,3/4)",
        "oeis": "A008977", "level": 2, "p1": [-24, -176, -384, -256], "p2": [0, 0, 0, 0],
        "formula": "(4n)!/(n!)^4",
        "points": [
            point("sqrt(-1/2)", "1/256", 2), point("(1+sqrt(-5))/2", "-1/1024", 5),
            point("sqrt(-3/2)", "1/2304", 6), point("(1+sqrt(-3))/2", "-1/144", 3),
            point("(1+3i)/2", "-1/17088", 1, fixed="-1/12288", note=TYPO), point("i", "1/648", 1),
            point("sqrt(-5/2)", "1/20736", 10), point("(1+sqrt(-7))/4", "1/81", 7),
            point("(1+sqrt(-13))/2", "-1/82944", 13), point("(1+sqrt(-7))/2", "-1/3969", 7),
            point("3sqrt(-2)/2", "1/614656", 2), point("sqrt(-11/2)", "1/9801", 22, fixed="1/2509056", note=TYPO),
            point("(1+5i)/2", "-1/6635520", 1), point("(1+sqrt(-37))/2", "-1/199148544", 37),
            point("sqrt(-29/2)", "1/24591257856", 58),
        ],
    },
    {
        "section": "9.3", "id": "hyp_1_3", "name": "3F2(1/3,1/2,2/3)",
        "oeis": "A151374", "level": 3, "p1": [-12, -78, -162, -108], "p2": [0, 0, 0, 0],
        "formula": "(3n)!(2n)!/(n!)^5",
        "points": [
            point("sqrt(-3)/3", "1/108", 3), point("(3+sqrt(-15))/6", "-1/27", 15),
            point("sqrt(-6)/3", "1/216", 6), point("(1+sqrt(-2))/3", "1/8", 2),
            point("(1+sqrt(-17/3))/2", "-1/1728", 51), point("(1+sqrt(-3))/2", "-1/192", 3),
            point("(1+sqrt(-11))/6", "1/1024", 11, fixed="1/64", note=TYPO), point("(1+5sqrt(-1/3))/2", "-1/8640", 3),
            point("sqrt(-5/3)", "1/3375", 15), point("(1+sqrt(-41/3))/2", "-1/110592", 123),
            point("(1+7sqrt(-1/3))/2", "-1/326592", 3),
        ],
    },
    {
        "section": "9.4", "id": "hyp_1_6", "name": "3F2(1/6,1/2,5/6)",
        "oeis": "A001421", "level": 1, "p1": [-120, -1104, -2592, -1728], "p2": [0, 0, 0, 0],
        "formula": "(6n)!/((3n)!(n!)^3)", "points": [],
    },
    {
        "section": "9.5", "id": "apery3", "name": "Apery zeta(3) / Almkvist gamma",
        "oeis": "A005259", "level": 6, "p1": [-5, -27, -51, -34], "p2": [1, 3, 3, 1],
        "formula": "sum C(n,k)^2 C(n+k,k)^2",
        "points": [point("(2+sqrt(-2))/6", "1", 2), point("(3+sqrt(-3))/6", "-1", 3)],
    },
    {
        "section": "9.6", "id": "almkvist_eta_tw", "name": "Almkvist eta (twisted)",
        "oeis": "A005810", "level": 5, "p1": [5, 21, 33, 22], "p2": [125, 375, 375, 125],
        "formula": "sum (-1)^(n-k) C(n,k)^3 C(4n-5k,3n)", "points": [],
    },
    {
        "section": "9.7", "id": "domb_tw", "name": "Domb / Almkvist alpha (twisted)",
        "oeis": "A002895", "level": 6, "p1": [4, 18, 30, 20], "p2": [64, 192, 192, 64],
        "formula": "(-1)^n sum C(n,k)^2 C(2k,k) C(2n-2k,n-k)",
        "points": [
            point("(9+sqrt(-15))/24", "1", 15, fixed="-1", note=SIGN, tau_shift="1/2"), point("(3+sqrt(-3))/6", "-1/2", 3, fixed="1/2", note=SIGN, tau_shift="1/2"),
            point("(3+sqrt(-3))/12", "1/4", 3, fixed="-1/4", note=SIGN, tau_shift="1/2"), point("(1+sqrt(-2))/6", "1/8", 2, fixed="-1/8", note=SIGN, tau_shift="1/2"),
            point("(1+sqrt(-2/3))/2", "-1/8", 6, fixed="1/8", note=SIGN, tau_shift="1/2"), point("sqrt(-3)/6", "1/16", 3, fixed="-1/16", note=SIGN, tau_shift="1/2"),
            point("(1+sqrt(-4/3))/2", "-1/32", 3, fixed="1/32", note=SIGN, tau_shift="1/2"), point("sqrt(-15)/30", "1/64", 15, fixed="-1/64", note=SIGN, tau_shift="1/2"),
        ],
    },
    {
        "section": "9.8", "id": "almkvist_zudilin_tw", "name": "Almkvist-Zudilin delta (twisted)",
        "oeis": "A125143", "level": 6, "p1": [3, 13, 21, 14], "p2": [81, 243, 243, 81],
        "formula": "sum_{k<=n/3} (-3)^(n-3k)(n+k)!/(k!^4(n-3k)!)",
        "points": [
            point("(3+sqrt(-2))/6", "-1", 2, fixed="1", note=SIGN, tau_shift="1/2"), point("(2+i)/6", "1/3", 1, fixed="-1/3", note=SIGN, tau_shift="1/2"),
            point("(3+sqrt(-3))/12", "1/9", 3, fixed="-1/9", note=SIGN, tau_shift="1/2"), point("(1+sqrt(-2/3))/2", "-1/9", 6, fixed="1/9", note=SIGN, tau_shift="1/2"),
            point("i/2", "1/27", 1, fixed="-1/27", note=SIGN, tau_shift="1/2"), point("(1+sqrt(-2))/2", "-1/81", 2, fixed="1/81", note=SIGN, tau_shift="1/2"),
        ],
    },
    {
        "section": "9.9", "id": "almkvist_epsilon", "name": "Almkvist epsilon",
        "oeis": "A290575", "level": 8, "p1": [-4, -20, -36, -24], "p2": [16, 48, 48, 16],
        "formula": "sum_{k>=ceil(n/2)} C(n,k)^2 C(2k,n)^2",
        "points": [
            point("(5+sqrt(-7))/16", "1", 7), point("(2+sqrt(-2))/4", "-1/4", 2),
            point("(1+i)/4", "1/4", 1), point("(1+sqrt(-7))/8", "1/16", 7),
        ],
    },
    {
        "section": "9.10", "id": "almkvist_zeta", "name": "Almkvist zeta",
        "oeis": "A290576", "level": 9, "p1": [-3, -15, -27, -18], "p2": [-27, -81, -81, -27],
        "formula": "sum_{k,l} C(n,k)^2 C(n,l) C(k,l) C(k+l,n)",
        "points": [
            point("(5+sqrt(-11))/18", "1", 11), point("(4+sqrt(-2))/9", "-1", 2),
            point("(1+sqrt(-1/3))/2", "-1/3", 3), point("(1+sqrt(-3))/6", "1/9", 3),
            point("sqrt(-2)/3", "1/27", 2), point("(3+sqrt(-11))/6", "-1/27", 11),
        ],
    },
    {
        "section": "9.11", "id": "cooper_s7", "name": "Cooper s7",
        "oeis": "A183204", "level": 7, "p1": [-4, -21, -39, -26], "p2": [-24, -78, -81, -27],
        "formula": "sum_{k>=ceil(n/2)} C(n,k)^2 C(n+k,k) C(2k,n)",
        "points": [
            point("(1+sqrt(-1/7))/2", "-1", 7), point("(2+sqrt(-3))/7", "1/2", 3),
            point("(3+sqrt(-19))/14", "1/8", 19), point("(1+sqrt(-5/7))/2", "-1/8", 35),
            point("(1+3sqrt(-3))/14", "1/24", 3), point("sqrt(-1/7)", "1/27", 7),
            point("(1+sqrt(-13/7))/2", "-1/64", 91), point("sqrt(-4/7)", "1/125", 7),
            point("(1+sqrt(-61/7))/2", "-1/10648", 427),
        ],
    },
    {
        "section": "9.12", "id": "cooper_s10", "name": "Yang-Zudilin / Cooper s10 / Franel-4",
        "oeis": "A005260", "level": 10, "p1": [-2, -10, -18, -12], "p2": [-60, -188, -192, -64],
        "formula": "sum C(n,k)^4",
        "points": [
            point("(5+sqrt(-15))/20", "1", 15), point("(2+i)/5", "-1/2", 1),
            point("(2+sqrt(-6))/10", "1/4", 6), point("(1+sqrt(-1/5))/2", "-1/4", 5),
            point("(1+sqrt(-3/5))/2", "-1/9", 15), point("(1+3i)/10", "1/12", 1),
            point("sqrt(-1/10)", "1/16", 10), point("(1+i)/2", "-1/20", 1),
            point("sqrt(-3/10)", "1/36", 30), point("(1+sqrt(-9/5))/2", "-1/64", 5),
            point("sqrt(-7/10)", "1/196", 70), point("(1+sqrt(-17/5))/2", "-1/324", 85),
            point("sqrt(-13/10)", "1/1296", 130), point("sqrt(-19/10)", "1/5776", 190),
        ],
    },
    {
        "section": "9.13", "id": "cooper_s18", "name": "Cooper s18",
        "oeis": "", "level": 18, "p1": [-6, -26, -42, -28], "p2": [180, 564, 576, 192],
        "formula": "Cooper s18 binomial sum",
        "points": [
            point("(3+sqrt(-5))/6", "-1/4", 5), point("(3+sqrt(-7))/12", "1/9", 7),
            point("(3+sqrt(-7))/6", "-1/9", 7), point("(1+i)/6", "1/12", 1),
            point("sqrt(-2)/6", "1/16", 2), point("i/3", "1/18", 1),
            point("sqrt(-10)/6", "1/36", 10), point("(3+sqrt(-13))/6", "-1/36", 13),
            point("sqrt(-22)/6", "1/144", 22), point("(3+5i)/6", "-1/180", 1),
            point("(3+sqrt(-37))/6", "-1/576", 37),
        ],
    },
    {
        "section": "9.14", "id": "apery2_central", "name": "Apery zeta(2) Hadamard central",
        "oeis": "", "level": 5, "p1": [-6, -34, -66, -44], "p2": [-12, -44, -48, -16],
        "formula": "C(2n,n) sum C(n,k)^2 C(n+k,k)",
        "points": [
            point("(3+sqrt(-11))/10", "1/4", 11), point("(1+sqrt(-3/5))/2", "-1/3", 15),
            point("(1+sqrt(-7/5))/2", "-1/28", 35), point("(1+sqrt(-19))/10", "1/36", 19),
            point("(-1+2i)/5", "1/18", 1), point("sqrt(-2/5)", "1/72", 10),
            point("sqrt(-3/5)", "1/588", 15, fixed="1/147", note=TYPO), point("(1+sqrt(-23/5))/2", "-1/828", 115),
            point("(1+sqrt(-47/5))/2", "-1/15228", 235),
        ],
    },
    {
        "section": "9.15", "id": "franel_central", "name": "Franel-3 Hadamard central",
        "oeis": "", "level": 6, "p1": [-4, -22, -42, -28], "p2": [-96, -352, -384, -128],
        "formula": "C(2n,n) sum C(n,k)^3",
        "points": [
            point("(1+sqrt(-1/3))/2", "-1/4", 3), point("(1+sqrt(-5))/6", "1/16", 5),
            point("(1+i)/2", "-1/16", 1), point("(3+sqrt(-15))/12", "1/5", 15),
            point("sqrt(-1/6)", "1/32", 6), point("sqrt(-1/2)", "1/96", 2),
            point("sqrt(-1/3)", "1/50", 3), point("(1+sqrt(-7/3))/2", "-1/116", 21, fixed="-1/112", note=TYPO),
            point("(1+sqrt(-5/3))/2", "-1/49", 15), point("sqrt(-5/6)", "1/320", 30),
            point("(1+sqrt(-11/3))/2", "-1/400", 33), point("sqrt(-7/6)", "1/896", 42),
            point("(1+sqrt(-19/3))/2", "-1/2704", 57), point("sqrt(-13/6)", "1/10400", 78),
            point("sqrt(-17/6)", "1/39200", 102),
        ],
    },
]


_ALTERNATIVE_PUBLISHED_CM_FIELDS = {("9.7", 4): 6, ("9.7", 5): 2, ("9.8", 3): 6, ("9.8", 4): 3}
for _family in FAMILIES:
    for _index, _point in enumerate(_family["points"], 1):
        _alternative = _ALTERNATIVE_PUBLISHED_CM_FIELDS.get((_family["section"], _index))
        if _alternative is not None:
            _point["published_cm_match_status"] = "identified_alternative_cm_pairing_numerical"
            _point["published_cm_d_alternative"] = _alternative
            _point["erratum"] += f"；原t本身另有CM配对d={_alternative}，不是非CM值"


COORDINATE_METADATA_FIELDS = (
    "alpha_published", "alpha_corrected", "tau_evaluated", "tau_shift",
    "t_published", "t_corrected", "correction_status",
    "published_cm_match_status", "published_cm_d_alternative",
)


def coordinate_metadata(pt: dict) -> dict:
    """Additive export contract: never rebind the legacy alpha field."""
    return {key: str(pt[key]) if isinstance(pt[key], Fraction) else pt[key]
            for key in COORDINATE_METADATA_FIELDS}


def poly(coefficients: list[int], x: int) -> int:
    return sum(coefficient * x**power for power, coefficient in enumerate(coefficients))


def generate_coefficients(family: dict, count: int) -> list[int]:
    if count <= 0:
        return []
    values = [1]
    for n in range(1, count):
        numerator = poly(family["p1"], n - 1) * values[n - 1]
        if n >= 2:
            numerator += poly(family["p2"], n - 2) * values[n - 2]
        numerator = -numerator
        denominator = n**3
        if numerator % denominator:
            raise ArithmeticError(f"{family['section']} recurrence nonintegral at n={n}")
        values.append(numerator // denominator)
    return values


def binom_integer(top: int, bottom: int) -> int:
    if bottom < 0:
        return 0
    if top >= 0:
        return math.comb(top, bottom) if bottom <= top else 0
    return (-1) ** bottom * math.comb(bottom - top - 1, bottom)


def direct_coefficient(section: str, n: int) -> int:
    c = math.comb
    f = math.factorial
    if section == "9.1": return c(2*n, n) ** 3
    if section == "9.2": return f(4*n) // f(n) ** 4
    if section == "9.3": return f(3*n) * f(2*n) // f(n) ** 5
    if section == "9.4": return f(6*n) // (f(3*n) * f(n) ** 3)
    if section == "9.5": return sum(c(n,k)**2 * c(n+k,k)**2 for k in range(n+1))
    if section == "9.6":
        return sum((-1)**(n-k) * c(n,k)**3 * binom_integer(4*n-5*k, 3*n) for k in range(n+1))
    if section == "9.7":
        return (-1)**n * sum(c(n,k)**2 * c(2*k,k) * c(2*n-2*k,n-k) for k in range(n+1))
    if section == "9.8":
        return sum((-3)**(n-3*k) * f(n+k) // (f(k)**4 * f(n-3*k)) for k in range(n//3+1))
    if section == "9.9":
        return sum(c(n,k)**2 * c(2*k,n)**2 for k in range((n+1)//2,n+1))
    if section == "9.10":
        return sum(c(n,k)**2*c(n,l)*c(k,l)*c(k+l,n) for k in range(n+1) for l in range(n+1))
    if section == "9.11":
        return sum(c(n,k)**2*c(n+k,k)*c(2*k,n) for k in range((n+1)//2,n+1))
    if section == "9.12": return sum(c(n,k)**4 for k in range(n+1))
    if section == "9.13":
        if n == 0: return 1
        return sum((-1)**k*c(n,k)*c(2*k,k)*c(2*n-2*k,n-k)*(binom_integer(2*n-3*k-1,n)+binom_integer(2*n-3*k,n)) for k in range(n//3+1))
    if section == "9.14": return c(2*n,n)*sum(c(n,k)**2*c(n+k,k) for k in range(n+1))
    if section == "9.15": return c(2*n,n)*sum(c(n,k)**3 for k in range(n+1))
    raise KeyError(section)


def self_check(limit: int = 9) -> list[dict]:
    results = []
    for family in FAMILIES:
        recurrence = generate_coefficients(family, limit)
        direct = [direct_coefficient(family["section"], n) for n in range(limit)]
        results.append({"section": family["section"], "pass": recurrence == direct,
                        "recurrence": recurrence, "direct": direct})
    return results


ERRATA = [(f["section"], i, p["t_published"], p["t"], p["erratum"])
          for f in FAMILIES for i, p in enumerate(f["points"], 1) if p["erratum"]]


def _kronecker(a: int, p: int) -> int:
    a %= p
    return 0 if a == 0 else (1 if pow(a, (p - 1) // 2, p) == 1 else -1)


def deuring_check(primes: Iterable[int] | None = None, use_published: bool = False) -> list[dict]:
    """Finite consistency screen, not a non-CM test or a growing-precision theorem.

    A mismatch calls for checking coordinates, CM labels, good reduction and the
    family-specific Hasse identification. A match does not prove those hypotheses.
    """
    if primes is None:
        primes = [q for q in range(5, 400) if all(q % r for r in range(2, math.isqrt(q) + 1))]
    primes = tuple(primes)
    out = []
    for family in FAMILIES:
        cache: dict[int, list[int]] = {}
        for index, pt in enumerate(family["points"], 1):
            t = pt["t_published"] if use_published else pt["t"]
            ok = n = 0
            for p in primes:
                if (6 * family["level"]) % p == 0 or t.denominator % p == 0 or t.numerator % p == 0:
                    continue
                chi = _kronecker(pt["delta"], p)
                if chi == 0:
                    continue
                if p not in cache:
                    a = [1]
                    for m in range(1, p):
                        num = -(poly(family["p1"], m - 1) * a[m - 1]
                                + (poly(family["p2"], m - 2) * a[m - 2] if m >= 2 else 0))
                        a.append(num * pow(m ** 3, -1, p) % p)
                    cache[p] = a
                x = t.numerator * pow(t.denominator, -1, p) % p
                h = sum(c * pow(x, k, p) for k, c in enumerate(cache[p])) % p
                n += 1
                ok += (h == 0) == (chi == -1)
            out.append({"point_id": f'{family["section"]}-t{index:02d}', "ok": ok, "n": n})
    return out
