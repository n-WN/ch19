from __future__ import annotations

from fractions import Fraction

from .bivar import Bivar, degree_x, degree_y, pow_bivar, shift_x, shift_y
from .elimination import resultant_in_x_by_interpolation
from .lll import lll_reduction

# 二元 Coppersmith（教学版，简化 Howgrave-Graham 思路）
# 目标：给定 F(x,y) ∈ Z[x,y]，模 N，若存在小根 |x0|<X, |y0|<Y 使 F(x0,y0) ≡ 0 (mod N)，
# 尝试恢复 (x0,y0)
# 仅面向教学和小规模参数；构造法与参数并非最优。


def construct_bivar_lattice(
    F: Bivar, N: int, X: int, Y: int, m: int, tx: int, ty: int
) -> tuple[list[list[int]], list[tuple[int, int]]]:
    # 基多项式：
    # 对 i=0..m-1：N^{m-i} * F(x,y)^i * x^ax * y^ay，0<=ax<=dx-1, 0<=ay<=dy-1
    # 以及 F(x,y)^m * x^ax * y^ay，0<=ax<tx, 0<=ay<ty
    dx = degree_x(F) + 1
    dy = degree_y(F) + 1

    polys: list[Bivar] = []
    for i in range(m):
        Fi = pow_bivar(F, i)
        # 标量放大
        FiN = {(ix, iy): coeff * pow(N, m - i) for (ix, iy), coeff in Fi.items()}
        polys.extend(
            shift_y(shift_x(FiN, ax), ay)
            for ax in range(max(dx - 1, 1))
            for ay in range(max(dy - 1, 1))
        )
    Fm = pow_bivar(F, m)
    polys.extend(shift_y(shift_x(Fm, ax), ay) for ax in range(tx) for ay in range(ty))

    # 列缩放：列 (ix,iy) 乘以 X^ix * Y^iy
    for P in polys:
        for (ix, iy), v in list(P.items()):
            P[(ix, iy)] = v * pow(X, ix) * pow(Y, iy)

    # 展平为整数矩阵
    # 需统一列顺序：按 (ix,iy) 字典序
    all_monos = set()
    for P in polys:
        all_monos.update(P.keys())
    cols = sorted(all_monos)
    col_index: dict[tuple[int, int], int] = {mon: i for i, mon in enumerate(cols)}

    B: list[list[int]] = []
    for P in polys:
        row = [0] * len(cols)
        for mon, v in P.items():
            row[col_index[mon]] = v
        B.append(row)

    return B, cols


def eval_unscaled_row_at(
    row: list[int], cols: list[tuple[int, int]], x0: int, y0: int, X: int, Y: int
) -> Fraction:
    # 行 row 是缩放后多项式的系数，列 (ix,iy) 经过了 X^ix Y^iy 缩放
    # 反缩放：sum row[idx]/(X^ix Y^iy) * x0^ix * y0^iy
    total = Fraction(0)
    for idx, (ix, iy) in enumerate(cols):
        ck = row[idx]
        if ck:
            total += Fraction(ck, pow(X, ix) * pow(Y, iy)) * pow(x0, ix) * pow(y0, iy)
    return total


def try_find_small_roots_bivar(
    F: Bivar, N: int, X: int, Y: int, m: int = 2, tx: int = 2, ty: int = 2
) -> list[tuple[int, int]]:
    B, cols = construct_bivar_lattice(F, N, X, Y, m, tx, ty)
    Bref = lll_reduction(B)

    # 使用最短的两条向量构造两个多项式 G1,G2（反缩放），再对 y 做结果式消元得到单变量 R(x)
    if len(Bref) < 2:
        return []
    g_rows = Bref[:2]

    # 构造分数系数的二元多项式（反缩放后）
    def row_to_bivarfrac(row: list[int]) -> dict[tuple[int, int], Fraction]:
        out: dict[tuple[int, int], Fraction] = {}
        for idx, (ix, iy) in enumerate(cols):
            ck = row[idx]
            if ck:
                out[(ix, iy)] = out.get((ix, iy), Fraction(0)) + Fraction(
                    ck, pow(X, ix) * pow(Y, iy)
                )
        return out

    G1 = row_to_bivarfrac(g_rows[0])
    G2 = row_to_bivarfrac(g_rows[1])

    # 计算关于 y 的结果式 R(x)
    R = resultant_in_x_by_interpolation(G1, G2, X, Y)

    # 在 |x|<X 范围寻找整数根（R(x)=0），并回代穷举少量 y 候选验证 F(x,y)≡0(mod N)
    candidates = set()
    for x0 in range(-X + 1, X):
        # R 可能仅差常数因子；以 R(x0)=0 判断
        valR = 0
        x0_pow = 1
        for ci in R:
            valR += ci * x0_pow
            x0_pow *= x0
        if valR != 0:
            continue
        # 小范围枚举 y
        for y0 in range(-Y + 1, Y):
            total = 0
            for (ix, iy), v in F.items():
                total += v * pow(x0, ix) * pow(y0, iy)
            if total % N == 0:
                candidates.add((x0, y0))
    return sorted(candidates)
