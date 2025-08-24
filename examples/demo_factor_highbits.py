#!/usr/bin/env python3
from __future__ import annotations

import random
from math import isqrt

# 因式分解：已知 p 的高位（近似值）
# 设 p≈p0，|p−p0|<X，令 q0=⌊N/p0⌋，构造 F(x,y)=(p0+x)(q0+y)−N。
# 若 |x|<X、|y| 也很小，则可用二元 Coppersmith + 结果式消元恢复 (x,y)。
from coppersmith.bivar import Bivar
from coppersmith.bivariate import try_find_small_roots_bivar


def gen_prime(bits: int) -> int:
    def is_prime(n: int) -> bool:
        if n < 2:
            return False
        if n % 2 == 0:
            return n == 2
        r = isqrt(n)
        f = 3
        while f <= r:
            if n % f == 0:
                return False
            f += 2
        return True

    while True:
        x = random.getrandbits(bits) | 1
        if is_prime(x):
            return x


def main() -> None:
    random.seed(131)
    b = 28  # 每个素数位数（示例较小，便于快速演示）
    p = gen_prime(b)
    q = gen_prime(b)
    N = p * q

    # 已知 p 的高位：取 k (> b/2) 位高位已知，构造 p0 为“保留高 k 位，其余置 0”
    k = (b * 3) // 5  # 约 0.6b，确保 |x|<2^{b-k} < 2^{b/2} ≈ N^{1/4}
    p0 = (p >> (b - k)) << (b - k)
    x_true = p - p0
    q0 = N // p0
    y_true = q - q0

    X = 1 << (b - k)
    Y = X  # 量级上 y_true 与 x_true 同阶

    # F(x,y) = (p0+x)(q0+y) - N = x y + q0 x + p0 y + p0 q0 - N
    F: Bivar = {
        (1, 1): 1,
        (1, 0): q0,
        (0, 1): p0,
        (0, 0): p0 * q0 - N,
    }

    roots = try_find_small_roots_bivar(F, N=N, X=X, Y=Y, m=2, tx=2, ty=2)

    recovers = []
    for x, y in roots:
        if 0 <= x < X and abs(y) < (Y * 4):  # 容许 y 稍放宽
            pp = p0 + x
            qq = q0 + y
            if pp * qq == N:
                recovers.append((pp, qq))
    print(
        {
            "N": N,
            "p": p,
            "q": q,
            "p0": p0,
            "q0": q0,
            "X": X,
            "Y": Y,
            "x_true": x_true,
            "y_true": y_true,
            "roots_found": roots[:10],
            "recovers": recovers[:3],
        }
    )

    assert any(pp == p and qq == q for (pp, qq) in recovers), (
        "未能恢复 (p,q)；可增大 m,tx,ty 或调整 b,k 再试"
    )
    print("OK: 通过已知高位成功因式分解 N")


if __name__ == "__main__":
    main()
