#!/usr/bin/env python3
from __future__ import annotations

import random
from math import isqrt

from coppersmith.univariate import find_small_roots_univariate

# 演示：构造一个小根 r，模数 N=p*q，二次多项式 f(x)=x^2 + c
# 使得 f(r) ≡ 0 (mod N)。
# 我们选取 |r| < X，期待算法找回 r。


def gen_prime(bits: int) -> int:
    # 简单试除法，教学演示足够
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
        n = random.getrandbits(bits) | 1
        if is_prime(n):
            return n


def main() -> None:
    random.seed(42)
    p = gen_prime(14)
    q = gen_prime(14)
    N = p * q

    # 选择小根 r，界 X
    X = 2**8
    r = random.randrange(-X // 2, X // 2)

    # 选择 c 使得 r^2 + c ≡ 0 (mod N) -> c ≡ -r^2 mod N
    c = (-(r * r)) % N

    # f(x) = x^2 + c
    f_coeffs = [c, 0, 1]

    roots = find_small_roots_univariate(f_coeffs, N=N, X=X, m=3, t=3)
    print(
        {
            "N": N,
            "X": X,
            "r_true": r,
            "roots_found": roots,
        }
    )

    assert r in roots, "未能找回小根 r（可调大 m,t 或减小 X 再试）"
    print("OK: 找回 r")


if __name__ == "__main__":
    main()
