#!/usr/bin/env python3
from __future__ import annotations

import random

# RSA 小指数 e=3 小消息教学演示：找回 m 使 m^3 ≡ c (mod N)，且 m < N^{1/3}
# 经典用法：f(x) = x^3 - c，X ≈ N^{1/3}
from coppersmith.univariate import find_small_roots_univariate


def gen_small_rsa_modulus(bits: int = 36) -> int:
    # 生成小规模 N=p*q，仅为演示（非安全）。
    def gen_prime(nbits: int) -> int:
        from math import isqrt

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
            x = random.getrandbits(nbits) | 1
            if is_prime(x):
                return x

    p = gen_prime(bits // 2)
    q = gen_prime(bits // 2)
    return p * q


def main() -> None:
    random.seed(2025)
    N = gen_small_rsa_modulus(36)
    e = 3

    # 选择小消息 m，保证 m^3 < N
    # 令 m ≈ floor(N^{1/3})/4 作为保守选择
    X = int(round(N ** (1 / 3)))
    m_true = max(2, X // 4)
    c = pow(m_true, e, N)

    # f(x) = x^3 - c
    f_coeffs = [-c, 0, 0, 1]

    roots = find_small_roots_univariate(f_coeffs, N=N, X=X, m=4, t=3)
    print(
        {
            "N": N,
            "e": e,
            "X": X,
            "m_true": m_true,
            "roots_found": roots,
        }
    )

    assert m_true in roots, "未能找回 m；可增大 m,t 或调小 X 再试"
    print("OK: 找回 m")


if __name__ == "__main__":
    main()
