#!/usr/bin/env python3
from __future__ import annotations

# 整数域小根（无模数情形）：给定 f(x)，若存在很小整数根，可直接在界内搜索（教学）。
# 在论文背景下，这一部分通常不单独强调，但有助于理解“X 界”与搜索思路。


def main() -> None:
    # 对比基线：无模数情形，直接搜索根
    # f(x) = (x - r)*(x + r) = x^2 - r^2
    r = 17
    X = 64

    # 直接搜索
    brute = [x for x in range(-X + 1, X) if x * x - r * r == 0]

    print(
        {
            "X": X,
            "r_true": r,
            "bruteforce_roots": brute,
        }
    )


if __name__ == "__main__":
    main()
