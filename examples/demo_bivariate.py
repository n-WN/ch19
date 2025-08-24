#!/usr/bin/env python3
from __future__ import annotations

import random

# 二元小根教学演示：构造 F(x,y) = x^2 + y + c，使 (r,s) 为模 N 的小根
# 注意：真实 Coppersmith 二元需要更细致的参数与消元，这里只做教学演示
from coppersmith.bivar import Bivar
from coppersmith.bivariate import try_find_small_roots_bivar


def main() -> None:
    random.seed(7)
    # 小规模 N
    p = 499
    q = 547
    N = p * q

    X = 24
    Y = 24
    r = random.randrange(-X // 2, X // 2)
    s = random.randrange(-Y // 2, Y // 2)

    # F(x,y) = x^2 + y + c，其中 c ≡ - (r^2 + s) (mod N)
    c = (-(r * r + s)) % N
    F: Bivar = {
        (2, 0): 1,
        (0, 1): 1,
        (0, 0): c,
    }

    roots = try_find_small_roots_bivar(F, N=N, X=X, Y=Y, m=2, tx=2, ty=2)
    print(
        {
            "N": N,
            "X": X,
            "Y": Y,
            "(r,s)_true": (r, s),
            "roots_found": roots[:10],
        }
    )

    # (r,s) 在模 N 意义下是一解（多解同余类），枚举会包含 (r,s)
    assert (r, s) in roots, "未能找回 (r,s)；可调大参数或缩小 X,Y 再试"
    print("OK: 找回 (r,s)")


if __name__ == "__main__":
    main()
