# Coppersmith 方法讲义：从直觉到实现（单变量、二元、消元与应用）

本讲义基于 Steven Galbraith《Mathematics of Public Key Cryptography》ch19（Coppersmith’s Method and Related Applications）的公开章节，提供一份针对“会一点数论（会看模运算、会 gcd），但不懂格”的读者的系统说明。配套代码为纯手工、零外部依赖的最小可运行实现：单变量小根、二元小根，以及二元消元（结果式）。

重要说明：本仓库仅用于教学演示，代码注重可读性而非性能或鲁棒性，请勿用于生产或安全关键用途。

---

## 0. 受众与前置知识
- 你需要知道：模运算（≡ 与 mod 的含义）、整除、最大公因数 gcd、基本多项式代数。
- 你暂时不需要知道：格与 LLL 的细节证明、代数数论、$p-1$ 素数判定等高级概念。

我们会用“尽量直觉”的方式解释为什么“构造一堆看似无关的多项式 + 一个神奇的短向量算法（LLL）”可以把“模方程的小根”变成“整数方程的真根”。

---

## 1. 问题设置（单变量小根）
给定：
- 模数 $N\in\mathbb{Z}_{>0}$
- 整系数多项式 $f(x)\in\mathbb{Z}[x]$，次数 $d=\deg f$
- 目标：若存在“足够小”的整数根 $x_0$（大小界 $X$）满足
  $$f(x_0) \equiv 0 \pmod{N}, \quad |x_0| < X,$$
  希望能在多项式时间内找回 $x_0$。

经典结论（非形式化）：若 $X\lesssim N^{1/d}$（更严格地说 $|x_0| < N^{1/d-\varepsilon}$），则可以用 Coppersmith 方法在多项式时间内恢复 $x_0$。

---

## 2. 关键思想：让“模 N 的零”变成“整数上的零”
直观目标：构造若干整系数多项式 $g_i(x)$，使它们在 $x_0$ 处“同时很小”且“模 $N^m$ 为 0”。如果某个整数多项式 $h(x)$ 同时满足
- $h(x_0) \equiv 0 \pmod{N^m}$（强整除），
- 且 $|h(x_0)| < N^m$（数值上比模数还小），

那么 $h(x_0)$ 只能等于 0（因为它是整数，既是 $N^m$ 的倍数，又严格小于 $N^m$ 的绝对值）。这就把“模方程的根”变成了“整数方程的真根”。

如何达成这两点？
- 让模 $N^m$ 为 0：若 $f(x_0) \equiv 0 \pmod N$，则 $f(x_0) = N\cdot t$。于是 $\forall\, i<m$，$N^{m-i}\, f(x)^i$ 在 $x_0$ 处必定被 $N^m$ 整除。再加上 $f(x)^m$ 自身也被 $N^m$ 整除。
- 让“数值上很小”：把变量缩放 $x\mapsto X\cdot u$，并把多项式的“第 $k$ 列系数”统一乘以 $X^k$。这样，一个多项式的“系数向量范数”就能很好地近似它在区间 $|x|\le X$ 上的大小。若 $|x_0|<X$，则 $|u_0|=|x_0|/X<1$。寻找“系数向量很短”的整数线性组合（这正是 LLL 擅长的事），常能得到在 $u_0$ 处数值很小的多项式。

这两点叠加，就得到前述“既是 $N^m$ 倍数又很小”从而“只能为 0”的结论。

---

## 3. 与论文一致的格构造（Howgrave–Graham 变体）
记 $d=\deg f$。我们构造如下多项式集合（均为整系数）：
$$
\begin{aligned}
\mathcal{G}
= {} & \big\{\, N^{m-i}\, x^j\, f(x)^i : 0\le i<m,\, 0\le j<d \,\big\} \\
&\cup\, \big\{\, x^j\, f(x)^m : 0\le j<t \,\big\}.
\end{aligned}
$$
- 对任意 $x_0$，若 $f(x_0)\equiv 0\pmod N$，则上集合中任意元素在 $x_0$ 处都被 $N^m$ 整除。
- 将变量替换为 $x\mapsto X\cdot u$，并对“系数矩阵”的第 $k$ 列乘以 $X^k$，得到一个整数矩阵 $B$（每行就是一个多项式的系数向量）。
- 使用 LLL 对 $B$ 做格基约化，得到若干“短向量”。这些短向量对应的线性组合在 $|u|\le 1$（等价于 $|x|\le X$）范围内通常很小。

设约化得到的某一行对应多项式为 $h(Xu)$。当 $u=u_0=x_0/X$ 时：
- $h(x_0)$ 被 $N^m$ 整除；
- 由于“短向量”导致系数小，进而 $|h(x_0)|$ 小；
于是 $h(x_0)=0$，即 $x_0$ 是 $h$ 的整数根。随后在区间 $(-X,X)$ 中直接搜索/验证即可。

> 我们的实现（coppersmith/univariate.py）按上述集合构造矩阵，LLL 后对前若干短向量反缩放评估，并在 $(-X,X)$ 内检查整数根且验证 $f(r)\equiv 0\pmod N$。

---

## 4. 参数如何选？
- 理论尺度：单变量 $d=\deg f$，通常需要 $X \leq N^{1/d}$ 才有希望成功。
- 起步建议：$m\approx d$，$t\approx d$。失败时可逐步增大 $m,t$（格维度提高，成功率提高但计算更慢），或调小 $X$。
- 维度估计：行数约 $m\cdot d + t$；列数约为“构造集中最高次数 + 1”。

在我们的教学工程中，示例采用小规模 N 与保守的 X，以保证几秒内完成演示。

---

## 5. 二元小根：$F(x,y)\equiv 0\pmod N$，$|x|<X,\,|y|<Y$
思路与单变量相同：构造在 $(x_0,y_0)$ 处被 $N^m$ 整除的一族多项式，并做双变量缩放：列 $(i,j)$ 统一乘以 $X^i Y^j$。我们采用 Howgrave–Graham 风格的集合（$d_x=\deg_x F+1,\ d_y=\deg_y F+1$）：
$$
\begin{aligned}
\mathcal{H}
= {} & \Big\{ N^{m-i}\, F(x,y)^i\, x^{a_x} y^{a_y}
: 0\le i<m,\ 0\le a_x<d_x,\ 0\le a_y<d_y \Big\} \\
&\cup\, \Big\{ F(x,y)^m\, x^{a_x} y^{a_y}
: 0\le a_x<t_x,\ 0\le a_y<t_y \Big\}.
\end{aligned}
$$
- 变量替换 $(x,y)\mapsto (X\,u, Y\,v)$，对列 $(i,j)$ 乘以 $X^i Y^j$，得到整数矩阵并 LLL。
- 取两条最短向量对应的 $G_1,G_2$，它们在 $(x_0,y_0)$ 处“既小又被 $N^m$ 整除”，于是满足 $G_1(x_0,y_0)=G_2(x_0,y_0)=0$ 的强条件（在整数上为 0）。

问题变为：从 $G_1,G_2$ 中消去 $y$（或 $x$），得到单变量多项式 $R(x)$，再在 $(-X,X)$ 搜索 $R(x)=0$ 的整数根并回代求 $y$。

---

## 6. 消元（结果式）：为什么以及怎么做（我们如何手写）
- 定义（思想）：对关于 $y$ 的两一元多项式 $P(y),Q(y)$，其结果式（关于 $y$）记为 $\operatorname{Res}_y(P,Q)$。性质：
  $$\operatorname{Res}_y(P,Q)=0\ \Longleftrightarrow\ \exists\ y\in\mathbb{C},\ P(y)=Q(y)=0.$$
  结果式可用 Sylvester 矩阵的行列式表示。
- 直接在 $\mathbb{Z}[x]$ 中构造 $\operatorname{Res}_y\big(G_1(x,y),G_2(x,y)\big)$ 会非常巨大。我们采用“多点专化 + 插值”的朴素而稳健的办法：
  1) 固定若干 $x_0$，把 $G_1,G_2$ 在 $x=x_0$ 处专化成关于 $y$ 的一元多项式（分数系数），清分母为整数；
  2) 构造 Sylvester 矩阵并用 Bareiss 无分式消元（整数算法）求行列式，得到一个整数值 $R'(x_0)$；
  3) 收集足够多的点 $\{(x_0, R'(x_0))\}$，用拉格朗日插值重建 $R'(x)$。此 $R'(x)$ 与真实 $R(x)$ 只差一个非零常数因子，不影响求解 $R(x)=0$ 的整数根。

> 这些都在 coppersmith/elimination.py 中手写实现：Bareiss 行列式、Sylvester 矩阵、插值与评估。

---

## 7. 我们的实现路线图（对应文件）
- coppersmith/poly.py：整数多项式基本运算。
- coppersmith/lll.py：Fraction 版 LLL，包含 Gram–Schmidt、size reduction、Lovász 条件检查。
- coppersmith/univariate.py：单变量小根（Howgrave–Graham 变体），列缩放与反缩放评估，区间搜索验证。
- coppersmith/bivar.py：二元多项式运算（加、乘、幂、移位、评估）。
- coppersmith/bivariate.py：二元小根（格构造、列缩放、LLL、两式消元、回代验证）。
- coppersmith/elimination.py：Bareiss 行列式、Sylvester 矩阵、插值求结果式。
- examples/：若干可运行的经典/教学案例。
- scripts/run_demos.sh：一键运行，采用安全 shell 规范（set -euo pipefail 等）。

---

## 8. 典型应用建模与复现实验
### 8.1 RSA 小指数（e=3）的小消息
- 模型：$c \equiv m^3 \pmod N$，且 $m < N^{1/3}$。
- 单变量建模：$f(x)=x^3-c$，选 $X\approx \lfloor N^{1/3} \rfloor$，按第 3 节流程求小根。
- 我们的演示：examples/demo_rsa_small_e.py（也演示了纯 CRT + 整数立方根的 Hastad 广播案例）。

### 8.2 已知素因子高位（Partial Key Exposure）因式分解
- 设 $N=pq$，已知 $p$ 的高 $k$ 位，写 $p=p_0+x$，其中 $|x|<2^{b-k}$（$b$ 为位长）。令 $q \approx N/p_0$，写 $q=q_0+y$。
- 二元建模：
  $$F(x,y) = (p_0+x)(q_0+y) - N = xy + q_0 x + p_0 y + (p_0q_0 - N).$$
  目标是在 $|x|<X,\ |y|<Y$ 内恢复 $(x,y)$，进而恢复 $(p,q)$。
- 要点：若 $k>b/2$，则 $X=2^{b-k} < 2^{b/2}\approx N^{1/4}$，通常在教学参数下可行。
- 我们的演示：examples/demo_factor_highbits.py。

### 8.3 二元教学例：$F(x,y)=x^2+y+c$（构造小根）
- 取 $c\equiv - (r^2+s)\pmod N$，则 $(r,s)$ 为小根。用二元流程 + 结果式消元可找回 $(r,s)$。
- 演示：examples/demo_bivariate.py。

### 8.4 Hastad 广播（e=3，无填充）
- 三个互素模数 $N_1,N_2,N_3$，同一明文 $m$，$c_i\equiv m^3\pmod{N_i}$。若 $m^3 < N_1N_2N_3$，则 CRT 合并得 $C=m^3$，整数立方根恢复 $m$（与 LLL 无关，但常与小根场景一起讲解）。
- 演示：examples/demo_hastad_broadcast.py。

---

## 9. 如何运行
- 一键脚本（推荐）：
```bash
bash scripts/run_demos.sh
```
- 单独运行任一示例：
```bash
python -m examples.demo_univar
python -m examples.demo_bivariate
python -m examples.demo_integer_smallroot
python -m examples.demo_rsa_small_e
python -m examples.demo_factor_highbits
python -m examples.demo_hastad_broadcast
```

---

## 10. 常见问题（FAQ）
- 找不到根？
  - 适当减小 $X$，或增大 $m,t$（二元中增大 $m,t_x,t_y$）。
  - 尝试提高 $N$ 规模但让根更“小”（更符合 $X\lesssim N^{1/d}$）。
  - 检查多项式是否按“升幂系数”传入；确认验证条件 $f(r)\bmod N=0$ 是否成立。
- 结果式插值失败/退化？
  - 特殊 $x_0$ 可能导致专化降阶或全零，跳过该点并取更多采样点。
  - 结果式只差一个常数因子，不影响找根；我们用整数插值并做最大公因子规约来稳定系数。
- 出现负根或对称根（如 $\pm r$）？
  - 由多项式结构决定，属正常现象。

---

## 11. 安全提示
- 本仓库仅为教学演示，严禁直接用于攻击真实系统。
- 现代密码协议会使用填充与防护（如 RSA-OAEP），避免本讲义中的“理想化弱设置”。

---

## 12. 参考
- Galbraith, Steven. Mathematics of Public Key Cryptography, Chapter 19.
- Coppersmith, D. Small solutions to polynomial equations, and low exponent RSA vulnerabilities.
- Howgrave-Graham, N. Finding small roots of univariate modular equations revisited.
- Hastad, J. On using RSA with low exponent.

---

## 13. 项目结构（供对照）
```
.
├── coppersmith/
│   ├── __init__.py
│   ├── poly.py                 # 整数多项式工具
│   ├── lll.py                  # Fraction 版 LLL
│   ├── univariate.py           # 单变量小根
│   ├── bivar.py                # 二元多项式运算
│   ├── bivariate.py            # 二元小根 + 结果式消元流程
│   └── elimination.py          # Bareiss 行列式 + Sylvester + 插值
├── examples/
│   ├── demo_univar.py
│   ├── demo_bivariate.py
│   ├── demo_integer_smallroot.py
│   ├── demo_rsa_small_e.py
│   ├── demo_factor_highbits.py
│   └── demo_hastad_broadcast.py
├── scripts/
│   └── run_demos.sh
├── README.md
└── pyproject.toml
```
