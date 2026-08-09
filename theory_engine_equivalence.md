# 引擎一致性的理论模型

> 配套脚本：`theory_engine_equivalence.py`（逐环节数值检验，日志 `theory_engine_equivalence_log.txt`）
> 与 `check_engine_vs_v45.py`（与 `lidar_histogram_sim_v45.ipynb` 的源码/比特级核对）。
>
> 本文回答：**为什么**「精确逐光子引擎」「快速引擎」「解析式」三者必然一致，
> 而不只是「实测数值碰巧对得上」。

---

## 0. 先把话说清楚：这里有四份代码，但只有一个数学对象

| 记号 | 代码 | 是什么 |
|---|---|---|
| **M** | —— | **连续时间数学模型**（下面第 1 节定义） |
| **E1** | `spad_binary_trace` | M 的**时间离散化**（步长 `DT_FINE`），逐光子抽样 |
| **E2** | `noise_macro_hist_fast` | M 的**精确采样器**（更新过程 + $H^{-1}$ 直查表） |
| **E3** | `p_bin_equilibrium` | M 的**平衡态闭式解** |
| **E4** | `binary_macro_stepping` | E1 的**向量化恒等变形**（同一离散模型，换个写法） |

要证的三件事：

1. **E2 与 E3 一致是恒等式，不是巧合**——它们是同一个概率律的「采样」与「解析」两种用法。
2. **E1 → M 的误差是 $O(\mathrm{d}t^2)$**，在生产参数 $\mathrm{d}t = 200\ \mathrm{ps}$ 下相对量级 $10^{-4}$。
3. **E4 ≡ E1**，靠一个初等恒等式。

第 1–7 节推导 M 与 E2、E3；第 8 节处理 E1/E4 的离散化；第 9 节处理抖动；
第 10 节把每个环节对应到数值证据。

---

## 1. 模型 M 的定义与记号

| 符号 | 含义 | 本项目取值 |
|---|---|---|
| $r_{\mathrm{amb}}$ | 单 SPAD 环境**光子到达率**（不含 PDE） | 由目标 bg 反解 |
| $\mathrm{PDE}$ | 满过电压时的探测效率 $=\mathrm{PDE}_{\max}$ | 0.30 |
| $r_{\det}$ | $\equiv r_{\mathrm{amb}}\cdot\mathrm{PDE}$ | —— |
| $\tau$ | RC 恢复时间常数 $R\,C_J$ | 8.7315 ns |
| $g(\cdot)$ | 响应函数，$g(0)=0,\ g(1)=1$ | `exp` 型，$g(x)=\frac{1-e^{-kx}}{1-e^{-k}}$，$k=3$ |
| $V_{\mathrm{th,frac}}$ | 计数阈值占满过电压比例 | 0.60 |
| $T_{\mathrm{OVER}}$ | 过阈窗宽 $=-\tau\ln(1-V_{\mathrm{th,frac}})$ | 8.0006 ns |
| $\sigma$ | IRF 抖动标准差 | 100 ps |

**M 的三条规则：**

1. 环境光子按**齐次泊松过程**到达，强度 $r_{\mathrm{amb}}$。
2. 设 $\Delta(t)=t-(\text{上一次雪崩时刻})$ 为**年龄**。过电压恢复为
   $$v(\Delta)=\frac{V_{ov}(\Delta)}{V_{ov,\max}}=1-e^{-\Delta/\tau}$$
   $t$ 时刻到达的光子**独立地**以概率 $\mathrm{PDE}\cdot g\big(v(\Delta(t))\big)$ 引发雪崩；
   雪崩发生则 $V_{ov}\leftarrow 0$，年龄归零。
3. 1 ns 采样点 $c$ 记 1 $\iff$ $(c-T_{\mathrm{OVER}},\,c]$ 内至少有一次雪崩。

> 注意规则 2 里的 $g$：**恢复期不是"完全不响应"，而是响应概率连续爬升**。
> 这正是用户关心的"电压在按 RC 恢复、响应能力也在恢复"。硬死时间模型对应的是
> $g=\mathbf{1}[\Delta>t_{\mathrm{dead}}]$ 这种阶跃，M 里没有这种阶跃。

---

## 2. 泊松稀释 ⟹ 条件强度只依赖"年龄"

对齐次泊松过程做**依赖历史的独立稀释**：每个点以概率 $p(t)$ 保留（$p$ 可依赖被保留点的历史），
则保留下来的点过程的条件强度为

$$\lambda\big(t \mid \mathcal F_{t^-}\big) = r_{\mathrm{amb}}\cdot p(t)$$

代入 $p(t)=\mathrm{PDE}\cdot g(v(\Delta(t)))$：

$$\boxed{\ \lambda\big(t\mid\mathcal F_{t^-}\big) = r_{\det}\; g\!\left(1-e^{-\Delta(t)/\tau}\right) \;\equiv\; h\big(\Delta(t)\big)\ }$$

**这一步是全部推导的枢纽**：条件强度对历史的依赖，**只**通过年龄 $\Delta(t)$ 进入，
不依赖更早的雪崩时刻、也不依赖绝对时间 $t$。原因是规则 2 里雪崩把 $V_{ov}$ 打回 **0**——
一个固定的状态，于是每次雪崩之后系统"失忆"。

---

## 3. 年龄依赖强度 ⟹ 更新过程，间隔律有闭式

条件强度只依赖年龄的点过程，按定义就是**更新过程（renewal process）**，
$h(\cdot)$ 就是间隔分布的**风险函数（hazard）**。于是相邻雪崩间隔 $X$ 满足

$$H(\Delta)=\int_0^{\Delta} h(s)\,\mathrm{d}s = \int_0^{\Delta} r_{\det}\,g\!\left(1-e^{-s/\tau}\right)\mathrm{d}s$$
$$S(\Delta)=\Pr(X>\Delta)=e^{-H(\Delta)},\qquad f_X(\Delta)=h(\Delta)\,S(\Delta)$$

且各间隔 **i.i.d.**（每次雪崩后状态完全相同）。

对应代码：`build_renewal_table(r_det)` 就是在数值积分这个 $H(\Delta)$。

小性质，后面要用：$h(0)=r_{\det}\,g(0)=0$。**刚雪崩完的瞬间，触发率恰好是 0**，
然后连续爬升。这条会在第 8 节把离散化误差压得极小。

---

## 4. 逆变换采样：快速引擎 E2 采的就是 $S$

若 $E\sim\mathrm{Exp}(1)$，令 $\Delta=H^{-1}(E)$，则

$$\Pr(\Delta>x)=\Pr\big(E>H(x)\big)=e^{-H(x)}=S(x)$$

所以 $\Delta$ 的分布**精确等于**第 3 节的间隔律。这就是 `noise_macro_hist_fast` 里那两行：

```python
E = rng.standard_exponential(m_, dtype=np.float32)     # E ~ Exp(1)
delta = inv[i0]*(1.0-fr) + inv[i0+1]*fr                # Δ = H⁻¹(E)，直查表 + 线性插值
```

`build_inv_table` 只是把 $H^{-1}$ 重采样到均匀 $E$ 网格，好让查表变成 $O(1)$——
是**实现优化**，不是模型近似。

> **T1 实测**（$n=2\times10^6$）：在 $x=0.5,1,2,3,5,8,12,20,30,50,80$ ns 十一个点上，
> 经验生存函数与 $e^{-H(x)}$ 的最大相对偏差 **0.110%**，而该样本量的 MC 相对误差量级为 0.071%。

---

## 5. 过阈窗的并集 = 不交并——差分数组技巧是恒等式

采样点 $c$ 点亮 $\iff$ $c$ 落在某个过阈窗内，即落在

$$U=\bigcup_k \big[a_k,\ a_k+T_{\mathrm{OVER}}\big)$$

其中 $a_1<a_2<\cdots$ 为雪崩时刻。窗口之间会重叠（这就是"顺延堆积"），
直接累加会重复计数。引擎的做法是把每个窗**在下一次雪崩处截断**：

$$U=\bigsqcup_k \Big[a_k,\ \min\big(a_k+T_{\mathrm{OVER}},\ a_{k+1}\big)\Big)$$

**证明**（两行）：记右端为 $U'$。
$U'\subseteq U$ 显然（每个区间都被对应的原窗包含）。
反之设 $x\in U$，取**最大**的 $k$ 使 $a_k\le x$，则 $x<a_{k+1}$；
又 $x\in U$ 意味着存在 $j$ 使 $a_j\le x<a_j+T_{\mathrm{OVER}}$，而 $a_j\le x$ 蕴含 $j\le k$，
于是 $a_k\ge a_j$ 给出 $x<a_j+T\le a_k+T$。故 $x\in[a_k,\min(a_k+T,a_{k+1}))\subseteq U'$。∎

区间两两不交后，"落在哪些 bin"就能用**差分数组 + `bincount` + `cumsum`** 一次算完，
这正是 `noise_macro_hist_fast` 末尾那几行。**它是恒等变形，不引入任何近似。**

> **T3 实测**：把同一批雪崩时刻分别用「暴力涂并集」和「截断+差分数组」处理，
> 无抖动 **3992/3992** 条轨迹逐 bin 完全相同；有抖动（$\sigma=100$ ps）**3987/3987** 完全相同。

**唯一的例外**：抖动可能把相邻雪崩的先后顺序颠倒（需要间隔 $\lesssim\sigma$），
此时 $b_{hi}\le b_{lo}$，该区间被 `b_hi > b_lo` 掩码丢弃。
其概率上界为 $\Pr(X<10\sigma=1\ \mathrm{ns})=5.15\times10^{-3}$（每个间隔），
且即便乱序，在 1 ns 采样格上通常也看不出差别——实测 0 次可见偏差。

---

## 6. 平衡态：更新-回报定理给出解析式 E3

定义**后向回溯时间** $B(c)=c-\max\{a_k\le c\}$。由规则 3，

$$\text{bin 点亮}\iff B(c)<T_{\mathrm{OVER}}$$

更新过程的后向回溯时间有标准的平稳（平衡态）密度

$$f_B(u)=\frac{S(u)}{\mu},\qquad \mu=\mathbb E[X]=\int_0^{\infty}S(u)\,\mathrm{d}u$$

于是

$$\boxed{\ p_{\mathrm{bin}}=\Pr\big(B<T_{\mathrm{OVER}}\big)=\frac{1}{\mu}\int_0^{T_{\mathrm{OVER}}}S(u)\,\mathrm{d}u\ }$$

这就是 `p_bin_equilibrium` 里的那行积分。

**等价的第二种写法（更新-回报定理）**：每个更新周期长度为 $X$，其中被点亮的时长是
$\min(X,\,T_{\mathrm{OVER}})$（若 $X<T$，下一次雪崩接着点亮，本周期全亮），故

$$p_{\mathrm{bin}}=\frac{\mathbb E\big[\min(X,T_{\mathrm{OVER}})\big]}{\mathbb E[X]}$$

两式相等，因为对非负随机变量有恒等式 $\mathbb E[\min(X,c)]=\int_0^c S(u)\,\mathrm{d}u$。

> **T2 实测**（$r_{\det}=3.2985\times10^7$ cps，即 noise27 = 6）：
>
> | 算法 | $p_{\mathrm{bin}}$ | 相对差 |
> |---|---|---|
> | ① `p_bin_equilibrium` 解析式 | 0.222266 | —— |
> | ② $(1/\mu)\int_0^T S$ 独立数值积分 | 0.222200 | −0.030% |
> | ③ $\mathbb E[\min(X,T)]/\mathbb E[X]$（更新-回报，采样） | 0.222229 | −0.017% |
> | ④ 快速引擎 MC（400,000 条） | 0.222332 | +0.030% |
>
> 平均间隔：解析 33.7494 ns vs 采样 33.7571 ns（+0.023%）。
>
> ②③④ 指向同一个数 ⟹ **快速引擎采样的就是解析式所描述的那个平衡态分布**。
> 剩下的 0.03% 是数值积分与 MC 噪声，不是模型差异。

---

## 7. 宏像元：为什么 $n_{\mathrm{pix}}$ 与 $N_{\mathrm{shots}}$ 可以折成 $n_{\mathrm{tr}}$

纯环境光下：不同 SPAD 看到的是**独立**的光子流、有**各自独立**的复位状态；
不同 shot 之间也彼此独立、且每发开始时都从同一初态出发。
因此 $n_{\mathrm{pix}}\times N_{\mathrm{shots}}$ 条二值轨迹是 **i.i.d.** 的，
它们进入结果的方式**只有"有多少条"这一个自由度**：

$$n_{\mathrm{tr}} \equiv n_{\mathrm{pix}}\cdot N_{\mathrm{shots}}$$

推论一（单 bin 边缘律，**精确**）：

$$\text{count}(c)=\sum_{k=1}^{n_{\mathrm{tr}}}\mathbf 1\{B_k(c)<T_{\mathrm{OVER}}\}\ \sim\ \mathrm{Binomial}\big(n_{\mathrm{tr}},\,p_{\mathrm{bin}}\big)$$

推论二：整条直方图的联合分布（因而 peak 分布、FAR 阈值）**只是 $n_{\mathrm{tr}}$ 与 $r_{\det}$ 的函数**。
这就是 `compare_macro_3x9_vs_3x6.py` 里
「3×6 跑 6 发 ≡ 3×9 跑 4 发（都是 $n_{\mathrm{tr}}=108$）」的**证明**，
而那边 24 档环境光上 `T@1%` 最大差 0 计数是它的**实证**。

> **T4 实测**（noise27 = 6，$p_{\mathrm{bin}}=0.22227$，各 120,000 条）：
>
> | $n_{\mathrm{tr}}$ | 均值（理论） | 方差（理论） | 与二项分布的总变差距离 |
> |---|---|---|---|
> | 27 | 6.0056（6.0012） | 4.6598（4.6673） | 0.0021 |
> | 36 | 7.9977（8.0016） | 6.2361（6.2231） | 0.0033 |
> | 108 | 23.9945（24.0047） | 18.6188（18.6693） | 0.0048 |

> ⚠️ 这只是**单个 bin 的边缘分布**。bin 之间因 $T_{\mathrm{OVER}}=8$ ns 而**正相关**
> （一次雪崩点亮约 8 个连续 bin），所以 peak **不能**按"152 个独立二项取最大"来算，
> 有效独立 bin 数 $M_{\mathrm{eff}}\approx 46\text{–}76$。详见 `theory_peak_bg_multishot.md`。

---

## 8. 精确引擎 E1 的离散化：把离散模型也精确解出来

E1 不是 M，而是 M 在步长 $\mathrm{d}t$（`DT_FINE` = 200 ps）上的离散化：
光子被归到网格点 `t_arr = np.repeat(tf, n_ph)`，年龄按网格点计算。

### 8.1 一步内至多一次雪崩，且概率有闭式

一个网格步内到达 $n\sim\mathrm{Poisson}(\mu)$ 个光子，$\mu=r_{\mathrm{amb}}\,\mathrm{d}t$，
各自以 $\varphi(a)=\mathrm{PDE}\cdot g(1-e^{-a/\tau})$ 触发（$a$ = 步初的年龄）。
**首个触发后年龄归零 ⟹ $g(0)=0$ ⟹ 同一步内后续光子不可能再触发**。于是

$$\Pr(\text{本步至少一次雪崩}\mid a)=1-\mathbb E\big[(1-\varphi)^n\big]=1-e^{-\mu\varphi(a)}\equiv q(a)$$

用到了泊松的概率母函数 $\mathbb E[z^n]=e^{-\mu(1-z)}$。

**这同时证明了 E4 ≡ E1**：`binary_macro_stepping` 里那行
`p = -np.expm1(-mu_all*phi[age])` 就是上式，只不过它按时间步向量化，
而 E1 是逐个光子抽。两者是同一个离散模型。

### 8.2 离散模型是离散时间更新过程，可精确求解

设一次雪崩后第 $m$ 步的年龄为 $m\,\mathrm{d}t$，则 $q(0)=0$ 且

$$\Pr(X_d>m\ \text{步})=\prod_{j=1}^{m}\big(1-q(j\,\mathrm{d}t)\big),\qquad
\mu_d=\mathrm{d}t\sum_{m\ge 0}\Pr(X_d>m)$$

过阈窗宽 $T_{\mathrm{OVER}}$ 是**连续量、与 $\mathrm{d}t$ 无关**（代码里就是按真实 8.0006 ns 涂窗），
故

$$p_{\mathrm{bin}}^{(d)}=\frac{1}{\mu_d}\int_0^{T_{\mathrm{OVER}}}\Pr(X_d\,\mathrm{d}t>u)\,\mathrm{d}u$$

（被积函数是阶梯函数，末段按 $T_{\mathrm{OVER}}-K\,\mathrm{d}t$ 补齐，$K=\lfloor T_{\mathrm{OVER}}/\mathrm{d}t\rfloor$。）

### 8.3 数值结果：误差是 $O(\mathrm{d}t^2)$，生产参数下可忽略

**T5a 实测**，noise27 = 6（连续模型 $p_{\mathrm{bin}}=0.222266$，$\mu=33.7494$ ns）：

| $\mathrm{d}t$ [ps] | 离散 $p_{\mathrm{bin}}$ | 相对差 | 离散 $\mu$ [ns] | $\mu$ 相对差 | 误差比上一行 |
|---|---|---|---|---|---|
| 3200 | 0.221882 | −0.1727% | 34.0925 | 1.0164% | — |
| 1600 | 0.222427 | +0.0725% | 33.8377 | 0.2615% | 3.89 |
| 800 | 0.222322 | +0.0251% | 33.7717 | 0.0659% | 3.97 |
| 400 | 0.222295 | +0.0130% | 33.7550 | 0.0165% | 3.99 |
| **200（生产值）** | **0.222288** | **+0.0100%** | **33.7508** | **0.0041%** | **4.00** |
| 100 | 0.222287 | +0.0092% | 33.7498 | 0.0010% | 4.00 |
| 50 | 0.222286 | +0.0090% | 33.7495 | 0.0003% | 4.00 |

noise27 = 12 完全同样的形态（$\mu$ 误差 2.5977% → 0.0007%，比值同样是 3.89/3.97/3.99/4.00/4.00/4.00）。

两个结论：

- **$\mathrm{d}t$ 每减半，误差降 4 倍 ⟹ 离散化误差是 $O(\mathrm{d}t^2)$**，不是 $O(\mathrm{d}t)$。
- 生产参数 $\mathrm{d}t=200$ ps 处，$p_{\mathrm{bin}}$ 偏差 $\approx 1\times10^{-4}$ 相对量级，
  **已低于 `p_bin_equilibrium` 自身数值积分的精度**（第 6 节 ② 与 ① 差 0.03%），
  所以 MC 根本不可能分辨出来。

**为什么误差这么小**：最大的离散化误差源本来应该是"一步内本该发生两次雪崩，却只算了一次"。
但 $h(0)=0$ 把它自动压掉了——刚雪崩完 $g\approx 0$，紧接着的 200 ps 内
$\varphi\approx 0.3\times g(0.0226)=0.021$，一步内触发概率
$q\approx 1-e^{-0.022\times 0.021}\approx 4.6\times10^{-4}$。
**RC 恢复模型自己保证了它的离散化是良态的**；换成硬死时间（$g$ 在 $t_{\mathrm{dead}}$ 处跳变）
反而会有 $O(\mathrm{d}t)$ 的边界误差。

> **T5b 交叉验证**（noise27 = 6，各 6000 条）：
>
> | $\mathrm{d}t$ [ps] | 逐光子 E1 | 步进 E4 | 离散模型解析 |
> |---|---|---|---|
> | 800 | 0.22182 ± 0.00114 | 0.21979 ± 0.00116 | 0.22232 |
> | 200 | 0.22444 ± 0.00115 | 0.22152 ± 0.00116 | 0.22229 |
>
> E1 与解析值一致（0.4σ / 1.9σ）。E4 在 $\mathrm{d}t=800$ ps 偏低 2.2σ、
> 在 200 ps 偏低 0.7σ，未达显著，但**方向一致偏低**，怀疑与它
> 「先出 bin、再处理本步雪崩」的半步对齐有关。
> 纯噪声扫描用的是 E2 不是 E4，不受影响；但 **PoD 的信号部分用的是 E4**，
> 这一点已记入待办。

---

## 9. 抖动为什么可以在解析式里忽略

抖动给每次雪崩时刻加 i.i.d. 的 $\varepsilon_k\sim\mathcal N(0,\sigma^2)$，
再涂宽度仍为 $T_{\mathrm{OVER}}$ 的窗。它**不改变**点过程的强度（平稳过程做 i.i.d. 平移仍平稳），
但会改变**窗口并集的测度**，因为重叠结构变了。

单个周期的点亮时长从 $\min(X,T)$ 变成 $\min(X+\delta,T)$，$\delta=\varepsilon_{k+1}-\varepsilon_k$，
$\mathbb E[\delta]=0$，$\mathrm{Var}[\delta]=2\sigma^2$。因为 $\min(\cdot,T)$ 是**凹函数**，
由 Jensen 不等式抖动**永远使点亮测度略微减小**。亏损只在 $X\approx T_{\mathrm{OVER}}$ 附近非零，
积分后的总亏损约为

$$\mathbb E[\text{亏损}]\ \approx\ f_X(T_{\mathrm{OVER}})\cdot\sigma^2$$

代入 noise27 = 6：$h(T_{\mathrm{OVER}})=r_{\det}\,g(0.600)=3.2985\times10^7\times0.8784=2.898\times10^{7}$，
$S(T_{\mathrm{OVER}})=0.8482$，故 $f_X(T_{\mathrm{OVER}})=2.458\times10^{7}\ \mathrm{s^{-1}}$，

$$\frac{\mathbb E[\text{亏损}]}{\mathbb E[\min(X,T)]}\approx\frac{2.458\times10^{7}\times(10^{-10})^2}{0.2223\times 33.75\times10^{-9}}\approx 3.3\times10^{-5}$$

即抖动使 $p_{\mathrm{bin}}$ 下降约 **0.003%**。所以 `p_bin_equilibrium` 不含抖动项是合理的，
而快速引擎里保留抖动（它影响 peak 的**空间结构**，不只是 $p_{\mathrm{bin}}$）也不矛盾。

---

## 10. 推导环节 ↔ 数值证据 一览

| 环节 | 命题 | 代码 | 检验 | 结果 |
|---|---|---|---|---|
| §2 | 泊松稀释 ⟹ 强度只依赖年龄 | —— | 解析 | 恒等 |
| §3–4 | $\Delta=H^{-1}(E)$ 的分布 $=e^{-H}$ | `build_inv_table` | **T1** | 最大偏差 0.110%（MC 量级 0.071%） |
| §5 | 并集 = 截断后的不交并 | `noise_macro_hist_fast` | **T3** | 3992/3992、3987/3987 逐 bin 相同 |
| §6 | $p_{\mathrm{bin}}=\frac1\mu\int_0^T S=\frac{\mathbb E[\min(X,T)]}{\mathbb E[X]}$ | `p_bin_equilibrium` | **T2** | 四法一致，最大差 0.030% |
| §7 | 单 bin 边缘 $=\mathrm{Binomial}(n_{\mathrm{tr}},p)$ | 折叠 $n_{\mathrm{tr}}$ | **T4** | 总变差距离 0.002–0.005 |
| §8 | E1 = M 的 $O(\mathrm{d}t^2)$ 离散化 | `spad_binary_trace` | **T5a** | 每减半降 4.00 倍；200 ps 处 $10^{-4}$ |
| §8.1 | E4 ≡ E1 | `binary_macro_stepping` | **T5b** | 一致（E4 偏低 0.7–2.2σ，见备注） |
| §9 | 抖动可忽略 | —— | 解析 | $3.3\times10^{-5}$ |

另有跨文件的端到端核对（`check_engine_vs_v45.py`）：
PoD 的 `spad_binary_trace` 与 v45 cell 32 同名函数**同种子下 60/60 条轨迹逐 bin 完全相同**。

---

## 11. 这套推导覆盖不到的地方（别误用）

1. **有信号时 M 不再是更新过程**。信号率 $r_{\mathrm{sig}}(t)$ 依赖**绝对时间**，
   条件强度变成 $h(t,\Delta)=\big(r_{\mathrm{sig}}(t)+r_{\mathrm{amb}}\big)\mathrm{PDE}\cdot g(v(\Delta))$，
   不再只依赖年龄 ⟹ §3–§6 全部失效，**没有闭式**，只能用 E1/E4 步进。
   这就是为什么信号部分用 `binary_macro_stepping_per_shot` 而不是快速引擎——
   不是懒得优化，是数学上不允许。
2. **§7 只给出单 bin 边缘律**。bin 间正相关没有闭式，$M_{\mathrm{eff}}$ 只能靠 MC。
3. **§8 的 $O(\mathrm{d}t^2)$ 依赖 $g(0)=0$**。若把响应函数换成硬死时间阶跃，
   结论不成立，需重新评估 `DT_FINE`。
4. **T5b 里 E4 系统性偏低**（虽未达显著）。E4 服务于 PoD 的信号支路，
   建议后续用更大样本单独复核其半步对齐。
