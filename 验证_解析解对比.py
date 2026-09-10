# -*- coding: utf-8 -*-
"""
A 题求解器的解析解交叉验证（独立于主程序，可单独复现）

包含三个算例：
  T1 常系数 + 第一类(Dirichlet)边界   → 圆柱 Bessel 级数（λ 为 J0 零点）
  T2 常系数 + 第三类(对流)边界 Bi=hR/k → 级数解（λ 满足 λJ1(λ)=Bi·J0(λ)）
  T3 常数解保持（T∞=T0、C∞=C0、Ṙ=0）
另附 C 场（水分）的同构 Dirichlet 检验。

用法：python 验证_解析解对比.py
"""
import numpy as np
from scipy.special import j0, j1, jn_zeros
from scipy.optimize import brentq
import 药材烘干求解 as M

R, kk, rho, cp = M.R0, M.Q1["k"], M.Q1["rho"], M.Q1["cp"]
alpha = kk / (rho * cp)


def analytic_dirichlet(lam, r, t, Ts, T0):
    Fo = alpha * t / R ** 2
    s = np.zeros_like(np.asarray(r, float))
    for L in lam:
        s = s + (2.0 / (L * j1(L))) * j0(L * np.asarray(r, float) / R) * np.exp(-(L ** 2) * Fo)
    return Ts + (T0 - Ts) * s


def analytic_convective(Bi, r, t, Ts, T0, n=30):
    """λ J1(λ) = Bi J0(λ) 的级数解（中心与一般 r）"""
    f = lambda L: L * j1(L) - Bi * j0(L)
    xs = np.linspace(1e-8, 60.0, 6000)
    v = f(xs)
    lam = []
    for i in range(len(xs) - 1):
        if v[i] * v[i + 1] < 0:
            lam.append(brentq(f, xs[i], xs[i + 1]))
        if len(lam) >= n:
            break
    lam = np.array(lam)
    An = 2.0 * Bi / ((lam ** 2 + Bi ** 2) * j0(lam))
    Fo = alpha * t / R ** 2
    s = np.zeros_like(np.asarray(r, float))
    for A, L in zip(An, lam):
        s = s + A * j0(L * np.asarray(r, float) / R) * np.exp(-(L ** 2) * Fo)
    return Ts + (T0 - Ts) * s


def props_const(rho_, cp_, k_, D_):
    def f(C, T):
        C = np.asarray(C, float)
        return (np.full_like(C, rho_), np.full_like(C, cp_),
                np.full_like(C, k_), np.full_like(C, D_))
    return f


def main():
    print("=" * 76)
    print("A 题求解器 解析解交叉验证    R=%.4f m  alpha=%.6e m^2/s" % (R, alpha))
    print("=" * 76)
    ok = True

    # ---------------- T1：Dirichlet 边界（h 极大） ----------------
    print("\n[T1] 常系数 + Dirichlet 边界 vs J0 零点级数解（中心温度）")
    lam0 = jn_zeros(0, 40)
    Ts = M.T0 + 20.0
    env = lambda t: (np.full(np.shape(t), Ts), np.full(np.shape(t), M.C0))
    print("  t(s)     numeric    analytic    rel.err(相对阶跃幅度)   rel.err(相对当前变化)")
    for t_end in (300.0, 600.0, 1800.0):
        cfg = M.SolverCfg(t_end=t_end, dt=1.0, N=80, props=props_const(rho, cp, kk, 0.0),
                          env_fun=env, out_every=t_end, tol=1e-12, h=1e5, hm=0.0)
        num = float(M.solve_1d(cfg)["T"][-1][0])
        ana = analytic_dirichlet(lam0, 0.0, t_end, Ts, M.T0) - M.K0
        rel = abs(num - ana) / abs(Ts - M.T0)              # 相对阶跃幅度（稳定指标）
        cur = abs(num - ana) / max(abs(ana - (M.T0 - M.K0)), 1e-12)
        print("  %6.0f  %9.4f  %9.4f  %19.4f%%  %22.4f%%" % (t_end, num, ana, rel * 100, cur * 100))
        ok &= rel < 0.01

    # ---------------- T2：对流边界（有限 Biot） ----------------
    Bi = M.Q1["h"] * R / kk
    print("\n[T2] 常系数 + 对流边界 Bi=hR/k=%.4f vs 级数解（中心温度）" % Bi)
    Ts2 = M.T0 + 13.5
    env2 = lambda t: (np.full(np.shape(t), Ts2), np.full(np.shape(t), M.C0))
    print("  t(s)     numeric    analytic    rel.err(相对阶跃幅度)")
    for t_end in (60.0, 300.0, 900.0, 1800.0):
        cfg = M.SolverCfg(t_end=t_end, dt=1.0, N=80, props=props_const(rho, cp, kk, 0.0),
                          env_fun=env2, out_every=t_end, tol=1e-12, h=M.Q1["h"], hm=0.0)
        num = float(M.solve_1d(cfg)["T"][-1][0])
        ana = analytic_convective(Bi, 0.0, t_end, Ts2, M.T0) - M.K0
        rel = abs(num - ana) / abs(Ts2 - M.T0)
        print("  %6.0f  %9.4f  %9.4f  %19.4f%%" % (t_end, num, ana, rel * 100))
        ok &= rel < 0.01

    # ---------------- C 场同构 Dirichlet 检验 ----------------
    print("\n[C]  常 D + Dirichlet 边界 vs 级数解（中心含水率）")
    Dc = 5e-9
    envC = lambda t: (np.full(np.shape(t), M.T0), np.zeros(np.shape(t)))
    lamC = jn_zeros(0, 40)
    for t_end in (600.0, 1800.0):
        cfg = M.SolverCfg(t_end=t_end, dt=1.0, N=80, props=props_const(rho, cp, kk, Dc),
                          env_fun=envC, out_every=t_end, tol=1e-12, h=0.0, hm=1e5)
        num = float(M.solve_1d(cfg)["C"][-1][0])
        Fo = Dc * t_end / R ** 2
        s = sum((2.0 / (L * j1(L))) * np.exp(-(L ** 2) * Fo) for L in lamC)
        ana = M.C0 * s
        rel = abs(num - ana) / abs(ana)
        print("  %6.0f  %9.4f  %9.4f  %9.4f%%" % (t_end, num, ana, rel * 100))
        ok &= rel < 0.01

    # ---------------- T3：常数解保持 ----------------
    cst = M.v_constant_solution(N=40, dt=1.0)
    print("\n[T3] 常数解保持最大误差 = %.3e  (%s)" % (cst, "通过" if cst < 1e-9 else "不通过"))
    ok &= cst < 1e-9

    print("\n" + "=" * 76)
    print("总判定：%s" % ("全部通过" if ok else "存在未通过项"))
    print("=" * 76)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
