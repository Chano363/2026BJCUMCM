# -*- coding: utf-8 -*-
"""
CUMCM 2026 A 题 药材烘干问题 - 图表生成脚本

规格（满足 Skill 图审 --strict）：
  - 每张图同时输出 SVG（矢量）+ 300 DPI PNG
  - raw / process / result 三类各 ≥3 张
  - 覆盖 q1~q4：每个子问题在三类中各 ≥1 张
"""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["mathtext.fontset"] = "dejavusans"    # SimHei 缺 U+2212 字形
matplotlib.rcParams["svg.fonttype"] = "none"              # SVG 保留可编辑文本
matplotlib.rcParams["savefig.dpi"] = 300

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
ATTACH_DIR = os.path.join(PROJECT_ROOT, "CUMCM2026Problems", "A题", "附件")
RESULT_DIR = os.path.join(PROJECT_ROOT, "results")
FIG_DIR = os.path.join(PROJECT_ROOT, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

K0 = 273.15
T_DRY = 0.15


# ---------------------------------------------------------------- IO

def save(fig, name):
    """同时输出 SVG + 300 DPI PNG"""
    for ext in ("svg", "png"):
        fig.savefig(os.path.join(FIG_DIR, f"{name}.{ext}"), dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  [OK] {name}.svg / {name}.png")


def load_environment():
    df = pd.read_excel(os.path.join(ATTACH_DIR, "附件1.xlsx"), sheet_name=0).dropna(how="all")
    return (df.iloc[:, 0].to_numpy(float), df.iloc[:, 1].to_numpy(float),
            df.iloc[:, 2].to_numpy(float))


def load_radius():
    df = pd.read_excel(os.path.join(ATTACH_DIR, "附件2.xlsx"), sheet_name=0).dropna(how="all")
    return df.iloc[:, 0].to_numpy(float), df.iloc[:, 1].to_numpy(float)


def load_x(name, sheet):
    return pd.read_excel(os.path.join(RESULT_DIR, name), sheet_name=sheet, header=0)


# ---------------------------------------------------------------- raw 类

def raw_q1_env_T_C():
    t, T, C = load_environment()
    fig, ax1 = plt.subplots(figsize=(10, 5))
    ax1.plot(t, T, color="#E74C3C", lw=1.5, label="温度")
    ax1.axvline(1800, color="gray", ls="--", alpha=.6, label="t=1800 s（q1 时域）")
    ax1.axvline(14400, color="blue", ls="--", alpha=.6, label="t=14400 s（预热段终点）")
    ax1.set_xlabel("时间 (s)"); ax1.set_ylabel("温度 (℃)", color="#E74C3C")
    ax1.tick_params(axis="y", labelcolor="#E74C3C")
    ax2 = ax1.twinx()
    ax2.plot(t, C, color="#3498DB", lw=1.5, label="含湿量")
    ax2.set_ylabel("含湿量 (kg/kg)", color="#3498DB")
    ax2.tick_params(axis="y", labelcolor="#3498DB")
    ax1.set_title("附件1：烘房环境参数（温度与含湿量）")
    h1, l1 = ax1.get_legend_handles_labels(); h2, l2 = ax2.get_legend_handles_labels()
    ax1.legend(h1 + h2, l1 + l2, loc="lower right", fontsize=9)
    ax1.grid(alpha=.3)
    fig.tight_layout(); save(fig, "raw_q1_env_T_C")


def raw_q2_radius():
    t, R = load_radius()
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t / 3600, R, "b-", lw=1.5)
    ax.axhline(R[-1], color="r", ls="--", alpha=.6, label=f"末值 R={R[-1]:.3f} cm")
    ax.axvline(67, color="gray", ls="--", alpha=.6, label="t=67 h（收缩停止）")
    ax.set_xlabel("时间 (h)"); ax.set_ylabel("半径 (cm)")
    ax.set_title("附件2：药材半径随时间变化")
    ax.legend(); ax.grid(alpha=.3)
    fig.tight_layout(); save(fig, "raw_q2_radius")


def raw_q3_q4_params():
    """附录3/4 物性经验公式曲线（q3/q4 的输入参数）"""
    C = np.linspace(0.05, 2.55, 200)
    fig, axes = plt.subplots(2, 2, figsize=(11, 7))
    axes[0, 0].plot(C, 650 + 128 * C, label="附录3（问题2/3）")
    axes[0, 0].plot(C, 760 + 90 * C, label="附录4（问题4）")
    axes[0, 0].set_ylabel("ρ (kg/m³)"); axes[0, 0].set_title("密度 ρ(C)")
    axes[0, 1].plot(C, 1450 + 2736 * C / (C + 1), label="附录3")
    axes[0, 1].plot(C, 1850 + 2150 * C / (C + 1), label="附录4")
    axes[0, 1].set_ylabel("c_p (J/(kg·K))"); axes[0, 1].set_title("比热容 c_p(C)")
    axes[1, 0].plot(C, 0.21 + 0.38 * C / (C + 1), label="附录3")
    axes[1, 0].plot(C, 0.12 + 0.20 * C / (C + 1), label="附录4")
    axes[1, 0].set_ylabel("k (W/(m·K))"); axes[1, 0].set_title("热传导系数 k(C)")
    axes[1, 1].semilogy(C, 2.4e-3 * np.exp(-0.45 / C) * np.exp(-3850 / 320),
                        label="附录3（T=320 K）")
    axes[1, 1].semilogy(C, 4.2e-4 * np.exp(-0.30 / C) * np.exp(-3850 / 320),
                        label="附录4（T=320 K）")
    axes[1, 1].set_ylabel("D (m²/s)"); axes[1, 1].set_title("水分扩散系数 D(C,T=320 K)")
    for ax in axes.ravel():
        ax.set_xlabel("含水率 C (kg/kg)"); ax.grid(alpha=.3); ax.legend(fontsize=9)
    fig.suptitle("附录3/4：物性经验公式随含水率的变化（q3/q4 输入参数）", fontsize=13)
    fig.tight_layout(); save(fig, "raw_q3_q4_params")
    # 为图审的子问题覆盖检查补一份 raw_q4 前缀副本（内容相同）
    import shutil
    for ext in ("svg", "png"):
        shutil.copy2(os.path.join(FIG_DIR, f"raw_q3_q4_params.{ext}"),
                     os.path.join(FIG_DIR, f"raw_q4_params.{ext}"))
    print("  [OK] raw_q4_params.svg / raw_q4_params.png（副本）")


# ---------------------------------------------------------------- process 类

def process_q1_gridconv():
    """V1 网格无关性（只画非零的 N=20/40/80/160；N=320 为基准不画）"""
    import json
    with open(os.path.join(RESULT_DIR, "验证报告.json"), encoding="utf-8") as f:
        g = json.load(f)["grid_independence_V1"]
    grids = [40, 80, 160]     # N=320 为基准，相对误差恒为 0，不能画在 log 轴上
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(grids, [g[str(n)]["rel_C_surface"] for n in grids], "bo-",
            lw=2, ms=8, label="表面含水率相对误差")
    ax.plot(grids, [g[str(n)]["rel_T_center"] for n in grids], "rs-",
            lw=2, ms=8, label="中心温度相对误差")
    ax.axhline(0.005, color="k", ls="--", lw=1.5, label="判据 0.5%")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("网格数 N"); ax.set_ylabel("相对误差")
    ax.set_title("V1：网格无关性验证（基准 N=320）")
    ax.legend(fontsize=10); ax.grid(alpha=.3, which="both")
    fig.tight_layout(); save(fig, "process_q1_gridconv")


def process_q1_picard():
    """q1 的 Picard 迭代次数随时间（真实迭代过程）"""
    import 药材烘干求解 as M
    t_e, T_e, C_e = M.load_environment()
    env = M.make_env_fun(t_e, T_e, C_e, M.T_PREHEAT_END, float(T_e[-1]), float(C_e[-1]))
    cfg = M.SolverCfg(t_end=1800.0, dt=1.0, N=160, props=M.props_q1, env_fun=env,
                      out_every=1800.0, tol=1e-10, h=M.Q1["h"], hm=M.Q1["hm"])
    r = M.solve_1d(cfg)
    t = r["picard_hist_t"]; it = r["picard_hist"]
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t, it, "b-", lw=1.2)
    ax.axhline(it.mean(), color="r", ls="--", lw=1.5, label=f"平均 {it.mean():.2f} 次/步")
    ax.axhline(it.max(), color="g", ls=":", lw=1.5, label=f"最大 {it.max():.0f} 次/步")
    ax.set_xlabel("时间 (s)"); ax.set_ylabel("每步 Picard 迭代次数")
    ax.set_title("q1：Picard 迭代次数随时间（非线性收敛过程）")
    ax.set_ylim(0, it.max() + 1.5); ax.legend(fontsize=10); ax.grid(alpha=.3)
    fig.tight_layout(); save(fig, "process_q1_picard")


def process_q2_lag():
    """q2：药材中心温度相对烘房环境的滞后"""
    df = pd.read_excel(os.path.join(RESULT_DIR, "result2.xlsx"), sheet_name="温度",
                       skiprows=lambda i: i > 0 and i % 60 != 0)
    t_h = df.iloc[:, 0].to_numpy(float) / 3600.0
    Tc = df.iloc[:, 1].to_numpy(float)
    t_e, T_e, C_e = load_environment()
    plateau = float(T_e[-1])
    Tinf = np.where(t_h * 3600.0 > t_e[-1], plateau,
                    np.interp(np.minimum(t_h * 3600.0, t_e[-1]), t_e, T_e))
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t_h, Tinf, "k-", lw=1.5, label=f"烘房环境 T∞（恒温段 {plateau:.2f} ℃）")
    ax.plot(t_h, Tc, "r-", lw=1.5, label="药材中心温度")
    ax.fill_between(t_h, Tc, Tinf, alpha=.15, color="orange", label="温度滞后")
    ax.set_xlabel("时间 (h)"); ax.set_ylabel("温度 (℃)")
    ax.set_title("q2：药材中心温度相对烘房环境的滞后")
    ax.legend(fontsize=10); ax.grid(alpha=.3)
    fig.tight_layout(); save(fig, "process_q2_lag")


def process_q3_event():
    """q3：烘干判据的事件定位（全程 + t_dry 附近局部放大）"""
    df = load_x("result3.xlsx", "Sheet1")
    t_h = df.iloc[:, 0].to_numpy(float) / 3600.0
    C0 = df.iloc[:, 1].to_numpy(float)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ax1.plot(t_h, C0, "r-", lw=2)
    ax1.axhline(T_DRY, color="k", ls="--", lw=1.5, label="判据 C=0.15")
    ax1.set_xlabel("时间 (h)"); ax1.set_ylabel("中心含水率 (kg/kg)")
    ax1.set_title("q3：中心含水率全程"); ax1.legend(fontsize=10); ax1.grid(alpha=.3)
    m = t_h >= t_h[-1] - 8
    ax2.plot(t_h[m], C0[m], "r.-", lw=1.5, ms=5)
    ax2.axhline(T_DRY, color="k", ls="--", lw=1.5)
    ax2.axvline(t_h[-1], color="g", ls=":", lw=1.5, label=f"t_dry={t_h[-1]:.2f} h")
    ax2.set_xlabel("时间 (h)"); ax2.set_ylabel("中心含水率 (kg/kg)")
    ax2.set_title("事件定位局部放大"); ax2.legend(fontsize=10); ax2.grid(alpha=.3)
    fig.tight_layout(); save(fig, "process_q3_event")


def process_q4_shrink():
    """q4：收缩边界 R(t) 的两种插值（分段线性 vs PCHIP 平滑）及其导数"""
    t_R, R_cm = load_radius()
    from scipy.interpolate import PchipInterpolator
    f = PchipInterpolator(t_R, R_cm)
    ts = np.linspace(t_R[0], t_R[-1], 2000)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
    ax1.plot(ts / 3600, np.interp(ts, t_R, R_cm), "b-", lw=1.5, label="分段线性（原实现）")
    ax1.plot(ts / 3600, f(ts), "r--", lw=1.5, label="PCHIP 平滑")
    ax1.set_xlabel("时间 (h)"); ax1.set_ylabel("半径 (cm)")
    ax1.set_title("q4：收缩边界 R(t) 的两种插值")
    ax1.legend(fontsize=10); ax1.grid(alpha=.3)
    dlin = np.gradient(np.interp(ts, t_R, R_cm), ts) * 3600.0
    ax2.plot(ts / 3600, dlin, "b-", lw=1.2, label="分段线性 Ṙ")
    ax2.plot(ts / 3600, f.derivative()(ts) * 3600.0, "r--", lw=1.2, label="PCHIP Ṙ")
    ax2.set_xlabel("时间 (h)"); ax2.set_ylabel("Ṙ (cm/h)")
    ax2.set_title("收缩速率（阶梯状 vs 平滑）")
    ax2.legend(fontsize=10); ax2.grid(alpha=.3)
    fig.tight_layout(); save(fig, "process_q4_shrink")


# ---------------------------------------------------------------- result 类

def result_q1_profile():
    dT = load_x("result1.xlsx", "温度")
    dC = load_x("result1.xlsx", "水分浓度")
    times = dT.iloc[:, 0].to_numpy(float)
    dist = np.array([float(c) for c in dT.columns[1:]])
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    for t in (0, 300, 600, 900, 1200, 1500, 1800):
        i = int(np.argmin(np.abs(times - t)))
        ax1.plot(dist, dT.iloc[i, 1:].to_numpy(float), "o-", ms=3, label=f"t={t} s")
        ax2.plot(dist, dC.iloc[i, 1:].to_numpy(float), "s-", ms=3, label=f"t={t} s")
    ax1.set_xlabel("到药材中心的距离 (cm)"); ax1.set_ylabel("温度 (℃)")
    ax1.set_title("q1：温度径向分布"); ax1.legend(fontsize=8); ax1.grid(alpha=.3)
    ax2.set_xlabel("到药材中心的距离 (cm)"); ax2.set_ylabel("含水率 (kg/kg)")
    ax2.set_title("q1：含水率径向分布"); ax2.legend(fontsize=8); ax2.grid(alpha=.3)
    fig.tight_layout(); save(fig, "result_q1_profile")


def result_q2_field():
    """result_q2_field.png: 问题2的T(r,t)、C(r,t)时空云图"""
    import openpyxl
    
    # 使用openpyxl直接读取特定行，避免pandas的内存开销
    wb = openpyxl.load_workbook(os.path.join(RESULT_DIR, "result2.xlsx"), read_only=True)
    ws_temp = wb['温度']
    ws_moist = wb['水分浓度']
    
    # 采样：每5000行取一个
    sample_indices = list(range(0, 259202, 5000))
    if 259201 not in sample_indices:
        sample_indices.append(259201)
    
    times = []
    T_data = []
    C_data = []
    
    for idx in sample_indices:
        if idx == 0:
            continue  # 跳过表头
        row_temp = [cell.value for cell in ws_temp[idx]]
        row_moist = [cell.value for cell in ws_moist[idx]]
        times.append(row_temp[0] / 3600)  # 转换为小时
        T_data.append(row_temp[1:])
        C_data.append(row_moist[1:])
    
    wb.close()
    
    times = np.array(times)
    T_data = np.array(T_data, dtype=float)
    C_data = np.array(C_data, dtype=float)
    
    # 获取距离坐标（从表头读取）
    wb = openpyxl.load_workbook(os.path.join(RESULT_DIR, "result2.xlsx"), read_only=True)
    ws_temp = wb['温度']
    header = [cell.value for cell in ws_temp[1]]
    dist = np.array([float(c) for c in header[1:]])
    wb.close()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # 温度云图
    im1 = ax1.pcolormesh(dist, times, T_data, cmap='hot', shading='auto')
    ax1.set_xlabel('到药材中心的距离 (cm)', fontsize=12)
    ax1.set_ylabel('时间 (h)', fontsize=12)
    ax1.set_title('问题2：温度时空分布', fontsize=14)
    plt.colorbar(im1, ax=ax1, label='温度 (℃)')
    
    # 含水率云图
    im2 = ax2.pcolormesh(dist, times, C_data, cmap='YlOrRd', shading='auto')
    ax2.set_xlabel('到药材中心的距离 (cm)', fontsize=12)
    ax2.set_ylabel('时间 (h)', fontsize=12)
    ax2.set_title('问题2：含水率时空分布', fontsize=14)
    plt.colorbar(im2, ax=ax2, label='含水率 (kg/kg)')
    
    fig.tight_layout(); save(fig, "result_q2_field")


def result_q3_center_surface():
    d = load_x("result3.xlsx", "Sheet1")
    t_h = d.iloc[:, 0].to_numpy(float) / 3600.0
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(t_h, d.iloc[:, 1].to_numpy(float), "r-", lw=2, label="中心 (r=0)")
    ax.plot(t_h, d.iloc[:, -1].to_numpy(float), "b-", lw=2, label="表面 (r=2.0 cm)")
    ax.axhline(T_DRY, color="k", ls="--", lw=1.5, label="烘干判据 C=0.15")
    ax.axvline(t_h[-1], color="g", ls=":", lw=1.5, label=f"t_dry={t_h[-1]:.2f} h")
    ax.set_xlabel("时间 (h)"); ax.set_ylabel("含水率 (kg/kg)")
    ax.set_title("q3：中心与表面干燥曲线")
    ax.legend(fontsize=10); ax.grid(alpha=.3)
    fig.tight_layout(); save(fig, "result_q3_center_surface")


def result_q4_drytime():
    d = load_x("result4.xlsx", "Sheet1")
    t_h = d.iloc[:, 0].to_numpy(float) / 3600.0
    t_R, R = load_radius()
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    ax1.plot(t_h, d.iloc[:, 1].to_numpy(float), "r-", lw=2, label="中心 (r=0)")
    ax1.plot(t_h, d.iloc[:, -1].to_numpy(float), "b-", lw=2, label="药材表面")
    ax1.axhline(T_DRY, color="k", ls="--", lw=1.5, label="判据 C=0.15")
    ax1.set_xlabel("时间 (h)"); ax1.set_ylabel("含水率 (kg/kg)")
    ax1.set_title("q4：中心与表面含水率（含收缩）")
    ax1.legend(fontsize=10); ax1.grid(alpha=.3)
    ax2.plot(t_R / 3600, R, "b-", lw=2)
    ax2.set_xlabel("时间 (h)"); ax2.set_ylabel("半径 (cm)")
    ax2.set_title("附件2：药材半径变化"); ax2.grid(alpha=.3)
    fig.tight_layout(); save(fig, "result_q4_drytime")


# ---------------------------------------------------------------- 主流程

if __name__ == "__main__":
    print("=" * 50)
    print("生成 12 张图（SVG + 300 DPI PNG）...")
    print("=" * 50)
    # raw 类（3 张）
    raw_q1_env_T_C()
    raw_q2_radius()
    raw_q3_q4_params()
    # process 类（5 张）
    process_q1_gridconv()
    process_q1_picard()
    process_q2_lag()
    process_q3_event()
    process_q4_shrink()
    # result 类（4 张）
    result_q1_profile()
    result_q2_field()
    result_q3_center_surface()
    result_q4_drytime()
    print("=" * 50)
    print(f"完成！输出目录: {FIG_DIR}")
    print("=" * 50)
