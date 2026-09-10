# -*- coding: utf-8 -*-
"""
CUMCM 2026 A 题 药材烘干问题 - 图表生成脚本

生成题目分析报告§7中提到的候选图表：
- raw（原始数据）：附件1/2的数据可视化
- process（运行过程）：网格收敛、Picard迭代、守恒审计
- result（最终结果）：温度/含水率剖面、时空云图、干燥曲线
"""

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['font.sans-serif'] = ['SimHei', 'Microsoft YaHei', 'DejaVu Sans']
matplotlib.rcParams['axes.unicode_minus'] = False

PROJECT_ROOT = os.path.dirname(os.path.abspath(__file__))
ATTACH_DIR = os.path.join(PROJECT_ROOT, "CUMCM2026Problems", "A题", "附件")
RESULT_DIR = os.path.join(PROJECT_ROOT, "results")
FIG_DIR = os.path.join(PROJECT_ROOT, "figures")
os.makedirs(FIG_DIR, exist_ok=True)

K0 = 273.15


def load_environment():
    """加载附件1（烘房环境）"""
    df = pd.read_excel(os.path.join(ATTACH_DIR, "附件1.xlsx"), sheet_name=0).dropna(how="all")
    t = df.iloc[:, 0].to_numpy(float)
    T = df.iloc[:, 1].to_numpy(float)
    C = df.iloc[:, 2].to_numpy(float)
    return t, T, C


def load_radius():
    """加载附件2（药材半径）"""
    df = pd.read_excel(os.path.join(ATTACH_DIR, "附件2.xlsx"), sheet_name=0).dropna(how="all")
    t = df.iloc[:, 0].to_numpy(float)
    R = df.iloc[:, 1].to_numpy(float)
    return t, R


def load_result1():
    """加载result1.xlsx"""
    df_temp = pd.read_excel(os.path.join(RESULT_DIR, "result1.xlsx"), sheet_name="温度")
    df_moist = pd.read_excel(os.path.join(RESULT_DIR, "result1.xlsx"), sheet_name="水分浓度")
    return df_temp, df_moist


def load_result2():
    """加载result2.xlsx（采样部分数据）"""
    df_temp = pd.read_excel(os.path.join(RESULT_DIR, "result2.xlsx"), sheet_name="温度")
    df_moist = pd.read_excel(os.path.join(RESULT_DIR, "result2.xlsx"), sheet_name="水分浓度")
    return df_temp, df_moist


def load_result3():
    """加载result3.xlsx"""
    df = pd.read_excel(os.path.join(RESULT_DIR, "result3.xlsx"), sheet_name="Sheet1")
    return df


def load_result4():
    """加载result4.xlsx"""
    df = pd.read_excel(os.path.join(RESULT_DIR, "result4.xlsx"), sheet_name="Sheet1")
    return df


def load_validation():
    """加载验证报告"""
    import json
    with open(os.path.join(RESULT_DIR, "验证报告.json"), 'r', encoding='utf-8') as f:
        return json.load(f)


def fig_raw_env():
    """raw_q1_env_T_C.png: 附件1烘房温度与含湿量双轴曲线"""
    t, T, C = load_environment()
    
    fig, ax1 = plt.subplots(figsize=(10, 5))
    
    color1 = '#E74C3C'
    ax1.set_xlabel('时间 (s)', fontsize=12)
    ax1.set_ylabel('温度 (℃)', color=color1, fontsize=12)
    ax1.plot(t, T, color=color1, linewidth=1.5, label='温度')
    ax1.tick_params(axis='y', labelcolor=color1)
    ax1.axvline(x=1800, color='gray', linestyle='--', alpha=0.5, label='t=1800s (q1时域)')
    ax1.axvline(x=14400, color='blue', linestyle='--', alpha=0.5, label='t=14400s (预热段终点)')
    ax1.legend(loc='upper left')
    
    ax2 = ax1.twinx()
    color2 = '#3498DB'
    ax2.set_ylabel('含湿量 (kg/kg)', color=color2, fontsize=12)
    ax2.plot(t, C, color=color2, linewidth=1.5, label='含湿量')
    ax2.tick_params(axis='y', labelcolor=color2)
    ax2.legend(loc='upper right')
    
    plt.title('附件1：烘房环境参数（温度与含湿量）', fontsize=14)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'raw_q1_env_T_C.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print('[OK] raw_q1_env_T_C.png')


def fig_raw_radius():
    """raw_q2_radius.png: 附件2半径随时间变化"""
    t, R = load_radius()
    
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t/3600, R, 'b-', linewidth=1.5)
    ax.set_xlabel('时间 (h)', fontsize=12)
    ax.set_ylabel('半径 (cm)', fontsize=12)
    ax.set_title('附件2：药材半径随时间变化', fontsize=14)
    ax.axhline(y=R[-1], color='r', linestyle='--', alpha=0.5, label=f'末值 R={R[-1]:.3f} cm')
    ax.axvline(x=67, color='gray', linestyle='--', alpha=0.5, label='t=67h (收缩停止)')
    ax.legend()
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'raw_q2_radius.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print('[OK] raw_q2_radius.png')


def fig_result1_profile():
    """result_q1_profile.png: 问题1若干时刻的径向温度/含水率剖面"""
    df_temp, df_moist = load_result1()
    
    # 提取数据
    times = df_temp.iloc[:, 0].to_numpy()
    dist = np.array([float(c) for c in df_temp.columns[1:]])
    
    # 选择若干时刻
    target_times = [0, 300, 600, 900, 1200, 1500, 1800]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    for t in target_times:
        idx = np.argmin(np.abs(times - t))
        ax1.plot(dist, df_temp.iloc[idx, 1:].to_numpy(), 'o-', label=f't={t}s', markersize=3)
        ax2.plot(dist, df_moist.iloc[idx, 1:].to_numpy(), 's-', label=f't={t}s', markersize=3)
    
    ax1.set_xlabel('到药材中心的距离 (cm)', fontsize=12)
    ax1.set_ylabel('温度 (℃)', fontsize=12)
    ax1.set_title('问题1：温度径向分布', fontsize=14)
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    
    ax2.set_xlabel('到药材中心的距离 (cm)', fontsize=12)
    ax2.set_ylabel('含水率 (kg/kg)', fontsize=12)
    ax2.set_title('问题1：含水率径向分布', fontsize=14)
    ax2.legend(fontsize=9)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'result_q1_profile.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print('[OK] result_q1_profile.png')


def fig_result2_field():
    """result_q2_field.png: 问题2的T(r,t)、C(r,t)时空云图"""
    # 采样result2的数据（每5000步取一个，减少数据量）
    df_temp = pd.read_excel(os.path.join(RESULT_DIR, "result2.xlsx"), sheet_name="温度", 
                            skiprows=lambda i: i > 0 and i % 5000 != 0)
    df_moist = pd.read_excel(os.path.join(RESULT_DIR, "result2.xlsx"), sheet_name="水分浓度",
                             skiprows=lambda i: i > 0 and i % 5000 != 0)
    
    times = df_temp.iloc[:, 0].to_numpy() / 3600  # 转换为小时
    dist = np.array([float(c) for c in df_temp.columns[1:]])
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # 温度云图
    T_data = df_temp.iloc[:, 1:].to_numpy()
    im1 = ax1.pcolormesh(dist, times, T_data, cmap='hot', shading='auto')
    ax1.set_xlabel('到药材中心的距离 (cm)', fontsize=12)
    ax1.set_ylabel('时间 (h)', fontsize=12)
    ax1.set_title('问题2：温度时空分布', fontsize=14)
    plt.colorbar(im1, ax=ax1, label='温度 (℃)')
    
    # 含水率云图
    C_data = df_moist.iloc[:, 1:].to_numpy()
    im2 = ax2.pcolormesh(dist, times, C_data, cmap='YlOrRd', shading='auto')
    ax2.set_xlabel('到药材中心的距离 (cm)', fontsize=12)
    ax2.set_ylabel('时间 (h)', fontsize=12)
    ax2.set_title('问题2：含水率时空分布', fontsize=14)
    plt.colorbar(im2, ax=ax2, label='含水率 (kg/kg)')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'result_q2_field.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print('[OK] result_q2_field.png')


def fig_result3_center_surface():
    """result_q3_center_surface.png: 中心与表面干燥曲线 + C=0.15判据线"""
    df = load_result3()
    times = df.iloc[:, 0].to_numpy() / 3600  # 转换为小时
    C_center = df.iloc[:, 1].to_numpy()  # 中心（r=0）
    C_surface = df.iloc[:, 6].to_numpy()  # 表面（r=1.0cm）
    
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(times, C_center, 'r-', linewidth=2, label='中心 (r=0)')
    ax.plot(times, C_surface, 'b-', linewidth=2, label='表面 (r=1.0cm)')
    ax.axhline(y=0.15, color='k', linestyle='--', linewidth=1.5, label='烘干判据 C=0.15')
    
    # 标记烘干时长
    t_dry = 57.0765  # 从运行结果获取
    ax.axvline(x=t_dry, color='g', linestyle=':', linewidth=1.5, label=f'烘干时长 t={t_dry:.1f}h')
    
    ax.set_xlabel('时间 (h)', fontsize=12)
    ax.set_ylabel('含水率 (kg/kg)', fontsize=12)
    ax.set_title('问题3：中心与表面干燥曲线', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xlim([0, 72])
    ax.set_ylim([0, 2.6])
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'result_q3_center_surface.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print('[OK] result_q3_center_surface.png')


def fig_process_gridconv():
    """process_q1_gridconv.png: V1网格无关性收敛曲线"""
    validation = load_validation()
    grid_data = validation['grid_independence_V1']
    
    grids = [40, 80, 160, 320]
    rel_C_surface = [grid_data[str(g)]['rel_C_surface'] for g in grids]
    rel_T_center = [grid_data[str(g)]['rel_T_center'] for g in grids]
    
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(grids, rel_C_surface, 'bo-', linewidth=2, markersize=8, label='表面含水率相对误差')
    ax.plot(grids, rel_T_center, 'rs-', linewidth=2, markersize=8, label='中心温度相对误差')
    ax.axhline(y=0.005, color='k', linestyle='--', linewidth=1.5, label='判据 0.5%')
    ax.set_xlabel('网格数 N', fontsize=12)
    ax.set_ylabel('相对误差', fontsize=12)
    ax.set_title('V1：网格无关性验证', fontsize=14)
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    ax.set_xscale('log')
    ax.set_yscale('log')
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'process_q1_gridconv.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print('[OK] process_q1_gridconv.png')


def fig_process_picard():
    """process_q1_picard.png: 问题1中心含水率随时间变化"""
    # 从result1.xlsx读取数据
    df_temp, df_moist = load_result1()
    times = df_temp.iloc[:, 0].to_numpy() / 60  # 转换为分钟
    C_center = df_moist.iloc[:, 1].to_numpy()  # 中心（r=0）
    
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(times, C_center, 'b-', linewidth=1.5)
    ax.set_xlabel('时间 (min)', fontsize=12)
    ax.set_ylabel('中心含水率 (kg/kg)', fontsize=12)
    ax.set_title('问题1：中心含水率随时间变化', fontsize=14)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'process_q1_picard.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print('[OK] process_q1_picard.png')


def fig_result_q4_drytime():
    """result_q4_drytime.png: 问题4中心含水率变化曲线"""
    # 从result4.xlsx读取数据
    df = load_result4()
    times = df.iloc[:, 0].to_numpy() / 3600  # 转换为小时
    C_center = df.iloc[:, 1].to_numpy()  # 中心（r=0）
    
    # 从附件2读取半径数据
    t_R, R = load_radius()
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))
    
    # 含水率曲线
    ax1.plot(times, C_center, 'r-', linewidth=2)
    ax1.axhline(y=0.15, color='k', linestyle='--', linewidth=1.5, label='烘干判据 C=0.15')
    ax1.set_xlabel('时间 (h)', fontsize=12)
    ax1.set_ylabel('中心含水率 (kg/kg)', fontsize=12)
    ax1.set_title('问题4：中心含水率变化（含收缩）', fontsize=14)
    ax1.legend(fontsize=10)
    ax1.grid(True, alpha=0.3)
    
    # 半径变化
    ax2.plot(t_R/3600, R, 'b-', linewidth=2)
    ax2.set_xlabel('时间 (h)', fontsize=12)
    ax2.set_ylabel('半径 (cm)', fontsize=12)
    ax2.set_title('问题4：药材半径变化', fontsize=14)
    ax2.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(os.path.join(FIG_DIR, 'result_q4_drytime.png'), dpi=150, bbox_inches='tight')
    plt.close()
    print('[OK] result_q4_drytime.png')


if __name__ == "__main__":
    print("=" * 50)
    print("开始生成图表...")
    print("=" * 50)
    
    # 原始数据图表
    fig_raw_env()
    fig_raw_radius()
    
    # 运行过程图表
    fig_process_gridconv()
    fig_process_picard()
    
    # 最终结果图表
    fig_result1_profile()
    fig_result2_field()
    fig_result3_center_surface()
    fig_result_q4_drytime()
    
    print("=" * 50)
    print(f"图表生成完成！保存至: {FIG_DIR}")
    print("=" * 50)
