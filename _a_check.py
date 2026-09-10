# -*- coding: utf-8 -*-
"""A 题附件：关键时点精确定位（只读）"""
import numpy as np
import pandas as pd

BASE = r"d:\Code\2026BJCUMCM\CUMCM2026Problems\A题\附件"
env = pd.read_excel(f"{BASE}/附件1.xlsx")
rad = pd.read_excel(f"{BASE}/附件2.xlsx")

t = env.iloc[:, 0].to_numpy(float)
T = env.iloc[:, 1].to_numpy(float)
C = env.iloc[:, 2].to_numpy(float)

print("T 峰值 %.4f @ t=%.0f s (%.2f h)" % (T.max(), t[T.argmax()], t[T.argmax()] / 3600))
print("C_air 峰值 %.5f @ t=%.0f s (%.2f h)" % (C.max(), t[C.argmax()], t[C.argmax()] / 3600))
print("T 末端 3 点:", [(int(tt), round(cc, 4)) for tt, cc in zip(t[-3:], T[-3:])])
print("C 末端 3 点:", [(int(tt), round(cc, 5)) for tt, cc in zip(t[-3:], C[-3:])])
print("T 非单调的处数 = %d / %d" % (int((np.diff(T) < 0).sum()), len(T) - 1))
print("C 非单调的处数 = %d / %d" % (int((np.diff(C) < 0).sum()), len(C) - 1))

tr = rad.iloc[:, 0].to_numpy(float)
R = rad.iloc[:, 1].to_numpy(float)
dR = np.diff(R)
nz = np.nonzero(dR)[0]
last = nz[-1]
print()
print("半径: 首次为 0 的增量下标 = %d, 最后一个非零增量下标 = %d" % (int(np.nonzero(dR == 0)[0][0]), int(last)))
print("半径最后一次变化到 t=%.0f s (%.2f h) 时 R=%.4f" % (tr[last + 1], tr[last + 1] / 3600, R[last + 1]))
print("此后到 t=%.0f s (%.2f h) 半径恒为 %.4f" % (tr[-1], tr[-1] / 3600, R[-1]))
print("R 抽样(每12点):", [round(v, 4) for v in R[::12]])

# 输出规模核算
print()
print("result2 全程 @1 s 行数(含 t=0) = %d" % (int(tr[-1]) + 1))
print("result3/4 全程 @60 s 行数(含 t=0) = %d" % (int(tr[-1] / 60) + 1))
print("result1 @1 s 0..1800 = %d" % (1801,))
