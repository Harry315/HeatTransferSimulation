"""
NT11 收敛时间自动计算器
author : you
"""

import numpy as np
from scipy.interpolate import Akima1DInterpolator
from scipy.ndimage import maximum_filter1d
import matplotlib.pyplot as plt

# -------------------------------------------------
# 1. 原始数据（直接复制你的两组 array 即可）
t_raw = np.array([
    0.0, 0.000692502304445952, 0.00100000004749745, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0,
    7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0,
    21.0, 22.0, 23.0, 24.0, 25.0, 26.0, 27.0, 28.0, 29.0, 30.0, 31.0, 32.0, 33.0, 34.0,
    35.0, 36.0, 37.0, 38.0, 39.0, 40.0, 41.0, 42.0, 43.0, 44.0, 45.0, 46.0, 47.0, 48.0,
    49.0, 50.0
])
y_raw = np.array([
    25.0, 25.0, 25.0, 25.0, 38.1890128742565, 48.8109729940241, 55.168213150718,
    58.7143534747037, 60.659498388117, 61.722357749939, 62.3026236620816, 62.6193591031161,
    62.7922397960316, 62.8866013613614, 62.9381050630049, 62.9662168676203, 62.9815604469993,
    62.9899355281483, 62.9945065758445, 62.9970018213445, 62.9983634081754, 62.9991070140492,
    62.999511978843, 62.9997329711914, 62.9998550415039, 62.9999198913574, 62.9999580383301,
    62.9999771118164, 62.9999885559082, 62.9999923706055, 62.9999961853027, 62.9999961853027,
    63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0,
    63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0, 63.0
])
# -------------------------------------------------

def convergence_time(t, y, tol=1e-4, window=1.0, dt_fine=1e-3):
    """
    计算 NT11 曲线收敛时间（可非整数）
    内部已处理：去重、排序、单调递增
    """
    # 1. 去重 + 排序
    t, uniq_idx = np.unique(t, return_index=True)
    y = y[uniq_idx]

    # 2. 稳态值（最后 10 % 均值）
    n_tail = max(1, int(0.1 * y.size))
    y_ss = y[-n_tail:].mean()

    # 3. 连续插值
    t_fine = np.arange(t[0], t[-1] + dt_fine, dt_fine)
    akima = Akima1DInterpolator(t, y)
    y_fine = akima(t_fine)

    # 4. 滑动窗口：往后看 window/2 秒，往前看 window/2 秒
    half_w = int(np.ceil(window / 2 / dt_fine))
    dev = np.abs(y_fine - y_ss)
    max_dev = maximum_filter1d(dev, size=2 * half_w + 1, mode='nearest')

    # 5. 找第一个满足“全程 <= tol”的时刻
    converged = max_dev <= tol
    if not converged.any():
        return np.nan
    t_conv = t_fine[np.where(converged)[0][0]]
    return float(t_conv)


# -------------------------------------------------
# 主程序
# -------------------------------------------------
if __name__ == "__main__":
    # 同步去重／排序
    t, uniq_idx = np.unique(t_raw, return_index=True)
    y = y_raw[uniq_idx]

    t_conv = convergence_time(t, y, tol=1e-4, window=1.0)
    print("收敛时间 = {:.3f} s".format(t_conv))

    # 可视化
    plt.figure(figsize=(6, 4))
    plt.plot(t, y, 'o', markersize=3, label='original')
    plt.axhline(y=y[-5:].mean(), ls='--', color='gray', label='steady state')
    if not np.isnan(t_conv):
        plt.axvline(t_conv, color='r', label=f'convergence @ {t_conv:.3f} s')
    plt.xlabel('Time (s)')
    plt.ylabel('Avg NT11')
    plt.title('NT11 Convergence Analysis')
    plt.grid(True)
    plt.legend()
    plt.tight_layout()
    plt.show()