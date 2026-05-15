"""
IMU误差标定与初始对准 -- 可视化画图脚本

使用方法：
  1. 在 C++ 端调用 runCalibrationAndAlignment() 获取结果，
     然后调用 saveCalibDataForPlot(out_dir, input, result) 导出数据
  2. python main/calib_plot.py --data-dir <导出目录>

依赖：matplotlib, numpy
"""

import argparse
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

# 设置中文字体（Windows 用 SimHei 或 Microsoft YaHei，Linux 可能需要其他字体）
plt.rcParams["font.sans-serif"] = ["SimHei", "Microsoft YaHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="绘制IMU标定与初始对准结果图")
    parser.add_argument(
        "--data-dir",
        type=str,
        default="data/calib_out",
        help="saveCalibDataForPlot 导出的数据目录，默认 data/calib_out",
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=140,
        help="输出图片DPI，默认140",
    )
    parser.add_argument(
        "--show",
        dest="show",
        action="store_true",
        help="绘图后弹窗展示（默认开启）",
    )
    parser.add_argument(
        "--no-show",
        dest="show",
        action="store_false",
        help="仅保存图片，不弹窗展示",
    )
    parser.set_defaults(show=True)
    return parser.parse_args()


def load_six_position(path: str) -> dict:
    """读取 calib_six_position.txt，返回六个位置的真值、测量值、残差等。"""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"文件不存在: {path}")

    metadata = {}
    data_lines = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                # 解析 # key=value 或 # columns: ...
                if "=" in line:
                    kv = line.lstrip("#").strip()
                    parts = kv.split("=", 1)
                    if len(parts) == 2:
                        key = parts[0].strip()
                        val = parts[1].strip()
                        try:
                            metadata[key] = float(val)
                        except ValueError:
                            metadata[key] = val
            else:
                data_lines.append(line)

    # columns: pos_label true_fx true_fy true_fz meas_fx meas_fy meas_fz resid_fx resid_fy resid_fz
    labels = []
    true_f = np.zeros((6, 3))
    meas_f = np.zeros((6, 3))
    resid = np.zeros((6, 3))
    for i, line in enumerate(data_lines):
        parts = line.split()
        if len(parts) < 10:
            continue
        labels.append(parts[0])
        true_f[i, 0] = float(parts[1])
        true_f[i, 1] = float(parts[2])
        true_f[i, 2] = float(parts[3])
        meas_f[i, 0] = float(parts[4])
        meas_f[i, 1] = float(parts[5])
        meas_f[i, 2] = float(parts[6])
        resid[i, 0] = float(parts[7])
        resid[i, 1] = float(parts[8])
        resid[i, 2] = float(parts[9])

    return {
        "metadata": metadata,
        "labels": labels,
        "true_f": true_f,
        "meas_f": meas_f,
        "resid": resid,
    }


def load_gyro_static(path: str) -> dict:
    """读取 calib_gyro_static.txt，返回静态陀螺原始数据序列。"""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"文件不存在: {path}")

    metadata = {}
    data = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith("#"):
                if "=" in line:
                    kv = line.lstrip("#").strip()
                    parts = kv.split("=", 1)
                    if len(parts) == 2:
                        key = parts[0].strip()
                        val = parts[1].strip()
                        try:
                            metadata[key] = float(val)
                        except ValueError:
                            metadata[key] = val
            else:
                parts = line.split()
                if len(parts) >= 4:
                    data.append([float(x) for x in parts[:4]])

    data = np.array(data)  # shape (N, 4): index, dtheta_x_deg, dtheta_y_deg, dtheta_z_deg
    dtheta_deg = data[:, 1:4] if data.shape[1] >= 4 else np.zeros((0, 3))
    dt = float(metadata.get("dt", 0.005))
    bias_deg_h = np.array([
        float(metadata.get("gyro_bias_deg_h", "0 0 0").split()[0]),
        float(metadata.get("gyro_bias_deg_h", "0 0 0").split()[1]),
        float(metadata.get("gyro_bias_deg_h", "0 0 0").split()[2]),
    ])
    scale_ppm = np.array([
        float(metadata.get("gyro_scale_ppm", "0 0 0").split()[0]),
        float(metadata.get("gyro_scale_ppm", "0 0 0").split()[1]),
        float(metadata.get("gyro_scale_ppm", "0 0 0").split()[2]),
    ])

    return {
        "metadata": metadata,
        "dt": dt,
        "dtheta_deg": dtheta_deg,
        "bias_deg_h": bias_deg_h,
        "scale_ppm": scale_ppm,
    }


def load_alignment(path: str) -> dict:
    """读取 calib_alignment.txt，返回对准结果。"""
    if not os.path.isfile(path):
        raise FileNotFoundError(f"文件不存在: {path}")

    data_line = None
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            data_line = line
            break

    if data_line is None:
        raise ValueError(f"{path} 中没有找到数据行")

    parts = data_line.split()
    return {
        "roll_deg": float(parts[0]),
        "pitch_deg": float(parts[1]),
        "yaw_deg": float(parts[2]),
        "valid": bool(int(parts[3])),
        "condition_number": float(parts[4]),
    }


def plot_six_position(data: dict, out_dir: str, dpi: int):
    """六位置加速度计标定对比图。"""
    labels = data["labels"]
    true_f = data["true_f"]   # (6, 3)
    meas_f = data["meas_f"]   # (6, 3)
    resid = data["resid"]     # (6, 3)
    metadata = data["metadata"]

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    axis_names = ["X", "Y", "Z"]
    colors_true = ["#2196F3", "#4CAF50", "#FF9800", "#2196F3", "#4CAF50", "#FF9800"]
    colors_meas = ["#90CAF9", "#A5D6A7", "#FFCC80", "#90CAF9", "#A5D6A7", "#FFCC80"]

    x_pos = np.arange(6)
    bar_width = 0.35

    for ax_idx in range(3):
        ax = axes[ax_idx]

        for i in range(6):
            ax.bar(i - bar_width / 2, true_f[i, ax_idx], bar_width,
                   color=colors_true[i], alpha=0.85, edgecolor="white", linewidth=0.5)
            ax.bar(i + bar_width / 2, meas_f[i, ax_idx], bar_width,
                   color=colors_meas[i], alpha=0.85, edgecolor="white", linewidth=0.5,
                   yerr=[[0], [abs(resid[i, ax_idx])]], capsize=3,
                   error_kw={"linewidth": 1.2, "color": "red"})

        ax.set_xticks(x_pos)
        ax.set_xticklabels(labels, fontsize=8)
        ax.set_ylabel("比力 (m/s^2)")
        ax.set_title(f"{axis_names[ax_idx]}轴", fontweight="bold")
        ax.grid(True, alpha=0.25, axis="y")

    # 图例（手动构造）
    legend_elements = [
        Patch(facecolor="#2196F3", alpha=0.85, label="真值 (左)"),
        Patch(facecolor="#90CAF9", alpha=0.85, label="测量值 (右)"),
    ]
    fig.legend(handles=legend_elements, loc="lower center", ncol=2, fontsize=9)

    bias = [metadata.get("bias_x", 0), metadata.get("bias_y", 0), metadata.get("bias_z", 0)]
    M = [[metadata.get(f"M_{r}{c}", 0) for c in range(3)] for r in range(3)]
    # 对角线: 比例因子误差 (ppm)，非对角线: 交轴耦合 (ppm)
    scale_err = [(M[i][i] - 1.0) * 1e6 for i in range(3)]
    cross_axis_vals = [abs(M[r][c]) * 1e6 for r in range(3) for c in range(3) if r != c]
    max_cross = max(cross_axis_vals) if cross_axis_vals else 0
    rms = metadata.get("residual_rms", 0)
    fig.suptitle(
        f"六位置加速度计标定\n"
        f"零偏 = [{bias[0]:.6f}, {bias[1]:.6f}, {bias[2]:.6f}] m/s^2  |  "
        f"比例因子误差 = [{scale_err[0]:.2f}, {scale_err[1]:.2f}, {scale_err[2]:.2f}] ppm  |  "
        f"最大交轴耦合 = {max_cross:.2f} ppm  |  "
        f"残差RMS = {rms:.6f} m/s^2",
        fontsize=10,
    )

    fig.tight_layout(rect=[0, 0.08, 1, 0.88])
    path = os.path.join(out_dir, "六位置标定对比图.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"  -> {path}")


def plot_gyro_static(data: dict, out_dir: str, dpi: int):
    """陀螺标定结果图：静态角速率时间序列 + 分布直方图。"""
    dtheta_deg = data["dtheta_deg"]  # (N, 3) 单位 deg（每帧角增量）
    dt = data["dt"]
    bias_deg_h = data["bias_deg_h"]  # (3,) 单位 deg/h
    scale_ppm = data.get("scale_ppm", np.zeros(3))  # (3,) 单位 ppm

    if dtheta_deg.shape[0] < 2:
        print("  [WARNING] 陀螺静态数据不足，跳过画图")
        return

    # 角增量(deg) -> 角速率(deg/h): deg / dt_s * 3600
    dtheta_degh = dtheta_deg / dt * 3600.0

    # 构造时间轴
    t = np.arange(dtheta_deg.shape[0]) * dt

    fig, axes = plt.subplots(3, 2, figsize=(14, 10))
    axis_names = ["X", "Y", "Z"]

    for i in range(3):
        # 左列：时间序列
        ax_ts = axes[i, 0]
        ax_ts.plot(t, dtheta_degh[:, i], linewidth=0.5, color="#1565C0", alpha=0.8)
        ax_ts.axhline(y=np.mean(dtheta_degh[:, i]), color="red",
                      linewidth=1.2, linestyle="--", label=f"均值 = {np.mean(dtheta_degh[:, i]):.2f} deg/h")
        ax_ts.set_xlabel("时间 (s)")
        ax_ts.set_ylabel("角速率 (deg/h)")
        ax_ts.set_title(f"{axis_names[i]}轴陀螺 -- 原始角速率序列")
        ax_ts.grid(True, alpha=0.25)
        ax_ts.legend(fontsize=8)

        # 右列：直方图
        ax_hist = axes[i, 1]
        vals = dtheta_degh[:, i]
        ax_hist.hist(vals, bins=60, density=True, color="#42A5F5", alpha=0.7, edgecolor="white", linewidth=0.3)

        # 叠加正态分布拟合
        mu, sigma = np.mean(vals), np.std(vals)
        if sigma > 0:
            x_fit = np.linspace(mu - 4 * sigma, mu + 4 * sigma, 200)
            y_fit = np.exp(-0.5 * ((x_fit - mu) / sigma) ** 2) / (sigma * np.sqrt(2 * np.pi))
            ax_hist.plot(x_fit, y_fit, color="red", linewidth=1.5, label=f"正态分布 N({mu:.2f}, {sigma:.2f}^2)")
        ax_hist.set_xlabel("角速率 (deg/h)")
        ax_hist.set_ylabel("概率密度")
        ax_hist.set_title(f"{axis_names[i]}轴陀螺 -- 分布")
        ax_hist.grid(True, alpha=0.25)
        ax_hist.legend(fontsize=8)

    fig.suptitle(
        f"陀螺标定结果（两位置法）\n"
        f"零偏 = [{bias_deg_h[0]:.4f}, "
        f"{bias_deg_h[1]:.4f}, "
        f"{bias_deg_h[2]:.4f}] deg/h  |  "
        f"比例因子误差 = [{scale_ppm[0]:.1f}, "
        f"{scale_ppm[1]:.1f}, "
        f"{scale_ppm[2]:.1f}] ppm\n"
        f"（注：比例因子由地球自转估计，MEMS噪声下不确定性较大，精确值需速率转台标定）",
        fontsize=10,
    )
    fig.tight_layout(rect=[0, 0.03, 1, 0.88])
    path = os.path.join(out_dir, "陀螺静态数据图.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"  -> {path}")


def plot_alignment(data: dict, out_dir: str, dpi: int):
    """对准结果可视化：姿态角柱状图 + 条件数指示。"""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # 左图：姿态角柱状图
    ax_bar = axes[0]
    angles = [data["roll_deg"], data["pitch_deg"], data["yaw_deg"]]
    names = ["Roll", "Pitch", "Yaw"]
    colors = ["#E53935", "#43A047", "#1E88E5"]
    bars = ax_bar.bar(names, angles, color=colors, alpha=0.85, edgecolor="white", width=0.5)
    ax_bar.axhline(y=0, color="black", linewidth=0.8)
    for bar_obj, val in zip(bars, angles):
        ax_bar.text(bar_obj.get_x() + bar_obj.get_width() / 2,
                    bar_obj.get_height() + (2 if val >= 0 else -4),
                    f"{val:.2f} °", ha="center", fontsize=10, fontweight="bold")
    ax_bar.set_ylabel("角度 (°)")
    ax_bar.set_title("TRIAD 初始对准 -- 姿态角")
    ax_bar.grid(True, alpha=0.25, axis="y")

    # 右图：信息面板
    ax_info = axes[1]
    ax_info.axis("off")
    valid_text = "是" if data["valid"] else "否"
    valid_color = "green" if data["valid"] else "red"
    cond = data["condition_number"]
    cond_warn = " （接近奇异，偏航角不可靠）" if cond > 100 else ""

    info_lines = [
        f"横滚角:  {data['roll_deg']:.4f} °",
        f"俯仰角: {data['pitch_deg']:.4f} °",
        f"偏航角:   {data['yaw_deg']:.4f} °",
        "",
        f"有效:          {valid_text}",
        f"条件数:    {cond:.2f}{cond_warn}",
    ]
    for i, line in enumerate(info_lines):
        color = valid_color if line.startswith("有效:") else "black"
        size = 14 if i < 3 else 12
        ax_info.text(0.1, 0.9 - i * 0.09, line, fontsize=size, color=color,
                     transform=ax_info.transAxes)

    # 说明文字
    explanations = [
        "说明:",
        "有效 -- 加速度计量程合理且陀螺信噪比足够，对准结果可靠",
        "条件数 -- 重力与地球自转轴夹角的正弦倒数，",
        "          值越大表示两矢量越接近共线，>100时偏航角不可靠",
    ]
    for i, line in enumerate(explanations):
        ax_info.text(0.1, 0.28 - i * 0.06, line, fontsize=8, color="#555555",
                     transform=ax_info.transAxes)

    ax_info.set_title("对准结果汇总", fontweight="bold")

    fig.tight_layout()
    path = os.path.join(out_dir, "对准结果图.png")
    fig.savefig(path, dpi=dpi, bbox_inches="tight")
    print(f"  -> {path}")


def main() -> None:
    args = parse_args()

    data_dir = args.data_dir
    if not os.path.isdir(data_dir):
        print(f"[ERROR] 数据目录不存在: {data_dir}")
        print(f"请先在 C++ 端调用 saveCalibDataForPlot(\"{data_dir}\", input, result) 导出数据")
        sys.exit(1)

    print(f"Reading data from: {data_dir}")

    # 六位置标定图
    six_path = os.path.join(data_dir, "calib_six_position.txt")
    if os.path.isfile(six_path):
        plot_six_position(load_six_position(six_path), data_dir, args.dpi)
    else:
        print(f"  [SKIP] {six_path} not found")

    # 静态陀螺数据图
    gyro_path = os.path.join(data_dir, "calib_gyro_static.txt")
    if os.path.isfile(gyro_path):
        plot_gyro_static(load_gyro_static(gyro_path), data_dir, args.dpi)
    else:
        print(f"  [SKIP] {gyro_path} not found")

    # 对准结果图
    align_path = os.path.join(data_dir, "calib_alignment.txt")
    if os.path.isfile(align_path):
        plot_alignment(load_alignment(align_path), data_dir, args.dpi)
    else:
        print(f"  [SKIP] {align_path} not found")

    print("\nDone.")

    if args.show:
        plt.show()
    else:
        plt.close("all")


if __name__ == "__main__":
    main()
