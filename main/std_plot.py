import argparse
import os

import matplotlib.pyplot as plt
import numpy as np


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(description="绘制组合导航标准差结果图")
	parser.add_argument(
		"--std",
		type=str,
		default="data/STD_result.txt",
		help="标准差文件路径（默认data/STD_result.txt）",
	)
	parser.add_argument(
		"--out-dir",
		type=str,
		default="data",
		help="图片输出目录（默认data）",
	)
	parser.add_argument(
		"--dpi",
		type=int,
		default=140,
		help="图片DPI",
	)
	parser.add_argument(
		"--show",
		dest="show",
		action="store_true",
		help="绘图后弹窗展示所有图片（默认开启）",
	)
	parser.add_argument(
		"--no-show",
		dest="show",
		action="store_false",
		help="仅保存图片，不弹窗展示",
	)
	parser.set_defaults(show=True)
	return parser.parse_args()


def load_std_file(path: str) -> np.ndarray:
	if not os.path.isfile(path):
		raise FileNotFoundError(f"标准差文件不存在: {path}")

	data = np.loadtxt(path, dtype=float)
	if data.ndim == 1:
		data = data.reshape(1, -1)

	if data.shape[1] < 22:
		raise ValueError(f"标准差文件列数不足22列: {path}")
	return data


def save_three_axis_plot(
	t: np.ndarray,
	x: np.ndarray,
	y: np.ndarray,
	z: np.ndarray,
	labels,
	ylabel: str,
	title: str,
	file_name: str,
	out_dir: str,
	dpi: int,
) -> None:
	fig, ax = plt.subplots(figsize=(12, 6))
	ax.plot(t, x, label=labels[0], linewidth=1.2)
	ax.plot(t, y, label=labels[1], linewidth=1.2)
	ax.plot(t, z, label=labels[2], linewidth=1.2)
	ax.set_xlabel("Time (s)")
	ax.set_ylabel(ylabel)
	ax.set_title(title)
	ax.grid(True, alpha=0.3)
	ax.legend()
	fig.tight_layout()
	fig.savefig(os.path.join(out_dir, file_name), dpi=dpi)


def main() -> None:
	args = parse_args()
	os.makedirs(args.out_dir, exist_ok=True)

	data = load_std_file(args.std)
	t = data[:, 0]

	# 单位转换：C++输出的是SI单位，需要转为工程单位
	deg_per_rad = 180.0 / np.pi
	mgal_per_mps2 = 100000.0
	ppm_per_ratio = 1e6

	# 位置 (m) - 无需转换
	save_three_axis_plot(
		t,
		data[:, 1],
		data[:, 2],
		data[:, 3],
		labels=("N", "E", "D"),
		ylabel="Position std (m)",
		title="Position Std",
		file_name="位置标准差图.png",
		out_dir=args.out_dir,
		dpi=args.dpi,
	)

	# 速度 (m/s) - 无需转换
	save_three_axis_plot(
		t,
		data[:, 4],
		data[:, 5],
		data[:, 6],
		labels=("Vn", "Ve", "Vd"),
		ylabel="Velocity std (m/s)",
		title="Velocity Std",
		file_name="速度标准差图.png",
		out_dir=args.out_dir,
		dpi=args.dpi,
	)

	# 姿态 (已在C++端转为deg，直接绘制)
	save_three_axis_plot(
		t,
		data[:, 7],
		data[:, 8],
		data[:, 9],
		labels=("Roll", "Pitch", "Yaw"),
		ylabel="Attitude std (deg)",
		title="Attitude Std",
		file_name="姿态标准差图.png",
		out_dir=args.out_dir,
		dpi=args.dpi,
	)

	# 陀螺零偏 (rad/s → deg/h)
	save_three_axis_plot(
		t,
		data[:, 10],
		data[:, 11],
		data[:, 12],
		labels=("Gyro bias X", "Gyro bias Y", "Gyro bias Z"),
		ylabel="Gyro bias std (deg/h)",
		title="Gyroscope Bias Std",
		file_name="陀螺零偏标准差图.png",
		out_dir=args.out_dir,
		dpi=args.dpi,
	)

	# 陀螺比例因子 (ratio → ppm)
	save_three_axis_plot(
		t,
		data[:, 16],
		data[:, 17],
		data[:, 18],
		labels=("Gyro scale X", "Gyro scale Y", "Gyro scale Z"),
		ylabel="Gyro scale std (ppm)",
		title="Gyroscope Scale Std",
		file_name="陀螺比例因子标准差图.png",
		out_dir=args.out_dir,
		dpi=args.dpi,
	)

	# 加速度计零偏 (m/s² → mGal)
	save_three_axis_plot(
		t,
		data[:, 13],
		data[:, 14],
		data[:, 15],
		labels=("Accel bias X", "Accel bias Y", "Accel bias Z"),
		ylabel="Accel bias std (mGal)",
		title="Accelerometer Bias Std",
		file_name="加速度计零偏标准差图.png",
		out_dir=args.out_dir,
		dpi=args.dpi,
	)

	# 加速度计比例因子 (ratio → ppm)
	save_three_axis_plot(
		t,
		data[:, 19],
		data[:, 20],
		data[:, 21],
		labels=("Accel scale X", "Accel scale Y", "Accel scale Z"),
		ylabel="Accel scale std (ppm)",
		title="Accelerometer Scale Std",
		file_name="加速度计比例因子标准差图.png",
		out_dir=args.out_dir,
		dpi=args.dpi,
	)

	print("标准差图生成完成")
	print(f"输出目录: {args.out_dir}")
	print(
		"已生成: 位置标准差图.png, 速度标准差图.png, 姿态标准差图.png, "
		"陀螺零偏标准差图.png, 陀螺比例因子标准差图.png, "
		"加速度计零偏标准差图.png, 加速度计比例因子标准差图.png"
	)

	if args.show:
		plt.show()
	else:
		plt.close("all")


if __name__ == "__main__":
	main()
