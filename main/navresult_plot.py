import argparse
import os
from dataclasses import dataclass

import matplotlib.pyplot as plt
import numpy as np

TRUTH_FILE_PATH = "data_ICM20602/truth.nav"


@dataclass
class PvaSeries:
	time: np.ndarray
	lat_deg: np.ndarray
	lon_deg: np.ndarray
	h_m: np.ndarray
	vn: np.ndarray
	ve: np.ndarray
	vd: np.ndarray
	roll_deg: np.ndarray
	pitch_deg: np.ndarray
	yaw_deg: np.ndarray


def parse_args() -> argparse.Namespace:
	parser = argparse.ArgumentParser(
		description="绘制组合导航解算结果与真值PVA对比图"
	)
	parser.add_argument(
		"--nav",
		type=str,
		default="data_ICM20602/Navresult.txt",
		help="解算结果文件路径，默认data_ICM20602/Navresult.txt",
	)
	parser.add_argument(
		"--out-dir",
		type=str,
		default="data_ICM20602",
		help="图片输出目录，默认data_ICM20602",
	)
	parser.add_argument(
		"--truth-angle-unit",
		type=str,
		choices=["deg", "rad"],
		default="deg",
		help="真值文件中的角度单位（lat/lon与roll/pitch/yaw），默认deg",
	)
	parser.add_argument(
		"--dpi",
		type=int,
		default=140,
		help="输出图片DPI",
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


def load_pva_file(path: str, angle_unit: str = "deg") -> PvaSeries:
	if not os.path.isfile(path):
		raise FileNotFoundError(f"文件不存在: {path}")

	data = np.loadtxt(path, dtype=float)
	if data.ndim == 1:
		data = data.reshape(1, -1)

	if data.shape[1] < 10:
		raise ValueError(f"文件列数不足10列: {path}")

	# 兼容两种格式:
	# 1) 10列: [time, lat, lon, h, vn, ve, vd, roll, pitch, yaw]
	# 2) 11列: [week, time, lat, lon, h, vn, ve, vd, roll, pitch, yaw]
	start_col = 0
	if data.shape[1] >= 11:
		col0_span = float(np.nanmax(data[:, 0]) - np.nanmin(data[:, 0]))
		col1_span = float(np.nanmax(data[:, 1]) - np.nanmin(data[:, 1]))
		if col0_span < 1e-9 and col1_span > 1e-9:
			start_col = 1

	time = data[:, start_col + 0]
	lat = data[:, start_col + 1]
	lon = data[:, start_col + 2]
	h = data[:, start_col + 3]
	vn = data[:, start_col + 4]
	ve = data[:, start_col + 5]
	vd = data[:, start_col + 6]
	roll = data[:, start_col + 7]
	pitch = data[:, start_col + 8]
	yaw = data[:, start_col + 9]

	if angle_unit == "rad":
		lat = np.rad2deg(lat)
		lon = np.rad2deg(lon)
		roll = np.rad2deg(roll)
		pitch = np.rad2deg(pitch)
		yaw = np.rad2deg(yaw)

	return PvaSeries(
		time=time,
		lat_deg=lat,
		lon_deg=lon,
		h_m=h,
		vn=vn,
		ve=ve,
		vd=vd,
		roll_deg=roll,
		pitch_deg=pitch,
		yaw_deg=yaw,
	)


def unwrap_deg(angle_deg: np.ndarray) -> np.ndarray:
	return np.rad2deg(np.unwrap(np.deg2rad(angle_deg)))


def wrap_to_180(angle_deg: np.ndarray) -> np.ndarray:
	return (angle_deg + 180.0) % 360.0 - 180.0


def interp_to_time(src_t: np.ndarray, src_y: np.ndarray, dst_t: np.ndarray) -> np.ndarray:
	return np.interp(dst_t, src_t, src_y)


def align_series(nav: PvaSeries, truth: PvaSeries):
	t_start = max(nav.time[0], truth.time[0])
	t_end = min(nav.time[-1], truth.time[-1])
	if t_end <= t_start:
		raise ValueError("解算结果与真值时间区间无重叠")

	mask = (nav.time >= t_start) & (nav.time <= t_end)
	t = nav.time[mask]
	if t.size < 2:
		raise ValueError("重叠区间样本数不足")

	nav_sel = {
		"lat": nav.lat_deg[mask],
		"lon": nav.lon_deg[mask],
		"h": nav.h_m[mask],
		"vn": nav.vn[mask],
		"ve": nav.ve[mask],
		"vd": nav.vd[mask],
		"roll": unwrap_deg(nav.roll_deg)[mask],
		"pitch": unwrap_deg(nav.pitch_deg)[mask],
		"yaw": unwrap_deg(nav.yaw_deg)[mask],
	}

	truth_interp = {
		"lat": interp_to_time(truth.time, truth.lat_deg, t),
		"lon": interp_to_time(truth.time, truth.lon_deg, t),
		"h": interp_to_time(truth.time, truth.h_m, t),
		"vn": interp_to_time(truth.time, truth.vn, t),
		"ve": interp_to_time(truth.time, truth.ve, t),
		"vd": interp_to_time(truth.time, truth.vd, t),
		"roll": interp_to_time(truth.time, unwrap_deg(truth.roll_deg), t),
		"pitch": interp_to_time(truth.time, unwrap_deg(truth.pitch_deg), t),
		"yaw": interp_to_time(truth.time, unwrap_deg(truth.yaw_deg), t),
	}

	return t, nav_sel, truth_interp


def geodetic_error_to_ned_m(nav, truth):
	lat_ref_rad = np.deg2rad(truth["lat"])
	h_ref = truth["h"]

	# WGS84参数。
	a = 6378137.0
	f = 1.0 / 298.257223563
	e2 = 2.0 * f - f * f

	sin_lat = np.sin(lat_ref_rad)
	denom = 1.0 - e2 * sin_lat * sin_lat
	sqrt_denom = np.sqrt(denom)
	rn = a / sqrt_denom
	rm = a * (1.0 - e2) / (denom * sqrt_denom)

	dlat_rad = np.deg2rad(nav["lat"] - truth["lat"])
	dlon_rad = np.deg2rad(nav["lon"] - truth["lon"])
	dh = nav["h"] - truth["h"]

	dn = dlat_rad * (rm + h_ref)
	de = dlon_rad * (rn + h_ref) * np.cos(lat_ref_rad)
	dd = -dh

	return dn, de, dd


def save_plot_position(t, nav, truth, out_dir, dpi):
	dn, de, dd = geodetic_error_to_ned_m(nav, truth)

	fig, ax = plt.subplots(figsize=(12, 6))
	ax.plot(t, dn, label="North error", linewidth=1.2)
	ax.plot(t, de, label="East error", linewidth=1.2)
	ax.plot(t, dd, label="Down error", linewidth=1.2)
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Distance error (m)")
	ax.set_title("Position Error")
	ax.grid(True, alpha=0.3)
	ax.legend()
	fig.tight_layout()
	fig.savefig(os.path.join(out_dir, "位置误差图.png"), dpi=dpi)


def save_plot_velocity(t, nav, truth, out_dir, dpi):
	evn = nav["vn"] - truth["vn"]
	eve = nav["ve"] - truth["ve"]
	evd = nav["vd"] - truth["vd"]

	fig, ax = plt.subplots(figsize=(12, 6))
	ax.plot(t, evn, label="Vn error", linewidth=1.2)
	ax.plot(t, eve, label="Ve error", linewidth=1.2)
	ax.plot(t, evd, label="Vd error", linewidth=1.2)
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Velocity error (m/s)")
	ax.set_title("Velocity Error")
	ax.grid(True, alpha=0.3)
	ax.legend()
	fig.tight_layout()
	fig.savefig(os.path.join(out_dir, "速度误差图.png"), dpi=dpi)


def save_plot_attitude(t, nav, truth, out_dir, dpi):
	eroll = wrap_to_180(nav["roll"] - truth["roll"])
	epitch = wrap_to_180(nav["pitch"] - truth["pitch"])
	eyaw = wrap_to_180(nav["yaw"] - truth["yaw"])

	fig, ax = plt.subplots(figsize=(12, 6))
	ax.plot(t, eroll, label="Roll error", linewidth=1.2)
	ax.plot(t, epitch, label="Pitch error", linewidth=1.2)
	ax.plot(t, eyaw, label="Yaw error", linewidth=1.2)
	ax.set_xlabel("Time (s)")
	ax.set_ylabel("Attitude error (deg)")
	ax.set_title("Attitude Error")
	ax.grid(True, alpha=0.3)
	ax.legend()
	fig.tight_layout()
	fig.savefig(os.path.join(out_dir, "姿态误差图.png"), dpi=dpi)


def save_plot_track(nav, truth, out_dir, dpi):
	fig, ax = plt.subplots(figsize=(8, 7))
	ax.plot(nav["lon"], nav["lat"], label="nav", linewidth=1.2)
	ax.plot(truth["lon"], truth["lat"], label="truth", linewidth=1.2)
	ax.set_xlabel("Longitude (deg)")
	ax.set_ylabel("Latitude (deg)")
	ax.set_title("2D Track Comparison")
	ax.grid(True, alpha=0.3)
	ax.legend()
	fig.tight_layout()
	fig.savefig(os.path.join(out_dir, "解算轨迹图.png"), dpi=dpi)


def main() -> None:
	args = parse_args()

	os.makedirs(args.out_dir, exist_ok=True)

	nav = load_pva_file(args.nav, angle_unit="deg")
	truth = load_pva_file(TRUTH_FILE_PATH, angle_unit=args.truth_angle_unit)

	t, nav_aligned, truth_aligned = align_series(nav, truth)

	save_plot_position(t, nav_aligned, truth_aligned, args.out_dir, args.dpi)
	save_plot_velocity(t, nav_aligned, truth_aligned, args.out_dir, args.dpi)
	save_plot_attitude(t, nav_aligned, truth_aligned, args.out_dir, args.dpi)
	save_plot_track(nav_aligned, truth_aligned, args.out_dir, args.dpi)

	print("对比图生成完成")
	print(f"输出目录: {args.out_dir}")
	print(f"真值文件: {TRUTH_FILE_PATH}")
	print("已生成: 解算轨迹图.png, 位置误差图.png, 速度误差图.png, 姿态误差图.png")

	if args.show:
		plt.show()
	else:
		plt.close("all")


if __name__ == "__main__":
	main()
