#include <algorithm>
#include <array>
#include <cstring>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include "../fileio/filesaver.h"
#include "../fileio/gnssfileloader.h"
#include "../fileio/imufileloader.h"
#include "fusion_workflow.h"
#include "is_update.h"
namespace ins {
// 作用：degPerHourToRadPerSec 函数。
inline double degPerHourToRadPerSec(double deg_per_hour) {
	// 作用：deg2rad 函数。
	return deg2rad(deg_per_hour) / 3600.0;
}
// 作用：degPerSqrtHourToRadPerSqrtSec 函数。
inline double degPerSqrtHourToRadPerSqrtSec(double deg_per_sqrt_hour) {
	// 作用：deg2rad 函数。
	return deg2rad(deg_per_sqrt_hour) / std::sqrt(3600.0);
}
// 作用：mGalToMps2 函数。
inline double mGalToMps2(double mgal) {
	return mgal * 1e-5;
}
// 作用：mpsPerSqrtHourToMps2PerSqrtHz 函数。
inline double mpsPerSqrtHourToMps2PerSqrtHz(double mps_per_sqrt_hour) {
	return mps_per_sqrt_hour / std::sqrt(3600.0);
}
// 作用：ppmToRatio 函数。
inline double ppmToRatio(double ppm) {
	return ppm * 1e-6;
}
// 作用：normalizeGnssAngleToRadIfNeeded 函数。
inline void normalizeGnssAngleToRadIfNeeded(GnssData& gnss) {
	const double kRadMax = 3.2;
	if (std::fabs(gnss.latitude) > kRadMax || std::fabs(gnss.longitude) > kRadMax) {
		gnss.latitude = deg2rad(gnss.latitude);
		gnss.longitude = deg2rad(gnss.longitude);
	}
}
// 作用：定义 InitialParameterConfig 数据结构。
struct InitialParameterConfig {
	bool use_first_gnss_position = false;
	std::array<double, 3> init_blh = {30.4604325443,114.4725046685,23.000};
	bool init_blh_in_degree = true;
	std::array<double, 3> init_vel_n = {0.0, 0.0, 0.0};
	std::array<double, 3> init_euler = {1.63086, -0.13135, 276.8827};
	bool init_euler_in_degree = true;
	ImuErrorData init_imu_error = [] {
		ImuErrorData v;
		v.gyro_random_walk = {0.2, 0.2, 0.2};
		v.accel_random_walk = {0.2, 0.2, 0.2};
		return v;
	}();
	std::array<double, 3> init_gyro_bias = {0.0, 0.0, 0.0};
	std::array<double, 3> init_accel_bias = {0.0, 0.0, 0.0};
	std::array<double, 3> init_gyro_scale = {0.0, 0.0, 0.0};
	std::array<double, 3> init_accel_scale = {0.0, 0.0, 0.0};
	std::array<double, 3> pos_std = {0.1, 0.1, 0.1};                          // m
	std::array<double, 3> vel_std = {0.1, 0.1, 0.1};                              // m/s
	std::array<double, 3> att_std = {deg2rad(0.1), deg2rad(0.1), deg2rad(0.2)};  // rad
	std::array<double, 3> gyro_bias_std = {600, 600, 600};
	std::array<double, 3> accel_bias_std = {3000, 3000, 3000};
	std::array<double, 3> gyro_scale_std = {600, 600, 600};
	std::array<double, 3> accel_scale_std = {3000, 3000, 3000};
	std::array<double, 3> gyro_bias_process_std = {200, 200, 200};
	std::array<double, 3> accel_bias_process_std = {1000, 1000, 1000};
	std::array<double, 3> gyro_scale_process_std = {200, 200, 200};
	std::array<double, 3> accel_scale_process_std = {1000, 1000, 1000};
};
// 作用：定义 MainProgramConfig 数据结构。
struct MainProgramConfig {
	std::string imu_file_path = "data/ICM20602.txt";
	std::string gnss_file_path = "data/GNSS_RTK.txt";
	std::string output_dir = "data";
	double corrtime_gb = 3600.0;
	double corrtime_ab = 3600.0;
	double corrtime_gs = 3600.0;
	double corrtime_as = 3600.0;
	std::array<double, 3> antenna_lever_arm_b = {-0.073, 0.302, 0.087};
	double process_start_s = 357473;
	double process_end_s = -1.0;
	double gnss_sync_tolerance = 1e-3;
	std::size_t progress_bar_width = 40;
	std::size_t progress_update_step = 5000;
	InitialParameterConfig init;
};
// 作用：fileExists 函数。
inline bool fileExists(const std::string& file_path) {
	// 作用：ifs 函数。
	std::ifstream ifs(file_path, std::ios::binary);
	return ifs.good();
}
// 作用：resolveImuFilePath 函数。
inline std::string resolveImuFilePath(const MainProgramConfig& cfg) {
	if (!cfg.imu_file_path.empty() && fileExists(cfg.imu_file_path)) {
		return cfg.imu_file_path;
	}
	throw std::runtime_error("没锟斤拷imu锟侥硷拷");
}
// 作用：countTextLinesFast 函数。
inline std::size_t countTextLinesFast(const std::string& file_path) {
	// 作用：ifs 函数。
	std::ifstream ifs(file_path);
	if (!ifs.is_open()) {
		return 0;
	}
	std::size_t lines = 0;
	std::string line;
	while (std::getline(ifs, line)) {
		if (!line.empty()) {
			++lines;
		}
	}
	return lines;
}
// 作用：estimateImuRecordCount 函数。
inline std::size_t estimateImuRecordCount(const std::string& imu_file_path) {
	ImuFileLoader loader;
	if (!loader.open(imu_file_path)) {
		return 0;
	}
	if (loader.isBinary()) {
		// 作用：ifs 函数。
		std::ifstream ifs(imu_file_path, std::ios::binary | std::ios::ate);
		if (!ifs.is_open()) {
			return 0;
		}
		const std::streamoff size = ifs.tellg();
		if (size <= 0) {
			return 0;
		}
		return static_cast<std::size_t>(size / static_cast<std::streamoff>(sizeof(double) * 7));
	}
	// 作用：countTextLinesFast 函数。
	return countTextLinesFast(imu_file_path);
}
// 作用：calcMarkovPhi 函数。
inline double calcMarkovPhi(double dt, double tau) {
	const double tau_safe = (tau > 1e-6) ? tau : 1e-6;
	return std::exp(-dt / tau_safe);
}
// 作用：skewMatrix 函数。
inline Matrix skewMatrix(const Vec3& v) {
	// 作用：s 函数。
	Matrix s(3, 3, 0.0);
	s(0, 1) = -v.z;
	s(0, 2) = v.y;
	s(1, 0) = v.z;
	s(1, 2) = -v.x;
	s(2, 0) = -v.y;
	s(2, 1) = v.x;
	return s;
}
// 作用：mat3ToMatrix 函数。
inline Matrix mat3ToMatrix(const Mat3& m) {
	// 作用：out 函数。
	Matrix out(3, 3, 0.0);
	for (std::size_t r = 0; r < 3; ++r) {
		for (std::size_t c = 0; c < 3; ++c) {
			out(r, c) = m.m[r][c];
		}
	}
	return out;
}
// 作用：diagFromVec3 函数。
inline Matrix diagFromVec3(const Vec3& v) {
	// 作用：d 函数。
	Matrix d(3, 3, 0.0);
	d(0, 0) = v.x;
	d(1, 1) = v.y;
	d(2, 2) = v.z;
	return d;
}
// 作用：addBlock3x3 函数。
inline void addBlock3x3(Matrix& dst,
						std::size_t row0,
						std::size_t col0,
						const Matrix& blk,
						double scale = 1.0) {
	for (std::size_t r = 0; r < 3; ++r) {
		for (std::size_t c = 0; c < 3; ++c) {
			dst(row0 + r, col0 + c) += scale * blk(r, c);
		}
	}
}
// 作用：buildAlignedTransitionF 函数。
inline Matrix buildAlignedTransitionF(double dt,
								const NavigationStatusData& nav_state,
								const ImuMeasureData& imu_measure,
								double t_gb,
								double t_ab,
								double t_gs,
								double t_as) {
	const double dt_safe = (dt > 1e-6) ? dt : 1e-6;
	const double tgb_safe = (t_gb > 1e-6) ? t_gb : 1e-6;
	const double tab_safe = (t_ab > 1e-6) ? t_ab : 1e-6;
	const double tgs_safe = (t_gs > 1e-6) ? t_gs : 1e-6;
	const double tas_safe = (t_as > 1e-6) ? t_as : 1e-6;
	// 修正Issue 1：使用pvapre_而不是pvacur_来进行线性化
	const PvaData& pva = nav_state.pvapre_;
	const Blh blh{pva.blh[0], pva.blh[1], pva.blh[2]};
	const Vec3 vel_n{pva.vel_n[0], pva.vel_n[1], pva.vel_n[2]};
	const LocalPosition pos_local{blh.lat, blh.lon, blh.h};
	const LocalVelocity vel_local{vel_n.x, vel_n.y, vel_n.z};
	const EarthParameters earth = updateEarthParameters(pos_local, vel_local);
	const double rm = earth.rm;
	const double rn = earth.rn;
	const double h = blh.h;
	const double rm_h = rm + h;
	const double rn_h = rn + h;
	const double v_n = vel_n.x;
	const double v_e = vel_n.y;
	const double v_d = vel_n.z;
	const double tan_lat = std::tan(blh.lat);
	const double cos_lat = std::cos(blh.lat);
	const double sin_lat = std::sin(blh.lat);
	const double sec2_lat = 1.0 / std::max(1e-12, cos_lat * cos_lat);
	const double safe_cos_lat = (std::fabs(cos_lat) < 1e-8)
								  ? ((cos_lat >= 0.0) ? 1e-8 : -1e-8)
								  : cos_lat;
	const Euler e{pva.euler[0], pva.euler[1], pva.euler[2]};
	const Quaternion q_nb = euler2quat(e);
	const Matrix c_nb = mat3ToMatrix(quat2dcm(q_nb));
	const ImuData imu_corr = compensateSingleImuByError(imu_measure.imucur_, nav_state.imuerror_, dt_safe);
	const Vec3 f_b{imu_corr.dvel_x / dt_safe, imu_corr.dvel_y / dt_safe, imu_corr.dvel_z / dt_safe};
	const Vec3 omega_ib_b{imu_corr.dtheta_x / dt_safe, imu_corr.dtheta_y / dt_safe, imu_corr.dtheta_z / dt_safe};
	const Vec3 omega_ie_n = iewn(blh);
	const Vec3 omega_en_n = enwn(blh, vel_n);
	const Vec3 omega_in_n = addVec3(omega_ie_n, omega_en_n);
	const double omega_e = earth.omega_ie;
	// 作用：fc 函数。
	Matrix fc(ErrorStateIndex21::kDim, ErrorStateIndex21::kDim, 0.0);
	// F_rr
	fc(ErrorStateIndex21::kPos + 0, ErrorStateIndex21::kPos + 0) = -v_d / rm_h;
	fc(ErrorStateIndex21::kPos + 0, ErrorStateIndex21::kPos + 2) = v_n / rm_h;
	fc(ErrorStateIndex21::kPos + 1, ErrorStateIndex21::kPos + 0) = v_e * tan_lat / rn_h;
	fc(ErrorStateIndex21::kPos + 1, ErrorStateIndex21::kPos + 1) = -(v_d + v_n * tan_lat) / rn_h;
	fc(ErrorStateIndex21::kPos + 1, ErrorStateIndex21::kPos + 2) = v_e / rn_h;
	fc(ErrorStateIndex21::kPos + 0, ErrorStateIndex21::kVel + 0) = 1.0;
	fc(ErrorStateIndex21::kPos + 1, ErrorStateIndex21::kVel + 1) = 1.0;
	fc(ErrorStateIndex21::kPos + 2, ErrorStateIndex21::kVel + 2) = 1.0;
	fc(ErrorStateIndex21::kVel + 0, ErrorStateIndex21::kPos + 0) =
		-2.0 * omega_e * v_e * safe_cos_lat / rm_h
		- (v_e * v_e * sec2_lat) / (rm_h * rn_h);
	fc(ErrorStateIndex21::kVel + 0, ErrorStateIndex21::kPos + 2) =
		(v_n * v_d) / (rm_h * rm_h)
		- (v_e * v_e * tan_lat) / (rn_h * rn_h);
	fc(ErrorStateIndex21::kVel + 1, ErrorStateIndex21::kPos + 0) =
		2.0 * omega_e * (v_n * safe_cos_lat - v_d * sin_lat) / rm_h
		+ (v_n * v_e * sec2_lat) / (rm_h * rn_h);
	fc(ErrorStateIndex21::kVel + 1, ErrorStateIndex21::kPos + 2) =
		(v_e * v_d) / (rn_h * rn_h)
		+ (v_n * v_e * tan_lat) / (rn_h * rn_h);
	fc(ErrorStateIndex21::kVel + 2, ErrorStateIndex21::kPos + 0) =
		2.0 * omega_e * v_e * sin_lat / rm_h;
	fc(ErrorStateIndex21::kVel + 2, ErrorStateIndex21::kPos + 2) =
		-(v_e * v_e) / (rn_h * rn_h)
		- (v_n * v_n) / (rm_h * rm_h)
		+ 2.0 * earth.gravity / (std::sqrt(rm * rn) + h);
	fc(ErrorStateIndex21::kVel + 0, ErrorStateIndex21::kVel + 0) = v_d / rm_h;
	fc(ErrorStateIndex21::kVel + 0, ErrorStateIndex21::kVel + 1) =
		-2.0 * (omega_e * sin_lat + v_e * tan_lat / rn_h);
	fc(ErrorStateIndex21::kVel + 0, ErrorStateIndex21::kVel + 2) = v_n / rm_h;
	fc(ErrorStateIndex21::kVel + 1, ErrorStateIndex21::kVel + 0) =
		2.0 * omega_e * sin_lat + v_e * tan_lat / rn_h;
	fc(ErrorStateIndex21::kVel + 1, ErrorStateIndex21::kVel + 1) =
		(v_d + v_n * tan_lat) / rn_h;
	fc(ErrorStateIndex21::kVel + 1, ErrorStateIndex21::kVel + 2) =
		2.0 * omega_e * safe_cos_lat + v_e / rn_h;
	fc(ErrorStateIndex21::kVel + 2, ErrorStateIndex21::kVel + 0) = -2.0 * v_n / rm_h;
	fc(ErrorStateIndex21::kVel + 2, ErrorStateIndex21::kVel + 1) =
		-2.0 * (omega_e * safe_cos_lat + v_e / rn_h);
	// F_vphi, F_vba, F_vsa
	addBlock3x3(fc,
				 ErrorStateIndex21::kVel,
				 ErrorStateIndex21::kAtt,
				 skewMatrix(matVec(quat2dcm(q_nb), f_b)), 1.0);
	addBlock3x3(fc,
				 ErrorStateIndex21::kVel,
				 ErrorStateIndex21::kAccelBias, c_nb, 1.0);
	addBlock3x3(fc,
				 ErrorStateIndex21::kVel,
				 // 作用：mul 函数。
				 ErrorStateIndex21::kAccelScale, mul(c_nb, diagFromVec3(f_b)), 1.0);
	fc(ErrorStateIndex21::kAtt + 0, ErrorStateIndex21::kPos + 0) = -omega_e * sin_lat / rm_h;
	fc(ErrorStateIndex21::kAtt + 0, ErrorStateIndex21::kPos + 2) = v_e / (rn_h * rn_h);
	fc(ErrorStateIndex21::kAtt + 1, ErrorStateIndex21::kPos + 2) = -v_n / (rm_h * rm_h);
	fc(ErrorStateIndex21::kAtt + 2, ErrorStateIndex21::kPos + 0) =
		omega_e * safe_cos_lat / rm_h - v_e * sec2_lat / (rm_h * rn_h);
	fc(ErrorStateIndex21::kAtt + 2, ErrorStateIndex21::kPos + 2) =
		-v_e * tan_lat / (rn_h * rn_h);
	// F_phiphi, F_phiv, F_phibg, F_phisg
	addBlock3x3(fc,
				 ErrorStateIndex21::kAtt,
				 ErrorStateIndex21::kAtt,
				 skewMatrix(omega_in_n),
				 -1.0);
	fc(ErrorStateIndex21::kAtt + 0, ErrorStateIndex21::kVel + 1) = 1.0 / rn_h;
	fc(ErrorStateIndex21::kAtt + 1, ErrorStateIndex21::kVel + 0) = -1.0 / rm_h;
	fc(ErrorStateIndex21::kAtt + 2, ErrorStateIndex21::kVel + 1) = -tan_lat / rn_h;
	addBlock3x3(fc,
				 ErrorStateIndex21::kAtt,
				 ErrorStateIndex21::kGyroBias,
				 c_nb,
				 -1.0);
	addBlock3x3(fc,
				 ErrorStateIndex21::kAtt,
				 ErrorStateIndex21::kGyroScale,
				 mul(c_nb, diagFromVec3(omega_ib_b)),
				 -1.0);
	for (std::size_t i = 0; i < 3; ++i) {
		fc(ErrorStateIndex21::kGyroBias + i, ErrorStateIndex21::kGyroBias + i) = -1.0 / tgb_safe;
		fc(ErrorStateIndex21::kAccelBias + i, ErrorStateIndex21::kAccelBias + i) = -1.0 / tab_safe;
		fc(ErrorStateIndex21::kGyroScale + i, ErrorStateIndex21::kGyroScale + i) = -1.0 / tgs_safe;
		fc(ErrorStateIndex21::kAccelScale + i, ErrorStateIndex21::kAccelScale + i) = -1.0 / tas_safe;
	}
	Matrix f = Matrix::identity(ErrorStateIndex21::kDim);
	for (std::size_t r = 0; r < ErrorStateIndex21::kDim; ++r) {
		for (std::size_t c = 0; c < ErrorStateIndex21::kDim; ++c) {
			f(r, c) += fc(r, c) * dt_safe;
		}
	}
	const double phi_gb = calcMarkovPhi(dt_safe, tgb_safe);
	const double phi_ab = calcMarkovPhi(dt_safe, tab_safe);
	const double phi_gs = calcMarkovPhi(dt_safe, tgs_safe);
	const double phi_as = calcMarkovPhi(dt_safe, tas_safe);
	for (std::size_t i = 0; i < 3; ++i) {
		f(ErrorStateIndex21::kGyroBias + i, ErrorStateIndex21::kGyroBias + i) = phi_gb;
		f(ErrorStateIndex21::kAccelBias + i, ErrorStateIndex21::kAccelBias + i) = phi_ab;
		f(ErrorStateIndex21::kGyroScale + i, ErrorStateIndex21::kGyroScale + i) = phi_gs;
		f(ErrorStateIndex21::kAccelScale + i, ErrorStateIndex21::kAccelScale + i) = phi_as;
	}
	return f;
}
// 作用：buildAlignedProcessNoiseQ 函数。
inline Matrix buildAlignedProcessNoiseQ(double dt,
												const Matrix& f,
												const ImuMeasureData& imu_measure,
												const NavigationStatusData& nav_state,
									const InitialParameterConfig& init_cfg,
									double t_gb,
									double t_ab,
									double t_gs,
									double t_as) {
	const double dt_safe = (dt > 1e-6) ? dt : 1e-6;
	const double tgb_safe = (t_gb > 1e-6) ? t_gb : 1e-6;
	const double tab_safe = (t_ab > 1e-6) ? t_ab : 1e-6;
	const double tgs_safe = (t_gs > 1e-6) ? t_gs : 1e-6;
	const double tas_safe = (t_as > 1e-6) ? t_as : 1e-6;

	// [wv, wphi, wgb, wab, wgs, was]
	// 作用：q_c 函数。
	Matrix q_c(18, 18, 0.0);
	for (std::size_t i = 0; i < 3; ++i) {
		const double vrw = mpsPerSqrtHourToMps2PerSqrtHz(
				std::fabs(init_cfg.init_imu_error.accel_random_walk[i]));
		const double arw = degPerSqrtHourToRadPerSqrtSec(
				std::fabs(init_cfg.init_imu_error.gyro_random_walk[i]));
		const double sigma_gb = degPerHourToRadPerSec(std::fabs(init_cfg.gyro_bias_process_std[i]));
		const double sigma_ab = mGalToMps2(std::fabs(init_cfg.accel_bias_process_std[i]));
		const double sigma_gs = ppmToRatio(std::fabs(init_cfg.gyro_scale_process_std[i]));
		const double sigma_as = ppmToRatio(std::fabs(init_cfg.accel_scale_process_std[i]));
		q_c(i, i) = vrw * vrw;
		q_c(3 + i, 3 + i) = arw * arw;
		q_c(6 + i, 6 + i) = 2.0 * sigma_gb * sigma_gb / tgb_safe;
		q_c(9 + i, 9 + i) = 2.0 * sigma_ab * sigma_ab / tab_safe;
		q_c(12 + i, 12 + i) = 2.0 * sigma_gs * sigma_gs / tgs_safe;
		q_c(15 + i, 15 + i) = 2.0 * sigma_as * sigma_as / tas_safe;
	}
	// 修正Issue 3：G(tk-1)使用k-1时刻姿态，G(tk)使用k时刻姿态
	const PvaData& pva_km1 = nav_state.pvapre_;
	const Euler e_km1{pva_km1.euler[0], pva_km1.euler[1], pva_km1.euler[2]};
	const Mat3 c_nb_km1_m3 = quat2dcm(euler2quat(e_km1));
	const Matrix c_nb_km1 = mat3ToMatrix(c_nb_km1_m3);

	const PvaData& pva_curr = nav_state.pvacur_;
	const Euler e_curr{pva_curr.euler[0], pva_curr.euler[1], pva_curr.euler[2]};
	const Mat3 c_nb_curr_m3 = quat2dcm(euler2quat(e_curr));
	const Matrix c_nb_curr = mat3ToMatrix(c_nb_curr_m3);

	// g_km1 使用 k-1 时刻状态 (c_nb_km1)
	Matrix g_km1(ErrorStateIndex21::kDim, 18, 0.0);
	addBlock3x3(g_km1, ErrorStateIndex21::kVel, 0, c_nb_km1, 1.0);
	addBlock3x3(g_km1, ErrorStateIndex21::kAtt, 3, c_nb_km1, 1.0);
	addBlock3x3(g_km1, ErrorStateIndex21::kGyroBias, 6, Matrix::identity(3), 1.0);
	addBlock3x3(g_km1, ErrorStateIndex21::kAccelBias, 9, Matrix::identity(3), 1.0);
	addBlock3x3(g_km1, ErrorStateIndex21::kGyroScale, 12, Matrix::identity(3), 1.0);
	addBlock3x3(g_km1, ErrorStateIndex21::kAccelScale, 15, Matrix::identity(3), 1.0);

	// g_k 使用 k 时刻状态 (c_nb_curr)
	Matrix g_k(ErrorStateIndex21::kDim, 18, 0.0);
	addBlock3x3(g_k, ErrorStateIndex21::kVel, 0, c_nb_curr, 1.0);
	addBlock3x3(g_k, ErrorStateIndex21::kAtt, 3, c_nb_curr, 1.0);
	addBlock3x3(g_k, ErrorStateIndex21::kGyroBias, 6, Matrix::identity(3), 1.0);
	addBlock3x3(g_k, ErrorStateIndex21::kAccelBias, 9, Matrix::identity(3), 1.0);
	addBlock3x3(g_k, ErrorStateIndex21::kGyroScale, 12, Matrix::identity(3), 1.0);
	addBlock3x3(g_k, ErrorStateIndex21::kAccelScale, 15, Matrix::identity(3), 1.0);
	// 作用：q_d 函数。
	Matrix q_d(ErrorStateIndex21::kDim, ErrorStateIndex21::kDim, 0.0);
	
	// G(t_{k-1}) * q(t_{k-1}) * G^T(t_{k-1})
	Matrix g_q_gt_km1 = mul(mul(g_km1, q_c), transpose(g_km1));
	// G(t_k) * q(t_k) * G^T(t_k)
	Matrix g_q_gt_k = mul(mul(g_k, q_c), transpose(g_k));
	
	// Phi * G(t_{k-1}) * q(t_{k-1}) * G^T(t_{k-1}) * Phi^T
	Matrix phi_g_q_gt_phit = mul(mul(f, g_q_gt_km1), transpose(f));
	
	q_d = add(phi_g_q_gt_phit, g_q_gt_k);
	for (std::size_t r = 0; r < q_d.rows(); ++r) {
		for (std::size_t c = 0; c < q_d.cols(); ++c) {
			q_d(r, c) *= (0.5 * dt_safe);
		}
	}
	return q_d;
}
// 作用：buildInitialCovarianceP0 函数。
inline Matrix buildInitialCovarianceP0(const InitialParameterConfig& init_cfg) {
	// 作用：p0 函数。
	Matrix p0(ErrorStateIndex21::kDim, ErrorStateIndex21::kDim, 0.0);
	for (std::size_t i = 0; i < 3; ++i) {
		p0(ErrorStateIndex21::kPos + i, ErrorStateIndex21::kPos + i) =
				init_cfg.pos_std[i] * init_cfg.pos_std[i];
		p0(ErrorStateIndex21::kVel + i, ErrorStateIndex21::kVel + i) =
				init_cfg.vel_std[i] * init_cfg.vel_std[i];
		p0(ErrorStateIndex21::kAtt + i, ErrorStateIndex21::kAtt + i) =
				init_cfg.att_std[i] * init_cfg.att_std[i];
		const double gb_std = degPerHourToRadPerSec(std::fabs(init_cfg.gyro_bias_std[i]));
		const double ab_std = mGalToMps2(std::fabs(init_cfg.accel_bias_std[i]));
		const double gs_std = ppmToRatio(std::fabs(init_cfg.gyro_scale_std[i]));
		const double as_std = ppmToRatio(std::fabs(init_cfg.accel_scale_std[i]));
		p0(ErrorStateIndex21::kGyroBias + i, ErrorStateIndex21::kGyroBias + i) =
				gb_std * gb_std;
		p0(ErrorStateIndex21::kAccelBias + i, ErrorStateIndex21::kAccelBias + i) =
				ab_std * ab_std;
		p0(ErrorStateIndex21::kGyroScale + i, ErrorStateIndex21::kGyroScale + i) =
				gs_std * gs_std;
		p0(ErrorStateIndex21::kAccelScale + i, ErrorStateIndex21::kAccelScale + i) =
				as_std * as_std;
	}
	return p0;
}
// 作用：buildInitialNavigationStatus 函数。
inline NavigationStatusData buildInitialNavigationStatus(
		const InitialParameterConfig& init_cfg,
		bool has_first_gnss,
		const GnssData& first_gnss) {
	NavigationStatusData nav_init;
	std::array<double, 3> init_blh_rad = init_cfg.init_blh;
	if (init_cfg.init_blh_in_degree) {
		init_blh_rad[0] = deg2rad(init_blh_rad[0]);
		init_blh_rad[1] = deg2rad(init_blh_rad[1]);
	}
	std::array<double, 3> init_euler_rad = init_cfg.init_euler;
	if (init_cfg.init_euler_in_degree) {
		init_euler_rad[0] = deg2rad(init_euler_rad[0]);
		init_euler_rad[1] = deg2rad(init_euler_rad[1]);
		init_euler_rad[2] = deg2rad(init_euler_rad[2]);
	}
	nav_init.pvapre_.blh = init_blh_rad;
	nav_init.pvapre_.vel_n = init_cfg.init_vel_n;
	nav_init.pvapre_.euler = init_euler_rad;
	if (init_cfg.use_first_gnss_position && has_first_gnss) {
		nav_init.pvapre_.blh = {first_gnss.latitude, first_gnss.longitude, first_gnss.altitude};
	}
	nav_init.pvacur_ = nav_init.pvapre_;
	nav_init.imuerror_ = init_cfg.init_imu_error;
	nav_init.imuerror_.gyro_bias = init_cfg.init_gyro_bias;
	nav_init.imuerror_.accel_bias = init_cfg.init_accel_bias;
	nav_init.imuerror_.gyro_scale = init_cfg.init_gyro_scale;
	nav_init.imuerror_.accel_scale = init_cfg.init_accel_scale;
	return nav_init;
}
// 作用：extractStd 函数。
inline std::array<double, 3> extractStd(const Matrix& p, std::size_t start_idx) {
	std::array<double, 3> out{0.0, 0.0, 0.0};
	for (std::size_t i = 0; i < 3; ++i) {
		const double v = p(start_idx + i, start_idx + i);
		out[i] = (v > 0.0) ? std::sqrt(v) : 0.0;
	}
	return out;
}
// 作用：printProgressBar 函数。
inline void printProgressBar(std::size_t processed, std::size_t total, std::size_t width) {
	if (width == 0) {
		width = 40;
	}
	if (total == 0) {
		const char spinner[4] = {'|', '/', '-', '\\'};
		const char c = spinner[processed % 4];
		std::cout << "\r[" << c << "] 锟窖达拷锟斤拷IMU锟斤拷锟斤拷: " << processed << std::flush;
		return;
	}
	if (processed > total) {
		processed = total;
	}
	const double ratio = static_cast<double>(processed) / static_cast<double>(total);
	const std::size_t filled = static_cast<std::size_t>(ratio * static_cast<double>(width));
	std::cout << "\r[";
	for (std::size_t i = 0; i < width; ++i) {
		std::cout << ((i < filled) ? '#' : '.');
	}
	std::cout << "] " << std::fixed << std::setprecision(1) << (ratio * 100.0)
			  << "% (" << processed << "/" << total << ")" << std::flush;
}
// 作用：RunGnssInsMain 函数。
inline int RunGnssInsMain(const MainProgramConfig& cfg = MainProgramConfig()) {
	const std::string imu_file_path = resolveImuFilePath(cfg);
	ImuFileLoader imu_loader;
	if (!imu_loader.open(imu_file_path)) {
		throw std::runtime_error("没锟斤拷imu锟侥硷拷");
	}
	GnssFileLoader gnss_loader;
	const bool gnss_open_ok = gnss_loader.open(cfg.gnss_file_path);
	FileSaver saver;
	if (!saver.openDefaultFiles(cfg.output_dir)) {
		throw std::runtime_error("锟津开斤拷锟斤拷锟斤拷锟侥硷拷失锟斤拷");
	}
	ImuData imu_pre;
	ImuData imu_cur;
	if (!imu_loader.readNext(imu_pre)) {
		throw std::runtime_error("没锟斤拷imu锟侥硷拷");
	}
	if (!imu_loader.readNext(imu_cur)) {
		throw std::runtime_error("没锟斤拷imu锟侥硷拷");
	}
	const double imu_file_t0 = imu_pre.time;
	const double start_time = (cfg.process_start_s < 0.0) ? imu_file_t0 : cfg.process_start_s;
	if (cfg.process_end_s >= 0.0 && cfg.process_end_s < start_time) {
		throw std::runtime_error("时锟戒窗锟斤拷锟斤拷锟斤拷?? end锟斤拷锟斤拷锟斤拷诘锟斤拷锟絪tart锟斤拷锟斤拷锟斤拷锟斤拷??1");
	}
	while (imu_cur.time < start_time) {
		imu_pre = imu_cur;
		if (!imu_loader.readNext(imu_cur)) {
			throw std::runtime_error("时锟戒窗锟斤拷愠拷锟絀MU锟侥硷拷锟斤拷围");
		}
	}
	GnssData gnss_next;
	bool has_next_gnss = false;
	if (gnss_open_ok) {
		has_next_gnss = gnss_loader.readNext(gnss_next);
		if (has_next_gnss) {
			normalizeGnssAngleToRadIfNeeded(gnss_next);
		}
	}
	const NavigationStatusData nav_init = buildInitialNavigationStatus(
			cfg.init,
			has_next_gnss,
			gnss_next);
	GnssInsFusionWorkflow workflow;
	// 作用：x0 函数。
	Matrix x0(ErrorStateIndex21::kDim, 1, 0.0);
	const Matrix p0 = buildInitialCovarianceP0(cfg.init);
	workflow.initialize(nav_init, x0, p0);
	const std::size_t total_records = estimateImuRecordCount(imu_file_path);
	std::size_t processed = 0;
	while (true) {
		ImuMeasureData imu_measure;
		imu_measure.imupre_ = imu_pre;
		imu_measure.imucur_ = imu_cur;
		if (cfg.process_end_s >= 0.0 && imu_cur.time > cfg.process_end_s) {
			break;
		}
		auto runProcessStep = [&](const ImuMeasureData& step_imu,
								bool has_gnss_step,
								const GnssData& gnss_step) {
			double step_dt = step_imu.imucur_.time - step_imu.imupre_.time;
			if (step_dt < 1e-6) {
				step_dt = 1e-6;
			}
			
			// 1. 获取k-1时刻状态（传播前的状态）
			const NavigationStatusData nav_km1 = workflow.navState();

			// 2. 使用k-1时刻状态构建F矩阵（误差状态方程在名义轨迹附近线性化）
			const Matrix f = buildAlignedTransitionF(
				step_dt,
				nav_km1,  // 使用k-1时刻状态
				step_imu,
				cfg.corrtime_gb,
				cfg.corrtime_ab,
				cfg.corrtime_gs,
				cfg.corrtime_as);

			// 3. 准备传播：将k-1时刻的当前状态设为传播的起始状态（即新的pvapre_）
			NavigationStatusData nav_after_mech = nav_km1;
			nav_after_mech.pvapre_ = nav_km1.pvacur_; // 修正：propagateIns使用pvapre_作为起始点，必须设为k-1

			// 4. 执行INS机械编排（速度/位置/姿态传播），获取k时刻状态
			propagateIns(step_imu, nav_after_mech);
			// 此时 nav_after_mech.pvapre_ 是 k-1, nav_after_mech.pvacur_ 是 k

			// 5. 使用k时刻状态构建Q矩阵（梯形积分需要k-1和k时刻的G矩阵）
			// nav_after_mech.pvapre_ 是 k-1, nav_after_mech.pvacur_ 是 k
			const Matrix q = buildAlignedProcessNoiseQ(
								step_dt,
								f,
								step_imu,
								nav_after_mech,  // 使用包含k-1和k时刻的状态
				cfg.init,
				cfg.corrtime_gb,
				cfg.corrtime_ab,
				cfg.corrtime_gs,
				cfg.corrtime_as);

			// 6. 更新workflow状态并执行EKF predict + update
			workflow.setNavState(nav_after_mech);
			workflow.processStep(
				step_imu,
				f,
				q,
				has_gnss_step,
				gnss_step,
				cfg.antenna_lever_arm_b);
		};
		GnssData gnss_use;
		bool has_gnss_for_step = false;
		while (has_next_gnss && gnss_next.time <= imu_cur.time + cfg.gnss_sync_tolerance) {
			gnss_use = gnss_next;
			has_gnss_for_step = true;
			has_next_gnss = gnss_loader.readNext(gnss_next);
			if (has_next_gnss) {
				normalizeGnssAngleToRadIfNeeded(gnss_next);
			}
		}
		if (has_gnss_for_step) {
			const fusion::GnssImuTimeDecision decision = fusion::decideGnssImuTimeRelation(
				gnss_use.time,
				imu_pre.time,
				imu_cur.time,
				cfg.gnss_sync_tolerance);
std::cout << "\nGNSS Time: " << gnss_use.time << " IMU cur: " << imu_cur.time << " Action: " << (int)decision.action << std::endl;
			if (decision.action == fusion::UpdateAction::kUpdateOldImu) {
				workflow.processGnssUpdate(gnss_use, cfg.antenna_lever_arm_b);
				runProcessStep(imu_measure, false, gnss_use);
			} else if (decision.action == fusion::UpdateAction::kUpdateNewImu) {
				runProcessStep(imu_measure, true, gnss_use);
			} else if (decision.action == fusion::UpdateAction::kInterpolateAndUpdate) {
				double ratio = (gnss_use.time - imu_pre.time) / (imu_cur.time - imu_pre.time);
				if (ratio < 0.0) ratio = 0.0;
				if (ratio > 1.0) ratio = 1.0;
				ImuData imu_mid = imu_cur;
				imu_mid.time = gnss_use.time;
				imu_mid.dtheta_x = imu_cur.dtheta_x * ratio;
				imu_mid.dtheta_y = imu_cur.dtheta_y * ratio;
				imu_mid.dtheta_z = imu_cur.dtheta_z * ratio;
				imu_mid.dvel_x   = imu_cur.dvel_x * ratio;
				imu_mid.dvel_y   = imu_cur.dvel_y * ratio;
				imu_mid.dvel_z   = imu_cur.dvel_z * ratio;
				ImuData imu_cur_rem = imu_cur;
				imu_cur_rem.dtheta_x = imu_cur.dtheta_x * (1.0 - ratio);
				imu_cur_rem.dtheta_y = imu_cur.dtheta_y * (1.0 - ratio);
				imu_cur_rem.dtheta_z = imu_cur.dtheta_z * (1.0 - ratio);
				imu_cur_rem.dvel_x   = imu_cur.dvel_x * (1.0 - ratio);
				imu_cur_rem.dvel_y   = imu_cur.dvel_y * (1.0 - ratio);
				imu_cur_rem.dvel_z   = imu_cur.dvel_z * (1.0 - ratio);
				ImuMeasureData step1_imu;
				step1_imu.imupre_ = imu_pre;
				step1_imu.imucur_ = imu_mid;
				runProcessStep(step1_imu, true, gnss_use);
				ImuMeasureData step2_imu;
				step2_imu.imupre_ = imu_mid;
				step2_imu.imucur_ = imu_cur_rem;
				runProcessStep(step2_imu, false, gnss_use);
			} else {
				runProcessStep(imu_measure, false, gnss_use);
			}
		} else {
			runProcessStep(imu_measure, false, gnss_use);
		}
		const NavigationStatusData& nav = workflow.navState();
		const Matrix& p = workflow.kalman().covariance();
		const std::array<double, 3> pos_std = extractStd(p, ErrorStateIndex21::kPos);
		const std::array<double, 3> vel_std = extractStd(p, ErrorStateIndex21::kVel);
		const std::array<double, 3> att_std = extractStd(p, ErrorStateIndex21::kAtt);
		const std::array<double, 3> gb_std = extractStd(p, ErrorStateIndex21::kGyroBias);
		const std::array<double, 3> ab_std = extractStd(p, ErrorStateIndex21::kAccelBias);
		const std::array<double, 3> gs_std = extractStd(p, ErrorStateIndex21::kGyroScale);
		const std::array<double, 3> as_std = extractStd(p, ErrorStateIndex21::kAccelScale);
		if (!saver.writeTruthLine(imu_cur.time, nav.pvacur_)) {
			throw std::runtime_error("写锟斤拷Navresult失锟斤拷");
		}
		if (!saver.writeImuErrorLine(imu_cur.time, nav.imuerror_)) {
			throw std::runtime_error("写锟斤拷Imu_Error失锟斤拷");
		}
		if (!saver.writeStateStdLine(imu_cur.time, pos_std, vel_std, att_std, gb_std, ab_std, gs_std, as_std)) {
			throw std::runtime_error("写锟斤拷STD_result失锟斤拷");
		}
		
		// 检查协方差对角线元素是否都为正
		auto checkCov = [&p]() {
			for (int i = 0; i < p.rows(); ++i) {
				if (p(i, i) < 0.0) {
					std::cerr << "警告: 协方差对角线元素 [" << i << ", " << i << "] 为负数: " << p(i, i) << std::endl;
					return false;
				}
			}
			return true;
		};
		checkCov();

		++processed;
		if (processed % cfg.progress_update_step == 0 || processed == 1) {
			printProgressBar(processed, total_records, cfg.progress_bar_width);
		}
		
		// 更新上一时刻的状态和 IMU 数据
		imu_pre = imu_cur;
		if (!imu_loader.readNext(imu_cur)) {
			break;
		}
	}
	printProgressBar((total_records > 0) ? total_records : processed,
					 total_records,
					 cfg.progress_bar_width);
	std::cout << "\n锟斤拷系锟斤拷锟斤拷锟斤拷锟斤拷锟缴ｏ拷锟斤拷锟斤拷驯锟斤拷锟? "
			  << FileSaver::kNavResultFileName << ", "
			  << FileSaver::kImuErrorFileName << ", "
			  << FileSaver::kStdResultFileName << std::endl;
	if (!gnss_open_ok) {
		std::cout << "锟斤拷示" << std::endl;
	}
	return 0;
}
// 作用：printMainUsage 函数。
inline void printMainUsage(const char* exe_name) {
	std::cout << "锟矫凤拷: " << exe_name
			  << " [--imu <imu_file>] [--gnss <gnss_file>] [--out <output_dir>]"
			  << " [--sync <seconds>]"
			  << " [--corrtime-gb <seconds>] [--corrtime-ab <seconds>]"
			  << " [--corrtime-gs <seconds>] [--corrtime-as <seconds>]"
			  << " [--start <imu_time>] [--end <imu_time|-1>] [--help]"
			  << std::endl;
}
// 作用：ParseMainConfigFromArgs 函数。
inline MainProgramConfig ParseMainConfigFromArgs(int argc, char** argv) {
	MainProgramConfig cfg;
	for (int i = 1; i < argc; ++i) {
		const std::string arg = argv[i];
		if (arg == "--help" || arg == "-h") {
			printMainUsage((argc > 0 && argv[0] != nullptr) ? argv[0] : "program");
			throw std::runtime_error("help requested");
		}
		auto needValue = [&](const char* opt_name) -> std::string {
			if (i + 1 >= argc) {
				throw std::runtime_error(std::string("锟斤拷锟斤拷缺锟斤拷取?? ") + opt_name);
			}
			++i;
			return std::string(argv[i]);
		};
		if (arg == "--imu") {
			cfg.imu_file_path = needValue("--imu");
		} else if (arg == "--gnss") {
			cfg.gnss_file_path = needValue("--gnss");
		} else if (arg == "--out") {
			cfg.output_dir = needValue("--out");
		} else if (arg == "--sync") {
			cfg.gnss_sync_tolerance = std::stod(needValue("--sync"));
		} else if (arg == "--corrtime") {
			const double tau = std::stod(needValue("--corrtime"));
			cfg.corrtime_gb = tau;
			cfg.corrtime_ab = tau;
			cfg.corrtime_gs = tau;
			cfg.corrtime_as = tau;
		} else if (arg == "--corrtime-gb") {
			cfg.corrtime_gb = std::stod(needValue("--corrtime-gb"));
		} else if (arg == "--corrtime-ab") {
			cfg.corrtime_ab = std::stod(needValue("--corrtime-ab"));
		} else if (arg == "--corrtime-gs") {
			cfg.corrtime_gs = std::stod(needValue("--corrtime-gs"));
		} else if (arg == "--corrtime-as") {
			cfg.corrtime_as = std::stod(needValue("--corrtime-as"));
		} else if (arg == "--start") {
			cfg.process_start_s = std::stod(needValue("--start"));
		} else if (arg == "--end") {
			cfg.process_end_s = std::stod(needValue("--end"));
		} else {
			throw std::runtime_error("未知错误: " + arg);
		}
	}
	return cfg;
}
}  // namespace ins

// 作用：main 函数。
int main(int argc, char** argv) {
	try {
		const ins::MainProgramConfig cfg = ins::ParseMainConfigFromArgs(argc, argv);
		return ins::RunGnssInsMain(cfg);
	} catch (const std::exception& e) {
		if (std::strcmp(e.what(), "help requested") == 0) {
			return 0;
		}
		std::cerr << "参数解析失败: " << e.what() << std::endl;
		ins::printMainUsage((argc > 0 && argv[0] != nullptr) ? argv[0] : "program");
		return 1;
	}
}
