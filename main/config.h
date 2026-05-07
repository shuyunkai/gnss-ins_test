#pragma once

#include <string>
#include <array>
#include <cmath>
#include "../common/Coordinate_2.h"
#include "../fileio/imufileloader.h"

namespace ins {

// =========================================================================
// 参数数据结构 (所有参数值从 gnss_ins.yaml 加载，此处仅定义结构)
// =========================================================================

inline double degPerHourToRadPerSec(double deg_per_hour) {
	return deg2rad(deg_per_hour) / 3600.0;
}
inline double degPerSqrtHourToRadPerSqrtSec(double val) {
	return deg2rad(val) / 60.0;
}
inline double mpsPerSqrtHourToMps2PerSqrtHz(double val) {
	return val / 60.0;
}
inline double mGalToMps2(double mgal) {
	return mgal * 1e-5;
}
inline double ppmToRatio(double ppm) {
	return ppm * 1e-6;
}

struct InitialParameterConfig {
	bool use_first_gnss_position;
	std::array<double, 3> init_blh;
	bool init_blh_in_degree;
	std::array<double, 3> init_vel_n;
	std::array<double, 3> init_euler;
	bool init_euler_in_degree;
	ImuErrorData init_imu_error;
	std::array<double, 3> init_gyro_bias;
	std::array<double, 3> init_accel_bias;
	std::array<double, 3> init_gyro_scale;
	std::array<double, 3> init_accel_scale;
	std::array<double, 3> pos_std;
	std::array<double, 3> vel_std;
	std::array<double, 3> att_std;
	std::array<double, 3> gyro_bias_std;
	std::array<double, 3> accel_bias_std;
	std::array<double, 3> gyro_scale_std;
	std::array<double, 3> accel_scale_std;
	std::array<double, 3> gyro_bias_process_std;
	std::array<double, 3> accel_bias_process_std;
	std::array<double, 3> gyro_scale_process_std;
	std::array<double, 3> accel_scale_process_std;
};

enum class FilterAlgorithm {
	ExtendedKalman,
	UnscentedKalman,
	GraphOptimization
};

struct MainProgramConfig {
	FilterAlgorithm algorithm;
	std::string imu_file_path;
	std::string gnss_file_path;
	std::string output_dir;
	double corrtime_gb;
	double corrtime_ab;
	double corrtime_gs;
	double corrtime_as;
	std::array<double, 3> antenna_lever_arm_b;
	double process_start_s;
	double process_end_s;
	double gnss_sync_tolerance;
	double ukf_alpha;
	double ukf_kappa;
	double ukf_beta;
	std::size_t progress_bar_width;
	std::size_t progress_update_step;
	InitialParameterConfig init;
};

} // namespace ins
