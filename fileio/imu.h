#pragma once

#include <array>
#include <cmath>

// 传感器输入：IMU 增量数据结构 (角速度与加速度增量)
struct ImuData {
	double time = 0.0;

		double dtheta_x = 0.0;
	double dtheta_y = 0.0;
	double dtheta_z = 0.0;

		double dvel_x = 0.0;
	double dvel_y = 0.0;
	double dvel_z = 0.0;

		bool isFinite() const {
		return std::isfinite(time) &&
			   std::isfinite(dtheta_x) &&
			   std::isfinite(dtheta_y) &&
			   std::isfinite(dtheta_z) &&
			   std::isfinite(dvel_x) &&
			   std::isfinite(dvel_y) &&
			   std::isfinite(dvel_z);
	}
};

	struct ImuMeasureData {
	ImuData imupre_;
	ImuData imucur_;

		void shiftAndSetCurrent(const ImuData& imu_cur) {
		imupre_ = imucur_;
		imucur_ = imu_cur;
	}
};

// 导航状态：瞬时位置、速度、姿态结构 (Pos/Vel/Att)
struct PvaData {
		std::array<double, 3> blh = {0.0, 0.0, 0.0};
		std::array<double, 3> vel_n = {0.0, 0.0, 0.0};
		std::array<double, 3> euler = {0.0, 0.0, 0.0};
};

struct ImuErrorData {
	std::array<double, 3> gyro_bias = {0.0, 0.0, 0.0};
	std::array<double, 3> accel_bias = {0.0, 0.0, 0.0};

	std::array<double, 3> gyro_random_walk = {0.0, 0.0, 0.0};
	std::array<double, 3> accel_random_walk = {0.0, 0.0, 0.0};

	std::array<double, 3> gyro_scale = {0.0, 0.0, 0.0};
	std::array<double, 3> accel_scale = {0.0, 0.0, 0.0};
};

struct NavigationStatusData {
	PvaData pvapre_;
	PvaData pvacur_;
	ImuErrorData imuerror_;
};
