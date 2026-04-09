#pragma once

#include <array>
#include <cmath>

// 作用：定义 ImuData 数据结构。
struct ImuData {
	double time = 0.0;

	// 作用：此处说明当前字段或步骤的用途。
	double dtheta_x = 0.0;
	double dtheta_y = 0.0;
	double dtheta_z = 0.0;

	// 作用：此处说明当前字段或步骤的用途。
	double dvel_x = 0.0;
	double dvel_y = 0.0;
	double dvel_z = 0.0;

	// 作用：isFinite 函数。
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

// 作用：此处说明当前字段或步骤的用途。
	// 作用：此处说明当前字段或步骤的用途。
// 作用：定义 ImuMeasureData 数据结构。
struct ImuMeasureData {
	ImuData imupre_;
	ImuData imucur_;

	// 作用：shiftAndSetCurrent 函数。
	void shiftAndSetCurrent(const ImuData& imu_cur) {
		imupre_ = imucur_;
		imucur_ = imu_cur;
	}
};

// 作用：定义 PvaData 数据结构。
struct PvaData {
	// 作用：此处说明当前字段或步骤的用途。
	std::array<double, 3> blh = {0.0, 0.0, 0.0};
	// 作用：此处说明当前字段或步骤的用途。
	std::array<double, 3> vel_n = {0.0, 0.0, 0.0};
	// 作用：此处说明当前字段或步骤的用途。
	std::array<double, 3> euler = {0.0, 0.0, 0.0};
};

// 作用：定义 ImuErrorData 数据结构。
struct ImuErrorData {
	// 作用：此处说明当前字段或步骤的用途。
	std::array<double, 3> gyro_bias = {0.0, 0.0, 0.0};
	std::array<double, 3> accel_bias = {0.0, 0.0, 0.0};

	// 作用：此处说明当前字段或步骤的用途。
	// 作用：此处说明当前字段或步骤的用途。
	std::array<double, 3> gyro_random_walk = {0.0, 0.0, 0.0};
	std::array<double, 3> accel_random_walk = {0.0, 0.0, 0.0};

	// 作用：此处说明当前字段或步骤的用途。
	std::array<double, 3> gyro_bias_instability = {0.0, 0.0, 0.0};
	std::array<double, 3> accel_bias_instability = {0.0, 0.0, 0.0};

	// 作用：此处说明当前字段或步骤的用途。
	std::array<double, 3> gyro_scale = {0.0, 0.0, 0.0};
	std::array<double, 3> accel_scale = {0.0, 0.0, 0.0};

	// 作用：此处说明当前字段或步骤的用途。
	std::array<double, 3> gyro_scale_instability = {0.0, 0.0, 0.0};
	std::array<double, 3> accel_scale_instability = {0.0, 0.0, 0.0};
};

// 作用：此处说明当前字段或步骤的用途。
// 作用：定义 NavigationStatusData 数据结构。
struct NavigationStatusData {
	PvaData pvapre_;
	PvaData pvacur_;
	ImuErrorData imuerror_;
};
