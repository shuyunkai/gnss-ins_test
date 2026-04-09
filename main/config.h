#pragma once

#include <string>
#include <array>
#include <cmath>
#include "../common/Coordinate_2.h"
#include "../fileio/imufileloader.h"

namespace ins {

// =========================================================================
// 核心参数配置文件 (InitialParameterConfig / MainProgramConfig)
// 所有可调参的过程噪声、初始零偏、解算起止时间等都在这里集中维护
// =========================================================================

// 作用：degPerHourToRadPerSec 函数 (提供参数转换)
inline double degPerHourToRadPerSec(double deg_per_hour) {
	return deg2rad(deg_per_hour) / 3600.0;
}
inline double mGalToMps2(double mgal) {
	return mgal * 1e-5;
}
inline double ppmToRatio(double ppm) {
	return ppm * 1e-6;
}

// 作用：定义 InitialParameterConfig 数据结构 (初始化先验参数集合)
struct InitialParameterConfig {
	bool use_first_gnss_position = false; // 是否强行使用GNSS的第一帧定位起步

	// 初始位置: 纬度(deg), 经度(deg), 高度(m)
	std::array<double, 3> init_blh = {30.4653881320, 114.4694083874, 23.821};
	bool init_blh_in_degree = true;

	// 初始速度: 北向、东向、地向 (m/s)
	std::array<double, 3> init_vel_n = {-0.04006, 10.91736, 0.05754};

	// 初始姿态: Roll, Pitch, Yaw (deg)
	std::array<double, 3> init_euler = {0.55606   ,-0.79210 , 91.21104};
	bool init_euler_in_degree = true;

	ImuErrorData init_imu_error = [] {
		ImuErrorData v;
		v.gyro_random_walk = {0.2, 0.2, 0.2};
		v.accel_random_walk = {0.2, 0.2, 0.2};
		return v;
	}();

	// --------------------- 初始零偏与比例因子收敛真值 ---------------------
	std::array<double, 3> init_gyro_bias = {-1372.74702497, -477.57869162, 665.15948298};  // deg/h
	std::array<double, 3> init_accel_bias = {3700.05593095, -12971.62640788, 2489.12090896}; // mGal
	std::array<double, 3> init_gyro_scale = {325.28151451, -145.75970550, -395.30985879};    // ppm
	std::array<double, 3> init_accel_scale = {1651.17392773, 614.16423855, 788.58869356};   // ppm

	// --------------------- 初始卡尔曼滤波初始误差置信度(P0) ---------------------
	std::array<double, 3> pos_std = {0.01, 0.01, 0.05};            // 位置初始标准差 (m)
	std::array<double, 3> vel_std = {0.01, 0.01, 0.01};            // 速度初始标准差 (m/s)
	std::array<double, 3> att_std = {0.01, 0.01, 0.01};            // 姿态角初始标准差 (deg)

	// --------------------- 初始零偏与比例因子置信度 ---------------------
	std::array<double, 3> gyro_bias_std = {5, 5, 5};               // deg/h
	std::array<double, 3> accel_bias_std = {20, 20, 20};           // mGal
	std::array<double, 3> gyro_scale_std = {50, 50, 50};           // ppm
	std::array<double, 3> accel_scale_std = {50, 50, 50};          // ppm

	// --------------------- 滤波过程噪声 (时间游走) ---------------------
	std::array<double, 3> gyro_bias_process_std = {200, 200, 200};          // 陀螺仪零偏不稳定性 (deg/h/sqrt(s))
	std::array<double, 3> accel_bias_process_std = {1000, 1000, 200};       // 加计零偏不稳定性 (mGal/sqrt(s))
	std::array<double, 3> gyro_scale_process_std = {1000, 1000, 1000};      // 陀螺仪比例过程噪声 (ppm/sqrt(s))
	std::array<double, 3> accel_scale_process_std = {1000, 1000, 1000};     // 加计比例过程噪声 (ppm/sqrt(s))
};


// 作用：定义 MainProgramConfig 数据结构 (全局运行相关配置)
struct MainProgramConfig {
	std::string imu_file_path = "data/ICM20602.txt";   // 输入的 IMU 二进制/文本路径
	std::string gnss_file_path = "data/GNSS_RTK.txt";  // 输入的 GNSS 测量信息
	std::string output_dir = "data";                   // 输出结算与标准差结果文件夹

	double corrtime_gb = 3600.0;                       // 陀螺一阶马尔科夫过程相关时间 
	double corrtime_ab = 3600.0;                       // 加计一阶马尔科夫过程相关时间
	double corrtime_gs = 3600.0;                       // 陀螺比例因子过程相关时间
	double corrtime_as = 3600.0;                       // 加计比例因子过程相关时间
	
	std::array<double, 3> antenna_lever_arm_b = {-0.073, 0.302, 0.087}; // 天线臂长杆补偿

	double process_start_s = 357599;                   // 【注意核心】GNSS 解算处理启动时间 (秒)
	double process_end_s = -1.0;                       // GNSS 解算处理终止时间 (秒)
	
	double gnss_sync_tolerance = 1e-3;                 // 传感器双时间戳容限 (ms)
	std::size_t progress_bar_width = 40;               // 进度条渲染长度
	std::size_t progress_update_step = 5000;           // UI渲染触发间隔步数
	
	InitialParameterConfig init;                       // 包含上一结构体的综合配置对象
};

} // namespace ins
