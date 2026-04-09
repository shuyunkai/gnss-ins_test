#pragma once

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

#include "imu.h"

namespace ins {

class FileSaver {
public:
	static constexpr const char* kNavResultFileName = "Navresult.txt";
	static constexpr const char* kImuErrorFileName = "Imu_Error.txt";
	static constexpr const char* kStdResultFileName = "STD_result.txt";

	FileSaver() = default;

	~FileSaver() {
		closeAll();
	}

		bool openTruthFile(const std::string& file_path) {
		truth_stream_.close();
		truth_stream_.open(file_path, std::ios::out | std::ios::trunc);
		return truth_stream_.is_open();
	}

		bool openImuErrorFile(const std::string& file_path) {
		error_stream_.close();
		error_stream_.open(file_path, std::ios::out | std::ios::trunc);
		return error_stream_.is_open();
	}

		bool openStateStdFile(const std::string& file_path) {
		std_stream_.close();
		std_stream_.open(file_path, std::ios::out | std::ios::trunc);
		return std_stream_.is_open();
	}

				bool openDefaultFiles(const std::string& output_dir = "") {
		const std::string nav_path = buildPath(output_dir, kNavResultFileName);
		const std::string imu_error_path = buildPath(output_dir, kImuErrorFileName);
		const std::string std_path = buildPath(output_dir, kStdResultFileName);

		const bool ok_nav = openTruthFile(nav_path);
		const bool ok_imu = openImuErrorFile(imu_error_path);
		const bool ok_std = openStateStdFile(std_path);
		return ok_nav && ok_imu && ok_std;
	}

		void closeTruthFile() {
		if (truth_stream_.is_open()) {
			truth_stream_.close();
		}
	}

		void closeImuErrorFile() {
		if (error_stream_.is_open()) {
			error_stream_.close();
		}
	}

		void closeStateStdFile() {
		if (std_stream_.is_open()) {
			std_stream_.close();
		}
	}

		void closeAll() {
		closeTruthFile();
		closeImuErrorFile();
		closeStateStdFile();
	}

		bool isTruthFileOpen() const {
		return truth_stream_.is_open();
	}

		bool isImuErrorFileOpen() const {
		return error_stream_.is_open();
	}

		bool isStateStdFileOpen() const {
		return std_stream_.is_open();
	}

		// time(s), lat(deg), lon(deg), h(m), vn(m/s), ve(m/s), vd(m/s), roll(deg), pitch(deg), yaw(deg)
		bool writeTruthLine(double time_s, const PvaData& pva) {
		if (!truth_stream_.is_open()) {
			return false;
		}

		truth_stream_ << std::fixed << std::setprecision(10)
					  << time_s << ' '
					  << rad2deg(pva.blh[0]) << ' '
					  << rad2deg(pva.blh[1]) << ' '
					  << pva.blh[2] << ' '
					  << pva.vel_n[0] << ' '
					  << pva.vel_n[1] << ' '
					  << pva.vel_n[2] << ' '
					  << rad2deg(pva.euler[0]) << ' '
					  << rad2deg(pva.euler[1]) << ' '
					  << rad2deg(pva.euler[2]) << '\n';

		return truth_stream_.good();
	}

		// time(s), gyro_bias_xyz(deg/h), accel_bias_xyz(mGal), gyro_scale_xyz(ppm), accel_scale_xyz(ppm)
		bool writeImuErrorLine(double time_s, const ImuErrorData& imu_error) {
		if (!error_stream_.is_open()) {
			return false;
		}

		error_stream_ << std::fixed << std::setprecision(10)
					  << time_s << ' '
					  << radps2degh(imu_error.gyro_bias[0]) << ' '
					  << radps2degh(imu_error.gyro_bias[1]) << ' '
					  << radps2degh(imu_error.gyro_bias[2]) << ' '
					  << mps22mgal(imu_error.accel_bias[0]) << ' '
					  << mps22mgal(imu_error.accel_bias[1]) << ' '
					  << mps22mgal(imu_error.accel_bias[2]) << ' '
					  << ratio2ppm(imu_error.gyro_scale[0]) << ' '
					  << ratio2ppm(imu_error.gyro_scale[1]) << ' '
					  << ratio2ppm(imu_error.gyro_scale[2]) << ' '
					  << ratio2ppm(imu_error.accel_scale[0]) << ' '
					  << ratio2ppm(imu_error.accel_scale[1]) << ' '
					  << ratio2ppm(imu_error.accel_scale[2]) << '\n';

		return error_stream_.good();
	}

		// time(s), pos_std_ned(m), vel_std_ned(m/s), att_std_rpy(deg),
	// gyro_bias_std_xyz(deg/h), accel_bias_std_xyz(mGal), gyro_scale_std_xyz(ppm), accel_scale_std_xyz(ppm)
		bool writeStateStdLine(
			double time_s,
			const std::array<double, 3>& pos_std,
			const std::array<double, 3>& vel_std,
			const std::array<double, 3>& att_std_rad,
			const std::array<double, 3>& gyro_bias_std_radps,
			const std::array<double, 3>& accel_bias_std_mps2,
			const std::array<double, 3>& gyro_scale_std_ratio,
			const std::array<double, 3>& accel_scale_std_ratio) {
		if (!std_stream_.is_open()) {
			return false;
		}

		std_stream_ << std::fixed << std::setprecision(10)
					<< time_s << ' '
					<< pos_std[0] << ' ' << pos_std[1] << ' ' << pos_std[2] << ' '
					<< vel_std[0] << ' ' << vel_std[1] << ' ' << vel_std[2] << ' '
					<< rad2deg(att_std_rad[0]) << ' '
					<< rad2deg(att_std_rad[1]) << ' '
					<< rad2deg(att_std_rad[2]) << ' '
					<< radps2degh(gyro_bias_std_radps[0]) << ' '
					<< radps2degh(gyro_bias_std_radps[1]) << ' '
					<< radps2degh(gyro_bias_std_radps[2]) << ' '
					<< mps22mgal(accel_bias_std_mps2[0]) << ' '
					<< mps22mgal(accel_bias_std_mps2[1]) << ' '
					<< mps22mgal(accel_bias_std_mps2[2]) << ' '
					<< ratio2ppm(gyro_scale_std_ratio[0]) << ' '
					<< ratio2ppm(gyro_scale_std_ratio[1]) << ' '
					<< ratio2ppm(gyro_scale_std_ratio[2]) << ' '
					<< ratio2ppm(accel_scale_std_ratio[0]) << ' '
					<< ratio2ppm(accel_scale_std_ratio[1]) << ' '
					<< ratio2ppm(accel_scale_std_ratio[2]) << '\n';

		return std_stream_.good();
	}

private:
		static double rad2deg(double rad) {
		return rad * 180.0 / 3.14159265358979323846;
	}

		static double radps2degh(double radps) {
				return rad2deg(radps) * 3600.0;
	}

		static double mps22mgal(double mps2) {
		return mps2 * 1.0e5;
	}

		static double ratio2ppm(double ratio) {
		return ratio * 1.0e6;
	}

		static std::string buildPath(const std::string& output_dir, const std::string& file_name) {
		if (output_dir.empty()) {
			return file_name;
		}
		const char last = output_dir.back();
		if (last == '/' || last == '\\') {
			return output_dir + file_name;
		}
		return output_dir + "/" + file_name;
	}

private:
	std::ofstream truth_stream_;
	std::ofstream error_stream_;
	std::ofstream std_stream_;
};

}  // namespace ins
