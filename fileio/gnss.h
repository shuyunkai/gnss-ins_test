#pragma once

#include <cmath>

struct GnssData {
	double time = 0.0;
	double longitude = 0.0;
	double latitude = 0.0;
	double altitude = 0.0;

		double std_x = 0.0;
	double std_y = 0.0;
	double std_z = 0.0;

		bool has_velocity = false;
	double vel_n = 0.0;
	double vel_e = 0.0;
	double vel_d = 0.0;

		double std_vn = 0.0;
	double std_ve = 0.0;
	double std_vd = 0.0;

		bool isFinite() const {
		return std::isfinite(time) &&
			   std::isfinite(longitude) &&
			   std::isfinite(latitude) &&
			   std::isfinite(altitude) &&
			   std::isfinite(std_x) &&
			   std::isfinite(std_y) &&
			   std::isfinite(std_z) &&
			   std::isfinite(vel_n) &&
			   std::isfinite(vel_e) &&
			   std::isfinite(vel_d) &&
			   std::isfinite(std_vn) &&
			   std::isfinite(std_ve) &&
			   std::isfinite(std_vd);
	}
};

struct GnssMeasureData {
	GnssData gnssdata_;
};
