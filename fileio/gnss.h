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

		bool isFinite() const {
			return std::isfinite(time) &&
				   std::isfinite(longitude) &&
				   std::isfinite(latitude) &&
				   std::isfinite(altitude) &&
				   std::isfinite(std_x) &&
				   std::isfinite(std_y) &&
				   std::isfinite(std_z);

	}
};
