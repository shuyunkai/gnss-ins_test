#pragma once

#include <cmath>

namespace ins {

struct LocalPosition {
		double latitude = 0.0;
	double longitude = 0.0;
	double height = 0.0;
};

struct LocalVelocity {
		double vn = 0.0;
	double ve = 0.0;
	double vd = 0.0;
};

struct EarthParameters {
		double rm = 0.0;
		double rn = 0.0;
		double gravity = 0.0;

		double omega_ie = 0.0;
	double omega_en_north = 0.0;
	double omega_en_east = 0.0;
	double omega_en_down = 0.0;
};

struct Wgs84Constants {
		double a = 6378137.0;
		double f = 1.0 / 298.257223563;
		double omega_ie = 7.292115e-5;
		double ge = 9.7803253359;
		double k = 0.00193185265241;
};

inline EarthParameters updateEarthParameters(
		const LocalPosition& pos,
		const LocalVelocity& vel,
		const Wgs84Constants& wgs84 = Wgs84Constants()) {
	EarthParameters out;

	const double e2 = 2.0 * wgs84.f - wgs84.f * wgs84.f;
	const double sin_lat = std::sin(pos.latitude);
	const double sin_lat2 = sin_lat * sin_lat;
	const double denom = 1.0 - e2 * sin_lat2;
	const double sqrt_denom = std::sqrt(denom);

		out.rn = wgs84.a / sqrt_denom;
	out.rm = wgs84.a * (1.0 - e2) / (denom * sqrt_denom);

		const double g0 = wgs84.ge * (1.0 + wgs84.k * sin_lat2) / std::sqrt(1.0 - e2 * sin_lat2);
	out.gravity = g0 - 3.086e-6 * pos.height;

	out.omega_ie = wgs84.omega_ie;

	const double rm_h = out.rm + pos.height;
	const double rn_h = out.rn + pos.height;
	const double cos_lat = std::cos(pos.latitude);

		out.omega_en_north = vel.ve / rn_h;
	out.omega_en_east = -vel.vn / rm_h;
	out.omega_en_down = -vel.ve * std::tan(pos.latitude) / rn_h;

		(void)cos_lat;
	return out;
}

}  // 鍛藉悕绌洪棿缁撴潫
