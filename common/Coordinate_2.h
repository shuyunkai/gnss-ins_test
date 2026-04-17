#pragma once

#include <algorithm>
#include <array>
#include <cmath>

namespace ins {

struct Vec3 {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

struct Mat3 {
	double m[3][3] = {{0.0, 0.0, 0.0},
						 {0.0, 0.0, 0.0},
						 {0.0, 0.0, 0.0}};
};

struct Quaternion {
		double w = 1.0;
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

struct Euler {
		double roll = 0.0;
	double pitch = 0.0;
	double yaw = 0.0;
};

struct Blh {
		double lat = 0.0;
	double lon = 0.0;
	double h = 0.0;
};

struct Wgs84 {
	double a = 6378137.0;
	double f = 1.0 / 298.257223563;
	double omega_ie = 7.292115e-5;
};

inline double deg2rad(double deg) {
	return deg * 3.14159265358979323846 / 180.0;
}

inline double rad2deg(double rad) {
	return rad * 180.0 / 3.14159265358979323846;
}

inline Mat3 identity3() {
	Mat3 r;
	r.m[0][0] = 1.0;
	r.m[1][1] = 1.0;
	r.m[2][2] = 1.0;
	return r;
}

inline Mat3 transpose(const Mat3& a) {
	Mat3 r;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			r.m[i][j] = a.m[j][i];
		}
	}
	return r;
}

inline Vec3 matVec(const Mat3& a, const Vec3& v) {
	Vec3 r;
	r.x = a.m[0][0] * v.x + a.m[0][1] * v.y + a.m[0][2] * v.z;
	r.y = a.m[1][0] * v.x + a.m[1][1] * v.y + a.m[1][2] * v.z;
	r.z = a.m[2][0] * v.x + a.m[2][1] * v.y + a.m[2][2] * v.z;
	return r;
}

inline Mat3 matMul(const Mat3& a, const Mat3& b) {
	Mat3 r;
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			r.m[i][j] = 0.0;
			for (int k = 0; k < 3; ++k) {
				r.m[i][j] += a.m[i][k] * b.m[k][j];
			}
		}
	}
	return r;
}

inline double vecNorm(const Vec3& v) {
	return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

inline Quaternion normalizeQuaternion(const Quaternion& q) {
	const double n = std::sqrt(q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z);
	if (n <= 0.0) {
		return Quaternion{};
	}
	Quaternion r;
	r.w = q.w / n;
	r.x = q.x / n;
	r.y = q.y / n;
	r.z = q.z / n;
	return r;
}

inline Mat3 euler2dcm(const Euler& e) {
	const double cr = std::cos(e.roll);
	const double sr = std::sin(e.roll);
	const double cp = std::cos(e.pitch);
	const double sp = std::sin(e.pitch);
	const double cy = std::cos(e.yaw);
	const double sy = std::sin(e.yaw);

		Mat3 c;
	c.m[0][0] = cy * cp;
	c.m[0][1] = cy * sp * sr - sy * cr;
	c.m[0][2] = cy * sp * cr + sy * sr;
	c.m[1][0] = sy * cp;
	c.m[1][1] = sy * sp * sr + cy * cr;
	c.m[1][2] = sy * sp * cr - cy * sr;
	c.m[2][0] = -sp;
	c.m[2][1] = cp * sr;
	c.m[2][2] = cp * cr;
	return c;
}

inline Euler dcm2euler(const Mat3& c) {
	Euler e;
	e.pitch = std::asin(std::max(-1.0, std::min(1.0, -c.m[2][0])));
	e.roll = std::atan2(c.m[2][1], c.m[2][2]);
	e.yaw = std::atan2(c.m[1][0], c.m[0][0]);
	return e;
}

inline Quaternion euler2quat(const Euler& e) {
	const double cr2 = std::cos(0.5 * e.roll);
	const double sr2 = std::sin(0.5 * e.roll);
	const double cp2 = std::cos(0.5 * e.pitch);
	const double sp2 = std::sin(0.5 * e.pitch);
	const double cy2 = std::cos(0.5 * e.yaw);
	const double sy2 = std::sin(0.5 * e.yaw);

	Quaternion q;
	q.w = cy2 * cp2 * cr2 + sy2 * sp2 * sr2;
	q.x = cy2 * cp2 * sr2 - sy2 * sp2 * cr2;
	q.y = cy2 * sp2 * cr2 + sy2 * cp2 * sr2;
	q.z = sy2 * cp2 * cr2 - cy2 * sp2 * sr2;
		return normalizeQuaternion(q);
}

inline Mat3 quat2dcm(const Quaternion& q_in) {
	const Quaternion q = normalizeQuaternion(q_in);

	const double ww = q.w * q.w;
	const double xx = q.x * q.x;
	const double yy = q.y * q.y;
	const double zz = q.z * q.z;

	Mat3 c;
	c.m[0][0] = ww + xx - yy - zz;
	c.m[0][1] = 2.0 * (q.x * q.y - q.w * q.z);
	c.m[0][2] = 2.0 * (q.x * q.z + q.w * q.y);

	c.m[1][0] = 2.0 * (q.x * q.y + q.w * q.z);
	c.m[1][1] = ww - xx + yy - zz;
	c.m[1][2] = 2.0 * (q.y * q.z - q.w * q.x);

	c.m[2][0] = 2.0 * (q.x * q.z - q.w * q.y);
	c.m[2][1] = 2.0 * (q.y * q.z + q.w * q.x);
	c.m[2][2] = ww - xx - yy + zz;
	return c;
}

inline Quaternion dcm2quat(const Mat3& c) {
	Quaternion q;
	const double tr = c.m[0][0] + c.m[1][1] + c.m[2][2];
	if (tr > 0.0) {
		const double s = std::sqrt(tr + 1.0) * 2.0;
		q.w = 0.25 * s;
		q.x = (c.m[2][1] - c.m[1][2]) / s;
		q.y = (c.m[0][2] - c.m[2][0]) / s;
		q.z = (c.m[1][0] - c.m[0][1]) / s;
	} else if ((c.m[0][0] > c.m[1][1]) && (c.m[0][0] > c.m[2][2])) {
		const double s = std::sqrt(1.0 + c.m[0][0] - c.m[1][1] - c.m[2][2]) * 2.0;
		q.w = (c.m[2][1] - c.m[1][2]) / s;
		q.x = 0.25 * s;
		q.y = (c.m[0][1] + c.m[1][0]) / s;
		q.z = (c.m[0][2] + c.m[2][0]) / s;
	} else if (c.m[1][1] > c.m[2][2]) {
		const double s = std::sqrt(1.0 + c.m[1][1] - c.m[0][0] - c.m[2][2]) * 2.0;
		q.w = (c.m[0][2] - c.m[2][0]) / s;
		q.x = (c.m[0][1] + c.m[1][0]) / s;
		q.y = 0.25 * s;
		q.z = (c.m[1][2] + c.m[2][1]) / s;
	} else {
		const double s = std::sqrt(1.0 + c.m[2][2] - c.m[0][0] - c.m[1][1]) * 2.0;
		q.w = (c.m[1][0] - c.m[0][1]) / s;
		q.x = (c.m[0][2] + c.m[2][0]) / s;
		q.y = (c.m[1][2] + c.m[2][1]) / s;
		q.z = 0.25 * s;
	}
		return normalizeQuaternion(q);
}

inline Euler quat2euler(const Quaternion& q) {
		return dcm2euler(quat2dcm(q));
}

inline Quaternion rotvec2quat(const Vec3& rv) {
	const double angle = vecNorm(rv);
	if (angle < 1e-12) {
		Quaternion q;
		q.w = 1.0;
		q.x = 0.5 * rv.x;
		q.y = 0.5 * rv.y;
		q.z = 0.5 * rv.z;
				return normalizeQuaternion(q);
	}

	const double half = 0.5 * angle;
	const double s = std::sin(half) / angle;
	Quaternion q;
	q.w = std::cos(half);
	q.x = rv.x * s;
	q.y = rv.y * s;
	q.z = rv.z * s;
		return normalizeQuaternion(q);
}

inline Vec3 quat2rotvec(const Quaternion& q_in) {
	const Quaternion q = normalizeQuaternion(q_in);
	const double sin_half = std::sqrt(q.x * q.x + q.y * q.y + q.z * q.z);

	if (sin_half < 1e-12) {
		return Vec3{2.0 * q.x, 2.0 * q.y, 2.0 * q.z};
	}

	const double half_angle = std::atan2(sin_half, q.w);
	const double scale = 2.0 * half_angle / sin_half;
	return Vec3{q.x * scale, q.y * scale, q.z * scale};
}

inline Mat3 rotvec2dcm(const Vec3& rv) {
		return quat2dcm(rotvec2quat(rv));
}

inline Vec3 dcm2rotvec(const Mat3& c) {
		return quat2rotvec(dcm2quat(c));
}

inline Mat3 cne(const Blh& blh) {
	const double slat = std::sin(blh.lat);
	const double clat = std::cos(blh.lat);
	const double slon = std::sin(blh.lon);
	const double clon = std::cos(blh.lon);

	Mat3 c;
	c.m[0][0] = -slat * clon;
	c.m[0][1] = -slon;
	c.m[0][2] = -clat * clon;

	c.m[1][0] = -slat * slon;
	c.m[1][1] = clon;
	c.m[1][2] = -clat * slon;

	c.m[2][0] = clat;
	c.m[2][1] = 0.0;
	c.m[2][2] = -slat;
	return c;
}

inline Quaternion qne(const Blh& blh) {
		return dcm2quat(cne(blh));
}

inline Blh blh(const Quaternion& q_ne) {
	const Mat3 c = quat2dcm(q_ne);
	Blh out;
	out.lat = std::atan2(-c.m[2][2], c.m[2][0]);
	out.lon = std::atan2(-c.m[0][1], c.m[1][1]);
	out.h = 0.0;
	return out;
}

inline Vec3 blh2ecef(const Blh& p, const Wgs84& wgs84 = Wgs84()) {
	const double e2 = 2.0 * wgs84.f - wgs84.f * wgs84.f;
	const double slat = std::sin(p.lat);
	const double clat = std::cos(p.lat);
	const double slon = std::sin(p.lon);
	const double clon = std::cos(p.lon);

	const double rn = wgs84.a / std::sqrt(1.0 - e2 * slat * slat);

	Vec3 e;
	e.x = (rn + p.h) * clat * clon;
	e.y = (rn + p.h) * clat * slon;
	e.z = (rn * (1.0 - e2) + p.h) * slat;
	return e;
}

inline Blh ecef2blh(const Vec3& ecef, const Wgs84& wgs84 = Wgs84()) {
	const double e2 = 2.0 * wgs84.f - wgs84.f * wgs84.f;
	const double b = wgs84.a * (1.0 - wgs84.f);
	const double ep2 = (wgs84.a * wgs84.a - b * b) / (b * b);

	const double x = ecef.x;
	const double y = ecef.y;
	const double z = ecef.z;
	const double p = std::sqrt(x * x + y * y);

	Blh out;
	out.lon = std::atan2(y, x);

	if (p < 1e-12) {
		out.lat = (z >= 0.0) ? 0.5 * 3.14159265358979323846 : -0.5 * 3.14159265358979323846;
		out.h = std::fabs(z) - b;
		return out;
	}

	const double theta = std::atan2(z * wgs84.a, p * b);
	const double st = std::sin(theta);
	const double ct = std::cos(theta);
	out.lat = std::atan2(z + ep2 * b * st * st * st, p - e2 * wgs84.a * ct * ct * ct);

	const double slat = std::sin(out.lat);
	const double rn = wgs84.a / std::sqrt(1.0 - e2 * slat * slat);
	out.h = p / std::cos(out.lat) - rn;
	return out;
}

inline Vec3 DRi(const Vec3& dr_n, const Blh& ref_blh) {
		return matVec(cne(ref_blh), dr_n);
}

inline Vec3 DR(const Vec3& dr_e, const Blh& ref_blh) {
		return matVec(transpose(cne(ref_blh)), dr_e);
}

inline Blh local2global(const Vec3& local_n,
						 const Blh& origin,
						 const Wgs84& wgs84 = Wgs84()) {
	const Vec3 e0 = blh2ecef(origin, wgs84);
	const Vec3 de = DRi(local_n, origin);
	Vec3 e;
	e.x = e0.x + de.x;
	e.y = e0.y + de.y;
	e.z = e0.z + de.z;
		return ecef2blh(e, wgs84);
}

inline Vec3 global2local(const Blh& global,
						 const Blh& origin,
						 const Wgs84& wgs84 = Wgs84()) {
	const Vec3 eg = blh2ecef(global, wgs84);
	const Vec3 e0 = blh2ecef(origin, wgs84);
	Vec3 de;
	de.x = eg.x - e0.x;
	de.y = eg.y - e0.y;
	de.z = eg.z - e0.z;
		return DR(de, origin);
}

inline Vec3 iewe(const Wgs84& wgs84 = Wgs84()) {
	return Vec3{0.0, 0.0, wgs84.omega_ie};
}

inline Vec3 iewn(const Blh& blh_pos, const Wgs84& wgs84 = Wgs84()) {
	const Vec3 omega_ie_e = iewe(wgs84);
		return matVec(transpose(cne(blh_pos)), omega_ie_e);
}

inline Vec3 enwn(const Blh& blh_pos,
				 const Vec3& vel_n,
				 const Wgs84& wgs84 = Wgs84()) {
	const double e2 = 2.0 * wgs84.f - wgs84.f * wgs84.f;
	const double slat = std::sin(blh_pos.lat);
	const double denom = 1.0 - e2 * slat * slat;
	const double sqrt_denom = std::sqrt(denom);

	const double rn = wgs84.a / sqrt_denom;
	const double rm = wgs84.a * (1.0 - e2) / (denom * sqrt_denom);

	const double rn_h = rn + blh_pos.h;
	const double rm_h = rm + blh_pos.h;

	Vec3 w_en_n;
	w_en_n.x = vel_n.y / rn_h;
	w_en_n.y = -vel_n.x / rm_h;
	w_en_n.z = -vel_n.y * std::tan(blh_pos.lat) / rn_h;
	return w_en_n;
}

}  // 鍛藉悕绌洪棿缁撴潫
