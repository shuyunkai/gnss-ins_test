#pragma once
#include <cmath>
#include "../common/Coordinate_2.h"
#include "../common/erathupdate.h"
#include "../fileio/imu.h"
namespace ins {
struct InsCompensatedIncrement {
	Vec3 dtheta_b;
	Vec3 dvel_b;
	Vec3 zeta_n;
	double dt = 0.0;
};
struct AttitudeUpdateResult {
	Quaternion q_nb_cur;
	Blh blh_half;
	Vec3 omega_in_n_half;
};
inline Quaternion quatInverseUnit(const Quaternion& q) {
	return Quaternion{q.w, -q.x, -q.y, -q.z};
}
inline Vec3 addVec3(const Vec3& a, const Vec3& b) {
	return Vec3{a.x + b.x, a.y + b.y, a.z + b.z};
}
inline Vec3 subVec3(const Vec3& a, const Vec3& b) {
	return Vec3{a.x - b.x, a.y - b.y, a.z - b.z};
}
inline Vec3 scaleVec3(const Vec3& v, double s) {
	return Vec3{v.x * s, v.y * s, v.z * s};
}
inline double dotVec3(const Vec3& a, const Vec3& b) {
	return a.x * b.x + a.y * b.y + a.z * b.z;
}
inline Vec3 crossVec3(const Vec3& a, const Vec3& b) {
	return Vec3{
		a.y * b.z - a.z * b.y,
		a.z * b.x - a.x * b.z,
		a.x * b.y - a.y * b.x};
}
// 注意：移除内部归一化，避免多次连乘累积精度损失
// 归一化应在最终结果统一进行
inline Quaternion quatMul(const Quaternion& p, const Quaternion& q) {
	Quaternion r;
	r.w = p.w * q.w - p.x * q.x - p.y * q.y - p.z * q.z;
	r.x = p.w * q.x + p.x * q.w + p.y * q.z - p.z * q.y;
	r.y = p.w * q.y - p.x * q.z + p.y * q.w + p.z * q.x;
	r.z = p.w * q.z + p.x * q.y - p.y * q.x + p.z * q.w;
	return r;  // 不归一化，由调用方统一处理
}
inline Quaternion normalizeQuaternionSafe(const Quaternion& q,
								const Quaternion& fallback = Quaternion{}) {
	const double n2 = q.w * q.w + q.x * q.x + q.y * q.y + q.z * q.z;
	if (!std::isfinite(n2) || n2 <= 1e-20) {
				return normalizeQuaternion(fallback);
	}
		return normalizeQuaternion(q);
}
inline Vec3 imuDeltaTheta(const ImuData& imu) {
	return Vec3{imu.dtheta_x, imu.dtheta_y, imu.dtheta_z};
}
inline Vec3 imuDeltaVel(const ImuData& imu) {
	return Vec3{imu.dvel_x, imu.dvel_y, imu.dvel_z};
}
inline double safeScaleDenominator(double scale_error) {
	const double denom = 1.0 + scale_error;
	if (std::fabs(denom) < 1e-8) {
		return (denom >= 0.0) ? 1e-8 : -1e-8;
	}
	return denom;
}
inline ImuData compensateSingleImuByError(
		const ImuData& imu_raw,
		const ImuErrorData& imu_error,
		double dt) {
	ImuData imu_corr = imu_raw;
	const double inv_gyro_scale_x = 1.0 / safeScaleDenominator(imu_error.gyro_scale[0]);
	const double inv_gyro_scale_y = 1.0 / safeScaleDenominator(imu_error.gyro_scale[1]);
	const double inv_gyro_scale_z = 1.0 / safeScaleDenominator(imu_error.gyro_scale[2]);
	imu_corr.dtheta_x = inv_gyro_scale_x * (imu_raw.dtheta_x - imu_error.gyro_bias[0] * dt);
	imu_corr.dtheta_y = inv_gyro_scale_y * (imu_raw.dtheta_y - imu_error.gyro_bias[1] * dt);
	imu_corr.dtheta_z = inv_gyro_scale_z * (imu_raw.dtheta_z - imu_error.gyro_bias[2] * dt);
	const double inv_accel_scale_x = 1.0 / safeScaleDenominator(imu_error.accel_scale[0]);
	const double inv_accel_scale_y = 1.0 / safeScaleDenominator(imu_error.accel_scale[1]);
	const double inv_accel_scale_z = 1.0 / safeScaleDenominator(imu_error.accel_scale[2]);
	imu_corr.dvel_x = inv_accel_scale_x * (imu_raw.dvel_x - imu_error.accel_bias[0] * dt);
	imu_corr.dvel_y = inv_accel_scale_y * (imu_raw.dvel_y - imu_error.accel_bias[1] * dt);
	imu_corr.dvel_z = inv_accel_scale_z * (imu_raw.dvel_z - imu_error.accel_bias[2] * dt);
	return imu_corr;
}
inline ImuMeasureData compensateImuMeasureByError(
		const ImuMeasureData& imu_measure_raw,
		const ImuErrorData& imu_error,
		double min_dt = 1e-6) {
	ImuMeasureData imu_measure_corr = imu_measure_raw;
	double dt = imu_measure_raw.imucur_.time - imu_measure_raw.imupre_.time;
	if (dt < min_dt) {
		dt = min_dt;
	}
	imu_measure_corr.imupre_ = compensateSingleImuByError(imu_measure_raw.imupre_, imu_error, dt);
	imu_measure_corr.imucur_ = compensateSingleImuByError(imu_measure_raw.imucur_, imu_error, dt);
	return imu_measure_corr;
}
inline Vec3 pvaVelToVec3(const PvaData& pva) {
	return Vec3{pva.vel_n[0], pva.vel_n[1], pva.vel_n[2]};
}
inline Blh pvaBlhToBlh(const PvaData& pva) {
	return Blh{pva.blh[0], pva.blh[1], pva.blh[2]};
}
inline Euler pvaEulerToEuler(const PvaData& pva) {
	return Euler{pva.euler[0], pva.euler[1], pva.euler[2]};
}
inline void vec3ToPvaVel(const Vec3& v, PvaData& pva) {
	pva.vel_n[0] = v.x;
	pva.vel_n[1] = v.y;
	pva.vel_n[2] = v.z;
}
inline void blhToPvaBlh(const Blh& blh_pos, PvaData& pva) {
	pva.blh[0] = blh_pos.lat;
	pva.blh[1] = blh_pos.lon;
	pva.blh[2] = blh_pos.h;
}
inline void eulerToPvaEuler(const Euler& e, PvaData& pva) {
	pva.euler[0] = e.roll;
	pva.euler[1] = e.pitch;
	pva.euler[2] = e.yaw;
}
inline Blh estimateMidpointBlh(const Blh& blh_pre,
							   const Vec3& vel_ref_n,
							   double dt) {
	const LocalPosition pos_pre{blh_pre.lat, blh_pre.lon, blh_pre.h};
	const LocalVelocity vel_ref_local{vel_ref_n.x, vel_ref_n.y, vel_ref_n.z};
	const EarthParameters earth_pre = updateEarthParameters(pos_pre, vel_ref_local);
	const double h_half = blh_pre.h - 0.5 * dt * vel_ref_n.z;
	const double lat_half = blh_pre.lat + 0.5 * dt * vel_ref_n.x / (earth_pre.rm + h_half);
	const double cos_lat_half = std::cos(lat_half);
	const double safe_cos_lat_half = (std::fabs(cos_lat_half) < 1e-8)
							  ? ((cos_lat_half >= 0.0) ? 1e-8 : -1e-8)
							  : cos_lat_half;
	const double lon_half = blh_pre.lon +
		0.5 * dt * vel_ref_n.y / ((earth_pre.rn + h_half) * safe_cos_lat_half);
	return Blh{lat_half, lon_half, h_half};
}
inline InsCompensatedIncrement compensateImuConingSculling(
		const ImuData& imu_pre,
		const ImuData& imu_cur,
		double min_dt = 1e-6) {
	const Vec3 dtheta_km1 = imuDeltaTheta(imu_pre);
	const Vec3 dtheta_k = imuDeltaTheta(imu_cur);
	const Vec3 dvel_km1 = imuDeltaVel(imu_pre);
	const Vec3 dvel_k = imuDeltaVel(imu_cur);
	InsCompensatedIncrement out;
	out.dt = imu_cur.time - imu_pre.time;
	if (out.dt < min_dt) {
		out.dt = min_dt;
	}
	out.dtheta_b = addVec3(dtheta_k, scaleVec3(crossVec3(dtheta_km1, dtheta_k), 1.0 / 12.0));
	const Vec3 term0 = dvel_k;
	const Vec3 term1 = scaleVec3(crossVec3(dtheta_k, dvel_k), 0.5);
	const Vec3 term2 = scaleVec3(
		addVec3(crossVec3(dtheta_km1, dvel_k), crossVec3(dvel_km1, dtheta_k)),
		1.0 / 12.0);
	out.dvel_b = addVec3(term0, addVec3(term1, term2));
	out.zeta_n = Vec3{0.0, 0.0, 0.0};
	return out;
}
// 修正Issue 4：通过速度积分计算中间位置，移除对 q_ne_cur 的依赖，符合先位置后姿态的顺序
inline Quaternion updateAttitude(
		const PvaData& pva_pre,
		const PvaData& pva_cur,
		const Vec3& dtheta_b,
		double dt) {
	const Vec3 vel_pre = pvaVelToVec3(pva_pre);
	const Vec3 vel_cur = pvaVelToVec3(pva_cur);
	const Vec3 vel_mid = scaleVec3(addVec3(vel_pre, vel_cur), 0.5);
	const Blh blh_pre = pvaBlhToBlh(pva_pre);

	// 1. 计算中间时刻位置 (通过速度积分，而非四元数插值)
	// 这样可以避免依赖 k 时刻位置 (q_ne_cur)，确保计算顺序正确
	const Vec3 dr_half_n = scaleVec3(vel_mid, 0.5 * dt);
	Blh blh_mid = local2global(dr_half_n, blh_pre);
	blh_mid.h = blh_pre.h - 0.5 * dt * vel_mid.z;

	// 2. 计算中间时刻的地球参数和导航系角速度
	const Vec3 omega_ie_n = iewn(blh_mid);
	const Vec3 omega_en_n = enwn(blh_mid, vel_mid);
	const Vec3 omega_in_n = addVec3(omega_ie_n, omega_en_n);

	// 3. 按 qnn * qnb_pre * qbb 更新姿态
	const Quaternion q_nb_pre = euler2quat(pvaEulerToEuler(pva_pre));
	const Quaternion q_nn = rotvec2quat(scaleVec3(omega_in_n, dt));
	const Quaternion q_bb = rotvec2quat(dtheta_b);
	const Quaternion q_nb_cur = quatMul(quatMul(q_nn, q_nb_pre), q_bb);
	
		return normalizeQuaternionSafe(q_nb_cur, q_nb_pre);
}
inline Vec3 updateVelocityCorrect(
		const PvaData& pva_pre,
		const Quaternion& q_nb_pre,
		const Vec3& dvel_b,
		double dt) {
	const Blh blh_pre = pvaBlhToBlh(pva_pre);
	const Vec3 vel_pre = pvaVelToVec3(pva_pre);
	const Mat3 c_nb_pre = quat2dcm(q_nb_pre);
	const Vec3 dv_fb = matVec(c_nb_pre, dvel_b);

	// 1) 在k-1时刻计算 n 系比力积分项与重力/哥氏积分项。
	const Vec3 omega_ie_n_pre = iewn(blh_pre);
	const Vec3 omega_en_n_pre = enwn(blh_pre, vel_pre);
	const Vec3 omega_in_half_pre = scaleVec3(addVec3(omega_ie_n_pre, omega_en_n_pre), 0.5 * dt);
	const Vec3 dv_fn_stage1 = subVec3(dv_fb, crossVec3(omega_in_half_pre, dv_fb));
	const Vec3 omega_sum_pre = addVec3(scaleVec3(omega_ie_n_pre, 2.0), omega_en_n_pre);
	const Vec3 coriolis_pre = crossVec3(omega_sum_pre, vel_pre);
	const LocalPosition pos_pre_local{blh_pre.lat, blh_pre.lon, blh_pre.h};
	const LocalVelocity vel_pre_local{vel_pre.x, vel_pre.y, vel_pre.z};
	const EarthParameters earth_pre = updateEarthParameters(pos_pre_local, vel_pre_local);
	const Vec3 gravity_pre{0.0, 0.0, earth_pre.gravity};
	const Vec3 dv_gn_stage1 = scaleVec3(subVec3(gravity_pre, coriolis_pre), dt);

	// 2) 得到k-1/2速度，并通过速度积分计算中间位置（修正Issue 4：避免四元数外推）
	const Vec3 vel_mid = addVec3(vel_pre, scaleVec3(addVec3(dv_fn_stage1, dv_gn_stage1), 0.5));
	
	// 通过速度积分计算中间位置，而非四元数外推
	const Vec3 dr_half_n = scaleVec3(vel_mid, 0.5 * dt);
	Blh blh_half = local2global(dr_half_n, blh_pre);
	blh_half.h = blh_pre.h - 0.5 * dt * vel_mid.z;

	// 3) 在k-1/2重算地球参数、n 系比力积分项和重力/哥氏积分项。
	const Vec3 omega_ie_n_half = iewn(blh_half);
	const Vec3 omega_en_n_half = enwn(blh_half, vel_mid);
	const Vec3 omega_in_half = scaleVec3(addVec3(omega_ie_n_half, omega_en_n_half), 0.5 * dt);
	const Vec3 dv_fn = subVec3(dv_fb, crossVec3(omega_in_half, dv_fb));
	const Vec3 omega_sum_half = addVec3(scaleVec3(omega_ie_n_half, 2.0), omega_en_n_half);
	const Vec3 coriolis_half = crossVec3(omega_sum_half, vel_mid);
	const LocalPosition pos_half_local{blh_half.lat, blh_half.lon, blh_half.h};
	const LocalVelocity vel_mid_local{vel_mid.x, vel_mid.y, vel_mid.z};
	const EarthParameters earth_half = updateEarthParameters(pos_half_local, vel_mid_local);
	const Vec3 gravity_half{0.0, 0.0, earth_half.gravity};
	const Vec3 dv_gn = scaleVec3(subVec3(gravity_half, coriolis_half), dt);

	// 4) 最终k时刻速度。
	return addVec3(vel_pre, addVec3(dv_fn, dv_gn));
}
inline Blh updatePosition(
		const PvaData& pva_pre,
		const Vec3& vel_cur,
		double dt) {
	const Blh blh_pre = pvaBlhToBlh(pva_pre);
	const Vec3 vel_pre = pvaVelToVec3(pva_pre);

	// 1) 先计算中间时刻速度和中间位置。
	const Vec3 vel_mid = scaleVec3(addVec3(vel_pre, vel_cur), 0.5);
	const Vec3 dr_half_n = scaleVec3(vel_mid, 0.5 * dt);
	Blh blh_mid = local2global(dr_half_n, blh_pre);
	blh_mid.h = blh_pre.h - 0.5 * dt * vel_mid.z;

	// 2) 在中间时刻重算导航系旋转矢量。
	const Vec3 omega_ie_n_mid = iewn(blh_mid);
	const Vec3 omega_en_n_mid = enwn(blh_mid, vel_mid);
	const Vec3 rot_n = scaleVec3(addVec3(omega_ie_n_mid, omega_en_n_mid), dt);
	const Quaternion q_nn = rotvec2quat(rot_n);

	// 3) 计算地球系旋转四元数，并用 qee*qne_pre*qnn 更新到当前时刻位置四元数。
	const Wgs84 wgs84;
	const Vec3 rot_e{0.0, 0.0, -wgs84.omega_ie * dt};
	const Quaternion q_ee = rotvec2quat(rot_e);
	const Quaternion q_ne_pre = qne(blh_pre);
	const Quaternion q_ne_cur = quatMul(quatMul(q_ee, q_ne_pre), q_nn);

	// 4) 由 q_ne_cur 反算当前经纬度，高程按中点速度积分更新。
	Blh blh_cur = blh(q_ne_cur);
	blh_cur.h = blh_pre.h - dt * vel_mid.z;
	return blh_cur;
}
inline void propagateIns(
		const ImuMeasureData& imu_measure,
		NavigationStatusData& nav_state,
		double min_dt = 1e-6) {
	const ImuMeasureData imu_measure_corr =
		compensateImuMeasureByError(imu_measure, nav_state.imuerror_, min_dt);
	const InsCompensatedIncrement comp =
		compensateImuConingSculling(imu_measure_corr.imupre_, imu_measure_corr.imucur_, min_dt);
	const PvaData& pva_pre = nav_state.pvapre_;
	PvaData pva_cur = pva_pre;
	const Quaternion q_nb_pre = euler2quat(pvaEulerToEuler(pva_pre));
	const Vec3 vel_cur = updateVelocityCorrect(
		pva_pre,
		q_nb_pre,
		comp.dvel_b,
		comp.dt);
	const Blh blh_cur = updatePosition(pva_pre, vel_cur, comp.dt);
	vec3ToPvaVel(vel_cur, pva_cur);
	blhToPvaBlh(blh_cur, pva_cur);
	const Quaternion q_nb_cur = updateAttitude(pva_pre, pva_cur, comp.dtheta_b, comp.dt);
	eulerToPvaEuler(quat2euler(q_nb_cur), pva_cur);
	nav_state.pvacur_ = pva_cur;
}
}  // 鍛藉悕绌洪棿缁撴潫
