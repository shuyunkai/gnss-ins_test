#pragma once

#include <cmath>

namespace fusion {

// 作用：定义同一时刻判定阈值（秒）。
constexpr double kDefaultTimeEqualThreshold = 1e-3;

// 作用：描述两个时刻的前后关系。
enum class TimeOrder {
	kBefore = -1,
	kSame = 0,
	kAfter = 1
};

// 作用：描述当前IMU步应执行的融合动作。
enum class UpdateAction {
	kPredictOnly = 0,
	kUpdateOldImu,
	kUpdateNewImu,
	kInterpolateAndUpdate
};

// 作用：保存GNSS与双子样IMU的时序关系及最终动作。
struct GnssImuTimeDecision {
	TimeOrder gnss_vs_old_imu = TimeOrder::kSame;
	TimeOrder gnss_vs_new_imu = TimeOrder::kSame;

	bool gnss_between_imus = false;
	bool old_imu_close_to_gnss = false;
	bool new_imu_close_to_gnss = false;

	UpdateAction action = UpdateAction::kPredictOnly;
};

// 作用：比较两个时刻，返回前后关系。
inline TimeOrder compareTime(double lhs_time,
							 double rhs_time,
							 double equal_threshold = kDefaultTimeEqualThreshold) {
	const double dt = lhs_time - rhs_time;
	if (std::fabs(dt) <= equal_threshold) {
		return TimeOrder::kSame;
	}
	return (dt < 0.0) ? TimeOrder::kBefore : TimeOrder::kAfter;
}

// 作用：根据GNSS与IMU时序关系输出融合动作。
inline GnssImuTimeDecision decideGnssImuTimeRelation(
	double gnss_time,
	double old_imu_time,
	double new_imu_time,
	double equal_threshold = kDefaultTimeEqualThreshold) {
	GnssImuTimeDecision decision;

	decision.gnss_vs_old_imu = compareTime(gnss_time, old_imu_time, equal_threshold);
	decision.gnss_vs_new_imu = compareTime(gnss_time, new_imu_time, equal_threshold);

	decision.old_imu_close_to_gnss =
		(decision.gnss_vs_old_imu == TimeOrder::kSame);
	decision.new_imu_close_to_gnss =
		(decision.gnss_vs_new_imu == TimeOrder::kSame);

	decision.gnss_between_imus =
		(decision.gnss_vs_old_imu == TimeOrder::kAfter &&
		 decision.gnss_vs_new_imu == TimeOrder::kBefore);

	// 作用：节点对齐优先，优先在旧IMU节点更新。
	if (decision.old_imu_close_to_gnss && !decision.new_imu_close_to_gnss) {
		decision.action = UpdateAction::kUpdateOldImu;
		return decision;
	}

	// 作用：节点对齐优先，在新IMU节点更新。
	if (decision.new_imu_close_to_gnss && !decision.old_imu_close_to_gnss) {
		decision.action = UpdateAction::kUpdateNewImu;
		return decision;
	}

	// 作用：若两端都在阈值内，选择更近的一端更新。
	if (decision.old_imu_close_to_gnss && decision.new_imu_close_to_gnss) {
		const double old_dt = std::fabs(gnss_time - old_imu_time);
		const double new_dt = std::fabs(gnss_time - new_imu_time);
		decision.action = (old_dt <= new_dt) ? UpdateAction::kUpdateOldImu
											 : UpdateAction::kUpdateNewImu;
		return decision;
	}

	// 作用：不在双子样区间内时仅做预测。
	if (!decision.gnss_between_imus) {
		decision.action = UpdateAction::kPredictOnly;
		return decision;
	}

	// 作用：在区间内但不靠近端点时，做插值更新。
	decision.action = UpdateAction::kInterpolateAndUpdate;
	return decision;
}

}  // namespace fusion
