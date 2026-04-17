#pragma once
#include <cmath>
#include <cstddef>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <vector>
#include "../common/Coordinate_2.h"
#include "../common/erathupdate.h"
#include "../fileio/imu.h"
namespace ins {
// 卡尔曼滤波：基础矩阵运算类 (支持加减乘转置等操作)
class Matrix {
public:
	Matrix() = default;
	Matrix(std::size_t rows, std::size_t cols, double init_value = 0.0)
		: rows_(rows), cols_(cols), data_(rows * cols, init_value) {}
		std::size_t rows() const { return rows_; }
		std::size_t cols() const { return cols_; }
		double& operator()(std::size_t r, std::size_t c) {
		return data_[r * cols_ + c];
	}
		double operator()(std::size_t r, std::size_t c) const {
		return data_[r * cols_ + c];
	}
		static Matrix identity(std::size_t n) {
				Matrix i(n, n, 0.0);
		for (std::size_t k = 0; k < n; ++k) {
			i(k, k) = 1.0;
		}
		return i;
	}
private:
	std::size_t rows_ = 0;
	std::size_t cols_ = 0;
	std::vector<double> data_;
};
inline void checkSameShape(const Matrix& a, const Matrix& b, const char* op) {
	if (a.rows() != b.rows() || a.cols() != b.cols()) {
		throw std::invalid_argument(std::string(op) + ": matrix shape mismatch");
	}
}
inline Matrix transpose(const Matrix& a) {
		Matrix t(a.cols(), a.rows(), 0.0);
	for (std::size_t r = 0; r < a.rows(); ++r) {
		for (std::size_t c = 0; c < a.cols(); ++c) {
			t(c, r) = a(r, c);
		}
	}
	return t;
}
inline Matrix add(const Matrix& a, const Matrix& b) {
	checkSameShape(a, b, "add");
		Matrix r(a.rows(), a.cols(), 0.0);
	for (std::size_t i = 0; i < a.rows(); ++i) {
		for (std::size_t j = 0; j < a.cols(); ++j) {
			r(i, j) = a(i, j) + b(i, j);
		}
	}
	return r;
}
inline Matrix sub(const Matrix& a, const Matrix& b) {
	checkSameShape(a, b, "sub");
		Matrix r(a.rows(), a.cols(), 0.0);
	for (std::size_t i = 0; i < a.rows(); ++i) {
		for (std::size_t j = 0; j < a.cols(); ++j) {
			r(i, j) = a(i, j) - b(i, j);
		}
	}
	return r;
}
inline Matrix mul(const Matrix& a, const Matrix& b) {
	if (a.cols() != b.rows()) {
		throw std::invalid_argument("mul: matrix shape mismatch");
	}
		Matrix r(a.rows(), b.cols(), 0.0);
	for (std::size_t i = 0; i < a.rows(); ++i) {
		for (std::size_t k = 0; k < a.cols(); ++k) {
			const double aik = a(i, k);
			for (std::size_t j = 0; j < b.cols(); ++j) {
				r(i, j) += aik * b(k, j);
			}
		}
	}
	return r;
}
inline Matrix inverse(const Matrix& a) {
	if (a.rows() != a.cols()) {
		throw std::invalid_argument("inverse: matrix must be square");
	}
	const std::size_t n = a.rows();
		Matrix aug(n, 2 * n, 0.0);
	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t j = 0; j < n; ++j) {
			aug(i, j) = a(i, j);
		}
		aug(i, n + i) = 1.0;
	}
	for (std::size_t col = 0; col < n; ++col) {
		std::size_t pivot = col;
		double pivot_abs = std::fabs(aug(col, col));
		for (std::size_t r = col + 1; r < n; ++r) {
			const double v = std::fabs(aug(r, col));
			if (v > pivot_abs) {
				pivot_abs = v;
				pivot = r;
			}
		}
		if (pivot_abs < 1e-15) {
			throw std::runtime_error("inverse: matrix is singular");
		}
		if (pivot != col) {
			for (std::size_t j = 0; j < 2 * n; ++j) {
				const double tmp = aug(col, j);
				aug(col, j) = aug(pivot, j);
				aug(pivot, j) = tmp;
			}
		}
		const double div = aug(col, col);
		for (std::size_t j = 0; j < 2 * n; ++j) {
			aug(col, j) /= div;
		}
		for (std::size_t r = 0; r < n; ++r) {
			if (r == col) {
				continue;
			}
			const double factor = aug(r, col);
			if (std::fabs(factor) < 1e-30) {
				continue;
			}
			for (std::size_t j = 0; j < 2 * n; ++j) {
				aug(r, j) -= factor * aug(col, j);
			}
		}
	}
		Matrix inv(n, n, 0.0);
	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t j = 0; j < n; ++j) {
			inv(i, j) = aug(i, n + j);
		}
	}
	return inv;
}
inline bool isFiniteMatrix(const Matrix& a) {
	for (std::size_t i = 0; i < a.rows(); ++i) {
		for (std::size_t j = 0; j < a.cols(); ++j) {
			if (!std::isfinite(a(i, j))) {
				return false;
			}
		}
	}
	return true;
}
inline Matrix symmetrizeMatrix(const Matrix& a) {
	if (a.rows() != a.cols()) {
		throw std::invalid_argument("symmetrizeMatrix");;
	}
		Matrix s(a.rows(), a.cols(), 0.0);
	for (std::size_t i = 0; i < a.rows(); ++i) {
		for (std::size_t j = 0; j < a.cols(); ++j) {
			s(i, j) = 0.5 * (a(i, j) + a(j, i));
		}
	}
	return s;
}
inline bool isSymmetricMatrix(const Matrix& a, double tol = 1e-10) {
	if (a.rows() != a.cols()) {
		return false;
	}
	for (std::size_t i = 0; i < a.rows(); ++i) {
		for (std::size_t j = i + 1; j < a.cols(); ++j) {
			if (std::fabs(a(i, j) - a(j, i)) > tol) {
				return false;
			}
		}
	}
	return true;
}
inline bool isPositiveDefiniteMatrix(const Matrix& a, double eps = 1e-12) {
	if (a.rows() != a.cols()) {
		return false;
	}
	const std::size_t n = a.rows();
		Matrix l(n, n, 0.0);
	for (std::size_t i = 0; i < n; ++i) {
		for (std::size_t j = 0; j <= i; ++j) {
			double sum = a(i, j);
			for (std::size_t k = 0; k < j; ++k) {
				sum -= l(i, k) * l(j, k);
			}
			if (i == j) {
				if (sum <= eps || !std::isfinite(sum)) {
					return false;
				}
				l(i, j) = std::sqrt(sum);
			} else {
				if (std::fabs(l(j, j)) <= eps) {
					return false;
				}
				l(i, j) = sum / l(j, j);
			}
		}
	}
	return true;
}
inline Matrix addDiagonalJitter(const Matrix& a, double jitter) {
	if (a.rows() != a.cols()) {
		throw std::invalid_argument("addDiagonalJitter");;
	}
	Matrix r = a;
	for (std::size_t i = 0; i < r.rows(); ++i) {
		r(i, i) += jitter;
	}
	return r;
}
inline Matrix inverseWithJitter(const Matrix& a,
								double init_jitter = 1e-12,
								int max_tries = 6) {
	Matrix s = symmetrizeMatrix(a);
	double jitter = init_jitter;
	for (int k = 0; k < max_tries; ++k) {
		try {
						return inverse(s);
		} catch (const std::exception&) {
			s = addDiagonalJitter(s, jitter);
			jitter *= 10.0;
		}
	}
	throw std::runtime_error("inverseWithJitter");;
}
inline double clampSymmetric(double v, double abs_limit) {
	if (abs_limit <= 0.0) {
		return v;
	}
	if (v > abs_limit) {
		return abs_limit;
	}
	if (v < -abs_limit) {
		return -abs_limit;
	}
	return v;
}
class ExtendedKalmanFilter {
public:
	ExtendedKalmanFilter() = default;
		explicit ExtendedKalmanFilter(std::size_t state_dim) {
		resize(state_dim);
	}
		void resize(std::size_t state_dim) {
		x_ = Matrix(state_dim, 1, 0.0);
		p_ = Matrix::identity(state_dim);
		initialized_ = true;
	}
		bool isInitialized() const {
		return initialized_;
	}
		void setState(const Matrix& x0) {
		if (x0.cols() != 1) {
			throw std::invalid_argument("setState: x must be column vector");
		}
		x_ = x0;
		if (!initialized_) {
			p_ = Matrix::identity(x0.rows());
			initialized_ = true;
		}
	}
		void setCovariance(const Matrix& p0) {
		if (p0.rows() != p0.cols()) {
			throw std::invalid_argument("setCovariance: P must be square");
		}
		if (initialized_ && p0.rows() != x_.rows()) {
			throw std::invalid_argument("setCovariance: dimension mismatch with state");
		}
		p_ = symmetrizeMatrix(p0);
		if (!isCovarianceValid()) {
			regularizeCovariance();
		}
		if (!initialized_) {
			x_ = Matrix(p0.rows(), 1, 0.0);
			initialized_ = true;
		}
	}
		const Matrix& state() const {
		return x_;
	}
		const Matrix& covariance() const {
		return p_;
	}
		const Matrix& kalmanGain() const {
		return k_;
	}
		const Matrix& innovation() const {
		return y_;
	}
		const Matrix& innovationCovariance() const {
		return s_;
	}
		double lastInnovationNis() const {
		return last_innovation_nis_;
	}
		bool lastUpdateAccepted() const {
		return last_update_accepted_;
	}
		void setInnovationTest(bool enabled, double gate_threshold) {
		innovation_test_enabled_ = enabled;
		innovation_gate_threshold_ = (gate_threshold > 0.0) ? gate_threshold : innovation_gate_threshold_;
	}
		std::size_t stateDim() const {
		return x_.rows();
	}
		bool isCovarianceValid(double sym_tol = 1e-10, double pd_eps = 1e-12) const {
		if (!initialized_) {
			return false;
		}
		if (p_.rows() != p_.cols()) {
			return false;
		}
		if (!isFiniteMatrix(p_)) {
			return false;
		}
		if (!isSymmetricMatrix(p_, sym_tol)) {
			return false;
		}
		for (std::size_t i = 0; i < p_.rows(); ++i) {
			if (!(p_(i, i) > pd_eps)) {
				return false;
			}
		}
				return isPositiveDefiniteMatrix(p_, pd_eps);
	}
		void regularizeCovariance(double diag_floor = 1e-12,
							double init_jitter = 1e-12,
							int max_tries = 8) {
		ensureInitialized();
		p_ = symmetrizeMatrix(p_);
		for (std::size_t i = 0; i < p_.rows(); ++i) {
			if (!std::isfinite(p_(i, i)) || p_(i, i) < diag_floor) {
				p_(i, i) = diag_floor;
			}
		}
		double jitter = init_jitter;
		for (int k = 0; k < max_tries; ++k) {
			if (isCovarianceValid()) {
				return;
			}
			p_ = addDiagonalJitter(symmetrizeMatrix(p_), jitter);
			jitter *= 10.0;
		}
	}
	// 作用：状态更新与预测 (Predict)
	void predict(const Matrix& f, const Matrix& q) {
		predict(f, q, Matrix::identity(x_.rows()));
	}
	// 作用：状态更新与预测 (Predict)
	void predict(const Matrix& f, const Matrix& q, const Matrix& gamma) {
		ensureInitialized();
		if (f.rows() != x_.rows() || f.cols() != x_.rows()) {
			throw std::invalid_argument("predict: F dimension mismatch");
		}
		if (q.rows() != x_.rows() || q.cols() != x_.rows()) {
			throw std::invalid_argument("predict: Q dimension mismatch");
		}
		if (gamma.rows() != x_.rows() || gamma.cols() != x_.rows()) {
			throw std::invalid_argument("predict: Gamma dimension mismatch");
		}
		x_ = mul(f, x_);
		const Matrix q_mapped = mul(mul(gamma, q), transpose(gamma));
		p_ = add(mul(mul(f, p_), transpose(f)), q_mapped);
		p_ = symmetrizeMatrix(p_);
		if (!isCovarianceValid()) {
			regularizeCovariance();
		}
	}
	// 作用：状态更新与预测 (Predict)
	void predict(const Matrix& f, const Matrix& q, const Matrix& b, const Matrix& u) {
		predict(f, q, Matrix::identity(x_.rows()), b, u);
	}
	// 作用：状态更新与预测 (Predict)
	void predict(const Matrix& f,
					 const Matrix& q,
					 const Matrix& gamma,
					 const Matrix& b,
					 const Matrix& u) {
		ensureInitialized();
		if (f.rows() != x_.rows() || f.cols() != x_.rows()) {
			throw std::invalid_argument("predict(control): F dimension mismatch");
		}
		if (b.rows() != x_.rows() || b.cols() != u.rows() || u.cols() != 1) {
			throw std::invalid_argument("predict(control): B/U dimension mismatch");
		}
		if (q.rows() != x_.rows() || q.cols() != x_.rows()) {
			throw std::invalid_argument("predict(control): Q dimension mismatch");
		}
		if (gamma.rows() != x_.rows() || gamma.cols() != x_.rows()) {
			throw std::invalid_argument("predict(control): Gamma dimension mismatch");
		}
		x_ = add(mul(f, x_), mul(b, u));
		const Matrix q_mapped = mul(mul(gamma, q), transpose(gamma));
		p_ = add(mul(mul(f, p_), transpose(f)), q_mapped);
		p_ = symmetrizeMatrix(p_);
		if (!isCovarianceValid()) {
			regularizeCovariance();
		}
	}
	// 作用：量测更新 (Update)
	void update(const Matrix& z, const Matrix& h, const Matrix& r) {
		ensureInitialized();
		if (z.cols() != 1) {
			throw std::invalid_argument("update: z must be column vector");
		}
		if (h.cols() != x_.rows() || h.rows() != z.rows()) {
			throw std::invalid_argument("update: H dimension mismatch");
		}
		if (r.rows() != z.rows() || r.cols() != z.rows()) {
			throw std::invalid_argument("update: R dimension mismatch");
		}
		if (!isFiniteMatrix(r)) {
			throw std::invalid_argument("update: R contains non-finite values");
		}
		const Matrix r_sym = symmetrizeMatrix(r);
		if (!isPositiveDefiniteMatrix(r_sym)) {
			throw std::invalid_argument("update: R must be positive definite");
		}
		if (!isCovarianceValid()) {
			regularizeCovariance();
		}
		const Matrix ht = transpose(h);
		y_ = sub(z, mul(h, x_));
		s_ = symmetrizeMatrix(add(mul(mul(h, p_), ht), r_sym));
		const Matrix s_inv = inverseWithJitter(s_);
		last_innovation_nis_ = 0.0;
		k_ = mul(mul(p_, ht), s_inv);
		
		x_ = add(x_, mul(k_, y_));
		
		const Matrix i = Matrix::identity(x_.rows());
		const Matrix ikh = sub(i, mul(k_, h));
		p_ = add(mul(mul(ikh, p_), transpose(ikh)), mul(mul(k_, r_sym), transpose(k_)));
		p_ = symmetrizeMatrix(p_);
		if (!isCovarianceValid()) {
			regularizeCovariance();
		}
		last_update_accepted_ = true;
	}
		void zeroStateSubrange(std::size_t start_row, std::size_t count) {
		ensureInitialized();
		if (start_row + count > x_.rows()) {
			throw std::invalid_argument("zeroStateSubrange: out of range");
		}
		for (std::size_t i = 0; i < count; ++i) {
			x_(start_row + i, 0) = 0.0;
		}
	}
		void applyCovarianceResetJacobian(const Matrix& j) {
		ensureInitialized();
		if (j.rows() != x_.rows() || j.cols() != x_.rows()) {
			throw std::invalid_argument("applyCovarianceResetJacobian: J dimension mismatch");
		}
		p_ = mul(mul(j, p_), transpose(j));
		p_ = symmetrizeMatrix(p_);
		if (!isCovarianceValid()) {
			regularizeCovariance();
		}
	}
private:
		void ensureInitialized() const {
		if (!initialized_) {
			throw std::runtime_error("kalman filter is not initialized");
		}
	}
private:
	bool initialized_ = false;
	Matrix x_;
	Matrix p_;
	Matrix k_;
	Matrix y_;
	Matrix s_;
	bool innovation_test_enabled_ = false;
	double innovation_gate_threshold_ = 11.344866730144373;  // chi2(3, 0.99)
	double last_innovation_nis_ = 0.0;
	bool last_update_accepted_ = true;
};
struct ErrorStateIndex21 {
	static constexpr std::size_t kPos = 0;
	static constexpr std::size_t kVel = 3;
	static constexpr std::size_t kAtt = 6;
	static constexpr std::size_t kGyroBias = 9;
	static constexpr std::size_t kAccelBias = 12;
	static constexpr std::size_t kGyroScale = 15;
	static constexpr std::size_t kAccelScale = 18;
	static constexpr std::size_t kDim = 21;
};
inline Matrix skewMatrixFromVec3(const Vec3& v) {
		Matrix s(3, 3, 0.0);
	s(0, 1) = -v.z;
	s(0, 2) = v.y;
	s(1, 0) = v.z;
	s(1, 2) = -v.x;
	s(2, 0) = -v.y;
	s(2, 1) = v.x;
	return s;
}
inline Matrix buildPositionObservationMatrix21(
		const NavigationStatusData& nav_state,
		const std::array<double, 3>& antenna_lever_arm_b = {0.0, 0.0, 0.0}) {
		Matrix h(3, ErrorStateIndex21::kDim, 0.0);
	h(0, ErrorStateIndex21::kPos + 0) = 1.0;
	h(1, ErrorStateIndex21::kPos + 1) = 1.0;
	h(2, ErrorStateIndex21::kPos + 2) = 1.0;
	const PvaData& pva = nav_state.pvacur_;
	const Euler e{pva.euler[0], pva.euler[1], pva.euler[2]};
	const Mat3 c_nb = quat2dcm(euler2quat(e));
	const Vec3 lever_b{antenna_lever_arm_b[0], antenna_lever_arm_b[1], antenna_lever_arm_b[2]};
	const Vec3 lever_n = matVec(c_nb, lever_b);
	const Matrix lever_cross = skewMatrixFromVec3(lever_n);
	for (std::size_t r = 0; r < 3; ++r) {
		for (std::size_t c = 0; c < 3; ++c) {
			h(r, ErrorStateIndex21::kAtt + c) = lever_cross(r, c);
		}
	}
	return h;
}
inline Matrix buildPositionMeasurementNoiseFromGnssStd(
		double std_n,
		double std_e,
		double std_d,
		double min_std = 1e-4) {
	const double s_n = std::max(std_n, min_std);
	const double s_e = std::max(std_e, min_std);
	const double s_d = std::max(std_d, min_std);
		Matrix r(3, 3, 0.0);
	r(0, 0) = s_n * s_n;
	r(1, 1) = s_e * s_e;
	r(2, 2) = s_d * s_d;
	return r;
}
inline Matrix buildPositionVelocityMeasurementNoiseFromGnssStd(
		double std_n,
		double std_e,
		double std_d,
		double std_vn,
		double std_ve,
		double std_vd,
		double min_pos_std = 1e-4,
		double min_vel_std = 0.05) {
	const double s_n = std::max(std::fabs(std_n), min_pos_std);
	const double s_e = std::max(std::fabs(std_e), min_pos_std);
	const double s_d = std::max(std::fabs(std_d), min_pos_std);
	const double s_vn = std::max(std::fabs(std_vn), min_vel_std);
	const double s_ve = std::max(std::fabs(std_ve), min_vel_std);
	const double s_vd = std::max(std::fabs(std_vd), min_vel_std);
		Matrix r(6, 6, 0.0);
	r(0, 0) = s_n * s_n;
	r(1, 1) = s_e * s_e;
	r(2, 2) = s_d * s_d;
	r(3, 3) = s_vn * s_vn;
	r(4, 4) = s_ve * s_ve;
	r(5, 5) = s_vd * s_vd;
	return r;
}
inline Matrix buildResetJacobianFrom21StateError(const Matrix& x_err) {
        if (x_err.rows() < ErrorStateIndex21::kDim || x_err.cols() != 1) {
                throw std::invalid_argument("buildResetJacobianFrom21StateError: x_err dimension mismatch");
        }
        Matrix j = Matrix::identity(ErrorStateIndex21::kDim);
        const double phi_n = x_err(ErrorStateIndex21::kAtt + 0, 0);
        const double phi_e = x_err(ErrorStateIndex21::kAtt + 1, 0);
        const double phi_d = x_err(ErrorStateIndex21::kAtt + 2, 0);
        j(ErrorStateIndex21::kAtt + 0, ErrorStateIndex21::kAtt + 1) = phi_d;
        j(ErrorStateIndex21::kAtt + 0, ErrorStateIndex21::kAtt + 2) = -phi_e;
        j(ErrorStateIndex21::kAtt + 1, ErrorStateIndex21::kAtt + 0) = -phi_d;
        j(ErrorStateIndex21::kAtt + 1, ErrorStateIndex21::kAtt + 2) = phi_n;
        j(ErrorStateIndex21::kAtt + 2, ErrorStateIndex21::kAtt + 0) = phi_e;
        j(ErrorStateIndex21::kAtt + 2, ErrorStateIndex21::kAtt + 1) = -phi_n;
        return j;
}
inline void feedbackInsClosedLoopFrom21State(
		const Matrix& x_err,
		NavigationStatusData& nav_state) {
	if (x_err.rows() < ErrorStateIndex21::kDim || x_err.cols() != 1) {
		throw std::invalid_argument("feedbackInsClosedLoopFrom21State: x_err dimension mismatch");
	}
	PvaData corrected = nav_state.pvacur_;
	const double d_n = x_err(ErrorStateIndex21::kPos + 0, 0);
	const double d_e = x_err(ErrorStateIndex21::kPos + 1, 0);
	const double d_d = x_err(ErrorStateIndex21::kPos + 2, 0);
	const LocalPosition pos_cur{corrected.blh[0], corrected.blh[1], corrected.blh[2]};
	const LocalVelocity vel_cur{corrected.vel_n[0], corrected.vel_n[1], corrected.vel_n[2]};
	const EarthParameters earth = updateEarthParameters(pos_cur, vel_cur);
	const double cos_lat = std::cos(corrected.blh[0]);
	const double safe_cos_lat = (std::fabs(cos_lat) < 1e-8)
								  ? ((cos_lat >= 0.0) ? 1e-8 : -1e-8)
								  : cos_lat;
	corrected.blh[0] -= d_n / (earth.rm + corrected.blh[2]);
	corrected.blh[1] -= d_e / ((earth.rn + corrected.blh[2]) * safe_cos_lat);
	corrected.blh[2] += d_d;
	corrected.vel_n[0] -= x_err(ErrorStateIndex21::kVel + 0, 0);
	corrected.vel_n[1] -= x_err(ErrorStateIndex21::kVel + 1, 0);
	corrected.vel_n[2] -= x_err(ErrorStateIndex21::kVel + 2, 0);
	const Vec3 dphi{
		x_err(ErrorStateIndex21::kAtt + 0, 0),
		x_err(ErrorStateIndex21::kAtt + 1, 0),
		x_err(ErrorStateIndex21::kAtt + 2, 0)};
	const Euler e_ins{corrected.euler[0], corrected.euler[1], corrected.euler[2]};
	const Mat3 c_ins = quat2dcm(euler2quat(e_ins));
	Mat3 i_plus_phi_cross = identity3();
	i_plus_phi_cross.m[0][1] = -dphi.z;
	i_plus_phi_cross.m[0][2] = dphi.y;
	i_plus_phi_cross.m[1][0] = dphi.z;
	i_plus_phi_cross.m[1][2] = -dphi.x;
	i_plus_phi_cross.m[2][0] = -dphi.y;
	i_plus_phi_cross.m[2][1] = dphi.x;
	const Mat3 c_corr = matMul(i_plus_phi_cross, c_ins);
	const Euler e_corr = dcm2euler(c_corr);
	corrected.euler[0] = e_corr.roll;
	corrected.euler[1] = e_corr.pitch;
	corrected.euler[2] = e_corr.yaw;
	constexpr double kFeedbackAlpha = 1.0;
	const double gb_limit = deg2rad(5000.0) / 3600.0;  // rad/s
	const double ab_limit = 2.0;                        // m/s^2 (00000 mGal)
	const double gs_limit = 5000.0e-6;                  // ratio
	const double as_limit = 20000.0e-6;                 // ratio
	nav_state.imuerror_.gyro_bias[0] = clampSymmetric(
		nav_state.imuerror_.gyro_bias[0] + kFeedbackAlpha * x_err(ErrorStateIndex21::kGyroBias + 0, 0),
		gb_limit);
	nav_state.imuerror_.gyro_bias[1] = clampSymmetric(
		nav_state.imuerror_.gyro_bias[1] + kFeedbackAlpha * x_err(ErrorStateIndex21::kGyroBias + 1, 0),
		gb_limit);
	nav_state.imuerror_.gyro_bias[2] = clampSymmetric(
		nav_state.imuerror_.gyro_bias[2] + kFeedbackAlpha * x_err(ErrorStateIndex21::kGyroBias + 2, 0),
		gb_limit);
	nav_state.imuerror_.accel_bias[0] = clampSymmetric(
		nav_state.imuerror_.accel_bias[0] + kFeedbackAlpha * x_err(ErrorStateIndex21::kAccelBias + 0, 0),
		ab_limit);
	nav_state.imuerror_.accel_bias[1] = clampSymmetric(
		nav_state.imuerror_.accel_bias[1] + kFeedbackAlpha * x_err(ErrorStateIndex21::kAccelBias + 1, 0),
		ab_limit);
	nav_state.imuerror_.accel_bias[2] = clampSymmetric(
		nav_state.imuerror_.accel_bias[2] + kFeedbackAlpha * x_err(ErrorStateIndex21::kAccelBias + 2, 0),
		ab_limit);
	nav_state.imuerror_.gyro_scale[0] = clampSymmetric(
		nav_state.imuerror_.gyro_scale[0] + kFeedbackAlpha * x_err(ErrorStateIndex21::kGyroScale + 0, 0),
		gs_limit);
	nav_state.imuerror_.gyro_scale[1] = clampSymmetric(
		nav_state.imuerror_.gyro_scale[1] + kFeedbackAlpha * x_err(ErrorStateIndex21::kGyroScale + 1, 0),
		gs_limit);
	nav_state.imuerror_.gyro_scale[2] = clampSymmetric(
		nav_state.imuerror_.gyro_scale[2] + kFeedbackAlpha * x_err(ErrorStateIndex21::kGyroScale + 2, 0),
		gs_limit);
	nav_state.imuerror_.accel_scale[0] = clampSymmetric(
		nav_state.imuerror_.accel_scale[0] + kFeedbackAlpha * x_err(ErrorStateIndex21::kAccelScale + 0, 0),
		as_limit);
	nav_state.imuerror_.accel_scale[1] = clampSymmetric(
		nav_state.imuerror_.accel_scale[1] + kFeedbackAlpha * x_err(ErrorStateIndex21::kAccelScale + 1, 0),
		as_limit);
	nav_state.imuerror_.accel_scale[2] = clampSymmetric(
		nav_state.imuerror_.accel_scale[2] + kFeedbackAlpha * x_err(ErrorStateIndex21::kAccelScale + 2, 0),
		as_limit);
	nav_state.pvacur_ = corrected;
	nav_state.pvapre_ = corrected;
}
inline void gnssPositionUpdateAndFeedback21(
		ExtendedKalmanFilter& kf,
		const Matrix& residual_z,
		const Matrix& r,
		NavigationStatusData& nav_state,
		const std::array<double, 3>& antenna_lever_arm_b = {0.0, 0.0, 0.0},
		bool reset_state_after_feedback = true) {
	if (kf.stateDim() != ErrorStateIndex21::kDim) {
		throw std::invalid_argument("gnssPositionUpdateAndFeedback21: KF must be 21-dim");
	}
	if (residual_z.rows() != 3 || residual_z.cols() != 1) {
		throw std::invalid_argument("gnssPositionUpdateAndFeedback21: residual_z must be 3x1");
	}
	if (r.rows() != 3 || r.cols() != 3) {
		throw std::invalid_argument("gnssPositionUpdateAndFeedback21: R must be 3x3");
	}
	if (!isFiniteMatrix(r) || !isPositiveDefiniteMatrix(symmetrizeMatrix(r))) {
		throw std::invalid_argument("gnssPositionUpdateAndFeedback21: R must be finite and positive definite");
	}
	if (!kf.isCovarianceValid()) {
		throw std::runtime_error("gnssPositionUpdateAndFeedback21: prior P is invalid");
	}
	const Matrix h = buildPositionObservationMatrix21(nav_state, antenna_lever_arm_b);
	kf.update(residual_z, h, r);
	const Matrix x_fb = kf.state();
	feedbackInsClosedLoopFrom21State(x_fb, nav_state);
	if (reset_state_after_feedback) {
		const Matrix j_reset = buildResetJacobianFrom21StateError(x_fb);
		kf.applyCovarianceResetJacobian(j_reset);
		kf.zeroStateSubrange(0, ErrorStateIndex21::kDim);
	}
	if (!kf.isCovarianceValid()) {
		throw std::runtime_error("gnssPositionUpdateAndFeedback21: posterior P is invalid");
	}
}
inline void gnssPositionVelocityUpdateAndFeedback21(
		ExtendedKalmanFilter& kf,
		const Matrix& residual_z,
		const Matrix& r,
		NavigationStatusData& nav_state,
		const std::array<double, 3>& antenna_lever_arm_b = {0.0, 0.0, 0.0},
		const Vec3& omega_ib_b_est = Vec3{0.0, 0.0, 0.0},
		bool reset_state_after_feedback = true) {
	if (kf.stateDim() != ErrorStateIndex21::kDim) {
		throw std::invalid_argument("gnssPositionVelocityUpdateAndFeedback21: KF must be 21-dim");
	}
	if (residual_z.rows() != 6 || residual_z.cols() != 1) {
		throw std::invalid_argument("gnssPositionVelocityUpdateAndFeedback21: residual_z must be 6x1");
	}
	if (r.rows() != 6 || r.cols() != 6) {
		throw std::invalid_argument("gnssPositionVelocityUpdateAndFeedback21: R must be 6x6");
	}
	if (!isFiniteMatrix(r) || !isPositiveDefiniteMatrix(symmetrizeMatrix(r))) {
		throw std::invalid_argument("gnssPositionVelocityUpdateAndFeedback21: R must be finite and positive definite");
	}
	if (!kf.isCovarianceValid()) {
		throw std::runtime_error("gnssPositionVelocityUpdateAndFeedback21: prior P is invalid");
	}
	throw std::runtime_error("gnssPositionVelocityUpdateAndFeedback21: function removed, GNSS signal has no velocity");
}
}  // 锟秸硷拷锟
