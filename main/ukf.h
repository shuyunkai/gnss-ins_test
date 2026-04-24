#pragma once

#include "../common/Coordinate_2.h"
#include "kalman.h"
#include <stdexcept>
#include <vector>
#include <functional>
#include <cmath>

namespace ins {

class UnscentedKalmanFilter {
public:
    UnscentedKalmanFilter() = default;

    // (1) 步骤1：初始化
    // 设定初始状态 x_0 与协方差 P_0
    explicit UnscentedKalmanFilter(std::size_t state_dim, double alpha = 1e-3, double kappa = 0.0, double beta = 2.0) {
        resize(state_dim, alpha, kappa, beta);
    }

    void resize(std::size_t state_dim, double alpha = 1e-3, double kappa = 0.0, double beta = 2.0) {
        state_dim_ = state_dim;
        alpha_ = alpha;
        kappa_ = kappa;
        beta_ = beta;
        
        x_ = Matrix(state_dim, 1, 0.0);
        p_ = Matrix::identity(state_dim);
        
        lambda_ = alpha_ * alpha_ * (state_dim_ + kappa_) - state_dim_;
        
        weights_m_.resize(2 * state_dim_ + 1);
        weights_c_.resize(2 * state_dim_ + 1);
        
        weights_m_[0] = lambda_ / (state_dim_ + lambda_);
        weights_c_[0] = weights_m_[0] + (1 - alpha_ * alpha_ + beta_);
        
        for (std::size_t i = 1; i <= 2 * state_dim_; ++i) {
            weights_m_[i] = 1.0 / (2.0 * (state_dim_ + lambda_));
            weights_c_[i] = weights_m_[i];
        }
        
        initialized_ = true;
    }

    bool isInitialized() const { return initialized_; }

    void setState(const Matrix& x0) {
        x_ = x0;
        initialized_ = true;
    }

    void setCovariance(const Matrix& p0) {
        p_ = symmetrizeMatrix(p0);
        initialized_ = true;
    }

    const Matrix& state() const { return x_; }
    const Matrix& covariance() const { return p_; }
    std::size_t stateDim() const { return state_dim_; }
    const Matrix& kalmanGain() const { return k_; }

    // 辅助函数: Cholesky 分解
    static Matrix cholesky(const Matrix& a) {
        std::size_t n = a.rows();
        Matrix l(n, n, 0.0);
        for (std::size_t i = 0; i < n; i++) {
            for (std::size_t j = 0; j <= i; j++) {
                double sum = 0;
                for (std::size_t k = 0; k < j; k++)
                    sum += l(i, k) * l(j, k);
                
                if (i == j) {
                    double val = a(i, i) - sum;
                    l(i, j) = std::sqrt(val > 1e-15 ? val : 1e-15); // 增加对底层奇异的下界保护，防止对角元直接置零
                } else {
                    double l_jj = l(j, j);
                    l(i, j) = (a(i, j) - sum) / (l_jj > 1e-12 ? l_jj : 1e-12);
                }
            }
        }
        return l;
    }

    // 2.1. 生成 sigma 点 
    // 从 x_{k-1} 和 P_{k-1} 中构造 2n + 1 个 sigma 点
    std::vector<Matrix> generateSigmaPoints() const {
        std::size_t n = state_dim_;
        std::vector<Matrix> sigmas(2 * n + 1, Matrix(n, 1, 0.0));
        sigmas[0] = x_;
        Matrix p_scaled = p_;
        for (std::size_t r = 0; r < n; ++r) {
            for (std::size_t c = 0; c < n; ++c) {
                p_scaled(r, c) *= (n + lambda_);
            }
        }
        Matrix L = cholesky(p_scaled);
        for (std::size_t i = 0; i < n; ++i) {
            Matrix l_col(n, 1, 0.0);
            for (std::size_t r = 0; r < n; ++r) l_col(r, 0) = L(r, i);
            sigmas[i + 1] = add(x_, l_col);
            sigmas[i + 1 + n] = sub(x_, l_col);
        }
        return sigmas;
    }

    // (2) 步骤2：预测阶段 (Time Update)
    void predict(const std::function<Matrix(const Matrix&)>& f_func, const Matrix& q) {
        if (!initialized_) throw std::runtime_error("UKF not initialized");
        
        std::size_t n = state_dim_;
        
        // 2.1. 生成 sigma 点
        std::vector<Matrix> sigmas = generateSigmaPoints();
        
        // 2.2. 状态预测: 将 sigma 点输入状态转移函数 f(\chi_i)
        std::vector<Matrix> predicted_sigmas(2 * n + 1, Matrix(n, 1, 0.0));
        Matrix x_pred(n, 1, 0.0);
        for (std::size_t i = 0; i < 2 * n + 1; ++i) {
            predicted_sigmas[i] = f_func(sigmas[i]);
            // 然后加权平均得出 x_{k|k-1}
            for (std::size_t r = 0; r < n; ++r) {
                x_pred(r, 0) += weights_m_[i] * predicted_sigmas[i](r, 0);
            }
        }
        
        // 2.3. 协方差预测: P_{k|k-1}
        Matrix p_pred(n, n, 0.0);
        for (std::size_t i = 0; i < 2 * n + 1; ++i) {
            Matrix diff = sub(predicted_sigmas[i], x_pred);
            Matrix diff_cov = mul(diff, transpose(diff));
            for (std::size_t r = 0; r < n; ++r) {
                for (std::size_t c = 0; c < n; ++c) {
                    p_pred(r, c) += weights_c_[i] * diff_cov(r, c);
                }
            }
        }
        p_pred = add(p_pred, q); // 加上系统噪声协方差 Q
        
        x_ = x_pred;
        p_ = symmetrizeMatrix(p_pred);
        if (!isCovarianceValid()) regularizeCovariance();
    }





    // 兼容原有的线性转移矩阵 (将其封装为Lambda非线性函数)
    // UKF时间更新(预测)步：将所有Sigma点分别代入非线性状态转移函数进行推演，最后利用均值和协方差权重重建出 k 时刻的先验状态与先验P阵

    void predict(const Matrix& phi, const Matrix& q) {
        auto linear_f = [&phi](const Matrix& chi_i) { return mul(phi, chi_i); };
        predict(linear_f, q);
    }

    void predict(const Matrix& f, const Matrix& q, const Matrix& gamma) {
        if (!initialized_) throw std::runtime_error("UKF not initialized");
        const Matrix q_mapped = mul(mul(gamma, q), transpose(gamma));
        predict(f, q_mapped);
    }

    // (3) 步骤3：更新阶段 (Measurement Update)
    // UKF量测更新步：将预测后的Sigma点映射到实际观测空间，加权求和得出预估观测值，并利用互协方差计算卡尔曼增益系数进而获取后验更新状态

    void update(const Matrix& z_measurement, const std::function<Matrix(const Matrix&)>& h_func, const Matrix& r) {
        if (!initialized_) throw std::runtime_error("UKF not initialized");
        
        std::size_t n = state_dim_;
        std::size_t m = z_measurement.rows(); 
        
        // 此处为了观测函数的传播，再次从预测结果生成 sigma 点
        std::vector<Matrix> sigmas = generateSigmaPoints();
        
        // 3.1. 将预测的 sigma 点传入观测函数: h(\chi_i^x)
        std::vector<Matrix> z_sigmas(2 * n + 1, Matrix(m, 1, 0.0));
        Matrix z_pred(m, 1, 0.0);
        for (std::size_t i = 0; i < 2 * n + 1; ++i) {
            z_sigmas[i] = h_func(sigmas[i]);
            // 3.2. 计算观测预测值 \hat{z}_k
            for (std::size_t r_idx = 0; r_idx < m; ++r_idx) {
                z_pred(r_idx, 0) += weights_m_[i] * z_sigmas[i](r_idx, 0);
            }
        }
        
        // 3.2续: 计算观测预测值的协方差 P_z
        Matrix p_zz(m, m, 0.0);
        for (std::size_t i = 0; i < 2 * n + 1; ++i) {
            Matrix z_diff = sub(z_sigmas[i], z_pred);
            Matrix z_cov = mul(z_diff, transpose(z_diff));
            for (std::size_t r_idx = 0; r_idx < m; ++r_idx) {
                for (std::size_t c_idx = 0; c_idx < m; ++c_idx) {
                    p_zz(r_idx, c_idx) += weights_c_[i] * z_cov(r_idx, c_idx);
                }
            }
        }
        p_zz = add(p_zz, symmetrizeMatrix(r)); // 加上测量噪声 R
        s_ = p_zz;
        
        // 3.3. 计算状态与观测之间的协方差: P_{xz}
        Matrix p_xz(n, m, 0.0);
        for (std::size_t i = 0; i < 2 * n + 1; ++i) {
            Matrix x_diff = sub(sigmas[i], x_);
            Matrix z_diff = sub(z_sigmas[i], z_pred);
            Matrix xz_cov = mul(x_diff, transpose(z_diff));
            for (std::size_t r_idx = 0; r_idx < n; ++r_idx) {
                for (std::size_t c_idx = 0; c_idx < m; ++c_idx) {
                    p_xz(r_idx, c_idx) += weights_c_[i] * xz_cov(r_idx, c_idx);
                }
            }
        }
        
        // 3.4. 计算卡尔曼增益: K_{k} = P_{xz} P_z^{-1}
        const Matrix p_zz_inv = inverseWithJitter(p_zz);
        k_ = mul(p_xz, p_zz_inv);
        
        // 3.5. 更新状态与协方差: 
        // \hat{x}_k = \hat{x}_{k|k-1} + K_k * (Z_k - \hat{Z}_k)
        y_ = sub(z_measurement, z_pred); // 观测方程残差
        x_ = add(x_, mul(k_, y_));
        
        // P_{k} = P_{k|k-1} - K_k P_z K_k^T
        Matrix kpzzkt = mul(mul(k_, p_zz), transpose(k_));
        p_ = sub(p_, kpzzkt);
        p_ = symmetrizeMatrix(p_);
        
        if (!isCovarianceValid()) {
            regularizeCovariance();
        }
    }





    // 兼容原有的线性观测矩阵 (将其封装为Lambda非线性函数)
    // UKF量测更新步：将预测后的Sigma点映射到实际观测空间，加权求和得出预估观测值，并利用互协方差计算卡尔曼增益系数进而获取后验更新状态

    void update(const Matrix& z_residual, const Matrix& h, const Matrix& r) {
        auto linear_h = [&h](const Matrix& chi_i) { return mul(h, chi_i); };
        update(z_residual, linear_h, r);
    }

    void applyCovarianceResetJacobian(const Matrix& j) {
        if (!initialized_) throw std::runtime_error("UKF not initialized");
        if (j.rows() != state_dim_ || j.cols() != state_dim_) {
            throw std::invalid_argument("UKF apply: dims mismatch");
        }
        p_ = symmetrizeMatrix(mul(mul(j, p_), transpose(j)));
        if (!isCovarianceValid()) regularizeCovariance();
    }

    void zeroStateSubrange(std::size_t start_row, std::size_t count) {
        if (!initialized_) throw std::runtime_error("UKF not initialized");
        for (std::size_t i = 0; i < count; ++i) {
            x_(start_row + i, 0) = 0.0;
        }
    }

    bool isCovarianceValid(double sym_tol = 1e-10, double pd_eps = 1e-12) const {
        if (!initialized_) return false;
        if (!isFiniteMatrix(p_) || !isSymmetricMatrix(p_, sym_tol)) return false;
        return isPositiveDefiniteMatrix(p_, pd_eps);
    }

    void regularizeCovariance(double diag_floor = 1e-10, double init_jitter = 1e-10, int max_tries = 10) {
        if (!initialized_) return;
        p_ = symmetrizeMatrix(p_);
        for (std::size_t i = 0; i < p_.rows(); ++i) {
            if (!std::isfinite(p_(i, i)) || p_(i, i) < diag_floor) {
                p_(i, i) = diag_floor;
            }
        }
        double jitter = init_jitter;
        for (int k = 0; k < max_tries; ++k) {
            if (isCovarianceValid()) return;
            p_ = addDiagonalJitter(symmetrizeMatrix(p_), jitter);
            jitter *= 10.0;
        }
    }

private:
    std::size_t state_dim_ = 0;
    bool initialized_ = false;
    double alpha_ = 1e-3;
    double kappa_ = 0.0;
    double beta_ = 2.0;
    double lambda_ = 0.0;
    
    std::vector<double> weights_m_;
    std::vector<double> weights_c_;

    Matrix x_;
    Matrix p_;
    Matrix k_;
    Matrix y_;
    Matrix s_;
};

inline void gnssPositionUpdateAndFeedback21(
        UnscentedKalmanFilter& kf,
        const Matrix& residual_z,
        const Matrix& r,
        NavigationStatusData& nav_state,
        const std::array<double, 3>& antenna_lever_arm_b = {0.0, 0.0, 0.0},
        bool reset_state_after_feedback = true) {
    if (kf.stateDim() != ErrorStateIndex21::kDim) {
        throw std::invalid_argument("UKF: KF must be 21-dim");
    }
    // 纯非线性观测方程 h(\chi_i)：计算真实误差状态 \chi_i 下预期产生的测距残差
    auto h_func = [&](const Matrix& chi_i) {
        // 直接提取位置和姿态误差进行前馈摄动，规避 feedback 函数内部零偏限幅器对非线性传递特性的跳变影响
        const double d_n = chi_i(ErrorStateIndex21::kPos + 0, 0);
        const double d_e = chi_i(ErrorStateIndex21::kPos + 1, 0);
        const double d_d = chi_i(ErrorStateIndex21::kPos + 2, 0);
        const Vec3 dphi{chi_i(ErrorStateIndex21::kAtt + 0, 0),
                        chi_i(ErrorStateIndex21::kAtt + 1, 0),
                        chi_i(ErrorStateIndex21::kAtt + 2, 0)};

        // 摄动后的姿态矩阵
        const Euler nom_euler{nav_state.pvacur_.euler[0], nav_state.pvacur_.euler[1], nav_state.pvacur_.euler[2]};
        const Mat3 c_ins = quat2dcm(euler2quat(nom_euler));
        Mat3 i_plus_phi_cross = identity3();
        i_plus_phi_cross.m[0][1] = -dphi.z; i_plus_phi_cross.m[0][2] = dphi.y;
        i_plus_phi_cross.m[1][0] = dphi.z;  i_plus_phi_cross.m[1][2] = -dphi.x;
        i_plus_phi_cross.m[2][0] = -dphi.y; i_plus_phi_cross.m[2][1] = dphi.x;
        const Mat3 c_nb = matMul(i_plus_phi_cross, c_ins);

        // 摄动后的天线位置
        const LocalPosition pos_cur{nav_state.pvacur_.blh[0], nav_state.pvacur_.blh[1], nav_state.pvacur_.blh[2]};
        const LocalVelocity vel_cur{nav_state.pvacur_.vel_n[0], nav_state.pvacur_.vel_n[1], nav_state.pvacur_.vel_n[2]};
        const EarthParameters earth = updateEarthParameters(pos_cur, vel_cur);
        const double cos_lat = std::cos(nav_state.pvacur_.blh[0]);
        const double safe_cos_lat = (std::fabs(cos_lat) < 1e-8) ? ((cos_lat >= 0.0) ? 1e-8 : -1e-8) : cos_lat;

        const Blh imu_blh = {
            nav_state.pvacur_.blh[0] - d_n / (earth.rm + nav_state.pvacur_.blh[2]),
            nav_state.pvacur_.blh[1] - d_e / ((earth.rn + nav_state.pvacur_.blh[2]) * safe_cos_lat),
            nav_state.pvacur_.blh[2] + d_d
        };

        const Vec3 lever_b{antenna_lever_arm_b[0], antenna_lever_arm_b[1], antenna_lever_arm_b[2]};
        const Vec3 lever_n = matVec(c_nb, lever_b);
        const Blh pert_ant_blh = local2global(lever_n, imu_blh);
        
        // 标称(nav_state)预估的天线位置
        const Blh nom_imu_blh{nav_state.pvacur_.blh[0], nav_state.pvacur_.blh[1], nav_state.pvacur_.blh[2]};
        const Euler nom_imu_euler{nav_state.pvacur_.euler[0], nav_state.pvacur_.euler[1], nav_state.pvacur_.euler[2]};
        const Quaternion nom_q_nb = euler2quat(nom_imu_euler);
        const Mat3 nom_c_nb = quat2dcm(nom_q_nb);
        const Vec3 nom_lever_n = matVec(nom_c_nb, lever_b);
        const Blh nom_ant_blh = local2global(nom_lever_n, nom_imu_blh);

        // 返回如果真实误差是 chi_i 时，标称 INS 位置减去 GNSS 位置会观察到的残差
        // global2local(target, origin) -> 从 origin 指向 target 的矢量在 NED 的坐标
        // 我们要计算 (INS_nom - GNSS_pert)
        const Vec3 dr_pert_minus_nom = global2local(pert_ant_blh, nom_ant_blh);
        Matrix z(3, 1, 0.0);
        z(0, 0) = -dr_pert_minus_nom.x;
        z(1, 0) = -dr_pert_minus_nom.y;
        z(2, 0) = -dr_pert_minus_nom.z;
        return z;
    };
    kf.update(residual_z, h_func, r);
    const Matrix x_fb = kf.state();
    feedbackInsClosedLoopFrom21State(x_fb, nav_state);
    if (reset_state_after_feedback) {
        kf.zeroStateSubrange(0, ErrorStateIndex21::kDim);
    }
}

} // namespace ins
