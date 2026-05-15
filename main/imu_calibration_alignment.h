#pragma once

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "../common/Coordinate_2.h"
#include "../common/erathupdate.h"
#include "../fileio/imu.h"

namespace ins {

// =========================================================================
// 1. 六位置加速度计标定（零偏 + 比例因子 + 交轴耦合矩阵 M）
// =========================================================================
//
// 误差模型：
//   f_meas = b + M * f_true
//
// 其中：b = [b_x, b_y, b_z]^T  （零偏，单位 m/s^2）
//        M = 3x3 标定矩阵
//           - 对角线 M_ii = 1 + s_i（比例因子，理想值为1）
//           - 非对角线 M_ij (i!=j) 为交轴耦合系数（理想值为0）
//
// 六个位置（依次使z-y-x轴朝上）：
//   位置1: z朝上 -> f_true = [ 0,  0,  g]   位置2: z朝下 -> f_true = [ 0,  0, -g]
//   位置3: y朝上 -> f_true = [ 0,  g,  0]   位置4: y朝下 -> f_true = [ 0, -g,  0]
//   位置5: x朝上 -> f_true = [ g,  0,  0]   位置6: x朝下 -> f_true = [-g,  0,  0]
//
// 每个轴独立最小二乘（12参数，每轴4参数）：
//   设计矩阵 A^T A = diag(6, 2g^2, 2g^2, 2g^2)，解析解如下：
//     b_j = y_j/6
//     M_j0 = (y_j|x_up - y_j|x_down) / (2g)
//     M_j1 = (y_j|y_up - y_j|y_down) / (2g)
//     M_j2 = (y_j|z_up - y_j|z_down) / (2g)

struct SixPositionData {
    Vec3 pos1_z_up;
    Vec3 pos2_z_down;
    Vec3 pos3_y_up;
    Vec3 pos4_y_down;
    Vec3 pos5_x_up;
    Vec3 pos6_x_down;
    double gravity = 9.80;
};

struct AccelCalibrationResult {
    std::array<double, 3> bias = {0.0, 0.0, 0.0};  // 零偏 b (m/s^2)
    Mat3 M;                                           // 标定矩阵 (3x3)
    double residual_rms = 0.0;                        // 残差均方根 (m/s^2)
};

// 从一组静态IMU帧中计算平均速度增量 dvel。
inline Vec3 averageDvelFromStaticFrames(const std::vector<ImuData>& static_frames) {
    if (static_frames.empty()) return Vec3{0.0, 0.0, 0.0};
    Vec3 sum{0.0, 0.0, 0.0};
    for (const auto& imu : static_frames) {
        sum.x += imu.dvel_x;
        sum.y += imu.dvel_y;
        sum.z += imu.dvel_z;
    }
    const double inv_n = 1.0 / static_cast<double>(static_frames.size());
    return Vec3{sum.x * inv_n, sum.y * inv_n, sum.z * inv_n};
}

// 从一组静态IMU帧中计算平均角增量 dtheta。
inline Vec3 averageDthetaFromStaticFrames(const std::vector<ImuData>& static_frames) {
    if (static_frames.empty()) return Vec3{0.0, 0.0, 0.0};
    Vec3 sum{0.0, 0.0, 0.0};
    for (const auto& imu : static_frames) {
        sum.x += imu.dtheta_x;
        sum.y += imu.dtheta_y;
        sum.z += imu.dtheta_z;
    }
    const double inv_n = 1.0 / static_cast<double>(static_frames.size());
    return Vec3{sum.x * inv_n, sum.y * inv_n, sum.z * inv_n};
}

inline AccelCalibrationResult calibrateAccelerometerSixPosition(
        const SixPositionData& data) {

    const double g = data.gravity;
    const double inv_2g = 1.0 / (2.0 * g);

    // 6个位置的真实比力（载体坐标系）
    const Vec3 true_f[6] = {
        { 0.0,  0.0,  g},    // z_up
        { 0.0,  0.0, -g},    // z_down
        { 0.0,  g,   0.0},   // y_up
        { 0.0, -g,   0.0},   // y_down
        { g,    0.0,  0.0},  // x_up
        {-g,    0.0,  0.0}   // x_down
    };

    const Vec3 meas_f[6] = {
        data.pos1_z_up, data.pos2_z_down,
        data.pos3_y_up, data.pos4_y_down,
        data.pos5_x_up, data.pos6_x_down
    };

    AccelCalibrationResult result;

    auto y = [&](int pos, int axis) -> double {
        if (axis == 0) return meas_f[pos].x;
        if (axis == 1) return meas_f[pos].y;
        return meas_f[pos].z;
    };

    for (int j = 0; j < 3; ++j) {
        const double b  = (y(0,j) + y(1,j) + y(2,j) + y(3,j) + y(4,j) + y(5,j)) / 6.0;
        const double mj0 = (y(4,j) - y(5,j)) * inv_2g;  // M_j0: X输入耦合
        const double mj1 = (y(2,j) - y(3,j)) * inv_2g;  // M_j1: Y输入耦合
        const double mj2 = (y(0,j) - y(1,j)) * inv_2g;  // M_j2: Z输入耦合

        result.bias[j] = b;
        result.M.m[j][0] = mj0;
        result.M.m[j][1] = mj1;
        result.M.m[j][2] = mj2;
    }

    // 总残差均方根
    double total_res = 0.0;
    for (int i = 0; i < 6; ++i) {
        const double pred_x = result.bias[0]
            + result.M.m[0][0] * true_f[i].x
            + result.M.m[0][1] * true_f[i].y
            + result.M.m[0][2] * true_f[i].z;
        const double pred_y = result.bias[1]
            + result.M.m[1][0] * true_f[i].x
            + result.M.m[1][1] * true_f[i].y
            + result.M.m[1][2] * true_f[i].z;
        const double pred_z = result.bias[2]
            + result.M.m[2][0] * true_f[i].x
            + result.M.m[2][1] * true_f[i].y
            + result.M.m[2][2] * true_f[i].z;
        const double ex = meas_f[i].x - pred_x;
        const double ey = meas_f[i].y - pred_y;
        const double ez = meas_f[i].z - pred_z;
        total_res += ex*ex + ey*ey + ez*ez;
    }
    result.residual_rms = std::sqrt(total_res / 18.0);

    return result;
}


// =========================================================================
// 2. 静态陀螺零偏标定
// =========================================================================
//
// 误差模型（每轴独立）：
//   omega_meas_i = b_i + (1 + s_i) * omega_true_i    (i = x, y, z)

struct GyroBiasCalibrationResult {
    std::array<double, 3> bias     = {0.0, 0.0, 0.0};  // 零偏 (rad/s)
    std::array<double, 3> bias_std = {0.0, 0.0, 0.0};  // 零偏标准差 (rad/s)
    std::array<double, 3> scale    = {0.0, 0.0, 0.0};  // 比例因子误差 (无量纲)
};

// 从静态数据标定陀螺零偏（输入为角增量 dtheta）。
inline GyroBiasCalibrationResult calibrateGyroBiasStatic(
        const std::vector<ImuData>& static_frames) {

    const std::size_t n = static_frames.size();
    GyroBiasCalibrationResult result;

    if (n < 2) {
        return result;
    }

    Vec3 sum_dtheta{0.0, 0.0, 0.0};
    for (const auto& imu : static_frames) {
        sum_dtheta.x += imu.dtheta_x;
        sum_dtheta.y += imu.dtheta_y;
        sum_dtheta.z += imu.dtheta_z;
    }

    const double inv_n = 1.0 / static_cast<double>(n);
    result.bias[0] = sum_dtheta.x * inv_n;
    result.bias[1] = sum_dtheta.y * inv_n;
    result.bias[2] = sum_dtheta.z * inv_n;

    Vec3 var{0.0, 0.0, 0.0};
    for (const auto& imu : static_frames) {
        const double dx = imu.dtheta_x - result.bias[0];
        const double dy = imu.dtheta_y - result.bias[1];
        const double dz = imu.dtheta_z - result.bias[2];
        var.x += dx * dx;
        var.y += dy * dy;
        var.z += dz * dz;
    }
    result.bias_std[0] = std::sqrt(var.x / static_cast<double>(n - 1));
    result.bias_std[1] = std::sqrt(var.y / static_cast<double>(n - 1));
    result.bias_std[2] = std::sqrt(var.z / static_cast<double>(n - 1));

    return result;
}

// 带已知采样间隔 dt 的陀螺零偏标定，输出单位为 rad/s。
inline GyroBiasCalibrationResult calibrateGyroBiasStaticWithDt(
        const std::vector<ImuData>& static_frames,
        double dt) {

    const double dt_safe = (dt > 1e-6) ? dt : 1e-6;
    auto result_inc = calibrateGyroBiasStatic(static_frames);
    GyroBiasCalibrationResult result;
    result.bias[0]     = result_inc.bias[0]     / dt_safe;
    result.bias[1]     = result_inc.bias[1]     / dt_safe;
    result.bias[2]     = result_inc.bias[2]     / dt_safe;
    result.bias_std[0] = result_inc.bias_std[0] / dt_safe;
    result.bias_std[1] = result_inc.bias_std[1] / dt_safe;
    result.bias_std[2] = result_inc.bias_std[2] / dt_safe;
    return result;
}


// =========================================================================
// 2b. 两位置陀螺标定（零偏 + 比例因子）
// =========================================================================

struct GyroTwoPositionInput {
    Vec3 gyro_up[3];     // 各轴朝上时的平均角速率 (rad/s): [x_up, y_up, z_up]
    Vec3 gyro_down[3];   // 各轴朝下时的平均角速率 (rad/s): [x_down, y_down, z_down]
};

inline GyroBiasCalibrationResult calibrateGyroTwoPosition(
        const GyroTwoPositionInput& data,
        double latitude_rad,
        double earth_rate = 7.292115e-5)
{
    const double omega_v = earth_rate * std::sin(latitude_rad);

    GyroBiasCalibrationResult result;

    auto get = [](const Vec3& v, int axis) -> double {
        if (axis == 0) return v.x;
        if (axis == 1) return v.y;
        return v.z;
    };

    for (int j = 0; j < 3; ++j) {
        const double up   = get(data.gyro_up[j],   j);
        const double down = get(data.gyro_down[j], j);

        result.bias[j] = (up + down) * 0.5;

        const double denom = 2.0 * omega_v;
        if (std::fabs(denom) > 1e-15) {
            const double k = (up - down) / denom;
            result.scale[j] = k - 1.0;
        } else {
            result.scale[j] = 0.0;
        }

        const double resid_up   = up   - (result.bias[j] + (1.0 + result.scale[j]) * omega_v);
        const double resid_down = down - (result.bias[j] + (1.0 + result.scale[j]) * (-omega_v));
        result.bias_std[j] = std::sqrt((resid_up*resid_up + resid_down*resid_down) / 2.0);
    }

    return result;
}

// 3x3 矩阵求逆（用于补偿时计算 M^{-1}）
inline Mat3 inverseMat3(const Mat3& a) {
    const double det = a.m[0][0] * (a.m[1][1] * a.m[2][2] - a.m[1][2] * a.m[2][1])
                     - a.m[0][1] * (a.m[1][0] * a.m[2][2] - a.m[1][2] * a.m[2][0])
                     + a.m[0][2] * (a.m[1][0] * a.m[2][1] - a.m[1][1] * a.m[2][0]);
    const double inv_det = 1.0 / ((std::fabs(det) > 1e-20) ? det : 1e-20);

    Mat3 r;
    r.m[0][0] =  (a.m[1][1] * a.m[2][2] - a.m[1][2] * a.m[2][1]) * inv_det;
    r.m[0][1] =  (a.m[0][2] * a.m[2][1] - a.m[0][1] * a.m[2][2]) * inv_det;
    r.m[0][2] =  (a.m[0][1] * a.m[1][2] - a.m[0][2] * a.m[1][1]) * inv_det;
    r.m[1][0] =  (a.m[1][2] * a.m[2][0] - a.m[1][0] * a.m[2][2]) * inv_det;
    r.m[1][1] =  (a.m[0][0] * a.m[2][2] - a.m[0][2] * a.m[2][0]) * inv_det;
    r.m[1][2] =  (a.m[0][2] * a.m[1][0] - a.m[0][0] * a.m[1][2]) * inv_det;
    r.m[2][0] =  (a.m[1][0] * a.m[2][1] - a.m[1][1] * a.m[2][0]) * inv_det;
    r.m[2][1] =  (a.m[0][1] * a.m[2][0] - a.m[0][0] * a.m[2][1]) * inv_det;
    r.m[2][2] =  (a.m[0][0] * a.m[1][1] - a.m[0][1] * a.m[1][0]) * inv_det;
    return r;
}

// 使用标定结果补偿原始IMU增量数据（输出仍为增量格式）。
inline ImuData compensateImuData(const ImuData& raw,
                                  const AccelCalibrationResult& accel_calib,
                                  const GyroBiasCalibrationResult& gyro_calib,
                                  double dt) {
    ImuData comp = raw;
    const Mat3 invM = inverseMat3(accel_calib.M);

    const double fx = raw.dvel_x / dt;
    const double fy = raw.dvel_y / dt;
    const double fz = raw.dvel_z / dt;
    const double fx_u = fx - accel_calib.bias[0];
    const double fy_u = fy - accel_calib.bias[1];
    const double fz_u = fz - accel_calib.bias[2];
    comp.dvel_x = (invM.m[0][0] * fx_u + invM.m[0][1] * fy_u + invM.m[0][2] * fz_u) * dt;
    comp.dvel_y = (invM.m[1][0] * fx_u + invM.m[1][1] * fy_u + invM.m[1][2] * fz_u) * dt;
    comp.dvel_z = (invM.m[2][0] * fx_u + invM.m[2][1] * fy_u + invM.m[2][2] * fz_u) * dt;

    comp.dtheta_x = (raw.dtheta_x - gyro_calib.bias[0] * dt) / (1.0 + gyro_calib.scale[0]);
    comp.dtheta_y = (raw.dtheta_y - gyro_calib.bias[1] * dt) / (1.0 + gyro_calib.scale[1]);
    comp.dtheta_z = (raw.dtheta_z - gyro_calib.bias[2] * dt) / (1.0 + gyro_calib.scale[2]);

    return comp;
}


// =========================================================================
// 3. 解析粗对准 -- TRIAD算法
// =========================================================================

struct AlignmentResult {
    std::array<double, 3> euler = {0.0, 0.0, 0.0};  // 横滚、俯仰、航向（rad）
    Mat3 c_nb;                                        // 姿态方向余弦矩阵
    bool valid = false;                               // 对准是否有效
    double condition_number = 0.0;                    // TRIAD条件数（奇异性指标）
};

inline AlignmentResult analyticCoarseAlignment(
        const Vec3& specific_force_b,
        const Vec3& angular_rate_b,
        double latitude_rad,
        double height_m = 0.0,
        double gyro_noise_std = 1e-12,
        double gravity_override = 0.0,
        double earth_rate_override = 0.0)
{
    AlignmentResult result;

    const LocalPosition pos{latitude_rad, 0.0, height_m};
    const LocalVelocity vel{0.0, 0.0, 0.0};
    const EarthParameters earth = updateEarthParameters(pos, vel);
    const double g   = (gravity_override > 0.0)   ? gravity_override   : earth.gravity;
    const double omega_e = (earth_rate_override > 0.0) ? earth_rate_override : earth.omega_ie;

    const Vec3 g_n = {0.0, 0.0, g};

    const double cos_lat = std::cos(latitude_rad);
    const double sin_lat = std::sin(latitude_rad);
    const Vec3 omega_ie_n = {
         omega_e * cos_lat,
         0.0,
        -omega_e * sin_lat
    };

    const Vec3 g_b = {-specific_force_b.x,
                      -specific_force_b.y,
                      -specific_force_b.z};

    const Vec3 omega_b = angular_rate_b;

    const double g_norm  = vecNorm(g_b);
    const double g_n_norm = vecNorm(g_n);

    if (g_norm < 0.1 * g || g_norm > 3.0 * g) {
        result.valid = false;
        return result;
    }

    const double omega_e_norm = omega_e;

    const double inv_g_n  = 1.0 / std::max(g_n_norm,  1e-20);
    const double inv_g_b  = 1.0 / std::max(g_norm,     1e-20);

    Vec3 n1 = {g_n.x * inv_g_n, g_n.y * inv_g_n, g_n.z * inv_g_n};

    Vec3 cross_n;
    cross_n.x = g_n.y * omega_ie_n.z - g_n.z * omega_ie_n.y;
    cross_n.y = g_n.z * omega_ie_n.x - g_n.x * omega_ie_n.z;
    cross_n.z = g_n.x * omega_ie_n.y - g_n.y * omega_ie_n.x;
    const double cross_n_norm = vecNorm(cross_n);
    const double inv_cross_n = 1.0 / std::max(cross_n_norm, 1e-20);
    Vec3 n2 = {cross_n.x * inv_cross_n, cross_n.y * inv_cross_n, cross_n.z * inv_cross_n};

    Vec3 n3;
    n3.x = n1.y * n2.z - n1.z * n2.y;
    n3.y = n1.z * n2.x - n1.x * n2.z;
    n3.z = n1.x * n2.y - n1.y * n2.x;

    Vec3 b1 = {g_b.x * inv_g_b, g_b.y * inv_g_b, g_b.z * inv_g_b};

    Vec3 cross_b;
    cross_b.x = g_b.y * omega_b.z - g_b.z * omega_b.y;
    cross_b.y = g_b.z * omega_b.x - g_b.x * omega_b.z;
    cross_b.z = g_b.x * omega_b.y - g_b.y * omega_b.x;
    const double cross_b_norm = vecNorm(cross_b);
    const double inv_cross_b = 1.0 / std::max(cross_b_norm, 1e-20);
    Vec3 b2 = {cross_b.x * inv_cross_b, cross_b.y * inv_cross_b, cross_b.z * inv_cross_b};

    Vec3 b3;
    b3.x = b1.y * b2.z - b1.z * b2.y;
    b3.y = b1.z * b2.x - b1.x * b2.z;
    b3.z = b1.x * b2.y - b1.y * b2.x;

    Mat3 c_nb;
    for (int r = 0; r < 3; ++r) {
        for (int c = 0; c < 3; ++c) {
            auto get_b = [&](int row, int col) -> double {
                if (col == 0) { if (row == 0) return b1.x; if (row == 1) return b1.y; return b1.z; }
                if (col == 1) { if (row == 0) return b2.x; if (row == 1) return b2.y; return b2.z; }
                if (row == 0) return b3.x; if (row == 1) return b3.y; return b3.z;
            };
            auto get_n = [&](int row, int col) -> double {
                if (col == 0) { if (row == 0) return n1.x; if (row == 1) return n1.y; return n1.z; }
                if (col == 1) { if (row == 0) return n2.x; if (row == 1) return n2.y; return n2.z; }
                if (row == 0) return n3.x; if (row == 1) return n3.y; return n3.z;
            };
            c_nb.m[r][c] = get_b(r, 0) * get_n(c, 0)
                         + get_b(r, 1) * get_n(c, 1)
                         + get_b(r, 2) * get_n(c, 2);
        }
    }

    result.c_nb = c_nb;
    const Mat3 c_bn = transpose(c_nb);
    const Euler euler = dcm2euler(c_bn);
    result.euler[0] = euler.roll;
    result.euler[1] = euler.pitch;
    result.euler[2] = euler.yaw;

    const double sin_alpha_n = cross_n_norm / (g_n_norm * vecNorm(omega_ie_n));
    result.condition_number = 1.0 / std::max(sin_alpha_n, 1e-10);

    result.valid = (gyro_noise_std < 3.0 * omega_e_norm)
                && (g_norm > 0.5 * g && g_norm < 2.0 * g);

    return result;
}


// =========================================================================
// 4. 组合标定与对准流程
// =========================================================================

struct CalibrationAndAlignmentInput {
    std::vector<ImuData> pos1_z_up;
    std::vector<ImuData> pos2_z_down;
    std::vector<ImuData> pos3_y_up;
    std::vector<ImuData> pos4_y_down;
    std::vector<ImuData> pos5_x_up;
    std::vector<ImuData> pos6_x_down;

    std::vector<ImuData> rot1_zp;
    std::vector<ImuData> rot2_zn;
    std::vector<ImuData> rot3_yp;
    std::vector<ImuData> rot4_yn;
    std::vector<ImuData> rot5_xp;
    std::vector<ImuData> rot6_xn;

    std::vector<ImuData> align_static;

    double imu_dt   = 1.0 / 200.0;
    double align_dt = 0.0;

    double latitude_rad  = 0.0;
    double longitude_rad = 0.0;
    double height_m      = 0.0;

    double gravity_override   = 0.0;
    double earth_rate_override = 0.0;

    double gyro_arw_rad_per_sqrt_sec = 0.001;
};

struct CalibrationAndAlignmentOutput {
    AccelCalibrationResult accel_calib;
    GyroBiasCalibrationResult gyro_calib;

    AlignmentResult alignment;

    bool has_accel_calib = false;
    bool has_gyro_calib  = false;
    bool has_alignment   = false;

    bool gyro_used_angular_position = false;
};

inline CalibrationAndAlignmentOutput runCalibrationAndAlignment(
        const CalibrationAndAlignmentInput& input) {

    CalibrationAndAlignmentOutput out;

    const double dt       = input.imu_dt;
    const double align_dt = (input.align_dt > 0.0) ? input.align_dt : dt;

    const bool has_sixpos = !input.pos1_z_up.empty()
                         && !input.pos2_z_down.empty()
                         && !input.pos3_y_up.empty()
                         && !input.pos4_y_down.empty()
                         && !input.pos5_x_up.empty()
                         && !input.pos6_x_down.empty();

    if (has_sixpos) {
        const EarthParameters earth = [&]() {
            const LocalPosition pos{input.latitude_rad, input.longitude_rad, input.height_m};
            const LocalVelocity vel{0.0, 0.0, 0.0};
            return updateEarthParameters(pos, vel);
        }();
        const double local_gravity  = (input.gravity_override > 0.0)
                                      ? input.gravity_override : earth.gravity;
        const double local_earth_rate = (input.earth_rate_override > 0.0)
                                        ? input.earth_rate_override : 7.292115e-5;

        SixPositionData six_pos;
        six_pos.gravity = local_gravity;

        six_pos.pos1_z_up   = averageDvelFromStaticFrames(input.pos1_z_up);
        six_pos.pos2_z_down = averageDvelFromStaticFrames(input.pos2_z_down);
        six_pos.pos3_y_up   = averageDvelFromStaticFrames(input.pos3_y_up);
        six_pos.pos4_y_down = averageDvelFromStaticFrames(input.pos4_y_down);
        six_pos.pos5_x_up   = averageDvelFromStaticFrames(input.pos5_x_up);
        six_pos.pos6_x_down = averageDvelFromStaticFrames(input.pos6_x_down);

        auto divByDt = [dt](Vec3& v) {
            v.x /= dt; v.y /= dt; v.z /= dt;
        };
        divByDt(six_pos.pos1_z_up);
        divByDt(six_pos.pos2_z_down);
        divByDt(six_pos.pos3_y_up);
        divByDt(six_pos.pos4_y_down);
        divByDt(six_pos.pos5_x_up);
        divByDt(six_pos.pos6_x_down);

        out.accel_calib = calibrateAccelerometerSixPosition(six_pos);
        out.has_accel_calib = true;

        const double omega_v = local_earth_rate * std::sin(input.latitude_rad);

        auto meanGyroRate = [dt](const std::vector<ImuData>& frames) -> Vec3 {
            if (frames.empty()) return {0, 0, 0};
            Vec3 sum{0, 0, 0};
            for (const auto& imu : frames) {
                sum.x += imu.dtheta_x;
                sum.y += imu.dtheta_y;
                sum.z += imu.dtheta_z;
            }
            double inv = 1.0 / static_cast<double>(frames.size());
            return {sum.x * inv / dt, sum.y * inv / dt, sum.z * inv / dt};
        };

        const Vec3 gyro_rates[6] = {
            meanGyroRate(input.pos1_z_up), meanGyroRate(input.pos2_z_down),
            meanGyroRate(input.pos3_y_up), meanGyroRate(input.pos4_y_down),
            meanGyroRate(input.pos5_x_up), meanGyroRate(input.pos6_x_down)
        };

        auto getRate = [](const Vec3& v, int axis) -> double {
            if (axis == 0) return v.x; if (axis == 1) return v.y; return v.z;
        };

        const int up_idx[3]   = {4, 2, 0};
        const int down_idx[3] = {5, 3, 1};

        for (int j = 0; j < 3; ++j) {
            out.gyro_calib.bias[j] = (getRate(gyro_rates[up_idx[j]], j)
                                   +  getRate(gyro_rates[down_idx[j]], j)) * 0.5;
        }

        for (int j = 0; j < 3; ++j) {
            double var = 0.0;
            for (int i = 0; i < 6; ++i) {
                double d = getRate(gyro_rates[i], j) - out.gyro_calib.bias[j];
                var += d * d;
            }
            out.gyro_calib.bias_std[j] = std::sqrt(var / 5.0);
        }

        const bool has_rotation_data = !input.rot1_zp.empty()
                                    && !input.rot2_zn.empty()
                                    && !input.rot3_yp.empty()
                                    && !input.rot4_yn.empty()
                                    && !input.rot5_xp.empty()
                                    && !input.rot6_xn.empty();

        if (has_rotation_data) {
            out.gyro_used_angular_position = true;

            struct {
                const std::vector<ImuData>* file_a;
                const std::vector<ImuData>* file_b;
                int axis;
            } const rot_pairs[3] = {
                {&input.rot5_xp, &input.rot6_xn, 0},
                {&input.rot3_yp, &input.rot4_yn, 1},
                {&input.rot1_zp, &input.rot2_zn, 2},
            };

            const double alpha = 2.0 * std::acos(-1.0);

            for (int i = 0; i < 3; ++i) {
                const auto& rp = rot_pairs[i];
                const int ax = rp.axis;

                const std::size_t n = std::min(rp.file_a->size(), rp.file_b->size());
                double sigma_a = 0, sigma_b = 0;
                for (std::size_t k = 0; k < n; ++k) {
                    const auto& da = (*rp.file_a)[k];
                    const auto& db = (*rp.file_b)[k];
                    if (ax == 0) {
                        sigma_a += da.dtheta_x;  sigma_b += db.dtheta_x;
                    } else if (ax == 1) {
                        sigma_a += da.dtheta_y;  sigma_b += db.dtheta_y;
                    } else {
                        sigma_a += da.dtheta_z;  sigma_b += db.dtheta_z;
                    }
                }

                const double sigma_pos = (sigma_a > sigma_b) ? sigma_a : sigma_b;
                const double sigma_neg = (sigma_a > sigma_b) ? sigma_b : sigma_a;

                const double k_val = (sigma_pos - sigma_neg) / (2.0 * alpha);
                out.gyro_calib.scale[ax] = k_val - 1.0;
            }

            {
                const char* axis_name[3] = {"X", "Y", "Z"};
                const double d2r = rad2deg(1.0);
                std::cout << "\n[Gyro Scale: Angular Position Method]  (1.3.2 Eq 8-11, equal-length)\n";
                std::cout << "  Axis | Sigma_pos(deg) | Sigma_neg(deg) | Sigma_p-Sigma_n(deg) | k       | scale(ppm)\n";
                for (int i = 0; i < 3; ++i) {
                    const auto& rp = rot_pairs[i];
                    const int ax = rp.axis;
                    const std::size_t n = std::min(rp.file_a->size(), rp.file_b->size());
                    double sa = 0, sb = 0;
                    for (std::size_t k = 0; k < n; ++k) {
                        if (ax == 0) {
                            sa += (*rp.file_a)[k].dtheta_x;  sb += (*rp.file_b)[k].dtheta_x;
                        } else if (ax == 1) {
                            sa += (*rp.file_a)[k].dtheta_y;  sb += (*rp.file_b)[k].dtheta_y;
                        } else {
                            sa += (*rp.file_a)[k].dtheta_z;  sb += (*rp.file_b)[k].dtheta_z;
                        }
                    }
                    double sp = (sa > sb) ? sa : sb;
                    double sn = (sa > sb) ? sb : sa;
                    std::cout << "    " << axis_name[ax]
                              << "   | " << std::setw(10) << std::fixed << std::setprecision(2) << sp * d2r
                              << " | " << std::setw(10) << sn * d2r
                              << " | " << std::setw(10) << (sp - sn) * d2r
                              << " | " << std::setw(7)  << std::setprecision(6) << (out.gyro_calib.scale[ax] + 1.0)
                              << " | " << std::setw(10) << std::setprecision(1) << out.gyro_calib.scale[ax] * 1e6 << "\n";
                }
                std::cout << "  (Sigma_pos-Sigma_neg ~= 720deg; k ~= 1+scale_err)\n";
                std::cout << "  Bias from static two-position, NOT from this angular position data.\n\n";
            }

        } else {
            out.gyro_used_angular_position = false;

            for (int j = 0; j < 3; ++j) {
                const double diff = getRate(gyro_rates[up_idx[j]], j)
                                  - getRate(gyro_rates[down_idx[j]], j);
                const double denom = 2.0 * omega_v;
                out.gyro_calib.scale[j] = (std::fabs(denom) > 1e-15)
                    ? (diff / denom - 1.0) : 0.0;
            }
        }

        {
            const char* axis_name[3] = {"X", "Y", "Z"};
            const double d2r = rad2deg(1.0);
            std::cout << "[Gyro Bias: Static Two-Position Method]  (1.3.2 Eq 12-14, 3min avg)\n";
            std::cout << "  omega_v = " << omega_v * d2r * 3600.0 << " deg/h\n";
            std::cout << "  Axis |  up(deg/h)  | down(deg/h) | bias(deg/h) | bias_std(deg/h)\n";
            for (int j = 0; j < 3; ++j) {
                const double up   = getRate(gyro_rates[up_idx[j]],   j) * d2r * 3600.0;
                const double down = getRate(gyro_rates[down_idx[j]], j) * d2r * 3600.0;
                std::cout << "    " << axis_name[j]
                          << "   | " << std::setw(10) << std::fixed << std::setprecision(4) << up
                          << " | " << std::setw(11) << down
                          << " | " << std::setw(10) << out.gyro_calib.bias[j] * d2r * 3600.0
                          << " | " << std::setw(13) << out.gyro_calib.bias_std[j] * d2r * 3600.0 << "\n";
            }
            std::cout << "  (bias = (up+down)/2, Earth rate cancels)\n\n";
        }

        out.has_gyro_calib = true;
    }

    if (!input.align_static.empty()) {
        const Vec3 dvel_mean   = averageDvelFromStaticFrames(input.align_static);
        const Vec3 dtheta_mean = averageDthetaFromStaticFrames(input.align_static);

        const Vec3 f_b = {dvel_mean.x   / align_dt, dvel_mean.y   / align_dt, dvel_mean.z   / align_dt};
        const Vec3 w_b = {dtheta_mean.x / align_dt, dtheta_mean.y / align_dt, dtheta_mean.z / align_dt};

        out.alignment = analyticCoarseAlignment(
            f_b, w_b,
            input.latitude_rad,
            input.height_m,
            input.gyro_arw_rad_per_sqrt_sec,
            input.gravity_override,
            input.earth_rate_override);
        out.has_alignment = true;
    }

    return out;
}


// =========================================================================
// 5. 诊断报告输出
// =========================================================================

inline void printCalibrationReport(const CalibrationAndAlignmentOutput& result) {
    std::cout << "\n";
    std::cout << "============================================================\n";
    std::cout << "       IMU Error Calibration & Initial Alignment Report\n";
    std::cout << "============================================================\n\n";

    if (!result.has_accel_calib && !result.has_gyro_calib && !result.has_alignment) {
        std::cout << "  No calibration or alignment was performed.\n";
        std::cout << "  (No data was provided for either function.)\n\n";
        std::cout << "============================================================\n\n";
        return;
    }

    if (result.has_accel_calib) {
    std::cout << "[Six-Position Accelerometer Calibration]\n";
    std::cout << "  Bias (m/s^2):  "
              << result.accel_calib.bias[0] << ", "
              << result.accel_calib.bias[1] << ", "
              << result.accel_calib.bias[2] << "\n";
    std::cout << "  Bias (mGal):   "
              << result.accel_calib.bias[0] * 1e5 << ", "
              << result.accel_calib.bias[1] * 1e5 << ", "
              << result.accel_calib.bias[2] * 1e5 << "\n";
    std::cout << "  Calibration Matrix M (f_meas = b + M * f_true):\n";
    for (int r = 0; r < 3; ++r) {
        std::cout << "    ["
                  << std::setw(10) << std::fixed << std::setprecision(6) << result.accel_calib.M.m[r][0] << "  "
                  << std::setw(10) << result.accel_calib.M.m[r][1] << "  "
                  << std::setw(10) << result.accel_calib.M.m[r][2] << " ]\n";
    }
    std::cout << "  Scale Error (ppm):       "
              << (result.accel_calib.M.m[0][0] - 1.0) * 1e6 << ", "
              << (result.accel_calib.M.m[1][1] - 1.0) * 1e6 << ", "
              << (result.accel_calib.M.m[2][2] - 1.0) * 1e6 << "\n";
    std::cout << "  Cross-axis (ppm, max):   "
              << std::max({std::fabs(result.accel_calib.M.m[0][1]),
                           std::fabs(result.accel_calib.M.m[0][2]),
                           std::fabs(result.accel_calib.M.m[1][0]),
                           std::fabs(result.accel_calib.M.m[1][2]),
                           std::fabs(result.accel_calib.M.m[2][0]),
                           std::fabs(result.accel_calib.M.m[2][1])}) * 1e6 << "\n";
    std::cout << "  Residual RMS (m/s^2):    " << result.accel_calib.residual_rms << "\n\n";
    }

    if (result.has_gyro_calib) {
    std::cout << (result.gyro_used_angular_position
                     ? "[Angular Position Gyro Calibration]  (±360° rotation)\n"
                     : "[Two-Position Gyro Calibration]\n");
    std::cout << "  Bias (deg/h):  "
              << result.gyro_calib.bias[0] * rad2deg(1.0) * 3600.0 << ", "
              << result.gyro_calib.bias[1] * rad2deg(1.0) * 3600.0 << ", "
              << result.gyro_calib.bias[2] * rad2deg(1.0) * 3600.0 << "\n";
    std::cout << "  Bias std (deg/h): "
              << result.gyro_calib.bias_std[0] * rad2deg(1.0) * 3600.0 << ", "
              << result.gyro_calib.bias_std[1] * rad2deg(1.0) * 3600.0 << ", "
              << result.gyro_calib.bias_std[2] * rad2deg(1.0) * 3600.0 << "\n";
    std::cout << "  Scale Error (ppm): "
              << result.gyro_calib.scale[0] * 1e6 << ", "
              << result.gyro_calib.scale[1] * 1e6 << ", "
              << result.gyro_calib.scale[2] * 1e6;
    if (std::fabs(result.gyro_calib.scale[0]) > 0.05 ||
        std::fabs(result.gyro_calib.scale[1]) > 0.05 ||
        std::fabs(result.gyro_calib.scale[2]) > 0.05) {
        std::cout << "  (note: MEMS noise may dominate, rate-table verification recommended)";
    }
    std::cout << "\n\n";
    }

    if (result.has_alignment) {
    std::cout << "[Analytic Coarse Alignment (TRIAD)]\n";
    std::cout << "  Roll  (deg):  " << rad2deg(result.alignment.euler[0]) << "\n";
    std::cout << "  Pitch (deg):  " << rad2deg(result.alignment.euler[1]) << "\n";
    std::cout << "  Yaw   (deg):  " << rad2deg(result.alignment.euler[2]) << "\n";
    std::cout << "  Valid:         " << (result.alignment.valid ? "YES" : "NO") << "\n";
    std::cout << "  Condition #:   " << result.alignment.condition_number
              << (result.alignment.condition_number > 100.0 ? " (near-singular, yaw unreliable)" : "") << "\n";
    }

    std::cout << "\n============================================================\n\n";
}


// =========================================================================
// 6. 标定数据导出（供 Python 画图使用）
// =========================================================================

inline void saveCalibDataForPlot(
        const std::string& out_dir,
        const CalibrationAndAlignmentInput& input,
        const CalibrationAndAlignmentOutput& result)
{
    {
        std::string cmd;
#ifdef _WIN32
        cmd = "mkdir \"" + out_dir + "\" 2>nul";
#else
        cmd = "mkdir -p \"" + out_dir + "\" 2>/dev/null";
#endif
        std::system(cmd.c_str());
    }

    const double dt       = input.imu_dt;
    const double align_dt = (input.align_dt > 0.0) ? input.align_dt : dt;

    int files_written = 0;

    if (result.has_accel_calib) {
        std::string path = out_dir + "/calib_six_position.txt";
        std::ofstream ofs(path);
        if (!ofs) {
            std::cerr << "[WARNING] saveCalibDataForPlot: cannot write " << path << "\n";
        } else {
            ofs << std::setprecision(8) << std::fixed;
            const auto avgDvel = [&](const std::vector<ImuData>& frames) -> Vec3 {
                if (frames.empty()) return {0, 0, 0};
                Vec3 s{0, 0, 0};
                for (const auto& imu : frames) { s.x += imu.dvel_x; s.y += imu.dvel_y; s.z += imu.dvel_z; }
                double inv = 1.0 / static_cast<double>(frames.size());
                return {s.x * inv / dt, s.y * inv / dt, s.z * inv / dt};
            };

            const double g = (input.gravity_override > 0.0) ? input.gravity_override : [&]() {
                LocalPosition pos{input.latitude_rad, input.longitude_rad, input.height_m};
                LocalVelocity vel{0, 0, 0};
                return updateEarthParameters(pos, vel).gravity;
            }();

            const Vec3 true_f[6] = {
                {0, 0, g}, {0, 0, -g}, {0, g, 0}, {0, -g, 0}, {g, 0, 0}, {-g, 0, 0}
            };
            const Vec3 meas_f[6] = {
                avgDvel(input.pos1_z_up), avgDvel(input.pos2_z_down),
                avgDvel(input.pos3_y_up), avgDvel(input.pos4_y_down),
                avgDvel(input.pos5_x_up), avgDvel(input.pos6_x_down)
            };

            const int W = 12;
            ofs << "# gravity      = " << g << "\n";
            ofs << "# bias (m/s^2) = "
                << std::setw(W) << result.accel_calib.bias[0]
                << std::setw(W) << result.accel_calib.bias[1]
                << std::setw(W) << result.accel_calib.bias[2] << "\n";
            ofs << "# M matrix:\n";
            for (int r = 0; r < 3; ++r) {
                ofs << "#        [" << std::setw(W) << result.accel_calib.M.m[r][0]
                              << std::setw(W) << result.accel_calib.M.m[r][1]
                              << std::setw(W) << result.accel_calib.M.m[r][2] << " ]\n";
            }
            ofs << "# residual_rms = " << result.accel_calib.residual_rms << " m/s^2\n";
            ofs << "# columns:\n";
            ofs << "# " << std::setw(7) << "pos"
                << std::setw(W) << "true_fx" << std::setw(W) << "true_fy" << std::setw(W) << "true_fz"
                << std::setw(W) << "meas_fx" << std::setw(W) << "meas_fy" << std::setw(W) << "meas_fz"
                << std::setw(W) << "resid_fx" << std::setw(W) << "resid_fy" << std::setw(W) << "resid_fz" << "\n";

            const auto& M = result.accel_calib.M.m;
            const double bx = result.accel_calib.bias[0];
            const double by = result.accel_calib.bias[1];
            const double bz = result.accel_calib.bias[2];
            const char* labels[6] = {"z_up", "z_down", "y_up", "y_down", "x_up", "x_down"};
            for (int i = 0; i < 6; ++i) {
                const double pred_x = bx + M[0][0]*true_f[i].x + M[0][1]*true_f[i].y + M[0][2]*true_f[i].z;
                const double pred_y = by + M[1][0]*true_f[i].x + M[1][1]*true_f[i].y + M[1][2]*true_f[i].z;
                const double pred_z = bz + M[2][0]*true_f[i].x + M[2][1]*true_f[i].y + M[2][2]*true_f[i].z;
                ofs << "  " << std::setw(5) << std::left << labels[i] << std::right
                    << std::setw(W) << true_f[i].x << std::setw(W) << true_f[i].y << std::setw(W) << true_f[i].z
                    << std::setw(W) << meas_f[i].x << std::setw(W) << meas_f[i].y << std::setw(W) << meas_f[i].z
                    << std::setw(W) << (meas_f[i].x - pred_x)
                    << std::setw(W) << (meas_f[i].y - pred_y)
                    << std::setw(W) << (meas_f[i].z - pred_z) << "\n";
            }
        }
        ++files_written;
    }

    if (!input.align_static.empty()) {
        std::string path = out_dir + "/calib_gyro_static.txt";
        std::ofstream ofs(path);
        if (!ofs) {
            std::cerr << "[WARNING] saveCalibDataForPlot: cannot write " << path << "\n";
        } else {
            const double deg_per_rad = rad2deg(1.0);
            const double deg_h_per_rad_s = deg_per_rad * 3600.0;
            const int W = 14;
            ofs << std::setprecision(6) << std::fixed;
            ofs << "# dt             = " << align_dt << " s\n";
            ofs << "# gyro_bias      = "
                << std::setw(W) << result.gyro_calib.bias[0] * deg_h_per_rad_s
                << std::setw(W) << result.gyro_calib.bias[1] * deg_h_per_rad_s
                << std::setw(W) << result.gyro_calib.bias[2] * deg_h_per_rad_s << "  deg/h\n";
            ofs << "# gyro_bias_std  = "
                << std::setw(W) << result.gyro_calib.bias_std[0] * deg_h_per_rad_s
                << std::setw(W) << result.gyro_calib.bias_std[1] * deg_h_per_rad_s
                << std::setw(W) << result.gyro_calib.bias_std[2] * deg_h_per_rad_s << "  deg/h\n";
            ofs << "# gyro_scale     = "
                << std::setw(W) << result.gyro_calib.scale[0] * 1e6
                << std::setw(W) << result.gyro_calib.scale[1] * 1e6
                << std::setw(W) << result.gyro_calib.scale[2] * 1e6 << "  ppm\n";
            ofs << "# columns:\n";
            ofs << "# " << std::setw(8) << "index"
                << std::setw(W) << "dtheta_x" << std::setw(W) << "dtheta_y"
                << std::setw(W) << "dtheta_z" << "  (deg)\n";
            ofs << std::setprecision(8);
            for (std::size_t i = 0; i < input.align_static.size(); ++i) {
                const auto& imu = input.align_static[i];
                ofs << "  " << std::setw(8) << i
                    << std::setw(W) << imu.dtheta_x * deg_per_rad
                    << std::setw(W) << imu.dtheta_y * deg_per_rad
                    << std::setw(W) << imu.dtheta_z * deg_per_rad << "\n";
            }
        }
        ++files_written;
    }

    if (result.has_alignment) {
        std::string path = out_dir + "/calib_alignment.txt";
        std::ofstream ofs(path);
        if (!ofs) {
            std::cerr << "[WARNING] saveCalibDataForPlot: cannot write " << path << "\n";
        } else {
            const int W = 10;
            ofs << std::setprecision(6) << std::fixed;
            ofs << "# columns:\n";
            ofs << "# " << std::setw(W) << "roll" << std::setw(W) << "pitch"
                << std::setw(W) << "yaw" << std::setw(8) << "valid  "
                << std::setw(W) << "cond_num  (deg, deg, deg, bool, -)\n";
            ofs << "  " << std::setw(W) << rad2deg(result.alignment.euler[0])
                << std::setw(W) << rad2deg(result.alignment.euler[1])
                << std::setw(W) << rad2deg(result.alignment.euler[2])
                << std::setw(8) << (result.alignment.valid ? 1 : 0)
                << std::setw(W) << result.alignment.condition_number << "\n";
        }
        ++files_written;
    }

    if (result.has_accel_calib || result.has_gyro_calib) {
        std::string path = out_dir + "/calib_compensated_sixpos.txt";
        std::ofstream ofs(path);
        if (!ofs) {
            std::cerr << "[WARNING] saveCalibDataForPlot: cannot write " << path << "\n";
        } else {
            const int WT = 13, WD = 13;
            ofs << std::setprecision(9) << std::fixed;
            ofs << "# Six-position data, every frame compensated by calibration results\n";
            ofs << "# dt=" << dt << " s,  dtheta=rad(increment), dvel=m/s(increment)\n";
            ofs << "# columns:\n";
            ofs << "# " << std::setw(6) << "pos"
                << std::setw(WT) << "time"
                << std::setw(WD) << "dtheta_x" << std::setw(WD) << "dtheta_y" << std::setw(WD) << "dtheta_z"
                << std::setw(WD) << "dvel_x"   << std::setw(WD) << "dvel_y"   << std::setw(WD) << "dvel_z" << "\n";
            ofs << std::setprecision(3) << std::fixed;

            const char* pos_labels[6] = {"z_up","z_down","y_up","y_down","x_up","x_down"};
            const std::vector<ImuData>* pos_data[6] = {
                &input.pos1_z_up, &input.pos2_z_down,
                &input.pos3_y_up, &input.pos4_y_down,
                &input.pos5_x_up, &input.pos6_x_down
            };

            for (int i = 0; i < 6; ++i) {
                const auto& frames = *pos_data[i];
                if (frames.empty()) continue;
                ofs << "# ---- " << pos_labels[i] << " (" << frames.size() << " frames) ----\n";
                for (const auto& imu : frames) {
                    ImuData comp = compensateImuData(imu, result.accel_calib,
                                                     result.gyro_calib, dt);
                    ofs << "  " << std::setw(5) << std::left << pos_labels[i] << std::right
                        << std::setw(WT) << std::setprecision(3) << imu.time
                        << std::setw(WD) << std::setprecision(9) << comp.dtheta_x
                        << std::setw(WD) << comp.dtheta_y
                        << std::setw(WD) << comp.dtheta_z
                        << std::setw(WD) << comp.dvel_x
                        << std::setw(WD) << comp.dvel_y
                        << std::setw(WD) << comp.dvel_z << "\n";
                }
            }
        }
        ++files_written;
    }

    if (result.has_gyro_calib && result.gyro_used_angular_position) {
        std::string path = out_dir + "/calib_compensated_rot.txt";
        std::ofstream ofs(path);
        if (!ofs) {
            std::cerr << "[WARNING] saveCalibDataForPlot: cannot write " << path << "\n";
        } else {
            const int WT = 13, WD = 13;
            ofs << std::setprecision(9) << std::fixed;
            ofs << "# Rotation data, every frame compensated by gyro calibration results\n";
            ofs << "# dt=" << dt << " s,  dtheta=rad(increment), dvel=m/s(increment)\n";
            ofs << "# columns:\n";
            ofs << "# " << std::setw(7) << "file"
                << std::setw(WT) << "time"
                << std::setw(WD) << "dtheta_x" << std::setw(WD) << "dtheta_y" << std::setw(WD) << "dtheta_z"
                << std::setw(WD) << "dvel_x"   << std::setw(WD) << "dvel_y"   << std::setw(WD) << "dvel_z" << "\n";
            ofs << std::setprecision(3) << std::fixed;

            struct {
                const std::vector<ImuData>* file_a;
                const std::vector<ImuData>* file_b;
                const char* name_a;
                const char* name_b;
            } const rot_pairs[3] = {
                {&input.rot5_xp, &input.rot6_xn, "x+360", "x-360"},
                {&input.rot3_yp, &input.rot4_yn, "y+360", "y-360"},
                {&input.rot1_zp, &input.rot2_zn, "z+360", "z-360"},
            };

            for (int i = 0; i < 3; ++i) {
                const auto& rp = rot_pairs[i];
                for (int j = 0; j < 2; ++j) {
                    const auto& frames = *(j == 0 ? rp.file_a : rp.file_b);
                    const char* fname = (j == 0 ? rp.name_a : rp.name_b);
                    if (frames.empty()) continue;
                    ofs << "# ---- " << fname << " (" << frames.size() << " frames) ----\n";
                    for (const auto& imu : frames) {
                        ImuData comp = compensateImuData(imu, result.accel_calib,
                                                         result.gyro_calib, dt);
                        ofs << "  " << std::setw(6) << std::left << fname << std::right
                            << std::setw(WT) << std::setprecision(3) << imu.time
                            << std::setw(WD) << std::setprecision(9) << comp.dtheta_x
                            << std::setw(WD) << comp.dtheta_y
                            << std::setw(WD) << comp.dtheta_z
                            << std::setw(WD) << comp.dvel_x
                            << std::setw(WD) << comp.dvel_y
                            << std::setw(WD) << comp.dvel_z << "\n";
                    }
                }
            }
        }
        ++files_written;
    }

    std::cout << "[saveCalibDataForPlot] " << files_written << " file(s) written to " << out_dir << "/\n";
}

}  // namespace ins
