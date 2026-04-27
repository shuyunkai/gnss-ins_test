#pragma once
#include <array>
#include <stdexcept>
#include "../common/Coordinate_2.h"
#include "../fileio/gnss.h"
#include "../fileio/imu.h"
#include "ins_main.h"
#include "kalman.h"
#include "ukf.h"
#include "graph_optimization.h"
namespace ins {

inline Mat3 transpose3(const Mat3& a) {
    Mat3 b;
    for(int i=0; i<3; ++i)
        for(int j=0; j<3; ++j)
            b.m[i][j] = a.m[j][i];
    return b;
}

inline Matrix extractErrorStateFrom21State(const NavigationStatusData& nominal, const NavigationStatusData& perturbed) {
    Matrix err(21, 1, 0.0);
    const LocalPosition pos_nom{nominal.pvacur_.blh[0], nominal.pvacur_.blh[1], nominal.pvacur_.blh[2]};
    const LocalVelocity vel_nom{nominal.pvacur_.vel_n[0], nominal.pvacur_.vel_n[1], nominal.pvacur_.vel_n[2]};
    const EarthParameters earth = updateEarthParameters(pos_nom, vel_nom);
    const double cos_lat = std::cos(nominal.pvacur_.blh[0]);
    const double safe_cos_lat = (std::fabs(cos_lat) < 1e-8) ? ((cos_lat >= 0.0) ? 1e-8 : -1e-8) : cos_lat;

    err(ErrorStateIndex21::kPos + 0, 0) = (nominal.pvacur_.blh[0] - perturbed.pvacur_.blh[0]) * (earth.rm + nominal.pvacur_.blh[2]);
    err(ErrorStateIndex21::kPos + 1, 0) = (nominal.pvacur_.blh[1] - perturbed.pvacur_.blh[1]) * ((earth.rn + nominal.pvacur_.blh[2]) * safe_cos_lat);
    err(ErrorStateIndex21::kPos + 2, 0) = perturbed.pvacur_.blh[2] - nominal.pvacur_.blh[2];

    err(ErrorStateIndex21::kVel + 0, 0) = nominal.pvacur_.vel_n[0] - perturbed.pvacur_.vel_n[0];
    err(ErrorStateIndex21::kVel + 1, 0) = nominal.pvacur_.vel_n[1] - perturbed.pvacur_.vel_n[1];
    err(ErrorStateIndex21::kVel + 2, 0) = nominal.pvacur_.vel_n[2] - perturbed.pvacur_.vel_n[2];

    Mat3 c_pert = quat2dcm(euler2quat(Euler{perturbed.pvacur_.euler[0], perturbed.pvacur_.euler[1], perturbed.pvacur_.euler[2]}));
    Mat3 c_nom = quat2dcm(euler2quat(Euler{nominal.pvacur_.euler[0], nominal.pvacur_.euler[1], nominal.pvacur_.euler[2]}));
    Mat3 phi_cross = matMul(c_pert, transpose3(c_nom));
    err(ErrorStateIndex21::kAtt + 0, 0) = phi_cross.m[2][1];
    err(ErrorStateIndex21::kAtt + 1, 0) = phi_cross.m[0][2];
    err(ErrorStateIndex21::kAtt + 2, 0) = phi_cross.m[1][0];

    err(ErrorStateIndex21::kGyroBias + 0, 0) = perturbed.imuerror_.gyro_bias[0] - nominal.imuerror_.gyro_bias[0];
    err(ErrorStateIndex21::kGyroBias + 1, 0) = perturbed.imuerror_.gyro_bias[1] - nominal.imuerror_.gyro_bias[1];
    err(ErrorStateIndex21::kGyroBias + 2, 0) = perturbed.imuerror_.gyro_bias[2] - nominal.imuerror_.gyro_bias[2];

    err(ErrorStateIndex21::kAccelBias + 0, 0) = perturbed.imuerror_.accel_bias[0] - nominal.imuerror_.accel_bias[0];
    err(ErrorStateIndex21::kAccelBias + 1, 0) = perturbed.imuerror_.accel_bias[1] - nominal.imuerror_.accel_bias[1];
    err(ErrorStateIndex21::kAccelBias + 2, 0) = perturbed.imuerror_.accel_bias[2] - nominal.imuerror_.accel_bias[2];

    err(ErrorStateIndex21::kGyroScale + 0, 0) = perturbed.imuerror_.gyro_scale[0] - nominal.imuerror_.gyro_scale[0];
    err(ErrorStateIndex21::kGyroScale + 1, 0) = perturbed.imuerror_.gyro_scale[1] - nominal.imuerror_.gyro_scale[1];
    err(ErrorStateIndex21::kGyroScale + 2, 0) = perturbed.imuerror_.gyro_scale[2] - nominal.imuerror_.gyro_scale[2];

    err(ErrorStateIndex21::kAccelScale + 0, 0) = perturbed.imuerror_.accel_scale[0] - nominal.imuerror_.accel_scale[0];
    err(ErrorStateIndex21::kAccelScale + 1, 0) = perturbed.imuerror_.accel_scale[1] - nominal.imuerror_.accel_scale[1];
    err(ErrorStateIndex21::kAccelScale + 2, 0) = perturbed.imuerror_.accel_scale[2] - nominal.imuerror_.accel_scale[2];
    return err;
}

template <typename FilterType>
class GnssInsFusionWorkflow {
public:
    GnssInsFusionWorkflow() : kf_(ErrorStateIndex21::kDim) {}
        void initialize(const NavigationStatusData& nav_init,
                    const Matrix& x0,
                    const Matrix& p0) {
        if (x0.rows() != ErrorStateIndex21::kDim || x0.cols() != 1) {
            throw std::invalid_argument("initialize: x0 must be 21x1");
        }
        if (p0.rows() != ErrorStateIndex21::kDim || p0.cols() != ErrorStateIndex21::kDim) {
            throw std::invalid_argument("initialize: p0 must be 21x21");
        }
        nav_state_ = nav_init;
        kf_.setState(x0);
        kf_.setCovariance(p0);
        initialized_ = true;
    }
        // 注意: 本函数只执行EKF predict/update,不负责INS机械编排
    // INS传播应由调用方在调用本函数前完成
    void processStep(
                     const Matrix& f,
                     const Matrix& q,
                     bool has_gnss,
                     const GnssData& gnss_data,
                     const std::array<double, 3>& antenna_lever_arm_b = {0.0, 0.0, 0.0}) {
        ensureInitialized();
        // 移除重复的propagateIns调用,避免状态被传播两次
        // propagateIns(imu_measure, nav_state_);  //  已删除,由main.cpp统一控制
        
        // C++17 'if constexpr' avoids SFINAE need for template method calls.
        if constexpr (std::is_same<FilterType, UnscentedKalmanFilter>::value) {
            this->doPredict(q);
        } else {
            this->doPredict(f, q);
        }
        if (!has_gnss) {
            nav_state_.pvapre_ = nav_state_.pvacur_;
            return;
        }
        const Matrix residual_z = buildPositionResidualInsMinusGnss(
            nav_state_.pvacur_,
            gnss_data,
            antenna_lever_arm_b);
        const Matrix r = buildPositionMeasurementNoiseFromGnssStd(
            gnss_data.std_x,
            gnss_data.std_y,
            gnss_data.std_z);
        std::cout << "\nGNSS RESIDUAL: " << residual_z(0,0) << " " << residual_z(1,0) << " " << residual_z(2,0) << std::endl; 
        gnssPositionUpdateAndFeedback21(kf_, residual_z, r, nav_state_, antenna_lever_arm_b, true);
        
        nav_state_.pvapre_ = nav_state_.pvacur_;
    }
        void processImuOnly(const ImuMeasureData& imu_measure,
                        const Matrix& f,
                        const Matrix& q) {
        ensureInitialized();
        propagateIns(imu_measure, nav_state_);
        doPredict(imu_measure, f, q);
        nav_state_.pvapre_ = nav_state_.pvacur_;
    }
        void processGnssUpdate(const GnssData& gnss_data,
                           const std::array<double, 3>& antenna_lever_arm_b = {0.0, 0.0, 0.0}) {
        ensureInitialized();
        const Matrix residual_z = buildPositionResidualInsMinusGnss(
            nav_state_.pvacur_,
            gnss_data,
            antenna_lever_arm_b);
        const Matrix r = buildPositionMeasurementNoiseFromGnssStd(
            gnss_data.std_x,
            gnss_data.std_y,
            gnss_data.std_z);
        std::cout << "\nGNSS RESIDUAL: " << residual_z(0,0) << " " << residual_z(1,0) << " " << residual_z(2,0) << std::endl; 
        gnssPositionUpdateAndFeedback21(kf_, residual_z, r, nav_state_, antenna_lever_arm_b, true);
        
        nav_state_.pvapre_ = nav_state_.pvacur_;
    }
        const NavigationStatusData& navState() const {
        return nav_state_;
    }
    // 作用：setNavState 函数（用于更新传播后的导航状态）
    void setNavState(const NavigationStatusData& new_state) {
        nav_state_ = new_state;
    }
        const FilterType& kalman() const {
        return kf_;
    }

    // ====== SFINAE Dispatch for predict ======
    template <typename T = FilterType>
    typename std::enable_if<!std::is_same<T, UnscentedKalmanFilter>::value>::type
    doPredict(const Matrix& f, const Matrix& q) {
        kf_.predict(f, q);
    }

    template <typename T = FilterType>
    typename std::enable_if<std::is_same<T, UnscentedKalmanFilter>::value>::type
    doPredict(const Matrix& q) {
        kf_.predict(Matrix(), q);
    }

private:
        static Blh estimateAntennaBlhFromIns(const PvaData& ins_pva,
                                         const std::array<double, 3>& antenna_lever_arm_b) {
        const Blh imu_blh{ins_pva.blh[0], ins_pva.blh[1], ins_pva.blh[2]};
        const Euler imu_euler{ins_pva.euler[0], ins_pva.euler[1], ins_pva.euler[2]};
        const Quaternion q_nb = euler2quat(imu_euler);
        const Mat3 c_nb = quat2dcm(q_nb);
        const Vec3 lever_b{antenna_lever_arm_b[0], antenna_lever_arm_b[1], antenna_lever_arm_b[2]};
        const Vec3 lever_n = matVec(c_nb, lever_b);
                return local2global(lever_n, imu_blh);
    }
        static Matrix buildPositionResidualInsMinusGnss(const PvaData& ins_pva,
                                                    const GnssData& gnss,
                                                    const std::array<double, 3>& antenna_lever_arm_b) {
        const Blh ins_antenna_blh = estimateAntennaBlhFromIns(ins_pva, antenna_lever_arm_b);
        const Blh gnss_blh{gnss.latitude, gnss.longitude, gnss.altitude};
        const Vec3 dr_gnss_minus_ins_n = global2local(gnss_blh, ins_antenna_blh);
                Matrix z(3, 1, 0.0);
        z(0, 0) = -dr_gnss_minus_ins_n.x;
        z(1, 0) = -dr_gnss_minus_ins_n.y;
        z(2, 0) = -dr_gnss_minus_ins_n.z;
        return z;
    }
        static Vec3 estimateOmegaIbBFromImu(const ImuMeasureData& imu_measure,
                                        const NavigationStatusData& nav_state) {
        const double dt = imu_measure.imucur_.time - imu_measure.imupre_.time;
        const double dt_safe = (dt > 1e-6) ? dt : 1e-6;
        const ImuData imu_corr = compensateSingleImuByError(imu_measure.imucur_, nav_state.imuerror_, dt_safe);
        return Vec3{imu_corr.dtheta_x / dt_safe, imu_corr.dtheta_y / dt_safe, imu_corr.dtheta_z / dt_safe};
    }
        void ensureInitialized() const {
        if (!initialized_) {
            throw std::runtime_error("GnssInsFusionWorkflow未初始化");
        }
    }
private:
    bool initialized_ = false;
    NavigationStatusData nav_state_;
    FilterType kf_;
};
}  // 命名空间结束
