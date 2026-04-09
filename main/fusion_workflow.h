#pragma once
#include <array>
#include <stdexcept>
#include "../common/Coordinate_2.h"
#include "../fileio/gnss.h"
#include "../fileio/imu.h"
#include "ins_main.h"
#include "kalman.h"
namespace ins {
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
    void processStep(const ImuMeasureData& imu_measure,
                     const Matrix& f,
                     const Matrix& q,
                     bool has_gnss,
                     const GnssData& gnss_data,
                     const std::array<double, 3>& antenna_lever_arm_b = {0.0, 0.0, 0.0}) {
        ensureInitialized();
        // 移除重复的propagateIns调用,避免状态被传播两次
        // propagateIns(imu_measure, nav_state_);  // ← 已删除,由main.cpp统一控制
        kf_.predict(f, q);
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
        kf_.predict(f, q);
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
        const StandardKalmanFilter& kalman() const {
        return kf_;
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
    StandardKalmanFilter kf_;
};
}  // 命名空间结束
