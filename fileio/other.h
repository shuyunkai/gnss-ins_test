#pragma once

#include <array>

struct OtherData {
        double gnss_imu_time_offset_s = 0.0;

        std::array<double, 3> antenna_lever_arm_b = {0.0, 0.0, 0.0};
};
