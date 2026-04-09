#pragma once

#include <array>

// 作用：定义 OtherData 数据结构。
struct OtherData {
    // 作用：此处说明当前字段或步骤的用途。
    double gnss_imu_time_offset_s = 0.0;

    // 作用：此处说明当前字段或步骤的用途。
    std::array<double, 3> antenna_lever_arm_b = {0.0, 0.0, 0.0};
};
