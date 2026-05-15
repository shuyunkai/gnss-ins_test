/**
 * IMU误差标定与初始对准 -- 独立可执行程序
 *
 * 编译（在项目根目录下，需先载入 MSVC 环境）：
 *   cl.exe /EHsc /std:c++17 main/run_calib.cpp /Fe:run_calib.exe
 *
 * 输入文件支持三种格式（自动检测）：
 *   1) 项目原生7列文本: time dtheta_x dtheta_y dtheta_z dvel_x dvel_y dvel_z
 *   2) NovAtel RAWIMUSA ASC文本: %RAWIMUSA,week,tow;...,6xIMUcounts*checksum
 *   3) NovAtel 二进制格式: 56字节定长记录（无需转码，直接读取）
 *
 * 六位置对应关系（z-y-x 轴顺序）：
 *   pos1: z轴朝上    pos2: z轴朝下
 *   pos3: y轴朝上    pos4: y轴朝下
 *   pos5: x轴朝上    pos6: x轴朝下
 */

#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "imu_calibration_alignment.h"
#include "../fileio/imufileloader.h"

// =========================================================================
// 【可配置参数】-- 根据实际采集环境和IMU硬件修改以下值
// =========================================================================

// ---- 六位置静态数据文件路径（必填） ----
// 支持 .ASC (RAWIMUSA)、项目原生文本、或二进制格式（自动检测）
#define CFG_POS1  "data_calib/Calibration/z_up_3min"         // 位置1: z轴朝上
#define CFG_POS2  "data_calib/Calibration/z_down_3min"       // 位置2: z轴朝下
#define CFG_POS3  "data_calib/Calibration/y_up_3min"         // 位置3: y轴朝上
#define CFG_POS4  "data_calib/Calibration/y_down_3min"       // 位置4: y轴朝下
#define CFG_POS5  "data_calib/Calibration/x_up_3min"         // 位置5: x轴朝上
#define CFG_POS6  "data_calib/Calibration/x_down_3min"       // 位置6: x轴朝下

// ---- 陀螺角位置法标定数据（绕轴360deg旋转） ----
#define CFG_ROT1_ZP  "data_calib/Calibration/z+360"          // 绕z轴正转360deg
#define CFG_ROT2_ZN  "data_calib/Calibration/z-360"          // 绕z轴反转360deg
#define CFG_ROT3_YP  "data_calib/Calibration/y+360"          // 绕y轴正转360deg
#define CFG_ROT4_YN  "data_calib/Calibration/y-360"          // 绕y轴反转360deg
#define CFG_ROT5_XP  "data_calib/Calibration/x+360"          // 绕x轴正转360deg
#define CFG_ROT6_XN  "data_calib/Calibration/x-360"          // 绕x轴反转360deg

// ---- 静态对准数据文件路径（必填） ----
#define CFG_ALIGN_STATIC  "data_calib/Align/Align_30min.txt"

// ---- IMU型号与比例因子 ----
// XW-GI7681 (100Hz):  acc=1.5258789063e-06, gyo=1.0850694444e-07
// NovAtel SPAN-100C (200Hz): acc=2.0e-08, gyo=1.0e-09
#define CFG_ACC_SCALE  1.5258789063e-06
#define CFG_GYO_SCALE  1.0850694444e-07

// ---- 对准IMU比例因子（若与六位置IMU不同，在此覆写） ----
#define CFG_ALIGN_ACC_SCALE  2.0e-08
#define CFG_ALIGN_GYO_SCALE  1.0e-09

// ---- 当地位置（必填） ----
#define CFG_LAT_DEG   30.531651244
#define CFG_LON_DEG  114.36
#define CFG_HEIGHT_M   28.2134

// ---- 重力加速度和地球自转角速率（可选） ----
#define CFG_GRAVITY    9.7936174
#define CFG_EARTH_RATE 7.292115e-5

// ---- IMU采样参数（必填） ----
#define CFG_IMU_DT     0.01
#define CFG_ALIGN_DT   0.005

// ---- IMU噪声参数（可选） ----
#define CFG_GYRO_ARW   0.0

// ---- 输出目录 ----
#define CFG_OUT_DIR    "data_calib/calib_out"

// =========================================================================
// 以下为程序逻辑，通常无需修改
// =========================================================================

namespace {

// ---- 二进制 IMU 记录结构（56字节定长） ----
struct RawImuBinaryRecord {
    unsigned char sync_header[16];
    double       gps_tow;
    unsigned     imu_status;
    int          accel_z;
    int          accel_y_neg;
    int          accel_x;
    int          gyro_z;
    int          gyro_y_neg;
    int          gyro_x;
    unsigned     crc32;
};
static_assert(sizeof(RawImuBinaryRecord) == 56, "RAWIMUS binary record must be 56 bytes");

inline ImuData rawCountsToImuData(double gps_tow,
                                   int accel_z, int accel_yn, int accel_x,
                                   int gyro_z,  int gyro_yn,  int gyro_x,
                                   double acc_scale, double gyo_scale) {
    ImuData imu;
    imu.time     = gps_tow;
    imu.dvel_x   = static_cast<double>(accel_x)  * acc_scale;
    imu.dvel_y   = static_cast<double>(-accel_yn) * acc_scale;
    imu.dvel_z   = static_cast<double>(accel_z)  * acc_scale;
    imu.dtheta_x = static_cast<double>(gyro_x)   * gyo_scale;
    imu.dtheta_y = static_cast<double>(-gyro_yn)  * gyo_scale;
    imu.dtheta_z = static_cast<double>(gyro_z)   * gyo_scale;
    return imu;
}

enum class DetectedFormat {
    kNativeText,
    kRawImuAsc,
    kRawImuBinary,
    kUnknown
};

std::streamoff findRawImuSync(const std::string& file_path, int max_search = 512) {
    std::ifstream f(file_path, std::ios::binary);
    if (!f.is_open()) return -1;
    std::vector<unsigned char> buf(max_search + 3);
    f.read(reinterpret_cast<char*>(buf.data()), buf.size());
    auto n = f.gcount();
    for (std::streamoff i = 0; i < n - 3; ++i) {
        if (buf[i] == 0xAA && buf[i+1] == 0x44 && buf[i+2] == 0x13) {
            return i;
        }
    }
    return -1;
}

DetectedFormat detectFileFormat(const std::string& file_path) {
    {
        std::streamoff sync_off = findRawImuSync(file_path);
        if (sync_off >= 0) {
            std::ifstream f(file_path, std::ios::binary);
            f.seekg(sync_off);
            RawImuBinaryRecord rec{};
            f.read(reinterpret_cast<char*>(&rec), sizeof(rec));
            if (f.gcount() == static_cast<std::streamsize>(sizeof(rec))) {
                double tow = rec.gps_tow;
                if (tow > 0 && tow < 604800.0) {
                    return DetectedFormat::kRawImuBinary;
                }
            }
        }
    }

    std::ifstream f(file_path);
    if (!f.is_open()) return DetectedFormat::kUnknown;

    std::string line;
    if (!std::getline(f, line)) return DetectedFormat::kUnknown;

    if (line.compare(0, 9, "%RAWIMUS") == 0) {
        return DetectedFormat::kRawImuAsc;
    }

    for (char c : line) {
        if (c == ' ' || c == '\t' || c == '\r') continue;
        if ((c >= '0' && c <= '9') || c == '-' || c == '+') {
            return DetectedFormat::kNativeText;
        }
        break;
    }
    return DetectedFormat::kUnknown;
}

std::vector<ImuData> loadRawImuAscFrames(const std::string& file_path,
                                          double acc_scale, double gyo_scale) {
    std::vector<ImuData> frames;
    std::ifstream f(file_path);
    if (!f.is_open()) {
        std::cerr << "[ERROR] Cannot open file: " << file_path << "\n";
        return frames;
    }
    frames.reserve(20000);

    std::string line;
    int skipped = 0;
    while (std::getline(f, line)) {
        if (line.empty()) continue;
        if (line[0] != '%') { ++skipped; continue; }

        auto star = line.find('*');
        if (star != std::string::npos) line = line.substr(0, star);

        std::vector<std::string> fields;
        {
            std::istringstream ss(line);
            std::string field;
            while (std::getline(ss, field, ',')) fields.push_back(field);
        }
        if (fields.size() < 11) { ++skipped; continue; }

        try {
            double tow = 0.0;
            {
                auto semi = fields[2].find(';');
                std::string ts = (semi != std::string::npos)
                    ? fields[2].substr(0, semi) : fields[2];
                auto comma = ts.find(',');
                if (comma != std::string::npos) ts = ts.substr(comma + 1);
                tow = std::stod(ts);
            }

            ImuData imu = rawCountsToImuData(
                tow,
                std::stoi(fields[5]),
                std::stoi(fields[6]),
                std::stoi(fields[7]),
                std::stoi(fields[8]),
                std::stoi(fields[9]),
                std::stoi(fields[10]),
                acc_scale, gyo_scale);

            if (imu.isFinite()) frames.push_back(imu);
            else ++skipped;
        } catch (...) { ++skipped; }
    }
    if (skipped > 0) std::cout << "  (skipped " << skipped << " lines)\n";
    return frames;
}

std::vector<ImuData> loadRawImuBinaryFrames(const std::string& file_path,
                                             double acc_scale, double gyo_scale) {
    std::vector<ImuData> frames;
    std::ifstream f(file_path, std::ios::binary);
    if (!f.is_open()) {
        std::cerr << "[ERROR] Cannot open file: " << file_path << "\n";
        return frames;
    }

    f.seekg(0, std::ios::end);
    std::streamoff size = f.tellg();
    if (size <= 0) return frames;
    f.seekg(0, std::ios::beg);
    std::vector<unsigned char> buf(static_cast<size_t>(size));
    f.read(reinterpret_cast<char*>(buf.data()), size);
    auto nbytes = static_cast<std::streamoff>(f.gcount());

    frames.reserve(static_cast<size_t>(nbytes / 56) + 1);

    const unsigned char SYNC0 = 0xAA, SYNC1 = 0x44, SYNC2 = 0x13;
    std::streamoff pos = 0;
    int skipped_sync = 0, skipped_gps = 0;

    while (pos + static_cast<std::streamoff>(sizeof(RawImuBinaryRecord)) <= nbytes) {
        if (buf[static_cast<size_t>(pos)] != SYNC0 ||
            buf[static_cast<size_t>(pos + 1)] != SYNC1 ||
            buf[static_cast<size_t>(pos + 2)] != SYNC2) {
            ++pos;
            ++skipped_sync;
            continue;
        }

        const RawImuBinaryRecord* rec =
            reinterpret_cast<const RawImuBinaryRecord*>(buf.data() + pos);
        double tow = rec->gps_tow;
        if (!std::isfinite(tow) || tow <= 0.0 || tow >= 604800.0) {
            ++pos;
            ++skipped_gps;
            continue;
        }

        ImuData imu = rawCountsToImuData(
            tow, rec->accel_z, rec->accel_y_neg, rec->accel_x,
            rec->gyro_z, rec->gyro_y_neg, rec->gyro_x,
            acc_scale, gyo_scale);

        if (imu.isFinite()) frames.push_back(imu);
        else ++skipped_gps;

        pos += static_cast<std::streamoff>(sizeof(RawImuBinaryRecord));
    }

    if (skipped_sync > 0) std::cout << "  (re-sync: skipped " << skipped_sync << " bytes)\n";
    if (skipped_gps > 0) std::cout << "  (invalid GPS: " << skipped_gps << " records)\n";
    return frames;
}

std::vector<ImuData> loadImuFramesAuto(const std::string& file_path,
                                        double acc_scale, double gyo_scale) {
    DetectedFormat fmt = detectFileFormat(file_path);
    std::vector<ImuData> frames;

    switch (fmt) {
    case DetectedFormat::kRawImuBinary:
        std::cout << "  [Binary RAWIMU] ";
        frames = loadRawImuBinaryFrames(file_path, acc_scale, gyo_scale);
        break;
    case DetectedFormat::kRawImuAsc:
        std::cout << "  [RAWIMUSA ASC]  ";
        frames = loadRawImuAscFrames(file_path, acc_scale, gyo_scale);
        break;
    case DetectedFormat::kNativeText:
        std::cout << "  [Native text]   ";
        {
            ImuFileLoader loader;
            if (loader.open(file_path)) {
                ImuData imu;
                while (loader.readNext(imu)) frames.push_back(imu);
            }
        }
        break;
    default:
        std::cerr << "  [ERROR] Unknown format: " << file_path << "\n";
        break;
    }

    std::cout << "  " << file_path << " -> " << frames.size() << " frames\n";
    return frames;
}

struct EstimatedScales {
    double acc_scale = 0.0;
    double gyo_scale = 0.0;
    bool acc_valid = false;
    bool gyo_valid = false;
};

EstimatedScales estimateScalesFromSixPosition(
        const std::string& z_up_path,
        const std::string& z_down_path,
        double g, double dt, double latitude_rad)
{
    EstimatedScales result;

    const int sample_frames = 500;

    auto loadMeanZ = [sample_frames](const std::string& path) -> std::pair<double, double> {
        auto frames = loadRawImuBinaryFrames(path, 1.0, 1.0);
        if (frames.empty()) return {0.0, 0.0};
        const int n = std::min(sample_frames, static_cast<int>(frames.size()));
        double sum_dvel_z = 0.0, sum_dtheta_z = 0.0;
        for (int i = 0; i < n; ++i) {
            sum_dvel_z   += frames[i].dvel_z;
            sum_dtheta_z += frames[i].dtheta_z;
        }
        return {sum_dvel_z / n, sum_dtheta_z / n};
    };

    auto [acc_z_up,   gyr_z_up]   = loadMeanZ(z_up_path);
    auto [acc_z_down, gyr_z_down] = loadMeanZ(z_down_path);

    if (acc_z_up == 0.0 && acc_z_down == 0.0) return result;

    {
        const double acc_diff = acc_z_up - acc_z_down;
        const double denom = 2.0 * g * dt;
        if (std::fabs(acc_diff) > 1.0) {
            result.acc_scale = denom / acc_diff;
            result.acc_valid = true;
        }
    }

    {
        const double omega_e = 7.292115e-5;
        const double omega_v = omega_e * std::sin(latitude_rad);
        const double gyr_diff = gyr_z_up - gyr_z_down;
        const double denom = 2.0 * omega_v * dt;
        if (std::fabs(gyr_diff) > 0.1) {
            result.gyo_scale = denom / gyr_diff;
            result.gyo_valid = true;
        }
    }

    return result;
}

void printUsage(const char* prog) {
    std::cout << "IMU Error Calibration & Initial Alignment Tool\n\n";
    std::cout << "Usage: " << prog << " [options]\n";
    std::cout << "  If no CLI options given, defaults in the source file will be used.\n\n";
    std::cout << "Options (override defaults):\n";
    std::cout << "  --pos1 .. --pos6 <file>   Six-position data files\n";
    std::cout << "  --rot1 .. --rot6 <file>   +/-360deg rotation data\n";
    std::cout << "  --align-static <file>     Static alignment data\n";
    std::cout << "  --lat  <deg>              Latitude (deg)\n";
    std::cout << "  --lon  <deg>              Longitude (deg)\n";
    std::cout << "  --height <m>              Height (m)\n";
    std::cout << "  --dt   <s>                IMU sampling interval (calibration)\n";
    std::cout << "  --align-dt <s>            IMU sampling interval (alignment)\n";
    std::cout << "  --acc-scale <s>           Accel scale factor (calibration)\n";
    std::cout << "  --gyo-scale <s>           Gyro scale factor (calibration)\n";
    std::cout << "  --align-acc-scale <s>     Accel scale factor (alignment)\n";
    std::cout << "  --align-gyo-scale <s>     Gyro scale factor (alignment)\n";
    std::cout << "  --gravity <m/s^2>       Local gravity (0=auto-compute)\n";
    std::cout << "  --earth-rate <rad/s>    Earth rotation rate (0=auto-compute)\n";
    std::cout << "  --arw  <rad/sqrt(s)>      Gyro ARW from datasheet\n";
    std::cout << "  --out  <dir>              Output directory\n";
    std::cout << "  --help                    Show this help\n";
    std::cout << "\nInput formats (auto-detected):\n";
    std::cout << "  1) Native text:  time dtheta_x dtheta_y dtheta_z dvel_x dvel_y dvel_z\n";
    std::cout << "  2) RAWIMUSA ASC: %RAWIMUSA,week,tow;...,6xIMUcounts*checksum\n";
    std::cout << "  3) Binary:        56-byte fixed-length records\n";
}

}  // namespace

int main(int argc, char** argv) {
    std::string pos_files[6] = {CFG_POS1, CFG_POS2, CFG_POS3, CFG_POS4, CFG_POS5, CFG_POS6};
    std::string rot_files[6] = {CFG_ROT1_ZP, CFG_ROT2_ZN, CFG_ROT3_YP, CFG_ROT4_YN, CFG_ROT5_XP, CFG_ROT6_XN};
    std::string align_static_file = CFG_ALIGN_STATIC;
    double latitude_deg  = CFG_LAT_DEG;
    double longitude_deg = CFG_LON_DEG;
    double height_m      = CFG_HEIGHT_M;
    double imu_dt        = CFG_IMU_DT;
    double align_dt      = CFG_ALIGN_DT;
    double gyro_arw      = CFG_GYRO_ARW;
    double acc_scale     = CFG_ACC_SCALE;
    double gyo_scale     = CFG_GYO_SCALE;
    double align_acc_scale = CFG_ALIGN_ACC_SCALE;
    double align_gyo_scale = CFG_ALIGN_GYO_SCALE;
    double gravity_override = CFG_GRAVITY;
    double earth_rate_override = CFG_EARTH_RATE;
    std::string out_dir  = CFG_OUT_DIR;
    bool has_arw = (CFG_GYRO_ARW > 0.0);

    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        auto needVal = [&](const char* name) -> std::string {
            if (i + 1 >= argc) {
                std::cerr << "[ERROR] " << name << " requires a value\n";
                std::exit(1);
            }
            return std::string(argv[++i]);
        };

        if (arg == "--help" || arg == "-h") {
            printUsage(argv[0]);
            return 0;
        } else if (arg == "--pos1") pos_files[0] = needVal("--pos1");
        else if (arg == "--pos2") pos_files[1] = needVal("--pos2");
        else if (arg == "--pos3") pos_files[2] = needVal("--pos3");
        else if (arg == "--pos4") pos_files[3] = needVal("--pos4");
        else if (arg == "--pos5") pos_files[4] = needVal("--pos5");
        else if (arg == "--pos6") pos_files[5] = needVal("--pos6");
        else if (arg == "--rot1") rot_files[0] = needVal("--rot1");
        else if (arg == "--rot2") rot_files[1] = needVal("--rot2");
        else if (arg == "--rot3") rot_files[2] = needVal("--rot3");
        else if (arg == "--rot4") rot_files[3] = needVal("--rot4");
        else if (arg == "--rot5") rot_files[4] = needVal("--rot5");
        else if (arg == "--rot6") rot_files[5] = needVal("--rot6");
        else if (arg == "--align-static") { align_static_file = needVal("--align-static"); }
        else if (arg == "--lat")    { latitude_deg  = std::stod(needVal("--lat")); }
        else if (arg == "--lon")    { longitude_deg = std::stod(needVal("--lon")); }
        else if (arg == "--height") { height_m = std::stod(needVal("--height")); }
        else if (arg == "--dt")     { imu_dt = std::stod(needVal("--dt")); }
        else if (arg == "--align-dt")  { align_dt = std::stod(needVal("--align-dt")); }
        else if (arg == "--acc-scale") { acc_scale = std::stod(needVal("--acc-scale")); }
        else if (arg == "--gyo-scale") { gyo_scale = std::stod(needVal("--gyo-scale")); }
        else if (arg == "--align-acc-scale") { align_acc_scale = std::stod(needVal("--align-acc-scale")); }
        else if (arg == "--align-gyo-scale") { align_gyo_scale = std::stod(needVal("--align-gyo-scale")); }
        else if (arg == "--gravity")    { gravity_override = std::stod(needVal("--gravity")); }
        else if (arg == "--earth-rate") { earth_rate_override = std::stod(needVal("--earth-rate")); }
        else if (arg == "--arw")    { gyro_arw = std::stod(needVal("--arw")); has_arw = true; }
        else if (arg == "--out")    { out_dir = needVal("--out"); }
        else {
            std::cerr << "[ERROR] Unknown option: " << arg << "\n";
            printUsage(argv[0]);
            return 1;
        }
    }

    if (latitude_deg < -90.0 || latitude_deg > 90.0) {
        std::cerr << "[ERROR] Valid latitude required (--lat or CFG_LAT_DEG)\n";
        return 1;
    }
    if (longitude_deg < -180.0 || longitude_deg > 180.0) {
        std::cerr << "[ERROR] Valid longitude required (--lon or CFG_LON_DEG)\n";
        return 1;
    }
    if (imu_dt <= 0.0) {
        std::cerr << "[ERROR] IMU sampling interval required (--dt or CFG_IMU_DT)\n";
        return 1;
    }

    std::cout << "\n========================================\n";
    std::cout << "  IMU Calibration & Alignment Tool\n";
    std::cout << "========================================\n";
    std::cout << "Location: lat=" << latitude_deg << " deg, lon=" << longitude_deg
              << " deg, height=" << height_m << " m\n";
    std::cout << "IMU dt:   " << imu_dt << " s (" << (1.0 / imu_dt) << " Hz)\n";
    std::cout << "Calib IMU: dt=" << imu_dt << "s, acc_scale=" << acc_scale
              << ", gyo_scale=" << gyo_scale << "\n";
    std::cout << "Align IMU: dt=" << (align_dt > 0 ? align_dt : imu_dt) << "s"
              << ", acc_scale=" << (align_acc_scale > 0 ? align_acc_scale : acc_scale)
              << ", gyo_scale=" << (align_gyo_scale > 0 ? align_gyo_scale : gyo_scale) << "\n";
    if (has_arw) {
        std::cout << "Gyro ARW: " << gyro_arw << " rad/sqrt(s)\n";
    } else {
        std::cout << "Gyro ARW: not specified, gyro SNR check will be skipped\n";
    }
    std::cout << "Gravity:  " << (gravity_override > 0 ? gravity_override : 0.0)
              << " m/s^2" << (gravity_override > 0 ? " (user)" : " (auto-compute)") << "\n";
    std::cout << "Earth rate: " << (earth_rate_override > 0 ? earth_rate_override : 7.292115e-5)
              << " rad/s" << (earth_rate_override > 0 ? " (user)" : " (default WGS84)") << "\n";
    std::cout << "Output:   " << out_dir << "/\n\n";

    if (acc_scale <= 0.0 || gyo_scale <= 0.0) {
        std::cout << "Auto-estimating IMU scale factors from data...\n";
        const double g_est = (gravity_override > 0.0) ? gravity_override : [&]() {
            ins::LocalPosition pos{ins::deg2rad(latitude_deg), ins::deg2rad(longitude_deg), height_m};
            ins::LocalVelocity vel{0, 0, 0};
            return ins::updateEarthParameters(pos, vel).gravity;
        }();
        auto est = estimateScalesFromSixPosition(
            pos_files[0], pos_files[1],
            g_est, imu_dt, ins::deg2rad(latitude_deg));

        if (est.acc_valid && acc_scale <= 0.0) {
            acc_scale = est.acc_scale;
            std::cout << "  Accel scale: auto-estimated = " << acc_scale
                      << " m/s per count\n";
        }
        if (est.gyo_valid && gyo_scale <= 0.0) {
            gyo_scale = est.gyo_scale;
            std::cout << "  Gyro scale:  auto-estimated = " << gyo_scale
                      << " rad per count";
            if (std::fabs(est.gyo_scale) < 1e-10 || std::fabs(est.gyo_scale) > 1e-5) {
                std::cout << "  (note: MEMS gyro scale from Earth rate is approximate)";
            }
            std::cout << "\n";
        }
        if (!est.acc_valid && !est.gyo_valid) {
            std::cout << "  [WARNING] Auto-estimation failed\n";
        }
    }

    std::cout << "Loading six-position data...\n";
    ins::CalibrationAndAlignmentInput input;
    input.pos1_z_up   = loadImuFramesAuto(pos_files[0], acc_scale, gyo_scale);
    input.pos2_z_down = loadImuFramesAuto(pos_files[1], acc_scale, gyo_scale);
    input.pos3_y_up   = loadImuFramesAuto(pos_files[2], acc_scale, gyo_scale);
    input.pos4_y_down = loadImuFramesAuto(pos_files[3], acc_scale, gyo_scale);
    input.pos5_x_up   = loadImuFramesAuto(pos_files[4], acc_scale, gyo_scale);
    input.pos6_x_down = loadImuFramesAuto(pos_files[5], acc_scale, gyo_scale);

    int empty_positions = 0;
    for (int i = 0; i < 6; ++i) {
        const auto& frames = (i == 0) ? input.pos1_z_up
                           : (i == 1) ? input.pos2_z_down
                           : (i == 2) ? input.pos3_y_up
                           : (i == 3) ? input.pos4_y_down
                           : (i == 4) ? input.pos5_x_up
                           : input.pos6_x_down;
        if (frames.empty()) ++empty_positions;
    }
    if (empty_positions == 6) {
        std::cout << "No six-position data provided, calibration will be skipped.\n";
    } else if (empty_positions > 0) {
        std::cerr << "[ERROR] Only " << (6 - empty_positions)
                  << "/6 positions have data. Six-position calibration requires all 6.\n";
        return 1;
    }

    std::cout << "\nLoading rotation data (+/-360deg)...\n";
    bool any_rotation = false;
    {
        input.rot1_zp = loadImuFramesAuto(rot_files[0], acc_scale, gyo_scale);
        input.rot2_zn = loadImuFramesAuto(rot_files[1], acc_scale, gyo_scale);
        input.rot3_yp = loadImuFramesAuto(rot_files[2], acc_scale, gyo_scale);
        input.rot4_yn = loadImuFramesAuto(rot_files[3], acc_scale, gyo_scale);
        input.rot5_xp = loadImuFramesAuto(rot_files[4], acc_scale, gyo_scale);
        input.rot6_xn = loadImuFramesAuto(rot_files[5], acc_scale, gyo_scale);
        any_rotation = !input.rot1_zp.empty() || !input.rot2_zn.empty()
                    || !input.rot3_yp.empty() || !input.rot4_yn.empty()
                    || !input.rot5_xp.empty() || !input.rot6_xn.empty();
    }
    if (any_rotation) {
        std::cout << "  Angular position method will be used for gyro calibration.\n";
    } else {
        std::cout << "  No rotation data found, will fall back to two-position method.\n";
    }

    std::cout << "\nLoading static alignment data...\n";
    {
        double a_asc = (align_acc_scale > 0.0) ? align_acc_scale : acc_scale;
        double a_gsc = (align_gyo_scale > 0.0) ? align_gyo_scale : gyo_scale;
        input.align_static = loadImuFramesAuto(align_static_file, a_asc, a_gsc);
    }
    if (input.align_static.empty()) {
        std::cout << "No static alignment data provided, TRIAD alignment will be skipped.\n";
    }

    if (empty_positions == 6 && input.align_static.empty()) {
        std::cerr << "[ERROR] No data provided. At least one of calibration or alignment data is required.\n";
        return 1;
    }

    input.latitude_rad  = ins::deg2rad(latitude_deg);
    input.longitude_rad = ins::deg2rad(longitude_deg);
    input.height_m      = height_m;
    input.imu_dt        = imu_dt;
    input.align_dt      = align_dt;
    input.gravity_override    = gravity_override;
    input.earth_rate_override = earth_rate_override;
    input.gyro_arw_rad_per_sqrt_sec = gyro_arw;

    std::cout << "\nRunning calibration and alignment...\n";
    const auto result = ins::runCalibrationAndAlignment(input);

    ins::printCalibrationReport(result);
    ins::saveCalibDataForPlot(out_dir, input, result);

    std::cout << "\nTo generate plots, run:\n";
    std::cout << "  python main/calib_plot.py --data-dir " << out_dir << "\n";

    return 0;
}
