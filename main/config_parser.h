#pragma once

#include <array>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

namespace ins {

// 简单的 YAML 风格配置文件解析器，不依赖任何外部库
class ConfigParser {
public:
    void load(const std::string& file_path) {
        std::ifstream ifs(file_path);
        if (!ifs.is_open()) {
            throw std::runtime_error("Cannot open config file: " + file_path);
        }
        std::string line;
        std::size_t lineno = 0;
        while (std::getline(ifs, line)) {
            ++lineno;
            // 去除首尾空白
            const std::string trimmed = trim(line);
            if (trimmed.empty() || trimmed[0] == '#') continue;

            const std::size_t colon = trimmed.find(':');
            if (colon == std::string::npos) {
                throw std::runtime_error("Config syntax error at line " +
                    std::to_string(lineno) + ": " + trimmed);
            }
            std::string key = trim(trimmed.substr(0, colon));
            std::string val = trim(trimmed.substr(colon + 1));
            map_[key] = val;
        }
    }

    bool has(const std::string& key) const {
        return map_.find(key) != map_.end();
    }

    std::string str(const std::string& key, const std::string& fallback = "") const {
        auto it = map_.find(key);
        return (it != map_.end()) ? it->second : fallback;
    }

    double num(const std::string& key, double fallback = 0.0) const {
        auto it = map_.find(key);
        if (it == map_.end()) return fallback;
        try { return std::stod(it->second); }
        catch (...) { return fallback; }
    }

    bool flag(const std::string& key, bool fallback = false) const {
        auto it = map_.find(key);
        if (it == map_.end()) return fallback;
        const std::string v = toLower(it->second);
        return v == "true" || v == "yes" || v == "1" || v == "on";
    }

    int integer(const std::string& key, int fallback = 0) const {
        auto it = map_.find(key);
        if (it == map_.end()) return fallback;
        try { return std::stoi(it->second); }
        catch (...) { return fallback; }
    }

    std::array<double, 3> arr3(const std::string& key,
                               const std::array<double, 3>& fallback = {0,0,0}) const {
        auto it = map_.find(key);
        if (it == map_.end()) return fallback;
        return parseArray3(it->second, fallback);
    }

private:
    std::unordered_map<std::string, std::string> map_;

    static std::string trim(const std::string& s) {
        std::size_t b = 0, e = s.size();
        while (b < e && (s[b] == ' ' || s[b] == '\t' || s[b] == '\r')) ++b;
        while (e > b && (s[e-1] == ' ' || s[e-1] == '\t' || s[e-1] == '\r')) --e;
        return s.substr(b, e - b);
    }

    static std::string toLower(std::string s) {
        for (auto& c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        return s;
    }

    static std::array<double, 3> parseArray3(const std::string& raw,
                                              const std::array<double, 3>& fallback) {
        // 去除 [ ]
        std::string s = raw;
        if (!s.empty() && s.front() == '[') s.erase(0, 1);
        if (!s.empty() && s.back()  == ']') s.pop_back();

        std::vector<double> vals;
        std::stringstream ss(s);
        std::string token;
        while (std::getline(ss, token, ',')) {
            token = trim(token);
            if (token.empty()) continue;
            try { vals.push_back(std::stod(token)); }
            catch (...) { return fallback; }
        }
        if (vals.size() < 3) return fallback;
        return {vals[0], vals[1], vals[2]};
    }
};

} // namespace ins
