#pragma once

#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <string>

#include "gnss.h"

class GnssFileLoader {
public:
	GnssFileLoader() = default;

		explicit GnssFileLoader(const std::string& file_path) {
		open(file_path);
	}

	~GnssFileLoader() {
		close();
	}

		bool open(const std::string& file_path) {
		close();
		file_path_ = file_path;
		line_number_ = 0;
		bad_line_count_ = 0;
		last_error_.clear();

		stream_.open(file_path_);
		if (!stream_.is_open()) {
			last_error_ = "failed to open file: " + file_path_;
			return false;
		}
		return true;
	}

		void close() {
		if (stream_.is_open()) {
			stream_.close();
		}
	}

		bool isOpen() const {
		return stream_.is_open();
	}

		bool readNext(GnssData& out_data) {
		if (!stream_.is_open()) {
			last_error_ = "file is not open";
			return false;
		}

		std::string line;
		while (std::getline(stream_, line)) {
			++line_number_;

			GnssData data;
			if (!parseLine(line, data)) {
				++bad_line_count_;
				last_error_ = "parse failed at line " + std::to_string(line_number_);
				continue;
			}

			if (!data.isFinite()) {
				++bad_line_count_;
				last_error_ = "non-finite value at line " + std::to_string(line_number_);
				continue;
			}

			out_data = data;
			return true;
		}

		if (!stream_.eof()) {
			last_error_ = "read error near line " + std::to_string(line_number_ + 1);
		}

		return false;
	}

		std::uint64_t lineNumber() const {
		return line_number_;
	}

		std::uint64_t badLineCount() const {
		return bad_line_count_;
	}

		const std::string& lastError() const {
		return last_error_;
	}

private:
			static bool isAsciiSpace(char c) {
		return std::isspace(static_cast<unsigned char>(c)) != 0;
	}

			static bool parseDoubleToken(const char*& p, double& value) {
		while (*p != '\0' && isAsciiSpace(*p)) {
			++p;
		}

		if (*p == '\0') {
			return false;
		}

		errno = 0;
		char* end = nullptr;
		value = std::strtod(p, &end);

		if (end == p || errno == ERANGE) {
			return false;
		}

		p = end;
		return true;
	}

					static bool parseLine(const std::string& line, GnssData& out_data) {
		const char* p = line.c_str();
		double first_angle = 0.0;
		double second_angle = 0.0;

		if (!parseDoubleToken(p, out_data.time)) {
			return false;
		}
		if (!parseDoubleToken(p, first_angle)) {
			return false;
		}
		if (!parseDoubleToken(p, second_angle)) {
			return false;
		}

				// 1) time, lon, lat, h, ...
		// 2) time, lat, lon, h, ...
				const double a1 = std::fabs(first_angle);
		const double a2 = std::fabs(second_angle);
		if (a1 <= 180.0 && a2 <= 90.0) {
			out_data.longitude = first_angle;
			out_data.latitude = second_angle;
		} else if (a1 <= 90.0 && a2 <= 180.0) {
			out_data.latitude = first_angle;
			out_data.longitude = second_angle;
		} else {
			out_data.longitude = first_angle;
			out_data.latitude = second_angle;
		}

		if (!parseDoubleToken(p, out_data.altitude)) {
			return false;
		}
		if (!parseDoubleToken(p, out_data.std_x)) {
			return false;
		}
		if (!parseDoubleToken(p, out_data.std_y)) {
			return false;
		}
		if (!parseDoubleToken(p, out_data.std_z)) {
			return false;
		}

		out_data.has_velocity = false;
		out_data.vel_n = 0.0;
		out_data.vel_e = 0.0;
		out_data.vel_d = 0.0;
		out_data.std_vn = 0.0;
		out_data.std_ve = 0.0;
		out_data.std_vd = 0.0;

				const char* p_backup = p;
		double vn = 0.0;
		double ve = 0.0;
		double vd = 0.0;
		double sv_n = 0.0;
		double sv_e = 0.0;
		double sv_d = 0.0;
		const bool has_vn = parseDoubleToken(p, vn);
		const bool has_ve = has_vn && parseDoubleToken(p, ve);
		const bool has_vd = has_ve && parseDoubleToken(p, vd);
		const bool has_svn = has_vd && parseDoubleToken(p, sv_n);
		const bool has_sve = has_svn && parseDoubleToken(p, sv_e);
		const bool has_svd = has_sve && parseDoubleToken(p, sv_d);

		if (has_svd) {
			out_data.has_velocity = true;
			out_data.vel_n = vn;
			out_data.vel_e = ve;
			out_data.vel_d = vd;
			out_data.std_vn = sv_n;
			out_data.std_ve = sv_e;
			out_data.std_vd = sv_d;
		} else {
						p = p_backup;
		}

		while (*p != '\0' && isAsciiSpace(*p)) {
			++p;
		}

		return *p == '\0';
	}

private:
	std::ifstream stream_;
	std::string file_path_;
	std::uint64_t line_number_ = 0;
	std::uint64_t bad_line_count_ = 0;
	std::string last_error_;
};
