#pragma once

#include <fstream>
#include <string>

#include "gnss.h"
#include "parser_utils.h"

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
		stream_.open(file_path);
		return stream_.is_open();
	}

	void close() {
		if (stream_.is_open()) {
			stream_.close();
		}
	}

	bool readNext(GnssData& out_data) {
		if (!stream_.is_open()) {
			return false;
		}

		std::string line;
		while (std::getline(stream_, line)) {
			GnssData data;
			if (!parseLine(line, data)) {
				continue;
			}

			if (!data.isFinite()) {
				continue;
			}

			out_data = data;
			return true;
		}

		return false;
	}

private:
	static bool parseLine(const std::string& line, GnssData& out_data) {
		const char* p = line.c_str();
		double first_angle = 0.0;
		double second_angle = 0.0;

		if (!ins::parseDoubleToken(p, out_data.time)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, first_angle)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, second_angle)) {
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

		if (!ins::parseDoubleToken(p, out_data.altitude)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, out_data.std_x)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, out_data.std_y)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, out_data.std_z)) {
			return false;
		}

		while (*p != '\0' && ins::isAsciiSpace(*p)) {
			++p;
		}

		return true; // We don't care about trailing data like velocity
	}

private:
	std::ifstream stream_;
};
