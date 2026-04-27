#pragma once

#include <array>
#include <fstream>
#include <string>

#include "imu.h"
#include "parser_utils.h"

class ImuFileLoader {
public:
	enum class FileFormat {
		kUnknown = 0,
		kText,
		kBinaryDouble7
	};

	ImuFileLoader() = default;

	explicit ImuFileLoader(const std::string& file_path) {
		open(file_path);
	}

	~ImuFileLoader() {
		close();
	}

	bool open(const std::string& file_path) {
		close();
		file_format_ = FileFormat::kUnknown;

		stream_.open(file_path, std::ios::in | std::ios::binary);
		if (!stream_.is_open()) {
			return false;
		}

		if (!detectFileFormat()) {
			return false;
		}
		return true;
	}

	void close() {
		if (stream_.is_open()) {
			stream_.close();
		}
		file_format_ = FileFormat::kUnknown;
	}

	bool readNext(ImuData& out_data) {
		if (!stream_.is_open()) {
			return false;
		}

		if (file_format_ == FileFormat::kBinaryDouble7) {
			return readNextBinary(out_data);
		}
		if (file_format_ == FileFormat::kText) {
			return readNextText(out_data);
		}

		return false;
	}

	bool isBinary() const {
		return file_format_ == FileFormat::kBinaryDouble7;
	}

private:
	bool readNextText(ImuData& out_data) {
		if (!stream_.is_open()) {
			return false;
		}

		std::string line;
		while (std::getline(stream_, line)) {
			ImuData data;
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

	bool readNextBinary(ImuData& out_data) {
		if (!stream_.is_open()) {
			return false;
		}

		while (true) {
			std::array<double, 7> values{};
			stream_.read(reinterpret_cast<char*>(values.data()),
						 static_cast<std::streamsize>(sizeof(values)));

			if (stream_.gcount() == 0) {
				return false;
			}

			if (stream_.gcount() != static_cast<std::streamsize>(sizeof(values))) {
				return false;
			}

			ImuData data;
			data.time = values[0];
			data.dtheta_x = values[1];
			data.dtheta_y = values[2];
			data.dtheta_z = values[3];
			data.dvel_x = values[4];
			data.dvel_y = values[5];
			data.dvel_z = values[6];

			if (!data.isFinite()) {
				continue;
			}

			out_data = data;
			return true;
		}
	}

	bool detectFileFormat() {
		constexpr std::streamsize kProbeBytes = 4096;
		char probe[kProbeBytes]{};

		stream_.clear();
		stream_.seekg(0, std::ios::beg);
		stream_.read(probe, kProbeBytes);
		const std::streamsize n = stream_.gcount();

		if (n <= 0) {
			return false;
		}

		bool has_binary_marker = false;
		for (std::streamsize i = 0; i < n; ++i) {
			const unsigned char c = static_cast<unsigned char>(probe[i]);
			if (c == 0) {
				has_binary_marker = true;
				break;
			}
			if (c < 32 && c != '\n' && c != '\r' && c != '\t' && c != ' ') {
				has_binary_marker = true;
				break;
			}
		}

		if (has_binary_marker) {
			stream_.clear();
			stream_.seekg(0, std::ios::end);
			const std::streamoff size = stream_.tellg();
			if (size <= 0) {
				return false;
			}

			if ((size % static_cast<std::streamoff>(sizeof(double) * 7)) != 0) {
				return false;
			}
			file_format_ = FileFormat::kBinaryDouble7;
		} else {
			file_format_ = FileFormat::kText;
		}

		stream_.clear();
		stream_.seekg(0, std::ios::beg);
		return true;
	}

	static bool parseLine(const std::string& line, ImuData& out_data) {
		const char* p = line.c_str();

		if (!ins::parseDoubleToken(p, out_data.time)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, out_data.dtheta_x)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, out_data.dtheta_y)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, out_data.dtheta_z)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, out_data.dvel_x)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, out_data.dvel_y)) {
			return false;
		}
		if (!ins::parseDoubleToken(p, out_data.dvel_z)) {
			return false;
		}

		while (*p != '\0' && ins::isAsciiSpace(*p)) {
			++p;
		}

		return *p == '\0';
	}

private:
	std::ifstream stream_;
	FileFormat file_format_ = FileFormat::kUnknown;
};