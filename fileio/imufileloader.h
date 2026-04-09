#pragma once

#include <array>
#include <cctype>
#include <cerrno>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <string>

#include "imu.h"

// 作用：定义 ImuFileLoader 类。
class ImuFileLoader {
public:
	enum class FileFormat {
		kUnknown = 0,
		kText,
		kBinaryDouble7
	};

	ImuFileLoader() = default;

	// 作用：ImuFileLoader 函数。
	explicit ImuFileLoader(const std::string& file_path) {
		open(file_path);
	}

	~ImuFileLoader() {
		close();
	}

	// 作用：open 函数。
	bool open(const std::string& file_path) {
		close();
		file_path_ = file_path;
		line_number_ = 0;
		bad_line_count_ = 0;
		last_error_.clear();
		file_format_ = FileFormat::kUnknown;

		stream_.open(file_path_, std::ios::in | std::ios::binary);
		if (!stream_.is_open()) {
			last_error_ = "failed to open file: " + file_path_;
			return false;
		}

		if (!detectFileFormat()) {
			return false;
		}
		return true;
	}

	// 作用：close 函数。
	void close() {
		if (stream_.is_open()) {
			stream_.close();
		}
		file_format_ = FileFormat::kUnknown;
	}

	// 作用：isOpen 函数。
	bool isOpen() const {
		return stream_.is_open();
	}

	// 作用：此处说明当前字段或步骤的用途。
	// 作用：readNext 函数。
	bool readNext(ImuData& out_data) {
		if (!stream_.is_open()) {
			last_error_ = "file is not open";
			return false;
		}

		if (file_format_ == FileFormat::kBinaryDouble7) {
			// 作用：readNextBinary 函数。
			return readNextBinary(out_data);
		}
		if (file_format_ == FileFormat::kText) {
			// 作用：readNextText 函数。
			return readNextText(out_data);
		}

		last_error_ = "unknown file format";
		return false;
	}

	// 作用：fileFormat 函数。
	FileFormat fileFormat() const {
		return file_format_;
	}

	// 作用：isBinary 函数。
	bool isBinary() const {
		return file_format_ == FileFormat::kBinaryDouble7;
	}

	// 作用：isText 函数。
	bool isText() const {
		return file_format_ == FileFormat::kText;
	}

	// 作用：此处说明当前字段或步骤的用途。
	// 作用：readNextText 函数。
	bool readNextText(ImuData& out_data) {
		if (!stream_.is_open()) {
			last_error_ = "file is not open";
			return false;
		}

		std::string line;
		while (std::getline(stream_, line)) {
			++line_number_;

			ImuData data;
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

	// 作用：readNextBinary 函数。
	bool readNextBinary(ImuData& out_data) {
		if (!stream_.is_open()) {
			last_error_ = "file is not open";
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
				++bad_line_count_;
				last_error_ = "incomplete binary imu record at index " + std::to_string(line_number_ + 1);
				return false;
			}

			++line_number_;

			ImuData data;
			data.time = values[0];
			data.dtheta_x = values[1];
			data.dtheta_y = values[2];
			data.dtheta_z = values[3];
			data.dvel_x = values[4];
			data.dvel_y = values[5];
			data.dvel_z = values[6];

			if (!data.isFinite()) {
				++bad_line_count_;
				last_error_ = "non-finite value at binary record " + std::to_string(line_number_);
				continue;
			}

			out_data = data;
			return true;
		}
	}

	// 作用：lineNumber 函数。
	std::uint64_t lineNumber() const {
		return line_number_;
	}

	// 作用：badLineCount 函数。
	std::uint64_t badLineCount() const {
		return bad_line_count_;
	}

	// 作用：lastError 函数。
	const std::string& lastError() const {
		return last_error_;
	}

private:
	// 作用：此处说明当前字段或步骤的用途。
	// 作用：detectFileFormat 函数。
	bool detectFileFormat() {
		constexpr std::streamsize kProbeBytes = 4096;
		char probe[kProbeBytes]{};

		stream_.clear();
		stream_.seekg(0, std::ios::beg);
		stream_.read(probe, kProbeBytes);
		const std::streamsize n = stream_.gcount();

		if (n <= 0) {
			last_error_ = "empty imu file: " + file_path_;
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
				last_error_ = "invalid binary imu file size";
				return false;
			}

			if ((size % static_cast<std::streamoff>(sizeof(double) * 7)) != 0) {
				last_error_ = "unsupported binary imu format (expect 7 doubles per record)";
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

	// 作用：此处说明当前字段或步骤的用途。
	// 作用：isAsciiSpace 函数。
	static bool isAsciiSpace(char c) {
		return std::isspace(static_cast<unsigned char>(c)) != 0;
	}

	// 作用：此处说明当前字段或步骤的用途。
	// 作用：parseDoubleToken 函数。
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

	// 作用：此处说明当前字段或步骤的用途。
	// 作用：parseLine 函数。
	static bool parseLine(const std::string& line, ImuData& out_data) {
		const char* p = line.c_str();

		if (!parseDoubleToken(p, out_data.time)) {
			return false;
		}
		if (!parseDoubleToken(p, out_data.dtheta_x)) {
			return false;
		}
		if (!parseDoubleToken(p, out_data.dtheta_y)) {
			return false;
		}
		if (!parseDoubleToken(p, out_data.dtheta_z)) {
			return false;
		}
		if (!parseDoubleToken(p, out_data.dvel_x)) {
			return false;
		}
		if (!parseDoubleToken(p, out_data.dvel_y)) {
			return false;
		}
		if (!parseDoubleToken(p, out_data.dvel_z)) {
			return false;
		}

		while (*p != '\0' && isAsciiSpace(*p)) {
			++p;
		}

		// 作用：此处说明当前字段或步骤的用途。
		return *p == '\0';
	}

private:
	std::ifstream stream_;
	std::string file_path_;
	FileFormat file_format_ = FileFormat::kUnknown;
	std::uint64_t line_number_ = 0;
	std::uint64_t bad_line_count_ = 0;
	std::string last_error_;
};
