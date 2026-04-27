#pragma once

#include <cctype>
#include <cerrno>
#include <cstdlib>

namespace ins {
//辅助读取文件函数
inline bool isAsciiSpace(char c) {
	return std::isspace(static_cast<unsigned char>(c)) != 0;
}

inline bool parseDoubleToken(const char*& p, double& value) {
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

} // namespace ins
