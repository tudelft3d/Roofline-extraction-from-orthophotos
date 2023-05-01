#pragma once
#include <string>
#include <iostream>

inline void trace(const int verbose_level, const int level, const std::string message) 
{
	if (level <= verbose_level) {
		std::cout << message << std::endl;
	}
}