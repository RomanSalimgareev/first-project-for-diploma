#include "error_handling.h"

void assertError(bool condition, const std::string& msg, const char* file)
{
	if (!(condition))									
	{												
		std::cerr << "Error : " << msg << " Location error: " << file << "\n";	
		std::exit(EXIT_FAILURE);					
	}
}

void log(LogLevel level, const std::string& message)
{
	std::ofstream logFile("error_log.txt", std::ios::app);
	if (logFile.is_open())
	{
		std::time_t now = std::time(nullptr);
		char buffer[26];
		std::string levelStr = "Nothing";
		switch (level)
		{
		case LogLevel::WARNING:
			levelStr = "WARNING";
			break;
		case LogLevel::ERROR:
			levelStr = "ERROR";
			break;
		}
		ctime_s(buffer, sizeof(buffer), &now);
		logFile << "[" << buffer << "] : " << levelStr << "" <<
			message << "\n\n";
	}
}

void error(std::string msg, const char* file, const long& line)
{
	msg += "Location error: ";
	msg += std::string(file) + "\t";

	msg += "Line error: ";
	msg += std::to_string(line) + "\n\n";

	log(LogLevel::ERROR, msg);
	throw std::runtime_error(msg);
}

void warning(std::string msg, const char* file, const long& line)
{
	msg += "Location error: ";
	msg += std::string(file) + "\t";

	msg += "Line error: ";
	msg += std::to_string(line) + "\n\n";
	log(LogLevel::WARNING, msg);
}

std::string messageOutOfRange()
{
	return "Out of range. ";
}

std::string messageDivideZero()
{
	return "Dividing by zero. ";
}
