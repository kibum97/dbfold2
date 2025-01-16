#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <cstdio>

enum LogLevel {
    INFO,
    WARNING,
    ERROR
};

class Logger {
public:
    static void initialize(const std::string &logFile);
    static void log(LogLevel level, const std::string &message);
    static void logf(LogLevel level, const char *format, ...);
    static void finalize();

private:
    static FILE *logFilePtr;
    static const char* getLogLevelString(LogLevel level);
};

#endif // LOGGER_H