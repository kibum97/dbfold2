#include "logger.h"
#include <iostream>
#include <cstdarg>
#include <ctime>

FILE* Logger::logFilePtr = nullptr;

void Logger::initialize(const std::string &logFile) {
    logFilePtr = fopen(logFile.c_str(), "w");
    if (!logFilePtr) {
        std::cerr << "Failed to open log file: " << logFile << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Logger::log(LogLevel level, const std::string &message) {
    if (!logFilePtr) {
        std::cerr << "Logger not initialized!" << std::endl;
        return;
    }

    std::time_t now = std::time(nullptr);
    char timeStr[20];
    std::strftime(timeStr, sizeof(timeStr), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    fprintf(logFilePtr, "[%s] [%s] %s\n", timeStr, getLogLevelString(level), message.c_str());
    fflush(logFilePtr);
}

void Logger::logf(LogLevel level, const char *format, ...) {
    if (!logFilePtr) {
        std::cerr << "Logger not initialized!" << std::endl;
        return;
    }

    std::time_t now = std::time(nullptr);
    char timeStr[20];
    std::strftime(timeStr, sizeof(timeStr), "%Y-%m-%d %H:%M:%S", std::localtime(&now));

    fprintf(logFilePtr, "[%s] [%s] ", timeStr, getLogLevelString(level));

    va_list args;
    va_start(args, format);
    vfprintf(logFilePtr, format, args);
    va_end(args);

    fprintf(logFilePtr, "\n");
    fflush(logFilePtr);
}

void Logger::finalize() {
    if (logFilePtr) {
        fclose(logFilePtr);
        logFilePtr = nullptr;
    }
}

const char* Logger::getLogLevelString(LogLevel level) {
    switch (level) {
        case INFO: return "INFO";
        case WARNING: return "WARNING";
        case ERROR: return "ERROR";
        default: return "UNKNOWN";
    }
}