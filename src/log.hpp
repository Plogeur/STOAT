#pragma once

#include <iostream>
#include <fstream>
#include <mutex>

enum class LogLevel {
    Error = 0,
    Warning = 1,
    Info = 2,
    Debug = 3,
    Trace = 4
};

class Logger {
public:
    static Logger& instance();

    void setLevel(LogLevel level);
    void log(LogLevel level, const std::string& message);

    void debug(const std::string& msg);
    void info(const std::string& msg);
    void warn(const std::string& msg);
    void error(const std::string& msg);

private:
    LogLevel logLevel = LogLevel::Info;
    std::mutex mutex;

    Logger() = default;

    std::string levelToString(LogLevel level) const {
        switch (level) {
            case LogLevel::Error: return "ERROR";
            case LogLevel::Warning: return "WARN";
            case LogLevel::Info: return "INFO";
            case LogLevel::Debug: return "DEBUG";
            case LogLevel::Trace: return "TRACE";
            default: return "UNKNOWN";
        }
    }
};
    