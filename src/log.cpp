#pragma once

#include "log.hpp"

#define LOG_ERROR(msg)   Logger::instance().error(msg)
#define LOG_WARN(msg)    Logger::instance().warn(msg)
#define LOG_INFO(msg)    Logger::instance().info(msg)
#define LOG_DEBUG(msg)   Logger::instance().debug(msg)
#define LOG_TRACE(msg)   Logger::instance().log(LogLevel::Trace, msg)

// EXAMPLE :
// LOG_INFO("Program started");
// LOG_DEBUG("Loaded " + std::to_string(node_count) + " nodes");
// LOG_WARN("Using fallback parameter");
// LOG_ERROR("Cannot open file");
// LOG_TRACE("Detailed trace info...");

static Logger& Logger::instance() {
    static Logger _instance;
    return _instance;
}

void Logger::setLevel(LogLevel level) {
    logLevel = level;
}

void Logger::log(LogLevel level, const std::string& message) {
    if (level <= logLevel) {
        std::lock_guard<std::mutex> lock(mutex);
        std::ostream& out = (level == LogLevel::Error) ? std::cerr : std::cout;
        out << levelToString(level) << " : " << message << std::endl;
    }
}

void Logger::debug(const std::string& msg) { log(LogLevel::Debug, msg); }
void Logger::info(const std::string& msg) { log(LogLevel::Info, msg); }
void Logger::warn(const std::string& msg) { log(LogLevel::Warning, msg); }
void Logger::error(const std::string& msg) { log(LogLevel::Error, msg); }
