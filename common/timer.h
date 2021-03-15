#pragma once
#include <chrono>
#ifdef USE_SPDLOG
#include <spdlog/spdlog.h>
#endif 
//#include <spdlog/fmt/ostr.h> // must be included
//#include <spdlog/sinks/stdout_sinks.h>

template <typename TimeT = std::chrono::milliseconds>
class Timer {
public:
    Timer() {
        start = std::chrono::system_clock::now();
        //logger = spdlog::stdout_logger_mt("console");
    }

    size_t value() const {
        auto now = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<TimeT>(now - start);
        return (size_t) duration.count();
    }

    size_t reset() {
        auto now = std::chrono::system_clock::now();
        auto duration = std::chrono::duration_cast<TimeT>(now - start);
        start = now;
        return (size_t) duration.count();
    }

    void beginStage(const std::string &name) {
        reset();

#ifdef USE_SPDLOG
        logger->info("{} ..", name);
		    logger->flush();
#endif
//        std::cout.flush();
    }

    void endStage(const std::string &str = "") {
#ifdef USE_SPDLOG
    	logger->info("done. (took {}ms, {})", value(), str);
#endif
    }
private:
    std::chrono::system_clock::time_point start;
#ifdef USE_SPDLOG
	  std::shared_ptr<spdlog::logger> logger;
#endif
};
