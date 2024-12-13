#ifndef LOGGER_HPP
#define LOGGER_HPP
#include <iostream>
#include <unordered_map>
#include <memory>
namespace frfl
{
    namespace logger
    {
        class Logger{
        public: 
            std::string logger_name;
            bool enabled = true;
            Logger(const std::string& logger_name){
                this->logger_name = logger_name;
            }
            virtual ~Logger() = default;

            virtual void _log_impl(const std::string& message) = 0;
            
            template <typename... Args>
            void log(Args&&... args) {
                if (!enabled) return;
                std::string message = concatenate(std::forward<Args>(args)...);
                _log_impl(message);
            }

            template <typename... Args>
            void log_line(Args&&... args){
                if(!enabled) return;
                std::string message = concatenate(std::forward<Args>(args)...) + "\n";
                log(message);
            }
        };

        class StdLogger : public Logger {
        public:
            StdLogger(const std::string& logger_name) : Logger(logger_name) {}

            void _log_impl(const std::string& message) override {
                std::cout << "[" << logger_name << "] " << message;
            }
        };

        class Loggers{
            public:
            static std::unordered_map<std::string, std::shared_ptr<Logger>> loggers;

            template<typename LoggerType>
            static std::shared_ptr<LoggerType> build_logger(const std::string& logger_name){
                if(loggers.find(logger_name) != loggers.end()){
                    std::shared_ptr<LoggerType> logger = 
                        std::make_shared<LoggerType>(logger_name);
                    loggers[logger_name] = logger;
                    return logger;
                }
                return nullptr;
            }

            // Enable a specific logger by name
            static void enable_logger(const std::string& logger_name) {
                loggers[logger_name]->enabled = true;
            }

            // Disable a specific logger by name
            static void disable_logger(const std::string& logger_name) {
                loggers[logger_name]->enabled = false;
            }

            // Log function with conditional logging based on the logger's state
            template <typename... Args>
            inline static void log(const std::string& logger_name, Args&&... args) {
                loggers[logger_name]->log(std::forward(args)...);
            }

            // Log function with a newline
            template <typename... Args>
            inline static void log_line(const std::string& logger_name, Args&&... args) {
                loggers[logger_name]->log_line(std::forward(args)...);
            }
        };
    } // namespace logger
    
    
}; // namespace name

#endif