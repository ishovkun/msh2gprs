#pragma once
#include "Color.hpp"  // colors::Color
#include <iostream>  // std::cout
#include <fstream>  // std::ofstream
#include <memory>   // std::unique_ptr

namespace logging {

enum class LogLevel : int
{
  Critical = 0,
  Silent   = 1,
  Warning  = 2,
  Brief    = 3,
  Message  = 4,
  Debug    = 5,
};

/**
 * This class implements a standard logger with colors.
 * It operates as a singleton.
 */
class Logger {
 public:
  /**
   * Return reference to the logger
   */
  static Logger & ref();
  /**
   * Set file to print secondary output into
   */
  void set_file(const std::string file_name);
  /**
   * Set debugging level of the subsequent messages
   */
  void set_verbosity(const LogLevel verbosity);
  /**
   * Set the level of the subsequent message
   */
  void set_message_level(const LogLevel level);
  /**
   * do the log
   */
  template <typename T>
   std::ostream & operator<<(T & data);
  /**
   * do the log
   */
   void print(const std::string & message, const LogLevel msg_level);
  /**
   * Destructor
   */
  ~Logger();

  Logger(Logger const&) = delete;
  void operator=(Logger const&) = delete;

 private:
  /**
   * Constructor
   */
  Logger();
  // print debug, brief, etc. + timestamp
  void print_prefix_(const LogLevel level);


  // ------------------ Variables ---------------------- //
  static Logger * _p_logger;  // pointer to the logger instannce
  std::ofstream _fout;  // output file stream
  LogLevel _verbosity;           // current level of debugging
  LogLevel _msg_level;
};

template <typename T>
std::ostream & Logger::operator<<(T & data)
{
  if (_verbosity >= _msg_level)
  {
    print_prefix_(_msg_level);
    std::cout << data << std::endl;
    if (_fout.is_open())
    {
      _fout << data << std::endl;
    }

    std::cout << colors::Default;
  }
  return std::cout;
}


}  // end namespace logger

// #define log logging::Logger::ref().set_message_level(logging::Debug); logging::Logger::ref()
// #define brief logging::Logger::ref().set_message_level(logging::Brief); logging::Logger::ref()
// #define warning logging::Logger::ref().set_message_level(logging::Brief); logging::Logger::ref()
// #define error logging::Logger::ref().set_message_level(logging::Critical); logging::Logger::ref()
