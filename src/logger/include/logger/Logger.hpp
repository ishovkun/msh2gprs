#pragma once
#include "Color.hpp"  // colors::Color
#include <iostream>  // std::cout
#include <fstream>  // std::ofstream
#include <vector>   // std::vector
#include <memory>   // std::unique_ptr

namespace logging {

enum class LogLevel : int
{
  Critical  = 0,
  Silent    = 1,
  Warning   = 2,
  Important = 3,
  Message   = 4,
  Debug     = 5,
};

/**
 * This class implements a standard logger with colors.
 * It operates as a singleton.
 */
class Logger {
 public:
  /** Return reference to the logger */
  static Logger & ref();
  /** Add output stream */
  void add_stream(std::ostream & out);
  /** Set file to print secondary output into */
  void set_file(const std::string file_name);
  /** Set debugging level of the subsequent messages */
  void set_verbosity(const LogLevel verbosity);
  /* Set the level of the subsequent message */
  void set_message_level(const LogLevel level);
  /* Get the current verbosity level */
  LogLevel verbosity() const { return _verbosity; }
  /* write the log */
  template <typename T>
  // output operator
  Logger & operator<<(const T & data);
  // Operator overloading for std::endl
  Logger & operator<<(std::ostream& (*os)(std::ostream&));
  /* * Destructor */
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
  void select_color_(const LogLevel level);
  void reset_color_();


  // ------------------ Variables ---------------------- //
  static Logger * _p_logger;  // pointer to the logger instannce
  std::ofstream _fout;  // output file stream
  static constexpr int _file_stream_undefined = -1;
  int _file_stream_index;
  std::vector<std::ostream*> _streams;
  LogLevel _verbosity;           // current level of debugging
  LogLevel _msg_level;
};

template <typename T>
Logger & Logger::operator<<(const T & msg)
{
  if (_verbosity >= _msg_level)
  {
    // print_prefix_(_msg_level);
    select_color_(_msg_level);
    for (auto * p_stream : _streams)
      (*p_stream) << msg;

    reset_color_();
  }
  return *this;
}

// use for Message loging level
Logger& log();
// use for debug messages
Logger& debug();
// use for warnings
Logger& warning();
// use for critical errors
Logger& critical();
// use for important messages
Logger& important();


}  // end namespace logger
