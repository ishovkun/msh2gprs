#include "Logger.hpp"

namespace logging {

Logger * Logger::_p_logger = NULL;

Logger::Logger()
    : _verbosity(LogLevel::Debug), _file_stream_index(-1)
{
  add_stream(std::cout);
}

void Logger::add_stream(std::ostream & out)
{
  _streams.emplace_back(&out);
}

Logger::~Logger()
{
  if (_fout.is_open())
    _fout.close();
}

Logger & Logger::ref()
{
  if (!_p_logger)  // create a logger
  {
    static Logger logger_instance;
    _p_logger = &logger_instance;
  }
  return *_p_logger;
}

void Logger::set_file(const std::string file_name)
{
  if (_fout.is_open())
    _fout.close();
  _fout.open(file_name.c_str());

  if (_file_stream_index == _file_stream_undefined)
  {
    _file_stream_index = _streams.size();
    add_stream(_fout);
  }
}

void Logger::set_verbosity(const LogLevel verbosity)
{
  _verbosity = verbosity;
}

void Logger::set_message_level(const LogLevel level)
{
  _msg_level = level;
}

void Logger::print_prefix_(const LogLevel level)
{
  if (!_fout.is_open())
    return;

  switch (level)
  {
    case LogLevel::Warning:
      std::cout << "Warning: ";
      _fout << "\nWarning   : ";
      break;
    case LogLevel::Critical:
      _fout << "\nCritical  : ";
      break;
     default:
       break;
  }
}

void Logger::select_color_(const LogLevel level)
{
  switch (level)
  {
    case LogLevel::Important:
      std::cout << colors::BrightBlue;
      break;
    case LogLevel::Warning:
      std::cout << colors::BrightRed;
      break;
    case LogLevel::Critical:
      std::cout << colors::BrightRed;
      break;
    default:
      break;
  }
}

void Logger::reset_color_()
{
  std::cout << colors::Default;
}

Logger& log()
{
  auto & logger = Logger::ref();
  logger.set_message_level(LogLevel::Message);
  return logger;
}

Logger& debug()
{
  auto & logger = Logger::ref();
  logger.set_message_level(LogLevel::Debug);
  return logger;
}

Logger& warning()
{
  auto & logger = Logger::ref();
  logger.set_message_level(LogLevel::Warning);
  return logger;
}

Logger& critical()
{
  auto & logger = Logger::ref();
  logger.set_message_level(LogLevel::Critical);
  return logger;
}

Logger& important()
{
  auto & logger = Logger::ref();
  logger.set_message_level(LogLevel::Important);
  return logger;
}

Logger & Logger::operator<<(std::ostream& (*os)(std::ostream&))
{
  if (_verbosity >= _msg_level)
  {
    // print_prefix_(_msg_level);
    std::cout << os;
    if (_fout.is_open())
      _fout << os;

    // reset_color_();
  }
  return *this;
}


}  // end namespace logger
