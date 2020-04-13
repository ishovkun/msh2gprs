#pragma once
#include <string>
#include <iostream>    // provides cout
// Unix-only headers to get terminal width
#include <sys/ioctl.h> //ioctl() and TIOCGWINSZ
#include <unistd.h> // for STDOUT_FILENO

namespace utils {

class ProgressBar
{
 public:
  ProgressBar(std::string prefix, const size_t n_items, const size_t width = 80)
      : _prefix(prefix), _nitems(n_items), _max_width(width)
  {
    _percentage = 0;
    print_();
  }

  size_t set_progress(const size_t item_index)
  {
    const size_t new_percentage = get_percentage(item_index);
    if (new_percentage != _percentage)
    {
      _percentage = new_percentage;
      print_();
    }
    return new_percentage;
  }

  void finalize()
  {
    _percentage = 100;
    print_();
    std::cout << std::endl;
  }

  size_t get_percentage(const size_t item_index) const
  {
    return size_t(100.0 * item_index / _nitems);
  }

 private:
  void print_() const
  {
    size_t width = get_window_width_();
    if (_max_width != 0)
      width = std::min( _max_width, width );

    // write prefix
    const std::string prefix = "\r"  + _prefix + " [";
    std::cout << prefix;
    width -= prefix.size();

    // make and remember suffix
    const std::string suffix = "] " + std::to_string(_percentage) + " %\r";
    width -= suffix.size();

    // print bar
    const size_t nfill = size_t(double(_percentage) / 100.0 * width);
    for (size_t i=0; i<nfill; ++i)
      std::cout << _fill;
    for (size_t i=nfill; i<width; ++i)
      std::cout << ".";

    // print suffix
    std::cout << suffix;
    std::cout << std::flush;
  }

  size_t get_window_width_() const
  {
    struct winsize win_size;
    ioctl(STDOUT_FILENO, TIOCGWINSZ, &win_size);
    return win_size.ws_col;
  }

  // ======================== Variables ===========================
  const std::string _prefix;
  const size_t _nitems;
  size_t _max_width;
  size_t _percentage;
  const std::string _fill = "█";
};

}  // end namespace util
