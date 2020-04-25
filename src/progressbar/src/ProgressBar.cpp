#include "ProgressBar.hpp"

// Unix-only headers to get terminal width
#include <sys/ioctl.h> //ioctl() and TIOCGWINSZ
#include <unistd.h> // for STDOUT_FILENO

namespace utils {

ProgressBar::ProgressBar(std::string prefix, const size_t n_items, const size_t width)
      : _prefix(prefix), _nitems(n_items), _max_width(width)
{
  _percentage = 0;
  print_();
}

int ProgressBar::set_progress(const size_t item_index)
{
  const int new_percentage = get_percentage(item_index);
  if (new_percentage != _percentage)
  {
    _percentage = new_percentage;
    print_();
  }
  return new_percentage;
}

void ProgressBar::finalize()
{
  _percentage = 100;
  print_();
  std::cout << std::endl;
}

int ProgressBar::get_percentage(const size_t item_index) const
{
  return int(100.0 * item_index / _nitems);
}

void ProgressBar::print_() const
{
  // get window width
  int width = get_window_width_();
  if (_max_width != 0)
    width = std::min( _max_width, width );

  // prefix and suffix
  const std::string prefix = "\r"  + _prefix + " [";
  const std::string suffix = "] " + std::to_string(_percentage) + " %\r";

  // don't print anything if the window is too small
  if (width <= prefix.size() + suffix.size())
    return;

  // write prefix
  std::cout << prefix << std::flush;
  width -= (prefix.size() + suffix.size());

  // print bar
  const int nfill = int(double(_percentage) / 100.0 * width);

  for (size_t i=0; i<nfill; ++i)
    std::cout << _fill;
  std::cout << std::flush;

  for (size_t i=nfill; i<width; ++i)
    std::cout << ".";
  std::cout << std::flush;

  // print suffix
  std::cout << suffix;
  std::cout << std::flush;
}

int ProgressBar::get_window_width_() const
{
  struct winsize win_size;
  ioctl(STDOUT_FILENO, TIOCGWINSZ, &win_size);
  return win_size.ws_col;
}

}  // end namespace utils
