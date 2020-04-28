#pragma once
#include <string>
#include <iostream>    // provides cout

namespace logging {

class ProgressBar
{
 public:
  ProgressBar(std::string prefix, const size_t n_items, const size_t width = 80);

  int set_progress(const size_t item_index);

  void finalize();

  int get_percentage(const size_t item_index) const;

 private:
  void print_() const;

  int get_window_width_() const;

  // ======================== Variables ===========================
  const std::string _prefix;
  const size_t _nitems;
  int _max_width;
  int _percentage;
  const std::string _fill = "█";
};

}  // end namespace util
