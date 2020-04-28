#pragma once

#include <string>
#include <ostream>

namespace colors
{

class Color
{
 public:
  Color(const std::string esc_sequence,
        const short win_color_attribute);
  friend std::ostream &operator<<(std::ostream &,
                                  const Color  &);
 protected:
  std::string seq;
  short attr;
};

// std::ostream &operator<<(std::ostream &, const Color  &);

extern const Color Default;
extern const Color Blue;
extern const Color Green;
extern const Color Cyan;
extern const Color Red;
extern const Color Magenta;
extern const Color Yellow;
extern const Color White;

extern const Color BrightBlue;
extern const Color BrightGreen;
extern const Color BrightCyan;
extern const Color BrightRed;
extern const Color BrightMagenta;
extern const Color BrightYellow;
extern const Color BrightWhite;

}  // end namespace colors
