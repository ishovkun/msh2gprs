#include "Color.hpp"
#include <iostream>  // std::cout on linux
#ifdef WIN32
#include <windows.h>   // WinApi header
#endif

namespace colors
{

Color::Color(const std::string esc_sequence,
             const short win_color_attribute)
    :
    seq(esc_sequence),
    attr(win_color_attribute)
{}



std::ostream &operator<<(std::ostream  & os,
                         const Color   & color)
{
#ifdef WIN32
  HANDLE hConsole = GetStdHandle(STD_OUTPUT_HANDLE);
  SetConsoleTextAttribute(hConsole, color.attr);
#elif __linux__
  std::cout << "\033[" << color.seq << "m";
#endif
  return os;
}

short get_default_color()
{
  short default_color_attribute = 0;
#ifdef WIN32  // windows platform
  CONSOLE_SCREEN_BUFFER_INFO info;
  GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
  default_color_attribute = info.wAttributes;
#endif
  return default_color_attribute;
}

const Color Default = Color("0;0",  get_default_color());
const Color Blue    = Color("0;34", 1);
const Color Green   = Color("0;32", 2);
const Color Cyan    = Color("0;36", 3);
const Color Red     = Color("0;31", 4);
const Color Magenta = Color("0;35", 5);
const Color Yellow  = Color("0;33", 6);
const Color White   = Color("0;37", 7);

const Color BrightBlue    = Color("1;34", 9);
const Color BrightGreen   = Color("1;32", 10);
const Color BrightCyan    = Color("1;36", 11);
const Color BrightRed     = Color("1;31", 12);
const Color BrightMagenta = Color("1;35", 13);
const Color BrightYellow  = Color("1;33", 14);
const Color BrightWhite   = Color("1;37", 15);

}
