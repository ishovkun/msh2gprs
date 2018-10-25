#pragma once
#include <exception>
#include <string>

namespace angem
{

struct NotImplemented : public std::exception
{
  NotImplemented(const std::string & msg)
      :
      message(msg)
  {
    // message = "this functionality is currently not implemented in angem: " + msg;
  }

	virtual const char * what () const throw ()
  {
    return message.c_str();
  }

 private:
  std::string message;
};

}
