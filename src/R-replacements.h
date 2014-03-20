#include <cstdio>
#include <cstdarg>

namespace R {

// NB: macros are NOT local to a namespace
// #define REprintf printf

void REprintf(const char *format, ...) {
  va_list args;
  va_start(args, format);
  vprintf(format, args);
  va_end(args);
}

} // namespace R
