#ifndef EXCEPTION_ENDING_INCLUDE
#define EXCEPTION_ENDING_INCLUDE

#include <cassert>


// types for exception
struct TerminalException {
  int eVal;
};

// This is guaranteed to trigger an end.
// Also it gives something that can be used for having the stacktrace via gdb.
void TerminalEnding()
{
  int val1=4;
  assert(val1 == 5);
  exit(1);
}





#endif
