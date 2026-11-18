#include <setjmp.h>

#if defined(_MSC_VER)
#define OPTIMOTU_TLS __declspec(thread)
#else
#define OPTIMOTU_TLS __thread
#endif

static OPTIMOTU_TLS jmp_buf optimotu_wfa_jmp;
static OPTIMOTU_TLS volatile int optimotu_wfa_jmp_ready;

int optimotu_wfa_call(void (*fn)(void *), void *arg) {
  if (setjmp(optimotu_wfa_jmp) != 0) {
    optimotu_wfa_jmp_ready = 0;
    return 1;
  }
  optimotu_wfa_jmp_ready = 1;
  fn(arg);
  optimotu_wfa_jmp_ready = 0;
  return 0;
}

__attribute__((noreturn)) void optimotu_wfa_fail(int code) {
  (void)code;
  if (optimotu_wfa_jmp_ready) {
    optimotu_wfa_jmp_ready = 0;
    longjmp(optimotu_wfa_jmp, 1);
  }
  /*
   * exit() does not return. A call with no handler installed has no
   * safe libc function left (exit, _exit, and abort are all banned).
   */
  for (;;) {
  }
}
