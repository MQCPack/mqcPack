#include <stdlib.h>
#include <stdio.h>

void mqc_abort_(void)
{
  fflush(stdout);
  fflush(stderr);
  abort();
}
