/* C Example */
#include <math.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
/*
#define DEBUG
*/

int main(int argc, char *argv[])
{
  char charBuf[1024];
  int system_failure;

  strcpy( charBuf, "chmod 755 */*.sh" );
#ifdef DEBUG
  printf( "execute %s\n", charBuf );
  fflush(stdout);
#endif
  system_failure = system(charBuf);
  if ( system_failure == -1 ) {
    printf( " Job failed. Command was:/n%s/n", charBuf );
    exit(EXIT_FAILURE);
  }

}

