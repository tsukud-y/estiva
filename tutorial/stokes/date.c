#include <stdio.h>
#include <time.h>

void date(char *str)
{
  time_t t;
  char *p;
  t = time( 0 );
  p = ctime(&t);
  p[19] = '\0';
  p += 11;
  fprintf(stderr,"%s %s\n",p,str);
  fflush(stderr);
}
