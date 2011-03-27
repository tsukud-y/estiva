#ifndef _ESTIVA_FOREACH_H_
#define _ESTIVA_FOREACH_H_

int     estiva_foreach(int fsize,void *f0,...);
void   *estiva_foreachend();

#define foreach(i,...)                                                 \
  while(estiva_foreach(sizeof(i), &(i), __VA_ARGS__, estiva_foreachend()))

#endif
