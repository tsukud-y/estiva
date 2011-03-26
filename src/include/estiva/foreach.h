#ifndef _ESTIVA_FOREACH_H_
#define _ESTIVA_FOREACH_H_

int     estiva_foreach(int fsize,void *f0,...);
void   *estiva_foreachend();

#define foreach(f) while(estiva_foreach(sizeof(f),&f, 
#define foreachend estiva_foreachend())) 

#endif
