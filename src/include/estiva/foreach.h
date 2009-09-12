#ifndef _ESTIVA_FOREACH_H_
#define _ESTIVA_FOREACH_H_

#define  foreach(f) while(estiva_foreach(sizeof(f),&f, 
#define  end estiva_end())) 

extern int estiva_foreach(int fsize,void *f0,...);
extern void *estiva_end();

#endif
