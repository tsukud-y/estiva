#ifndef FOREACH_H
#define FOREACH_H
extern int estiva_foreach(int fsize,void *f0,...);
#define  foreach(f) while(estiva_foreach(sizeof(f),&f, 
extern void *estiva_end();
#define  end estiva_end())) 
#endif
