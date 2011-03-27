#ifndef _ESTIVA_QUE_H_
#define _ESTIVA_QUE_H_

typedef struct {
  void *elem;
  void *next;
} que;

void  estiva_enq(  que  *q, void   *e, size_t n);
void  estiva_initq(que **q);
void *estiva_lup(  que  *q, void **s2, size_t n);
void  estiva_push( que **q, void   *e, size_t n);
void *estiva_pop(  que **q, void   *e, size_t n);
void  estiva_forq(que *q, void **e);
int   estiva_forq_loop(void **e);

#define enq(q,e)      estiva_enq(   q, (void*)&e, sizeof(e))
#define initq(q)      estiva_initq(&q)
#define lup(q,e)      estiva_lup(   q, (void*)&e, sizeof(e))
#define pop(q,e)      estiva_pop(  &q, (void*)&e, sizeof(e))
#define push(q,e)     estiva_push( &q, (void*)&e, sizeof(e))
#define forq(q,e)							\
  for( estiva_forq(q,(void*)&e); estiva_forq_loop((void*)&e);)

#endif
