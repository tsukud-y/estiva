#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/ary.h"
#include "estiva/que.h"


void estiva_initq(que **q)
{
  if ( *q == NULL ) ary1(*q,1); //*q = calloc(1,sizeof(que));
}


void estiva_enq(que *q, void *e, size_t n)
{
  void **addr;
  
  while (1) {
    if (q->elem == NULL) {
      addr = q->elem = malloc(sizeof(void *));
      *addr = malloc(n);
      memcpy(*addr,e,n);
      q->next = calloc(1,sizeof(que));
      return ;
    }
    q = q->next;
  }
}


void estiva_push(que **q, void *e, size_t n)
{
  void **addr;
  que *p;

  p = *q;
  *q = calloc(1,sizeof(que));
  addr = (*q)->elem = malloc(sizeof(void *));
  *addr = malloc(n);
  memcpy(*addr,e,n);
  (*q)->next = p;
}


void *estiva_pop(que **q, void *e, size_t n)
{
  que *p;
  p = *q;
  if (p->elem==NULL) return NULL;
  memcpy(e,*(void**)p->elem,n);
  *q = p->next;
  free(*(void**)p->elem);
  free(p->elem);
  free(p);
  return e;
}


void *estiva_lup(que *q, void **s2, size_t n)
{
  void **addr, **s1;

  while( q->elem != NULL ) {
    addr = q->elem;
    s1 = *addr;
    if(!strcmp(*s1,*s2)) return s1;
    q = q->next;
  }
  return NULL;
}
