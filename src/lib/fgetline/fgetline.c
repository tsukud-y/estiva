#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "estiva/op.h"
#include "estiva/ary.h"

typedef struct {
  int c;
  void *next;
  void *tail;
} list;


static void addlist(list *elem, int c)
{
  list *top = elem;

  if ( top->tail == NULL )
    ;
  else
    elem = top->tail;

  while (1) {
    if ( elem->next == NULL ) {
      elem->c = c;
      elem->next = calloc(1,sizeof(list));
      top->tail = elem;
      return ;
    }
    elem = elem->next;
  }
}


static unsigned long cntlist(list *top)
{
  unsigned long n=0;
  while (1) {
    n++;
    top = top->next;
    if ( top == NULL ) return n;
  }
}


static char *list2str(list *top, char *str)
{
  unsigned long i=0;
  for ( ;top->next != NULL;i++) {
    str[i] = top->c;
    top = top->next;
  }
  str[i] = 0;
  return str;
}


static void rmlist(list *p)
{
  list *tmp;
  for ( p = p->next; p; p = tmp ) {
    tmp = p->next;
    free(p);
  }
}


char *estiva_fgetline(FILE *fp)
{
  static list *top;
  static char *buf, *append;
  unsigned long i, j, n=1000;
  int c; 

  ary1(buf,n+1);

  ary1(top,1);
  rmlist(top);
  top->next = NULL;
  top->tail = NULL;

  for ( i=0, j=0; EOF != (c = fgetc(fp)) && c != '\n'; ){
    if ( i < n )
      buf[i++] = c;
    else 
      addlist(top,c);
  }
  buf[i] = 0;

  ary1(append, cntlist(top)+1);
  return strcat(buf,list2str(top,append));
}
