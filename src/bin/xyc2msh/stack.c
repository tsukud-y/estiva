#include "estiva.h"
#include "stack.h"
#include "confary.h"

static int top=0, *stack;

void estiva_push(int x)
     /* push(n) estiva_push(n) */
{ static int sup=1;
  ary1(stack,sup);
  if(sup <= top){
    static int i, *tmp;
    sup *= 2;
    ary1(tmp,sup);   forall(0,i,dim1(stack)) tmp[i] = stack[i];
    ary1(stack,sup); forall(0,i,dim1(stack)) stack[i] = tmp[i];
    ary1(tmp,0);
  }
  stack[top++]=x;
}
int estiva_pop(int *p)
     /* pop(n) estiva_pop(&(n)) */
{ 
  if(top!=0){ *p=stack[--top]; return 1;}
  return 0;
}

