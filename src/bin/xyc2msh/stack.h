#ifndef STACK_H
#define STACK_H
extern void estiva_push(int x);
#define  push(n) estiva_push(n) 
extern int estiva_pop(int *p);
#define  pop(n) estiva_pop(&(n)) 
#endif
