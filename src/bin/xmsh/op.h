#ifndef _OP_H_
#define _OP_H_
extern void  op_initop(int pargc, char **pargv);
extern int   op_defop(char *str);
extern int   op_cmpop(char *opt, char *str);
extern char *op_getop(char *str);
extern void *argf(int argc, char **argv);

#define initop(pargc,pargv) op_initop(pargc,pargv)
#define defop(str) op_defop(str) 
#define cmpop(opt,str) op_cmpop(opt,str) 
#define getop(str) op_getop(str) 
#endif
