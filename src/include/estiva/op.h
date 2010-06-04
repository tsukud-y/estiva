#ifndef _ESTIVA_OP_H_
#define _ESTIVA_OP_H_

#define initop(pargc,pargv) estiva_initop(pargc,pargv)
#define defop(str)          estiva_defop(str) 
#define getop(str)          estiva_getop(str) 
#define stdfp()             estiva_stdfp()
#define ofp()               estiva_ofp()

extern void  estiva_initop(int pargc, char **pargv);
extern int   estiva_defop(char *str);
extern char *estiva_getop(char *str);
extern FILE *estiva_stdfp(void);
extern FILE *estiva_ofp(void);

#endif
