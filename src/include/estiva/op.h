#ifndef _ESTIVA_OP_H_
#define _ESTIVA_OP_H_

int     estiva_defop(char *str);
char   *estiva_getop(char *str);
void    estiva_initop(int pargc, char **pargv);
void   *estiva_ofp(void);
void   *estiva_stdfp(void);

#define defop(str)             estiva_defop(str) 
#define getop(str)             estiva_getop(str) 
#define initop(pargc,pargv)    estiva_initop(pargc,pargv)
#define ofp()                  estiva_ofp()
#define stdfp()                estiva_stdfp()

#endif
