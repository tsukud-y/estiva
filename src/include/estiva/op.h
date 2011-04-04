#ifndef _ESTIVA_OP_H_
#define _ESTIVA_OP_H_

extern int    *estiva_pargc;
extern char ***estiva_pargv;

int     estiva_defop(char *str);
char   *estiva_getop(char *str);
void    estiva_initop(int *pargc, char ***pargv);
void   *estiva_ofp(void);
void   *estiva_stdfp(void);

#define defop(str)             estiva_defop(str) 
#define getop(str)             estiva_getop(str) 
#define initop(argc,argv)      estiva_initop(&argc,&argv)
#define ofp()                  estiva_ofp()
#define stdfp()                estiva_stdfp()

#endif
