#ifndef FILE_H
#define FILE_H
extern int estiva_forFILE(FILE *fp);
#define  forFILE(fp) while(estiva_forFILE(fp)) 
extern char *estiva_S(int n);
#define  S(n) estiva_S(n) 
extern void estiva_FILE_cp(FILE *fp,FILE *out);
#define  FILE_cp(fp,out) estiva_FILE_cp(fp,out) 
#endif
