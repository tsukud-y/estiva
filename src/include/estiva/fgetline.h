#ifndef _ESTIVA_FGETLINE_H_
#define _ESTIVA_FGETLINE_H_

char   *estiva_fgetline(FILE *fp);
char   *estiva_chomp(char *str);
long    estiva_fsize(void *fp);
long    estiva_flines(void *fp);

#define chomp(str)    estiva_chomp(str)
#define fgetline(fp)  estiva_fgetline(fp)
#define flines(fp)    estiva_flines(fp)
#define fsize(fp)     estiva_fsize(fp)

#endif
