#ifndef _ESTIVA_FGETLINE_H_
#define _ESTIVA_FGETLINE_H_

#define fgetline(fp)  estiva_fgetline(fp)
#define chomp(str)    estiva_chomp(str)

char *estiva_fgetline(FILE *fp);
char *estiva_chomp(char *str);

#endif
