#ifndef _ESTIVA_TMPFILE_H_
#define _ESTIVA_TMPFILE_H_

#define tmpopen()    estiva_tmpopen()
#define tmpname(fp)  estiva_tmpname(fp)
#define tmpclose(fp) estiva_tmpclose(fp)

FILE *estiva_tmpopen(void);
char *estiva_tmpname(FILE *fp);
void estiva_tmpclose(FILE *fp);

#endif
