#ifndef _ESTIVA_MATVECMPI2_H_
#define _ESTIVA_MATVECMPI2_H_

void    estiva_distributemx(void *A);
void    estiva_efence(int flag);
int     estiva_matvecmpi2();

#define efence(flag)    estiva_efence(flag)
#define distributemx(A) estiva_distributemx(A)

#endif
