#ifndef _MATUTIL_H_
#define _MATUTIL_H_
extern int matutil_matf(void *argv0,char *fmt,...);

static int (*matf)(void *argv0,char *fmt,...) = matutil_matf;
#endif
