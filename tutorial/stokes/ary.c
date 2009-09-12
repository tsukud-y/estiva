#include <stdio.h>
#include <stdlib.h>

typedef struct{
  size_t dim0;
  long dim1, dim2;
  long padding;
} dim;

static dim* r;

static void check(void* v)
{
  if(v != NULL){
    r = v;
    r--;
    return;
  }
  fprintf(stderr,"dim(): You need alloc!\n");
  exit(1);
}

size_t dim0(void* v)
{
  check(v);
  return r->dim0;
}

long dim1(void* v)
{
  check(v);
  return r->dim1;
}

long dim2(void* v)
{
  check(v);
  return r->dim2;
}

static void* alloc(size_t n)
{
  void *p;
  p = calloc(1,n);
  if(p != NULL) return p;
  fprintf(stderr,"ary(): Can't alloc memory!\n");
  exit(1);
}

static void new_ary1(void** v, long n_1, size_t o)
{
  long n;
  dim* r;
  
  n = n_1-1;
  
  r = alloc(sizeof(dim)+n_1*o);
  
  r->dim2 = 0;   r->dim1 = n;   r->dim0 = o;
  
  r++;
  
  *v = r;
}

static void del_ary1(void** v)
{
  dim* r;

  r = *v;
  r--;
  free(r);

  *v == NULL;
}

void femlib_ary1(void** v, long n_1, size_t o)
{
  if(*v == NULL){ new_ary1(v,n_1,o); return;}

  if(n_1 == dim1(*v)+1) return;
  
  del_ary1(v);

  if(n_1 != 0) new_ary1(v,n_1,o);
}


static void new_ary2(void** v, long m_1, long n_1, size_t o)
{
  long i, m, n;
  dim* r;
  char** a;

  m = m_1-1;    n = n_1-1;

  r = alloc(sizeof(dim) + m_1*sizeof(void *));
  
  r->dim2 = m;    r->dim1 = n;    r->dim0 = o;
  
  r++;
  
  *v = r;
  a = *v;
  
  a[0] = alloc(m_1*n_1*o);
  for(i=1;i<=m;i++) a[i] = &a[i-1][n_1*o];
}

static void del_ary2(void** v)
{
  dim* r;
  char** a;
  long i, m;

  a = *v;
  m = dim2(a);
  for(i=0;i<=m;i++) free(a[i]);

  r = *v;
  r--;
  free(r);

  *v = NULL;
}

void femlib_ary2(void** v, long m_1, long n_1, size_t o)
{
  if(*v == NULL){ new_ary2(v,m_1,n_1,o); return;}
  
  if(m_1==dim2(*v)+1 && n_1==dim1(*v)+1) return;
  
  del_ary2(v); 

  if(m_1!=0 && n_1!=0) new_ary2(v,m_1,n_1,o); 
}
