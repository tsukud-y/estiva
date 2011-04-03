#include <stdio.h>
#include <stdlib.h>
#include <estiva/ary.h>

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
  abort();
}

size_t estiva_dim0(void* v)
{
  check(v);
  return r->dim0;
}

long estiva_dim1(void* v)
{
  check(v);
  return r->dim1;
}

long estiva_dim2(void* v)
{
  check(v);
  return r->dim2;
}

static void *pointer_array[1024];
static long array_index;

static void freefunc(void)
{
  long i;
  for (i=0; i<array_index; i++) free(pointer_array[i]);
}

static void atexit_free(void)
{
  static int init;
  if ( init == 0 ) {
    init = 1;    
    atexit(freefunc);
  }
  if ( array_index > 1000 ) { printf("Too many ary\n"); abort(); }
}

static void* alloc(size_t n)
{
  void *p;
  p = calloc(1,n);
  if(p != NULL) {
    pointer_array[array_index++] = p;
    atexit_free();
    return p;
  }
  else{
    fprintf(stderr,"ary(): Can't alloc memory!\n");
    abort();
  }
  return NULL;
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

  *v = NULL;
}

void estiva_ary1(void** v, long n_1, size_t o)
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

void estiva_ary2(void** v, long m_1, long n_1, size_t o)
{
  if(*v == NULL){ new_ary2(v,m_1,n_1,o); return;}
  
  if(m_1==dim2(*v)+1 && n_1==dim1(*v)+1) return;
  
  del_ary2(v); 

  if(m_1!=0 && n_1!=0) new_ary2(v,m_1,n_1,o); 
}
