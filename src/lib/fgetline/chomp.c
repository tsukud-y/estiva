char *estiva_chomp(char *str) {
  char *ret;
  ret = str;
  while(*str != 0 ){
    if( *str == '\n' ) *str = 0;
    str++;
  }
  return ret;
}

