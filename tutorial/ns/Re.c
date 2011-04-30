#include "ns.h"

double estiva_Re(void)
{
  static double Reynolds = 1000.0;
  if ( defop("-Re") ) Reynolds = atof(getop("-Re"));
  return Reynolds;
}

double estiva_tau(void)
{
  static double delta_t = 0.001;
  if ( defop("-tau") ) delta_t = atof(getop("-tau"));
  return delta_t;
}
