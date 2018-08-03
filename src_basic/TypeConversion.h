#ifndef TYPE_CONVERSION_INCLUDE
#define TYPE_CONVERSION_INCLUDE


template<typename T1, typename T2>
T1 UniversalTypeConversion(T2 const& a)
{
  T1 ret;
  TYPE_CONVERSION(a, ret);
  return ret;
}

#endif
