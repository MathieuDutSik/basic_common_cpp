#ifndef TYPE_CONVERSION_INCLUDE
#define TYPE_CONVERSION_INCLUDE


template<typename T1, typename T2>
T1 UniversalTypeConversion(T2 const& a)
{
  T1 ret;
  TYPE_CONVERSION(a, ret);
  return ret;
}


template<typename T>
T T_abs(T const& eVal)
{
  if (eVal > 0)
    return eVal;
  T fVal= - eVal;
  return fVal;
}

void NearestInteger_double_int(double const& xI, int & xO)
{
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d=round(xI);
  int xRnd_z=int(xRnd_d);
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr=[&](int const& u) -> double {
    double diff = double(u) - xI;
    return T_abs(diff);
  };
  double err=GetErr(xRnd_z);
  //  std::cerr << "err=" << err << "\n";
  while(true) {
    bool IsOK=true;
    for (int i=0; i<2; i++) {
      int shift=2*i -1;
      int xTest = xRnd_z + shift;
      double TheErr=GetErr(xTest);
      //      std::cerr << "i=" << i << " shift=" << shift << " xTest=" << xTest << " TheErr=" << TheErr << "\n";
      if (TheErr < err) {
	IsOK=false;
	xRnd_z=xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO=xRnd_z;
}





template<typename To>
void NearestInteger_double_To(double const& xI, To & xO)
{
  //  std::cerr << "Temp_common : NearestInteger\n";
  double xRnd_d=round(xI);
  int xRnd_i=int(xRnd_d);
  To xRnd_To=xRnd_i;
  //  std::cerr << "xI=" << xI << "\n";
  auto GetErr=[&](To const& u) -> double {
    double u_doubl;
    GET_DOUBLE(u, u_doubl);
    double diff = u_doubl - xI;
    return T_abs(diff);
  };
  double err=GetErr(xRnd_To);
  //  std::cerr << "err=" << err << "\n";
  while(true) {
    bool IsOK=true;
    for (int i=0; i<2; i++) {
      int shift=2*i -1;
      To xTest = xRnd_To + shift;
      double TheErr=GetErr(xTest);
      //      std::cerr << "i=" << i << " shift=" << shift << " xTest=" << xTest << " TheErr=" << TheErr << "\n";
      if (TheErr < err) {
	IsOK=false;
	xRnd_To=xTest;
      }
    }
    if (IsOK)
      break;
  }
  xO=xRnd_To;
}





template<typename To, typename Ti>
inline typename std::enable_if<(not is_mpreal<Ti>::value) && (not is_double_type<Ti>::value),To>::type UniversalNearestInteger(Ti const& a)
{
  To ret;
  NearestInteger(a, ret);
  //  std::cerr << "Temp_common a=" << a << " ret=" << ret << "\n";
  return ret;
}


template<typename To, typename Ti>
inline typename std::enable_if<is_double_type<Ti>::value,To>::type UniversalNearestInteger(Ti const& a)
{
  To ret;
  NearestInteger_double_To<To>(a, ret);
  //  std::cerr << "Temp_common(I) a=" << a << " ret=" << ret << "\n";
  return ret;
}





#endif
