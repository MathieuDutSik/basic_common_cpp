#ifndef PLUS_INFINITY_INCLUDE
#define PLUS_INFINITY_INCLUDE



template<typename T>
struct Tplusinfinity {
  Tplusinfinity(bool const& val1, T const& val2)
  {
    IsInfinity = val1;
    value = val2;
  }
  Tplusinfinity(T const& val)
  {
    IsInfinity=false;
    value = val;
  }
  void SetToInfinity()
  {
    IsInfinity=true;
  }
  Tplusinfinity<T> operator=(Tplusinfinity<T> const& x)
  {
    IsInfinity = x.IsInfinity;
    value = x.value;
  }
  Tplusinfinity<T> operator=(T const& val)
  {
    IsInfinity = false;
    value = val;
  }
  bool GetInfinity() const
  {
    return IsInfinity;
  }
  T GetValue() const
  {
    return value;
  }
private:
  bool IsInfinity;
  T value;
};


template<typename T>
bool operator==(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() && y.GetInfinity())
    return true;
  if (x.GetInfinity() != y.GetInfinity())
    return false;
  return x.GetValue() == y.GetValue();
}


template<typename T>
bool operator<(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() && y.GetInfinity())
    return false;
  if (!x.GetInfinity() && y.GetInfinity())
    return true;
  if (x.GetInfinity() && !y.GetInfinity())
    return false;
  return x.GetValue() < y.GetValue();
}



template<typename T>
bool operator<=(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() && y.GetInfinity())
    return true;
  if (!x.GetInfinity() && y.GetInfinity())
    return true;
  if (x.GetInfinity() && !y.GetInfinity())
    return false;
  return x.GetValue() <= y.GetValue();
}


template<typename T>
bool operator>(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() && y.GetInfinity())
    return false;
  if (!x.GetInfinity() && y.GetInfinity())
    return false;
  if (x.GetInfinity() && !y.GetInfinity())
    return true;
  return x.GetValue() > y.GetValue();
}


template<typename T>
Tplusinfinity<T> operator+(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() || y.GetInfinity())
    return Tplusinfinity<T>(true, 0);
  return Tplusinfinity<T>(false, x.GetValue() + y.GetValue());
}

template<typename T>
Tplusinfinity<T> operator+(Tplusinfinity<T> const& x, T const& y)
{
  return Tplusinfinity<T>(x.GetInfinity(), x.GetValue() + y);
}

template<typename T>
Tplusinfinity<T> operator+(T const& x, Tplusinfinity<T> const& y)
{
  return Tplusinfinity<T>(y.GetInfinity(), x + y.GetValue());
}

template<typename T>
Tplusinfinity<T> operator*(Tplusinfinity<T> const& x, Tplusinfinity<T> const& y)
{
  if (x.GetInfinity() || y.GetInfinity())
    return Tplusinfinity<T>(true, 0);
  return Tplusinfinity<T>(false, x.GetValue() * y.GetValue());
}

template<typename T>
Tplusinfinity<T> operator*(Tplusinfinity<T> const& x, T const& y)
{
  return Tplusinfinity<T>(x.GetInfinity(), x.GetValue() * y);
}

template<typename T>
Tplusinfinity<T> operator*(T const& x, Tplusinfinity<T> const& y)
{
  return Tplusinfinity<T>(y.GetInfinity(), x * y.GetValue());
}




#endif
