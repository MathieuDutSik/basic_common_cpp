Arithmetics
===========

Many of the programs of this repository allow to take different arithmetic.
This is both from a functional aspect and a speed aspect. There is no systematic
way this is done but the following rules are applied. The types are always
template parameter.

Functional types:
  * **Qsqrt2**, **Qsqrt3** and **Qsqrt5**. Those types are used for the fields **Q(sqrt(2))** and similar. For each field we need a type. So, if you need say **Q(sqrt(6))**, you need to modify the code.
  * The function **RealAlgebraic=FileDesc** is for working with real algebraic numbers with the description in **FileDesc**. That file should contain the continous fraction approximant up to some chosen precision. If insufficient, a clean error of failure will be reported.

For both those types, the coefficient are written in entries like **1**, **1+x**, **x/4**, **-3-3x/4**, **3+(3*4)*x^3** or such. No space in the entry.

Speed types:
  * The **mpq_class** is from gmp and is the standard rational type being used.
  * The **mpz_class** is from gmp and is the standard integer type being used.
  * The **boost::multiprecision::cpp_rational** is the pure boost based rational type. Not as fast as gmp, but purely heaer based.
  * The **boost::multiprecision::cpp_int** same as above but for the integers.
  * The **boost::multiprecision::mpq_rational** is boost based encapsulation of gmp library.
  * The **boost::multiprecision::mpq_int** same as above but for the integers.
  * The **Rational<T>** For the rational implementation from an integer type.
  * The **SafeInt64** for integer computation in int64. If the computation go over the limit, then an exception is raised and typically the computation stops. But no wrong result is reported.
  * The **Rational<SafeInt64>** the same as above but for rational fields.

The choice is usually for **mpq_class** which has the least issues.
