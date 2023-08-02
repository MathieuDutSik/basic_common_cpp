# Basic Common Cpp

This is a set of code used in C++. All of this is basic but used in many projects.
Author is Mathieu Dutour Sikiric (e-mail: mathieu.dutour@gmail.com), but anyone
can contribute.


# Access to the source code

Since this repository uses submodules, the cloning command is

```sh
$ git clone https://github.com/MathieuDutSik/basic_common_cpp.git --recursive
```

In order to update the submodule the command is
```sh
$ git submodule update --remote
```

# Design

## Namelist

Forran provides a quite powerful system for doing input files. It is named Namelist
and the extension of those files are *nml*. We provide this functionality in C++ with
a more rigorous way: An entry can be defined only one time for example.

## Scalar types

By scalar types we mean the one used for matrices, vectors, etc. Example could be integer,
**mpz_class**, **mpq_class**, etc.

We want to use the C++ template capabilities for having multiple types. This requires
some strict decisions about how things are done.

The conversion is done via **UniversalScalarConversion<To,Ti>** and similar with **To** the
output type and **Ti** the input type. The types are always strict, there is never any
implicit conversions.

The operation are of the kind **x > y** with **x** and **y** of the same type or **x** of
a constructed type and **y** an integer. The integer type is always second. The code
looks like
* **friend bool operator>(T const& x, T const& y)**
* **friend bool operator>(T const& x, int const& y)**

This goes the same way for **==**, **!=**, **<=**, **>=**, **<** and **>**.

The **operator<<**, **operator>>**, **hash()**, **load/save** (for boost serialization) have
to be defined for all the types.

The rules for the constructors are not as well defined. But it looks like the constructor
from integer is forbidden because it leads to inconsistencies with **x > 0** because
either **operator(T const& x, T const& y)** or **operator(T const& x, int const& y)** could be used.

