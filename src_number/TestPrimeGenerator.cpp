// Copyright (C) 2022 Mathieu Dutour Sikiric <mathieu.dutour@gmail.com>
// clang-format off
#include "NumberTheory.h"
#include "NumberTheoryRealField.h"
#include "NumberTheorySafeInt.h"
#include "NumberTheoryQuadField.h"
#include "factorizations.h"
#include <iostream>
// clang-format on

template <typename T>
void test_prime_generator(int num_primes) {
  PrimeGenerator<T> gen;

  std::cerr << "Generating " << num_primes << " primes starting from 10001:\n";

  for (int i = 0; i < num_primes; i++) {
    T prime = gen.get_prime();
    std::cout << prime;
    if (i < num_primes - 1) {
      std::cout << " ";
    }
  }
  std::cout << "\n";
}

void process(std::string const &arith, int num_primes) {
  if (arith == "safe_integer") {
    using T = SafeInt64;
    return test_prime_generator<T>(num_primes);
  }
  if (arith == "int") {
    using T = int;
    return test_prime_generator<T>(num_primes);
  }
  if (arith == "mpz") {
    using T = mpz_class;
    return test_prime_generator<T>(num_primes);
  }
  std::cerr << "Failed to find a matching entry for arith\n";
  throw TerminalException{1};
}

int main(int argc, char *argv[]) {
  HumanTime time;
  try {
    if (argc != 3) {
      std::cerr << "This program is used as\n";
      std::cerr << "TestPrimeGenerator [arith] [num_primes]\n";
      std::cerr << "\n";
      std::cerr << "    where\n";
      std::cerr << "arith: The arithmetic type (safe_integer, int, long, mpz)\n";
      std::cerr << "num_primes: The number of primes to generate\n";
      std::cerr << "\n";
      std::cerr << "Example: TestPrimeGenerator safe_integer 10\n";
      return -1;
    }
    //
    std::string arith = argv[1];
    int num_primes = std::stoi(argv[2]);

    if (num_primes <= 0) {
      std::cerr << "Number of primes must be positive\n";
      return -1;
    }

    process(arith, num_primes);
    //
    std::cerr << "Normal termination of the program\n";
  } catch (TerminalException const &e) {
    std::cerr << "Something wrong happened in the computation\n";
    exit(e.eVal);
  } catch (std::exception const &e) {
    std::cerr << "Exception: " << e.what() << "\n";
    return -1;
  }
  runtime(time);
}
