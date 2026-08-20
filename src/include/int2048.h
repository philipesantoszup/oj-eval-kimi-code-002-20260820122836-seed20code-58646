#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

// Integer 1:
// Implement a signed big integer class that only needs to support simple addition and subtraction

// Integer 2:
// Implement a signed big integer class that supports addition, subtraction, multiplication, and division, and overload related operators

// Do not use any header files other than the following
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>

// Do not use "using namespace std;"

namespace sjtu {

// Constants for digit compression
const long long BASE = 1000000000;
const int BASE_DIGITS = 9;

class int2048 {
private:
    std::vector<long long> digits;
    bool is_negative;

    // Helper functions
    static int compare_abs(const int2048 &a, const int2048 &b);
    static void add_abs(int2048 &result, const int2048 &a, const int2048 &b);
    static void subtract_abs(int2048 &result, const int2048 &a, const int2048 &b); // a >= b
    static void multiply_abs(int2048 &result, const int2048 &a, const int2048 &b);
    static void divide_abs(int2048 &quotient, int2048 &remainder, const int2048 &dividend, const int2048 &divisor);
    static int2048 multiply_by_int(const int2048 &a, long long x);
    static void fft(std::vector<std::complex<double>> & a, bool invert);
    static std::vector<long long> multiply_fft(const std::vector<long long> &a, const std::vector<long long> &b);

public:
    // Constructors
    int2048();
    int2048(long long);
    int2048(const std::string &);
    int2048(const int2048 &);
    ~int2048();

    // The parameter types of the following functions are for reference only, you can choose to use constant references or not
    // If needed, you can add other required functions yourself
    // ===================================
    // Integer1
    // ===================================

    // Read a big integer
    void read(const std::string &);
    // Output the stored big integer, no need for newline
    void print();

    // Add a big integer
    int2048 &add(const int2048 &);
    // Return the sum of two big integers
    friend int2048 add(int2048, const int2048 &);

    // Subtract a big integer
    int2048 &minus(const int2048 &);
    // Return the difference of two big integers
    friend int2048 minus(int2048, const int2048 &);

    // ===================================
    // Integer2
    // ===================================

    int2048 operator+() const;
    int2048 operator-() const;

    int2048 &operator=(const int2048 &);

    int2048 &operator+=(const int2048 &);
    friend int2048 operator+(int2048, const int2048 &);

    int2048 &operator-=(const int2048 &);
    friend int2048 operator-(int2048, const int2048 &);

    int2048 &operator*=(const int2048 &);
    friend int2048 operator*(int2048, const int2048 &);

    int2048 &operator/=(const int2048 &);
    friend int2048 operator/(int2048, const int2048 &);

    int2048 &operator%=(const int2048 &);
    friend int2048 operator%(int2048, const int2048 &);

    friend std::istream &operator>>(std::istream &, int2048 &);
    friend std::ostream &operator<<(std::ostream &, const int2048 &);

    friend bool operator==(const int2048 &, const int2048 &);
    friend bool operator!=(const int2048 &, const int2048 &);
    friend bool operator<(const int2048 &, const int2048 &);
    friend bool operator>(const int2048 &, const int2048 &);
    friend bool operator<=(const int2048 &, const int2048 &);
    friend bool operator>=(const int2048 &, const int2048 &);
};
} // namespace sjtu

#endif
