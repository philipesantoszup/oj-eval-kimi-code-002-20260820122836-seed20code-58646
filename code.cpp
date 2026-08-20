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
#include <cmath>
#include <algorithm>

// Do not use "using namespace std;"

namespace sjtu {

// Constants for digit compression
const long long BASE = 1000000000;
const int BASE_DIGITS = 9;

using cd = std::complex<double>;
const double PI = acos(-1);

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

// FFT implementation
void int2048::fft(std::vector<cd> & a, bool invert) {
    int n = a.size();
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;
        for (; j & bit; bit >>= 1)
            j ^= bit;
        j ^= bit;
        if (i < j)
            swap(a[i], a[j]);
    }

    for (int len = 2; len <= n; len <<= 1) {
        double ang = 2 * PI / len * (invert ? -1 : 1);
        cd wlen(cos(ang), sin(ang));
        for (int i = 0; i < n; i += len) {
            cd w(1);
            for (int j = 0; j < len / 2; j++) {
                cd u = a[i+j], v = a[i+j+len/2] * w;
                a[i+j] = u + v;
                a[i+j+len/2] = u - v;
                w *= wlen;
            }
        }
    }

    if (invert) {
        for (cd & x : a)
            x /= n;
    }
}

std::vector<long long> int2048::multiply_fft(const std::vector<long long> &a, const std::vector<long long> &b) {
    std::vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
    int n = 1;
    while (n < a.size() + b.size()) 
        n <<= 1;
    fa.resize(n);
    fb.resize(n);

    fft(fa, false);
    fft(fb, false);
    for (int i = 0; i < n; i++)
        fa[i] *= fb[i];
    fft(fa, true);

    std::vector<long long> result(n);
    for (int i = 0; i < n; i++)
        result[i] = std::round(fa[i].real());
    
    long long carry = 0;
    for (int i = 0; i < n; i++) {
        long long total = result[i] + carry;
        result[i] = total % BASE;
        carry = total / BASE;
    }
    
    while (result.size() > 1 && result.back() == 0)
        result.pop_back();
    
    return result;
}

// Helper functions
int int2048::compare_abs(const int2048 &a, const int2048 &b) {
    if (a.digits.size() != b.digits.size()) {
        return a.digits.size() > b.digits.size() ? 1 : -1;
    }
    for (int i = a.digits.size() - 1; i >= 0; i--) {
        if (a.digits[i] != b.digits[i]) {
            return a.digits[i] > b.digits[i] ? 1 : -1;
        }
    }
    return 0;
}

void int2048::add_abs(int2048 &result, const int2048 &a, const int2048 &b) {
    result.digits.clear();
    long long carry = 0;
    int max_len = std::max(a.digits.size(), b.digits.size());
    
    for (int i = 0; i < max_len || carry; i++) {
        long long sum = carry;
        if (i < a.digits.size()) sum += a.digits[i];
        if (i < b.digits.size()) sum += b.digits[i];
        
        result.digits.push_back(sum % BASE);
        carry = sum / BASE;
    }
}

void int2048::subtract_abs(int2048 &result, const int2048 &a, const int2048 &b) {
    // a >= b
    result.digits.clear();
    long long borrow = 0;
    
    for (int i = 0; i < a.digits.size(); i++) {
        long long digit_a = a.digits[i] - borrow;
        long long digit_b = (i < b.digits.size()) ? b.digits[i] : 0;
        
        if (digit_a < digit_b) {
            digit_a += BASE;
            borrow = 1;
        } else {
            borrow = 0;
        }
        
        result.digits.push_back(digit_a - digit_b);
    }
    
    // Remove trailing zeros
    while (result.digits.size() > 1 && result.digits.back() == 0) {
        result.digits.pop_back();
    }
}

void int2048::multiply_abs(int2048 &result, const int2048 &a, const int2048 &b) {
    if ((a.digits.size() == 1 && a.digits[0] == 0) || 
        (b.digits.size() == 1 && b.digits[0] == 0)) {
        result.digits = {0};
        return;
    }
    
    // Use FFT for large numbers, simple multiplication for small ones
    if (a.digits.size() > 100 || b.digits.size() > 100) {
        result.digits = multiply_fft(a.digits, b.digits);
    } else {
        result.digits.assign(a.digits.size() + b.digits.size(), 0);
        
        for (int i = 0; i < a.digits.size(); i++) {
            long long carry = 0;
            for (int j = 0; j < b.digits.size() || carry; j++) {
                long long product = result.digits[i + j] + carry;
                if (j < b.digits.size()) {
                    product += a.digits[i] * b.digits[j];
                }
                result.digits[i + j] = product % BASE;
                carry = product / BASE;
            }
        }
        
        while (result.digits.size() > 1 && result.digits.back() == 0) {
            result.digits.pop_back();
        }
    }
}

int2048 int2048::multiply_by_int(const int2048 &a, long long x) {
    if (x == 0) {
        return int2048(0);
    }
    
    int2048 result;
    result.digits.clear();
    long long carry = 0;
    
    for (int i = 0; i < a.digits.size() || carry; i++) {
        long long product = carry;
        if (i < a.digits.size()) {
            product += a.digits[i] * x;
        }
        result.digits.push_back(product % BASE);
        carry = product / BASE;
    }
    
    while (result.digits.size() > 1 && result.digits.back() == 0) {
        result.digits.pop_back();
    }
    
    return result;
}

void int2048::divide_abs(int2048 &quotient, int2048 &remainder, const int2048 &dividend, const int2048 &divisor) {
    // Uses binary search for each digit
    if (compare_abs(dividend, divisor) < 0) {
        quotient = int2048(0);
        remainder = dividend;
        return;
    }

    remainder.digits.clear();
    quotient.digits.resize(dividend.digits.size() - divisor.digits.size() + 1, 0);

    for (int i = dividend.digits.size() - 1; i >= 0; i--) {
        remainder.digits.insert(remainder.digits.begin(), dividend.digits[i]);
        while (remainder.digits.size() > 1 && remainder.digits.back() == 0) {
            remainder.digits.pop_back();
        }

        long long q = 0;
        long long l = 0, r = BASE - 1;
        
        while (l <= r) {
            long long mid = (l + r) / 2;
            int2048 temp = multiply_by_int(divisor, mid);
            if (compare_abs(temp, remainder) <= 0) {
                q = mid;
                l = mid + 1;
            } else {
                r = mid - 1;
            }
        }

        if (i <= (int)quotient.digits.size() - 1) {
            quotient.digits[i] = q;
        }

        int2048 temp = multiply_by_int(divisor, q);
        int2048 new_remainder;
        subtract_abs(new_remainder, remainder, temp);
        remainder = new_remainder;
    }

    // Clean up quotient
    while (quotient.digits.size() > 1 && quotient.digits.back() == 0) {
        quotient.digits.pop_back();
    }
}

// Constructors
int2048::int2048() : is_negative(false) {
    digits.push_back(0);
}

int2048::int2048(long long n) {
    if (n < 0) {
        is_negative = true;
        n = -n;
    } else {
        is_negative = false;
    }
    if (n == 0) {
        digits.push_back(0);
    } else {
        while (n > 0) {
            digits.push_back(n % BASE);
            n /= BASE;
        }
    }
}

int2048::int2048(const std::string &s) : is_negative(false) {
    digits.push_back(0);
    read(s);
}

int2048::int2048(const int2048 &other) : digits(other.digits), is_negative(other.is_negative) {}

int2048::~int2048() {}

// Assignment operator
int2048 &int2048::operator=(const int2048 &other) {
    if (this != &other) {
        digits = other.digits;
        is_negative = other.is_negative;
    }
    return *this;
}

// Read and print
void int2048::read(const std::string &s) {
    digits.clear();
    is_negative = false;
    
    int start = 0;
    if (s[0] == '-') {
        is_negative = true;
        start = 1;
    } else if (s[0] == '+') {
        start = 1;
    }

    // Skip leading zeros
    while (start < (int)s.size() && s[start] == '0') {
        start++;
    }

    if (start == (int)s.size()) {
        digits.push_back(0);
        is_negative = false;
        return;
    }

    std::string num_str = s.substr(start);
    
    // Process from the end
    for (int i = num_str.size(); i > 0; i -= BASE_DIGITS) {
        int begin = std::max(0, i - BASE_DIGITS);
        std::string part = num_str.substr(begin, i - begin);
        digits.push_back(std::stoll(part));
    }
}

void int2048::print() {
    if (is_negative && !(digits.size() == 1 && digits[0] == 0)) {
        std::cout << '-';
    }
    
    // Print first digit without leading zeros
    std::cout << digits.back();
    
    // Print remaining digits with leading zeros to make BASE_DIGITS
    char prev_fill = std::cout.fill('0');
    for (int i = digits.size() - 2; i >= 0; i--) {
        std::cout.width(BASE_DIGITS);
        std::cout << digits[i];
    }
    std::cout.fill(prev_fill);
}

// Unary operators
int2048 int2048::operator+() const {
    return *this;
}

int2048 int2048::operator-() const {
    int2048 result = *this;
    if (!(result.digits.size() == 1 && result.digits[0] == 0)) {
        result.is_negative = !result.is_negative;
    }
    return result;
}

// Add operations
int2048 &int2048::add(const int2048 &other) {
    *this += other;
    return *this;
}

int2048 add(int2048 a, const int2048 &b) {
    return a + b;
}

int2048 &int2048::operator+=(const int2048 &other) {
    int2048 result;

    if (is_negative == other.is_negative) {
        add_abs(result, *this, other);
        result.is_negative = is_negative;
    } else {
        int cmp = compare_abs(*this, other);
        if (cmp >= 0) {
            subtract_abs(result, *this, other);
            result.is_negative = is_negative;
        } else {
            subtract_abs(result, other, *this);
            result.is_negative = other.is_negative;
        }
    }

    // Handle zero case
    if (result.digits.size() == 1 && result.digits[0] == 0) {
        result.is_negative = false;
    }

    *this = result;
    return *this;
}

int2048 operator+(int2048 a, const int2048 &b) {
    a += b;
    return a;
}

// Subtract operations
int2048 &int2048::minus(const int2048 &other) {
    *this -= other;
    return *this;
}

int2048 minus(int2048 a, const int2048 &b) {
    return a - b;
}

int2048 &int2048::operator-=(const int2048 &other) {
    *this += -other;
    return *this;
}

int2048 operator-(int2048 a, const int2048 &b) {
    a -= b;
    return a;
}

// Multiply operations
int2048 &int2048::operator*=(const int2048 &other) {
    if ((digits.size() == 1 && digits[0] == 0) || 
        (other.digits.size() == 1 && other.digits[0] == 0)) {
        *this = int2048(0);
        return *this;
    }

    int2048 result;
    multiply_abs(result, *this, other);
    result.is_negative = is_negative != other.is_negative;

    *this = result;
    return *this;
}

int2048 operator*(int2048 a, const int2048 &b) {
    a *= b;
    return a;
}

// Divide operations
int2048 &int2048::operator/=(const int2048 &other) {
    if (other.digits.size() == 1 && other.digits[0] == 0) {
        // Division by zero - undefined behavior per problem statement
        return *this;
    }

    int2048 quotient, remainder;
    divide_abs(quotient, remainder, *this, other);

    // Handle signs for floor division
    bool result_negative = is_negative != other.is_negative;
    
    // If there's a remainder and signs differ, we need to adjust for floor division
    if (!(remainder.digits.size() == 1 && remainder.digits[0] == 0) && result_negative) {
        // Add 1 to quotient (since we're rounding down)
        int2048 one(1);
        int2048 new_quotient;
        add_abs(new_quotient, quotient, one);
        quotient = new_quotient;
    }

    quotient.is_negative = result_negative;
    
    // Handle zero case
    if (quotient.digits.size() == 1 && quotient.digits[0] == 0) {
        quotient.is_negative = false;
    }

    *this = quotient;
    return *this;
}

int2048 operator/(int2048 a, const int2048 &b) {
    a /= b;
    return a;
}

// Modulo operations
int2048 &int2048::operator%=(const int2048 &other) {
    // x % y = x - (x / y) * y
    int2048 temp = (*this / other) * other;
    *this -= temp;
    return *this;
}

int2048 operator%(int2048 a, const int2048 &b) {
    a %= b;
    return a;
}

// Stream operators
std::istream &operator>>(std::istream &is, int2048 &num) {
    std::string s;
    is >> s;
    num.read(s);
    return is;
}

std::ostream &operator<<(std::ostream &os, const int2048 &num) {
    if (num.is_negative && !(num.digits.size() == 1 && num.digits[0] == 0)) {
        os << '-';
    }
    
    os << num.digits.back();
    
    char prev_fill = os.fill('0');
    for (int i = num.digits.size() - 2; i >= 0; i--) {
        os.width(BASE_DIGITS);
        os << num.digits[i];
    }
    os.fill(prev_fill);
    
    return os;
}

// Comparison operators
bool operator==(const int2048 &a, const int2048 &b) {
    if (a.is_negative != b.is_negative) {
        // Check if both are zero
        return (a.digits.size() == 1 && a.digits[0] == 0) &&
               (b.digits.size() == 1 && b.digits[0] == 0);
    }
    return int2048::compare_abs(a, b) == 0;
}

bool operator!=(const int2048 &a, const int2048 &b) {
    return !(a == b);
}

bool operator<(const int2048 &a, const int2048 &b) {
    if (a.is_negative != b.is_negative) {
        return a.is_negative;
    }
    
    int cmp = int2048::compare_abs(a, b);
    if (a.is_negative) {
        return cmp > 0;
    } else {
        return cmp < 0;
    }
}

bool operator>(const int2048 &a, const int2048 &b) {
    return b < a;
}

bool operator<=(const int2048 &a, const int2048 &b) {
    return !(a > b);
}

bool operator>=(const int2048 &a, const int2048 &b) {
    return !(a < b);
}

} // namespace sjtu

#endif
