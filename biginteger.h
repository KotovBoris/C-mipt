#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>
#include <vector>

enum class Sign : int8_t { Plus = 1, Minus = -1, Zero = 0 };

void Invert(Sign& sign) {
  if (sign == Sign::Plus) {
    sign = Sign::Minus;
    return;
  }
  if (sign == Sign::Zero) {
    sign = Sign::Zero;
    return;
  }
  sign = Sign::Plus;
}

class BigInteger {
 public:
  BigInteger() : sign(Sign::Zero) {}

  BigInteger(long long number);

  BigInteger(std::string string);

  bool operator<(const BigInteger& second) const;

  bool operator>(const BigInteger& second) const { return second < *this; }

  bool operator==(const BigInteger& second) const {
    return !(*this < second || *this > second);
  }

  bool operator!=(const BigInteger& second) const {
    return (*this < second || *this > second);
  }

  bool operator>=(const BigInteger& second) const { return !(*this < second); }

  bool operator<=(const BigInteger& second) const { return !(*this > second); }

  BigInteger operator-() const;

  BigInteger& operator+=(const BigInteger& second);

  BigInteger& operator-=(const BigInteger& second);

  BigInteger& operator++();

  BigInteger operator++(int);

  BigInteger& operator--();

  BigInteger operator--(int);

  std::string toString() const;

  BigInteger& operator<<=(size_t shift);

  BigInteger operator<<(size_t shift) const;

  BigInteger& operator>>=(size_t shift);

  BigInteger operator>>(size_t shift) const;

  BigInteger& operator*=(int factor);

  BigInteger& operator*=(const BigInteger& factor);

  BigInteger& operator/=(const BigInteger& divisor);

  explicit operator bool() const { return (*this != 0); }

  static const size_t kDigitLength = 9;

 private:
  static const size_t kMod = 1'000'000'000;
  Sign sign;
  std::vector<int> data;

  BigInteger abs() const;

  void deleteZeros();

  int BinSearchDigit(const BigInteger& subtracted);
};

BigInteger BigInteger::abs() const {
  BigInteger copy = *this;
  copy.sign = Sign::Plus;
  return copy;
}

void BigInteger::deleteZeros() {
  long long index = data.size() - 1;
  while (index >= 0 && data[index] == 0) {
    data.pop_back();
    --index;
  }

  if (data.size() == 0) {
    sign = Sign::Zero;
  }
}

int BigInteger::BinSearchDigit(const BigInteger& subtracted) {
  int left = 0;
  if (!data.empty()) {
    left = data[data.size() - 1] /
           (subtracted.data[subtracted.data.size() - 1] + 1);
  }
  int right = kMod;
  while (right - left > 1) {
    int mid = (left + right) / 2;
    BigInteger copy = subtracted;
    copy *= mid;
    if (copy <= *this) {
      left = mid;
    } else {
      right = mid;
    }
  }
  BigInteger copy = subtracted;
  copy *= left;
  *this -= copy;

  return left;
}

BigInteger::BigInteger(long long number) {
  if (number == 0) {
    sign = Sign::Zero;
    return;
  }
  sign = Sign::Plus;
  if (number < 0) {
    sign = Sign::Minus;
    number *= -1;
  }
  while (number > 0) {
    data.push_back(number % kMod);
    number /= kMod;
  }
}

BigInteger::BigInteger(std::string string) {
  if (string == "0") {
    sign = Sign::Zero;
    return;
  }
  sign = Sign::Plus;
  bool minus = false;
  if (string[0] == '-') {
    sign = Sign::Minus;
    minus = true;
  }
  std::string digit(9, '0');
  long long index = string.size() - kDigitLength;
  for (; index >= minus; index -= kDigitLength) {
    for (size_t i = 0; i < kDigitLength; ++i) {
      digit[i] = string[index + i];
    }
    data.push_back(std::stoi(digit));
  }
  if (index + kDigitLength != minus) {
    std::string last_digit(index + kDigitLength - minus, '0');
    for (size_t i = minus; i < index + kDigitLength; ++i) {
      last_digit[i] = string[i];
    }
    data.push_back(std::stoi(last_digit));
  }
}

bool BigInteger::operator<(const BigInteger& second) const {
  if (sign != second.sign) {
    return sign < second.sign;
  }
  if (data.size() != second.data.size()) {
    return data.size() < second.data.size();
  }

  for (long long i = data.size() - 1; i >= 0; --i) {
    if (data[i] != second.data[i]) {
      return data[i] < second.data[i];
    }
  }

  return false;
}

BigInteger BigInteger::operator-() const {
  BigInteger copy = *this;
  Invert(copy.sign);
  return copy;
}

BigInteger& BigInteger::operator+=(const BigInteger& second) {
  if (second == 0) {
    return *this;
  }
  if (*this == 0) {
    *this = second;
    return *this;
  }

  size_t for_next = 0;
  size_t max_size = std::max(data.size(), second.data.size());
  data.resize(max_size + 1, 0);
  BigInteger copy = second;
  copy.data.resize(max_size + 1, 0);

  if (sign == second.sign) {
    for (size_t i = 0; i < max_size + 1; ++i) {
      data[i] += copy.data[i] + for_next;
      for_next = data[i] / kMod;
      data[i] %= kMod;
    }
    deleteZeros();
    return *this;
  }

  int invert = 1;
  if (abs() < copy.abs()) {
    sign = copy.sign;
    invert = -1;
  }

  for_next = 0;
  for (size_t i = 0; i < max_size; ++i) {
    data[i] = invert * (data[i] - copy.data[i]) - for_next;
    for_next = 0;
    if (data[i] < 0) {
      data[i] += kMod;
      for_next = 1;
    }
  }

  deleteZeros();
  return *this;
}

BigInteger& BigInteger::operator-=(const BigInteger& second) {
  *this += -second;
  return *this;
}

BigInteger& BigInteger::operator++() {
  *this += 1;
  return *this;
}

BigInteger BigInteger::operator++(int) {
  BigInteger copy = *this;
  *this += 1;
  return copy;
}

BigInteger& BigInteger::operator--() {
  *this -= 1;
  return *this;
}

BigInteger BigInteger::operator--(int) {
  BigInteger copy = *this;
  *this -= 1;
  return copy;
}

std::string BigInteger::toString() const {
  if (*this == 0) {
    return "0";
  }
  std::string string = std::string((sign == Sign::Minus), '-') +
                       std::to_string(data[data.size() - 1]);

  for (long long i = data.size() - 2; i >= 0; --i) {
    std::string digit = std::to_string(data[i]);
    digit = std::string(kDigitLength - digit.size(), '0') + digit;
    string += digit;
  }
  return string;
}

BigInteger& BigInteger::operator<<=(size_t shift) {
  data.insert(data.begin(), shift, 0);
  return *this;
}

BigInteger BigInteger::operator<<(size_t shift) const {
  BigInteger copy = *this;
  copy <<= shift;
  return copy;
}

BigInteger& BigInteger::operator>>=(size_t shift) {
  data.erase(data.begin(), data.begin() + shift);
  return *this;
}

BigInteger BigInteger::operator>>(size_t shift) const {
  BigInteger copy = *this;
  copy.data.erase(copy.data.begin(), copy.data.begin() + shift);
  return copy;
}

BigInteger& BigInteger::operator*=(int factor) {
  if (factor == 0 || *this == 0) {
    *this = BigInteger();
    return *this;
  }

  if (factor < 0) {
    Invert(sign);
    factor *= -1;
  }

  size_t for_next = 0;
  for (size_t i = 0; i < data.size(); ++i) {
    long long digit = static_cast<long long>(data[i]) * factor + for_next;
    for_next = digit / kMod;
    data[i] = digit % kMod;
  }
  if (for_next > 0) {
    data.push_back(for_next);
  }
  return *this;
}

BigInteger& BigInteger::operator*=(const BigInteger& factor) {
  if (factor == 0 || *this == 0) {
    *this = BigInteger();
    return *this;
  }
  BigInteger product = 0;
  for (size_t i = 0; i < factor.data.size(); ++i) {
    BigInteger summand = *this;
    summand *= factor.data[i];
    product += summand << i;
  }
  *this = product;
  if (factor.sign == Sign::Minus) {
    Invert(sign);
  }
  return *this;
}

BigInteger& BigInteger::operator/=(const BigInteger& divisor) {
  if (divisor.sign == Sign::Zero) {
    return *this;
  }
  long long shift = data.size() - divisor.data.size();
  Sign final_sign = sign;
  if (divisor.sign == Sign::Minus) {
    Invert(final_sign);
  }

  if (shift < 0) {
    *this = 0;
    return *this;
  }
  BigInteger subtracted = divisor << shift;
  subtracted.sign = Sign::Plus;
  sign = Sign::Plus;
  std::vector<int> quotient;
  while (shift >= 0) {
    quotient.push_back(BinSearchDigit(subtracted));
    subtracted >>= 1;
    --shift;
    deleteZeros();
  }

  std::reverse(quotient.begin(), quotient.end());
  data = quotient;
  sign = final_sign;
  deleteZeros();
  return *this;
}

bool operator<(long long first, const BigInteger& second) {
  return (BigInteger(first) < second);
}

bool operator>(long long first, const BigInteger& second) {
  return (BigInteger(first) > second);
}

bool operator==(long long first, const BigInteger& second) {
  return (BigInteger(first) == second);
}

bool operator!=(long long first, const BigInteger& second) {
  return (BigInteger(first) != second);
}

bool operator<=(long long first, const BigInteger& second) {
  return (BigInteger(first) <= second);
}

bool operator>=(long long first, const BigInteger& second) {
  return (BigInteger(first) >= second);
}

std::ostream& operator<<(std::ostream& out, const BigInteger& number) {
  out << number.toString();
  return out;
}

std::istream& operator>>(std::istream& in, BigInteger& number) {
  std::string string;
  in >> string;
  number = BigInteger(string);
  return in;
}

BigInteger operator+(const BigInteger& first, const BigInteger& second) {
  BigInteger copy = first;
  copy += second;
  return copy;
}

BigInteger operator-(const BigInteger& first, const BigInteger& second) {
  BigInteger copy = first;
  copy -= second;
  return copy;
}

BigInteger operator*(const BigInteger& first, const BigInteger& second) {
  BigInteger copy = first;
  copy *= second;
  return copy;
}

BigInteger operator/(const BigInteger& divisible, const BigInteger& divisor) {
  BigInteger copy = divisible;
  copy /= divisor;
  return copy;
}

BigInteger& operator%=(BigInteger& divisible, const BigInteger& divisor) {
  divisible = divisible - (divisible / divisor) * divisor;
  return divisible;
}

BigInteger operator%(const BigInteger& divisible, const BigInteger& divisor) {
  BigInteger copy = divisible;
  copy %= divisor;
  return copy;
}

BigInteger operator"" _bi(unsigned long long number) {
  return BigInteger(number);
}

BigInteger operator"" _bi(const char* cstring, size_t) {
  std::string string = cstring;
  return BigInteger(string);
}

BigInteger GCD(BigInteger first, BigInteger second) {
  return second ? GCD(second, first % second) : first;
}

class Rational {
 public:
  Rational() : numerator(0) {}

  Rational(long long numerator) : numerator(numerator) {}

  Rational(BigInteger numerator) : numerator(numerator) {}

  bool operator<(const Rational& second) const {
    return numerator * second.denominator < second.numerator * denominator;
  }

  bool operator>(const Rational& second) const { return second < *this; }

  bool operator==(const Rational& second) const {
    return (numerator == second.numerator) &&
           (denominator == second.denominator);
  }

  bool operator!=(const Rational& second) const {
    return (*this < second || *this > second);
  }

  bool operator>=(const Rational& second) const { return !(*this < second); }

  bool operator<=(const Rational& second) const { return !(*this > second); }

  Rational operator-() const;

  std::string toString() const;

  Rational& operator+=(const Rational& second);

  Rational& operator-=(const Rational& second);

  Rational& operator*=(const Rational& factor);

  Rational& operator/=(const Rational& divisor);

  std::string asDecimal(size_t presition);

  explicit operator double() {
    return std::stod(asDecimal(BigInteger::kDigitLength * 3));
  }

 private:
  BigInteger numerator;
  BigInteger denominator = 1;

  void shorten();
};

void Rational::shorten() {
  BigInteger gcd = GCD(numerator, denominator);
  if (gcd < 0) {
    gcd = -gcd;
  }
  numerator /= gcd;
  denominator /= gcd;
}

Rational Rational::operator-() const {
  Rational copy = *this;
  copy.numerator = -copy.numerator;
  return copy;
}

std::string Rational::toString() const {
  std::string string = numerator.toString();
  if (denominator != 1) {
    string += "/" + denominator.toString();
  }
  return string;
}

Rational& Rational::operator+=(const Rational& second) {
  numerator = numerator * second.denominator + denominator * second.numerator;
  denominator *= second.denominator;
  shorten();
  return *this;
}

Rational& Rational::operator-=(const Rational& second) {
  *this += -second;
  return *this;
}

Rational& Rational::operator*=(const Rational& factor) {
  numerator *= factor.numerator;
  denominator *= factor.denominator;
  shorten();
  return *this;
}

Rational& Rational::operator/=(const Rational& divisor) {
  numerator *= divisor.denominator;
  denominator *= divisor.numerator;
  shorten();
  if (denominator < 0) {
    numerator = -numerator;
    denominator = -denominator;
  }
  return *this;
}

std::string Rational::asDecimal(size_t presition) {
  size_t real_presition = presition / BigInteger::kDigitLength +
                          (presition % BigInteger::kDigitLength != 0);
  numerator <<= real_presition;
  std::string string = (numerator / denominator).toString();
  size_t after_point = real_presition * BigInteger::kDigitLength;
  if (*this < 0 && string.size() - 1 <= after_point) {
    string.erase(string.begin());
  }
  if (string.size() <= after_point) {
    string = "0." + std::string(after_point - (string.size()), '0') + string;
    if (*this < 0) {
      string = "-" + string;
    }
  } else {
    string.insert(string.end() - after_point, '.');
  }
  string.resize(string.size() - (after_point - presition));
  numerator >>= real_presition;
  return string;
}

std::ostream& operator<<(std::ostream& out, const Rational& number) {
  out << number.toString();
  return out;
}

Rational operator+(const Rational& first, const Rational& second) {
  std::cerr << "rational: " << first << " + " << second << std::endl;
  Rational copy = first;
  copy += second;
  return copy;
}

Rational operator-(const Rational& first, const Rational& second) {
  std::cerr << "rational: " << first << " - " << second << std::endl;
  Rational copy = first;
  copy -= second;
  return copy;
}

Rational operator*(const Rational& first, const Rational& second) {
  std::cerr << "rational: " << first << " * " << second << std::endl;
  Rational copy = first;
  copy *= second;
  return copy;
}

Rational operator/(const Rational& divisible, const Rational& divisor) {
  std::cerr << "rational: " << divisible << " / " << divisor << std::endl;
  Rational copy = divisible;
  copy /= divisor;
  return copy;
}
