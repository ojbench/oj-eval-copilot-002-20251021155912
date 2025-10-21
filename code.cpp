// Combined submission file for int2048
#include <complex>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <vector>
#include <string>

#pragma once
#ifndef SJTU_BIGINTEGER
#define SJTU_BIGINTEGER

namespace sjtu {
class int2048 {
private:
  static const int base = 10000; // 1e4
  static const int base_digits = 4;
  std::vector<int> d; // little-endian digits
  bool neg;
  void trim();
  static int absCmp(const int2048 &, const int2048 &);
  static int2048 absAdd(const int2048 &, const int2048 &);
  static int2048 absSub(const int2048 &, const int2048 &); // assume |a|>=|b|
  static int2048 mulSchool(const std::vector<int> &, const std::vector<int> &);
  static std::pair<int2048, int2048> divmod(const int2048 &, const int2048 &);
public:
  // 构造函数
  int2048();
  int2048(long long);
  int2048(const std::string &);
  int2048(const int2048 &);

  // 读入一个大整数
  void read(const std::string &);
  // 输出储存的大整数，无需换行
  void print();

  // 加上一个大整数
  int2048 &add(const int2048 &);
  // 返回两个大整数之和
  friend int2048 add(int2048, const int2048 &);

  // 减去一个大整数
  int2048 &minus(const int2048 &);
  // 返回两个大整数之差
  friend int2048 minus(int2048, const int2048 &);

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

// Implementation
namespace sjtu {
namespace {
using cd = std::complex<double>;
const double PI = 3.1415926535897932384626433832795;

void fft(std::vector<cd> &a, bool invert) {
  int n = (int)a.size();
  for (int i = 1, j = 0; i < n; ++i) {
    int bit = n >> 1;
    for (; j & bit; bit >>= 1) j ^= bit;
    j ^= bit;
    if (i < j) std::swap(a[i], a[j]);
  }
  for (int len = 2; len <= n; len <<= 1) {
    double ang = 2 * PI / len * (invert ? -1 : 1);
    cd wlen = std::polar(1.0, ang);
    for (int i = 0; i < n; i += len) {
      cd w(1.0, 0.0);
      for (int j = 0; j < len / 2; ++j) {
        cd u = a[i + j];
        cd v = a[i + j + len / 2] * w;
        a[i + j] = u + v;
        a[i + j + len / 2] = u - v;
        w *= wlen;
      }
    }
  }
  if (invert) {
    for (int i = 0; i < n; ++i) a[i] /= n;
  }
}

static void multiplyVectors(const std::vector<int> &a, const std::vector<int> &b, std::vector<int> &res) {
  if (a.empty() || b.empty()) { res.clear(); return; }
  const int BASE = 10000;
  std::vector<cd> fa(a.begin(), a.end()), fb(b.begin(), b.end());
  int n = 1;
  while (n < (int)a.size() + (int)b.size()) n <<= 1;
  fa.resize(n); fb.resize(n);
  fft(fa, false); fft(fb, false);
  for (int i = 0; i < n; i++) fa[i] *= fb[i];
  fft(fa, true);
  res.assign(n, 0);
  long long carry = 0;
  for (int i = 0; i < n; i++) {
    double realv = fa[i].real();
    long long t = (long long)(realv + (realv >= 0 ? 0.5 : -0.5)) + carry;
    res[i] = int(t % BASE);
    carry = t / BASE;
  }
  while (carry) { res.push_back(int(carry % BASE)); carry /= BASE; }
  while (!res.empty() && res.back() == 0) res.pop_back();
}
}

void int2048::trim() {
  while (!d.empty() && d.back() == 0) d.pop_back();
  if (d.empty()) neg = false;
}

int int2048::absCmp(const int2048 &a, const int2048 &b) {
  if (a.d.size() != b.d.size()) return a.d.size() < b.d.size() ? -1 : 1;
  for (int i = (int)a.d.size() - 1; i >= 0; --i) if (a.d[i] != b.d[i]) return a.d[i] < b.d[i] ? -1 : 1;
  return 0;
}

int2048 int2048::absAdd(const int2048 &a, const int2048 &b) {
  int2048 r; r.neg = false; r.d.resize(std::max(a.d.size(), b.d.size()) + 1, 0);
  int carry = 0; size_t i = 0;
  for (; i < a.d.size() || i < b.d.size() || carry; ++i) {
    long long sum = carry;
    if (i < a.d.size()) sum += a.d[i];
    if (i < b.d.size()) sum += b.d[i];
    r.d[i] = int(sum % base);
    carry = int(sum / base);
  }
  r.d.resize(i);
  r.trim();
  return r;
}

int2048 int2048::absSub(const int2048 &a, const int2048 &b) {
  int2048 r; r.neg = false; r.d.resize(a.d.size(), 0);
  int carry = 0; size_t i = 0;
  for (; i < a.d.size(); ++i) {
    long long cur = a.d[i] - carry - (i < b.d.size() ? b.d[i] : 0);
    if (cur < 0) cur += base, carry = 1; else carry = 0;
    r.d[i] = int(cur);
  }
  r.trim();
  return r;
}

int2048 int2048::mulSchool(const std::vector<int> &a, const std::vector<int> &b) {
  int2048 r; if (a.empty() || b.empty()) return r;
  r.d.assign(a.size() + b.size(), 0); r.neg = false;
  for (size_t i = 0; i < a.size(); ++i) {
    long long carry = 0;
    for (size_t j = 0; j < b.size() || carry; ++j) {
      long long cur = r.d[i + j] + carry + 1LL * a[i] * (j < b.size() ? b[j] : 0);
      r.d[i + j] = int(cur % base);
      carry = cur / base;
    }
  }
  r.trim();
  return r;
}

std::pair<int2048, int2048> int2048::divmod(const int2048 &A, const int2048 &B) {
  int2048 zero;
  if (B.d.empty()) return {zero, zero};
  int2048 a = A, b = B; a.neg = b.neg = false;
  if (absCmp(a, b) < 0) return {zero, a};
  int m = (int)b.d.size();
  int n = (int)a.d.size();
  int norm = base / (b.d.back() + 1);
  if (norm > 1) {
    long long carry = 0; for (size_t i = 0; i < a.d.size() || carry; ++i) {
      if (i == a.d.size()) a.d.push_back(0);
      long long cur = 1LL * (i < a.d.size() ? a.d[i] : 0) * norm + carry;
      a.d[i] = int(cur % base); carry = cur / base;
    }
    carry = 0; for (size_t i = 0; i < b.d.size() || carry; ++i) {
      if (i == b.d.size()) b.d.push_back(0);
      long long cur = 1LL * (i < b.d.size() ? b.d[i] : 0) * norm + carry;
      b.d[i] = int(cur % base); carry = cur / base;
    }
    a.trim(); b.trim();
  }
  n = (int)a.d.size(); m = (int)b.d.size();
  int2048 q; q.d.assign(n - m + 1, 0); q.neg = false;
  std::vector<int> u = a.d; u.push_back(0);
  const std::vector<int> &v = b.d;
  for (int k = n - m; k >= 0; --k) {
    long long u2 = 1LL * u[k + m] * base + u[k + m - 1];
    long long qhat = u2 / v[m - 1];
    long long rhat = u2 % v[m - 1];
    if (qhat >= base) { qhat = base - 1; rhat = u2 - qhat * v[m - 1]; }
    while (m > 1 && qhat * v[m - 2] > base * rhat + u[k + m - 2]) {
      --qhat; rhat += v[m - 1]; if (rhat >= base) break;
    }
    long long borrow = 0;
    for (int j = 0; j < m; ++j) {
      long long cur = u[k + j] - qhat * 1LL * v[j] - borrow;
      borrow = 0;
      if (cur < 0) { long long brr = (-cur + base - 1) / base; cur += brr * base; borrow = brr; }
      u[k + j] = int(cur);
    }
    long long curk = u[k + m] - borrow;
    if (curk < 0) {
      --qhat;
      long long c = 0;
      for (int j = 0; j < m; ++j) {
        long long cur2 = u[k + j] + v[j] + c;
        if (cur2 >= base) { cur2 -= base; c = 1; } else c = 0;
        u[k + j] = int(cur2);
      }
      u[k + m] += int(c);
    } else {
      u[k + m] = int(curk);
    }
    q.d[k] = int(qhat);
  }
  int2048 r; r.d.assign(m, 0);
  for (int i = 0; i < m; ++i) r.d[i] = u[i];
  if (norm > 1) {
    long long carry = 0;
    for (int i = m - 1; i >= 0; --i) {
      long long cur = r.d[i] + carry * base;
      r.d[i] = int(cur / norm);
      carry = cur % norm;
    }
  }
  q.trim(); r.trim();
  return {q, r};
}

int2048::int2048() : d(), neg(false) {}

int2048::int2048(long long v) : d(), neg(false) {
  if (v < 0) { neg = true; v = -v; }
  while (v) { d.push_back(int(v % base)); v /= base; }
}

int2048::int2048(const std::string &s) : d(), neg(false) { read(s); }

int2048::int2048(const int2048 &o) = default;

void int2048::read(const std::string &s) {
  d.clear(); neg = false;
  int i = 0; while (i < (int)s.size() && isspace((unsigned char)s[i])) ++i;
  bool nneg = false;
  if (i < (int)s.size() && (s[i] == '+' || s[i] == '-')) { nneg = (s[i] == '-'); ++i; }
  while (i < (int)s.size() && s[i] == '0') ++i;
  for (int j = (int)s.size() - 1; j >= i; j -= base_digits) {
    int x = 0;
    int l = std::max(i, j - base_digits + 1);
    for (int k = l; k <= j; ++k) x = x * 10 + (s[k] - '0');
    d.push_back(x);
  }
  neg = nneg && !d.empty();
  trim();
}

void int2048::print() {
  if (d.empty()) { std::cout << 0; return; }
  if (neg) std::cout << '-';
  int nsz = (int)d.size();
  std::cout << d.back();
  char buf[16];
  for (int i = nsz - 2; i >= 0; --i) {
    std::snprintf(buf, sizeof(buf), "%0*d", base_digits, d[i]);
    std::cout << buf;
  }
}

int2048 &int2048::add(const int2048 &rhs) { return (*this) += rhs; }
int2048 add(int2048 lhs, const int2048 &rhs) { lhs += rhs; return lhs; }
int2048 &int2048::minus(const int2048 &rhs) { return (*this) -= rhs; }
int2048 minus(int2048 lhs, const int2048 &rhs) { lhs -= rhs; return lhs; }

int2048 int2048::operator+() const { return *this; }
int2048 int2048::operator-() const { int2048 r(*this); if (!r.d.empty()) r.neg = !r.neg; return r; }

int2048 &int2048::operator=(const int2048 &o) = default;

int2048 &int2048::operator+=(const int2048 &rhs) {
  if (neg == rhs.neg) {
    *this = absAdd(*this, rhs);
    this->neg = rhs.neg;
  } else {
    int cmp = absCmp(*this, rhs);
    if (cmp >= 0) {
      *this = absSub(*this, rhs);
    } else {
      int2048 r = absSub(rhs, *this);
      r.neg = rhs.neg;
      *this = r;
    }
  }
  trim();
  return *this;
}

int2048 operator+(int2048 lhs, const int2048 &rhs) { lhs += rhs; return lhs; }

int2048 &int2048::operator-=(const int2048 &rhs) {
  int2048 t = rhs; if (!t.d.empty()) t.neg = !t.neg; return (*this) += t;
}

int2048 operator-(int2048 lhs, const int2048 &rhs) { lhs -= rhs; return lhs; }

int2048 &int2048::operator*=(const int2048 &rhs) {
  bool res_neg = neg ^ rhs.neg;
  int n = (int)d.size(), m = (int)rhs.d.size();
  int2048 r;
  if ((long long)n * m <= 20000) r = mulSchool(d, rhs.d);
  else { std::vector<int> prod; multiplyVectors(d, rhs.d, prod); r.d.swap(prod); }
  r.neg = res_neg && !r.d.empty();
  *this = r;
  return *this;
}

int2048 operator*(int2048 lhs, const int2048 &rhs) { lhs *= rhs; return lhs; }

int2048 &int2048::operator/=(const int2048 &rhs) {
  bool res_neg = neg ^ rhs.neg;
  auto qr = divmod(*this, rhs);
  int2048 q = qr.first; int2048 r = qr.second; // q,r are magnitudes for |a|/|b|
  if (res_neg) {
    if (!r.d.empty()) {
      int2048 one(1);
      q = absAdd(q, one);
      q.neg = true;
    } else {
      q.neg = true;
    }
  } else {
    q.neg = false;
  }
  q.trim();
  *this = q;
  return *this;
}

int2048 operator/(int2048 lhs, const int2048 &rhs) { lhs /= rhs; return lhs; }

int2048 &int2048::operator%=(const int2048 &rhs) {
  bool diff = neg ^ rhs.neg;
  auto qr = divmod(*this, rhs);
  int2048 r = qr.second; // magnitude remainder for |a|/|b|
  if (!r.d.empty() && diff) {
    int2048 bb = rhs; bb.neg = false;
    r = absSub(bb, r);
  }
  r.neg = (!r.d.empty()) && rhs.neg;
  r.trim();
  *this = r;
  return *this;
}

int2048 operator%(int2048 lhs, const int2048 &rhs) { lhs %= rhs; return lhs; }

std::istream &operator>>(std::istream &is, int2048 &x) {
  std::string s; is >> s; x.read(s); return is;
}
std::ostream &operator<<(std::ostream &os, const int2048 &x) {
  if (x.d.empty()) { os << 0; return os; }
  if (x.neg) os << '-';
  os << x.d.back();
  char buf[16];
  for (int i = (int)x.d.size() - 2; i >= 0; --i) {
    std::snprintf(buf, sizeof(buf), "%0*d", int2048::base_digits, x.d[i]);
    os << buf;
  }
  return os;
}

bool operator==(const int2048 &lhs, const int2048 &rhs) {
  return lhs.neg == rhs.neg && lhs.d == rhs.d;
}
bool operator!=(const int2048 &lhs, const int2048 &rhs) { return !(lhs == rhs); }

bool operator<(const int2048 &lhs, const int2048 &rhs) {
  if (lhs.neg != rhs.neg) return lhs.neg;
  int cmp = int2048::absCmp(lhs, rhs);
  return lhs.neg ? (cmp > 0) : (cmp < 0);
}
bool operator>(const int2048 &lhs, const int2048 &rhs) { return rhs < lhs; }
bool operator<=(const int2048 &lhs, const int2048 &rhs) { return !(rhs < lhs); }
bool operator>=(const int2048 &lhs, const int2048 &rhs) { return !(lhs < rhs); }

} // namespace sjtu
