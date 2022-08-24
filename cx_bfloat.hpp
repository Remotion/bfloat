#pragma once
#ifndef HUGEFLOAT_RE_HPP_
#define HUGEFLOAT_RE_HPP_
// =====================================================================================================================
//  cx_bfloat.hpp
//  #include "cx_bfloat.hpp"
//  Remotion (C) 2022 - All Rights Reserved
// =====================================================================================================================
// https://bellard.org/libbf/
// =====================================================================================================================

#include "libbf.hpp"
#include <initializer_list>

namespace cx {

// number of bits per base 10 digit 
constexpr auto BITS_PER_DIGIT = 3.32192809488736234786;

namespace detail {

struct bf_ctx_cleanup 
{
	 bf_ctx_cleanup() { }
	~bf_ctx_cleanup() { bf_clear_cache(&g_bf_ctx_default); }
	void force_instantiate() const {}
};

static const bf_ctx_cleanup g_bf_ctx_cleanup;

} // namespace detail

namespace detail {

template <class T>
constexpr T ipow(const T& t0, int32_t n) {
	T a = t0;
	int32_t b = n;
	bool const recip = b < 0;
	T r{static_cast<T>(1)};
	while (1) {
		if (b & 1) { r *= a; }
		b /= 2;
		if (b == 0) { break; }
		a *= a;
	}
	return recip ? T(1) / r : r;
}

template <class T, int32_t N>
constexpr T ipow(const T& t0) {
	T a = t0;
	constexpr bool recip = N < 0;
	int32_t b = N;
	T r {1};
	while (1) {
		if (b & 1) { r *= a; }
		b /= 2;
		if (b == 0) { break; }
		a *= a;
	}
	return recip ? T(1) / r : r;
}


template<typename T, int32_t exp>
constexpr T powi_fast(T x) {
	if constexpr (exp > 1) {
		if constexpr (exp % 2) { const T a = powi_fast<T,exp / 2>(x); return x * a * a; } 
		else { const T a = powi_fast<T,exp / 2>(x); return a * a; }
	} else 
	if constexpr (exp == 0) {
		return T(1.0);  // x^0 = 1
	} else 
	if constexpr (exp < 0) { // exp < 0
		return T(1.0) / powi_fast<T,-exp>(x); // x^(-n) = 1/x^n 
	}
	return x; // unused ?
}

/// initial approximations for nth roots generally for IEEE 754 doubles 
template <int32_t N>
__forceinline constexpr double nth_rootd_aprox(double x) {
	constexpr int32_t ebits = 11;
	constexpr int32_t fbits = 52;
	constexpr int64_t bias = (1 << (ebits - 1)) - 1;

	int64_t i = std::bit_cast<int64_t>(x);
	i = (i - (bias << fbits)) / N + (bias << fbits);

	double ax = std::bit_cast<double>(i);

	/// now do 4 Newton�Raphson (NR) iterations to improve accuracy!
	constexpr int32_t n_1 = N - 1;
	constexpr double inv_n = 1 / double(N);
	ax = (n_1 * inv_n) * ax + (x * inv_n) / powi_fast<double, n_1>(ax);
	ax = (n_1 * inv_n) * ax + (x * inv_n) / powi_fast<double, n_1>(ax);
	ax = (n_1 * inv_n) * ax + (x * inv_n) / powi_fast<double, n_1>(ax);
	ax = (n_1 * inv_n) * ax + (x * inv_n) / powi_fast<double, n_1>(ax);

	return ax;
}

} // namespace detail


// =====================================================================================================================
/// huge, big float type
//  PREC -> maximal precision.
//  FLAGS -> default flags.
// =====================================================================================================================
template<uint64_t PREC = 8192, uint32_t FLAGS = BF_RNDN>
class bf
{
public:
	bf_t data;
public:
	constexpr static auto precision = PREC;
public:
	constexpr bf() noexcept { bf_init(&g_bf_ctx_default, &data); }
	
	bf(double val) { bf_init(&g_bf_ctx_default, &data); bf_set_float64(&data, val); }

	bf(int64_t val) { bf_init(&g_bf_ctx_default, &data); bf_set_si(&data, val); }
	bf(int32_t val) { bf_init(&g_bf_ctx_default, &data); bf_set_si(&data, val); }
	
	bf(uint64_t val) { bf_init(&g_bf_ctx_default, &data); bf_set_ui(&data, val); }
	bf(uint32_t val) { bf_init(&g_bf_ctx_default, &data); bf_set_ui(&data, val); }

	bf(limb_t limb0, slimb_t exp) { bf_init(&g_bf_ctx_default, &data); bf_resize(&data, 1); data.expn = exp; data.tab[0] = limb0; }

	bf(std::initializer_list<limb_t> il, slimb_t exp) { 
		bf_init(&g_bf_ctx_default, &data); 
		bf_resize(&data, il.size());
		data.expn = exp;
		int32_t i = il.size()-1;
		for(const auto limb: il){ data.tab[i--] = limb; }
	}

	// Create a real number from an ASCII representation, radix = 0 is for automatic detection, 10 for decimal and 16 for hexadecimal.
	bf(char const* str, int32_t radix = 0) {
		bf_init(&g_bf_ctx_default, &data);
		if (str!=nullptr) {	bf_atof(&data, str, nullptr, radix, PREC/*BF_PREC_INF*/, FLAGS); }
	}

	bf(const bf& other) { bf_init(&g_bf_ctx_default, &data); bf_set(&data,&other.data); }

	template<uint64_t OPREC, uint32_t OFLAGS>
	explicit bf(const bf<OPREC, OFLAGS>& other) { bf_init(&g_bf_ctx_default, &data); bf_set(&data, &other.data); }

	bf& operator =(const bf& other) {
		if (this == &other) { return *this; }
		bf_delete(&data);
		bf_set(&data,&other.data);
		return *this;
	}

	template<uint64_t OPREC, uint32_t OFLAGS>
	bf& operator =(const bf<OPREC, OFLAGS>& other) {
		if constexpr (PREC==OPREC && FLAGS==OFLAGS) {
			if (this == &other) { return *this; }
		}
		bf_delete(&data);
		bf_set(&data, &other.data);
		return *this;
	}

	constexpr bf(bf&& other) noexcept {
		data = std::move(other.data);
		other.data.len = 0;
		other.data.tab = nullptr;
	}

	template<uint64_t OPREC, uint32_t OFLAGS>
	explicit constexpr bf(bf<OPREC, OFLAGS>&& other) noexcept {
		data = std::move(other.data);
		other.data.len = 0;
		other.data.tab = nullptr;
	}

	bf& operator =(bf&& other) noexcept {
		if (this == &other) { return *this; }
		bf_delete(&data);
		data = std::move(other.data);
		other.data.len = 0;
		other.data.tab = nullptr;
		return *this;
	}
	
	~bf() noexcept { bf_delete(&data); detail::g_bf_ctx_cleanup.force_instantiate(); }


	explicit operator double() const {
		double res;	bf_get_float64(&data, &res, bf_rnd_t(FLAGS));	return res;
	}

	explicit constexpr operator int32_t() const {
		int32_t res; bf_get_int32(&res, &data, 0);	return res;
	}

	explicit constexpr operator int64_t() const {
		int64_t res; bf_get_int64(&res, &data, 0);	return res;
	}

	friend constexpr void convert_to(double& res, const bf& x)  { bf_get_float64(&res, &x.data, 0);	}
	friend constexpr void convert_to(int32_t& res, const bf& x) { bf_get_int64(&res, &x.data, 0);	}
	friend constexpr void convert_to(int64_t &res, const bf& x) { bf_get_int64(&res, &x.data, 0);	}

	// get_sign return +1 if a > 0, 0 if a = 0, and -1 if a < 0. 
	friend constexpr int get_sign(const bf& x) noexcept { return bf_sgn(&x.data); }
	friend constexpr bool is_zero(const bf& x) noexcept { return bf_is_zero(&x.data); }

	constexpr bool is_negative() const noexcept { return data.sign != 0; }

	constexpr bool is_zero() const noexcept { return bf_is_zero(&data); } // expn == BF_EXP_ZERO
	constexpr bool is_one() const noexcept { return bf_is_one(&data); }

	constexpr bool is_finite() const noexcept { return bf_is_finite(&data); } // expn < BF_EXP_INF
	constexpr bool is_nan() const noexcept { return bf_is_nan(&data); } // expn == BF_EXP_NAN
	constexpr bool is_inf() const noexcept { return !bf_is_finite(&data); } // expn == BF_EXP_INF

	constexpr bool is_even() const noexcept { return bf_is_even(&data); }
	constexpr bool is_odd() const noexcept { return bf_is_odd(&data); }

	constexpr operator bool() const noexcept { return !is_zero() && !is_nan();	}

	constexpr int64_t exponent() const noexcept { return data.expn; }
	constexpr void set_exponent(slimb_t new_exp) noexcept { data.expn = new_exp; }
	constexpr auto size() const noexcept { return data.len; }

	void negate() {
		if (data.expn != BF_EXP_NAN) { data.sign ^= 1; }
	}

	// bf operator +() const { return *this; } 
	bf operator - () const noexcept {
		auto ret = *this;
		ret.data.sign ^= 1;
		return ret;
	}

	friend constexpr uint64_t precision(const bf& a) { return PREC; }

	constexpr void swap(bf& rhs) { std::swap(data, rhs.data); }
	friend constexpr void swap(bf& lhs, bf& rhs) { std::swap(lhs.data, rhs.data); }

	/// add +
	friend void add(bf& r, const bf& lhs, const bf& rhs, limb_t prec = PREC) { bf_add(&r.data, &lhs.data, &rhs.data, prec, FLAGS); }
	friend void add(bf& r, const bf& lhs, int64_t rhs, limb_t prec = PREC) { bf_add_si(&r.data, &lhs.data, rhs, prec, FLAGS); }

	friend bf operator + (const bf& lhs, const bf& rhs) {
		bf r; bf_add(&r.data, &lhs.data, &rhs.data, PREC, FLAGS);	return r;
	}
	friend bf operator + (const bf& lhs, int64_t rhs) {
		bf r; bf_add_si(&r.data, &lhs.data, rhs, PREC, FLAGS);		return r;
	}
	friend bf operator + (int64_t lhs, const bf& rhs) {
		bf r; bf_add_si(&r.data, &rhs.data, lhs, PREC, FLAGS);		return r;
	}
	friend bf operator + (const bf& lhs, double rhs) {
		bt_static tmp; bf_static_init(&g_bf_ctx_default, &tmp);
		bf_set_float64(&tmp, rhs);
		bf r; bf_add(&r.data, &lhs.data, &tmp, PREC, FLAGS);	return r;
	}

	/// sub -
	friend void sub(bf& r, const bf& lhs, const bf& rhs, limb_t prec = PREC) { bf_sub(&r.data, &lhs.data, &rhs.data, prec, FLAGS); }

	friend bf operator - (const bf& lhs, const bf& rhs) {
		bf r; bf_sub(&r.data, &lhs.data, &rhs.data, PREC, FLAGS);	return r;
	}
	friend bf operator - (const bf& lhs, double rhs) {
		bt_static tmp; bf_static_init(&g_bf_ctx_default, &tmp);
		bf_set_float64(&tmp, rhs);
		bf r; bf_sub(&r.data, &lhs.data, &tmp, PREC, FLAGS);	return r;
	}

	/// div /
	friend void div(bf& r, const bf& lhs, const bf& rhs, limb_t prec = PREC) { bf_div(&r.data, &lhs.data, &rhs.data, prec, FLAGS); }
	friend void div(bf& r, const bf& lhs, int64_t rhs, limb_t prec = PREC) { bf_div_si(&r.data, &lhs.data, rhs, prec, FLAGS); }
	friend void div(bf& r_lhs, int64_t rhs, limb_t prec = PREC) {
		if (is_pow2(rhs)) { bf_div_2exp(&r_lhs.data, ceil_log2(rhs), prec, FLAGS); }  
		else { bf_div_si(&r_lhs.data, &r_lhs.data, rhs, prec, FLAGS); }
	}

	friend bf operator / (const bf& lhs, const bf& rhs) {
		bf r; bf_div(&r.data, &lhs.data, &rhs.data, PREC, FLAGS);	return r;
	}
	friend bf operator / (int64_t lhs, const bf& rhs) {
		bt_static tmp; bf_static_init(&g_bf_ctx_default, &tmp);
		bf_set_float64(&tmp, lhs);
		bf r; bf_div(&r.data, &tmp, &rhs.data, PREC, FLAGS);	return r;
	}
	friend bf operator / (const bf& lhs, int64_t rhs) {
		bt_static tmp; bf_static_init(&g_bf_ctx_default, &tmp);
		bf_set_float64(&tmp, rhs);
		bf r; bf_div(&r.data, &lhs.data, &tmp, PREC, FLAGS);	return r;
	}

	//TODO: modulo, remainder?  mod, rem !
	friend bf operator % (const bf& lhs, const bf& rhs) {
		bf r; bf_mod(&r.data, &lhs.data, &rhs.data); return r;
	}

	// integer devision, do it only work for PREC = BF_PREC_INF?
	friend bf div(const bf& lhs, const bf& rhs) {
		bf r; bf t; bf_divrem(&r.data, &t.data, &lhs.data, &rhs.data, PREC, FLAGS, BF_DIVREM_EUCLIDIAN);	return r;
	}
	

	/// mul *
	friend void mul(bf& r, const bf& lhs, const bf& rhs, limb_t prec = PREC) { bf_mul(&r.data, &lhs.data, &rhs.data, prec, FLAGS); }
	friend void mul(bf& r, const bf& lhs, int64_t rhs, limb_t prec = PREC) { bf_mul_si(&r.data, &lhs.data, rhs, prec, FLAGS); }
	friend void mul(bf& r_lhs, int64_t rhs, limb_t prec = PREC) {
		if (is_pow2(rhs)) { bf_mul_2exp(&r_lhs.data, ceil_log2(rhs), prec, FLAGS); }  
		else { bf_mul_si(&r_lhs.data, &r_lhs.data, rhs, prec, FLAGS); }
	}

	friend bf operator * (const bf& lhs, const bf& rhs) {
		bf r; bf_mul(&r.data, &lhs.data, &rhs.data, PREC, FLAGS);	return r;
	}
	friend bf operator * (const bf& lhs, uint64_t rhs) {
		bf r; bf_mul_ui(&r.data, &lhs.data, rhs, PREC, FLAGS);		return r;
	}
	friend bf operator * (const bf& lhs, int64_t rhs) {
		bf r; bf_mul_si(&r.data, &lhs.data, rhs, PREC, FLAGS);	return r;
	}
	friend bf operator * (int64_t lhs, const bf& rhs) {
		bf r; bf_mul_si(&r.data, &rhs.data, lhs, PREC, FLAGS);	return r;
	}

	friend bf operator * (const bf& lhs, int32_t rhs) {
		bf r; bf_mul_si(&r.data, &lhs.data, rhs, PREC, FLAGS);	return r;
	}
	friend bf operator * (const bf& lhs, double rhs) {
		bt_static tmp; bf_static_init(&g_bf_ctx_default, &tmp);
		bf_set_float64(&tmp, rhs);
		bf r; bf_mul(&r.data, &lhs.data, &tmp, PREC, FLAGS);	return r;
	}


	bf const& operator += (const bf& x) {
		bf_add(&data, &data, &x.data, PREC, FLAGS); return *this;
	}
	bf const& operator += (double d) {
		bt_static tmp; bf_static_init(&g_bf_ctx_default, &tmp);
		bf_set_float64(&tmp, d);
		bf_add(&data, &data, &tmp, PREC, FLAGS); return *this;
	}
	bf const& operator += (const int64_t x) {
		bf_add_si(&data, &data, x, PREC, FLAGS); return *this;
	}
	bf const& operator += (const int32_t x) {
		bf_add_si(&data, &data, x, PREC, FLAGS); return *this;
	}

	bf const& operator -= (const bf& x) {
		bf_sub(&data, &data, &x.data, PREC, FLAGS); return *this;
	}
	bf const& operator -= (double d) {
		bt_static tmp; bf_static_init(&g_bf_ctx_default, &tmp);
		bf_set_float64(&tmp, d);
		bf_sub(&data, &data, &tmp, PREC, FLAGS); return *this;
	}
	bf const& operator -= (const int64_t x) {
		bf_add_si(&data, &data, -x, PREC, FLAGS); return *this;
	}
	bf const& operator -= (const int32_t x) {
		bf_add_si(&data, &data, -x, PREC, FLAGS); return *this;
	}

	bf const& operator *= (const bf& x) {
		bf_mul(&data, &data, &x.data, PREC, FLAGS); return *this;
	}
	bf const& operator *= (double d) {
		bt_static tmp; bf_static_init(&g_bf_ctx_default, &tmp);
		bf_set_float64(&tmp, d);
		bf_mul(&data, &data, &tmp, PREC, FLAGS); return *this;
	}
	bf const& operator *= (const int64_t x) {
		if (is_pow2(x)) { bf_mul_2exp(&data, ceil_log2(x), PREC, FLAGS); }
		else { bf_mul_si(&data, &data, x, PREC, FLAGS); }
		return *this;
	}
	bf const& operator *= (const int32_t x) {
		if (is_pow2(x)) { bf_mul_2exp(&data, ceil_log2(x), PREC, FLAGS); }
		else { bf_mul_si(&data, &data, x, PREC, FLAGS); }
		return *this;
	}

	bf const& operator /= (const bf& x) {
		bf_div(&data, &data, &x.data, PREC, FLAGS); return *this;
	}
	bf const& operator /= (const int64_t x) {
		if (is_pow2(x)) { bf_div_2exp(&data, ceil_log2(x), PREC, FLAGS); }
		else { bf_div_si(&data, &data, x, PREC, FLAGS); }
		return *this;
	}
	bf const& operator /= (const int32_t x) {
		if (is_pow2(x)) { bf_div_2exp(&data, ceil_log2(x), PREC, FLAGS); }
		else { bf_div_si(&data, &data, x, PREC, FLAGS); }
		return *this;
	}

	// left-shift multiplies * this by 2^bits
	bf const& operator <<= (int64_t bits) {
		bf_mul_2exp(&data, bits, PREC, FLAGS);
		return *this;
	}

	// right-shift divides / this by 2^bits
	bf const& operator >>= (int64_t bits) {
		bf_div_2exp(&data, bits, PREC, FLAGS);
		return *this;
	}
#if 0
	// multiplies by 2^bits
	friend bf operator << (const bf& lhs,int64_t bits) {
		bf r{lhs}; bf_mul_2exp(&r.data, bits, PREC, FLAGS); return r;
	}
	// divides by 2^bits
	friend bf operator >> (const bf& lhs,int64_t bits) {
		bf r{lhs}; bf_div_2exp(&r.data, bits, PREC, FLAGS); return r;
	}
#endif
	friend bf operator | (const bf& lhs, const bf& rhs) {
		bf r; bf_logic_or(&r.data, &lhs.data, &rhs.data);	return r;
	}
	friend bf operator ^ (const bf& lhs, const bf& rhs) {
		bf r; bf_logic_xor(&r.data, &lhs.data, &rhs.data);	return r;
	}
	friend bf operator & (const bf& lhs, const bf& rhs) {
		bf r; bf_logic_and(&r.data, &lhs.data, &rhs.data);	return r;
	}

	friend constexpr bool operator == (const bf& lhs, const bf& rhs) noexcept { return bf_cmp(&lhs.data, &rhs.data) == 0; }
	friend constexpr bool operator <= (const bf& lhs, const bf& rhs) noexcept { return bf_cmp(&lhs.data, &rhs.data) <= 0; }
	friend constexpr bool operator <  (const bf& lhs, const bf& rhs) noexcept { return bf_cmp(&lhs.data, &rhs.data) < 0; }
	friend constexpr bool operator >= (const bf& lhs, const bf& rhs) noexcept { return bf_cmp(&lhs.data, &rhs.data) >= 0; }
	friend constexpr bool operator >  (const bf& lhs, const bf& rhs) noexcept { return bf_cmp(&lhs.data, &rhs.data) > 0; }
	friend constexpr bool operator != (const bf& lhs, const bf& rhs) noexcept { return bf_cmp(&lhs.data, &rhs.data) != 0; }

	friend constexpr bool operator == (const bf& lhs, int64_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) == 0; }
	friend constexpr bool operator <= (const bf& lhs, int64_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) <= 0; }
	friend constexpr bool operator <  (const bf& lhs, int64_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) < 0; }
	friend constexpr bool operator >= (const bf& lhs, int64_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) >= 0; }
	friend constexpr bool operator >  (const bf& lhs, int64_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) > 0; }
	friend constexpr bool operator != (const bf& lhs, int64_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) != 0; }

	friend constexpr bool operator == (const bf& lhs, int32_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) == 0; }
	friend constexpr bool operator <= (const bf& lhs, int32_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) <= 0; }
	friend constexpr bool operator <  (const bf& lhs, int32_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) < 0; }
	friend constexpr bool operator >= (const bf& lhs, int32_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) >= 0; }
	friend constexpr bool operator >  (const bf& lhs, int32_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) > 0; }
	friend constexpr bool operator != (const bf& lhs, int32_t rhs) noexcept { return bf_cmp_si(&lhs.data, rhs) != 0; }

	// decomposes x such as  x = ret * 2^exp
	friend bf frexp(const bf& x, slimb_t* exp) { //TODO: possible wrong rounding !
		if (!x) { *exp = 0;	return x; }
		*exp = x.data.expn;
		auto ret = x;
		ret.data.expn = 0;
		return ret;
	}

	// ret = x * 2^exp, reverse of frexp(ret, exp)
	friend bf ldexp(const bf& x, slimb_t exp) { //TODO: possible wrong rounding !
		auto ret = x;
		if (ret) ret.data.expn += exp;
		return ret;
	}

	friend void ldexp(bf& result, slimb_t e) {
		if (e > 0) {
			bf_mul_2exp(&result.data,  e, PREC, BF_RNDN);
		}
		else if (e < 0) {
			bf_div_2exp(&result.data, -e, PREC, BF_RNDN);
		}
	}

	bf& abs() noexcept { data.sign = 0; return *this; }

	//TODO: not classical round, find better name!
	bf& round() {
		bf_round(&data, PREC, FLAGS);
		return *this;
	}

	// round to integer.
	bf& rint() {
		bf_rint(&data, FLAGS);
		return *this;
	}

	// this = floor(this)
	bf& floor() {
		if (data.sign) bf_rint(&data, BF_RNDA);
		else bf_rint(&data, BF_RNDZ);
		return *this;
	}

	// this = ceil(this)
	bf& ceil() {
		if (data.sign) bf_rint(&data, BF_RNDZ);
		else bf_rint(&data, BF_RNDA);
		return *this;
	}

	// this = trunc(this)
	bf& trunc() {
		bf_rint(&data, BF_RNDZ);
		return *this;
	}

	// set this to pi ~ 3.1415926535897...
	bf& set_pi()    { bf_const_pi(&data, PREC, FLAGS); return *this; }

	// set this to log(2) ~ 0.6931471805599453...
	bf& set_log2()  { bf_const_log2(&data, PREC, FLAGS); return *this; }

	// friend bf calc_pi() { bf r; bf_const_pi(&r.data, PREC, FLAGS); return r; }
	// friend bf calc_log2() { bf r; bf_const_log2(&r.data, PREC, FLAGS); return r; }

	friend constexpr bf min(const bf& a, const bf& b) { return (a < b) ? a : b; }
	friend constexpr bf max(const bf& a, const bf& b) { return (a > b) ? a : b; }

	friend bf sqrt(const bf& v) {
		bf r; bf_sqrt(&r.data, &v.data, PREC, FLAGS);	return r;
	}
	friend bf isqrt(const bf& v) {
		bf r; bf_sqrtrem(&r.data, nullptr, &v.data);	return r;
	}
	
	friend bf pow(const bf& x, const bf& y) {
		bf r; bf_pow(&r.data, &x.data, &y.data, PREC, FLAGS);	return r;
	}
	friend bf ipow(const bf& x, uint64_t y) {
		bf r; bf_pow_ui(&r.data, &x.data, y, PREC, FLAGS);	return r;
	}

	friend bf exp(const bf& v) {
		bf r; bf_exp(&r.data, &v.data, PREC, FLAGS);	return r;
	}

	friend bf log(const bf& v, limb_t prec = PREC) {
		bf r; bf_log(&r.data, &v.data, prec, FLAGS);	return r;
	}

	/// Trigonometric functions

	friend bf cos(const bf& v) {
		bf r; bf_cos(&r.data, &v.data, PREC, FLAGS);	return r;
	}
	friend bf sin(const bf& v) {
		bf r; bf_sin(&r.data, &v.data, PREC, FLAGS);	return r;
	}
	friend bf tan(const bf& v) {
		bf r; bf_tan(&r.data, &v.data, PREC, FLAGS);	return r;
	}

	friend bf atan(const bf& v) {
		bf r; bf_atan(&r.data, &v.data, PREC, FLAGS);	return r;
	}
	// atan(y/x)
	friend bf atan2(const bf& y, const bf& x) {
		bf r; bf_atan2(&r.data, &y.data, &x.data, PREC, FLAGS);	return r;
	}
	friend bf asin(const bf& v) {
		bf r; bf_asin(&r.data, &v.data, PREC, FLAGS);	return r;
	}
	friend bf acos(const bf& v) {
		bf r; bf_acos(&r.data, &v.data, PREC, FLAGS);	return r;
	}

	/// t0 to the integer power N.
	template <int32_t N>
	friend bf ipow(const bf& t0, limb_t prec) {
		bf a = t0;
		constexpr bool recip = N < 0;
		int32_t b = N;
		bf r{ 1 };
		while (1) {
			if (b & 1) { mul(r, r, a, prec); } //  r *= a;
			b /= 2;
			if (b == 0) { break; }
			mul(a, a, a, prec); // a *= a;
		}
		if (recip) div(r, bf{1}, r, prec); // r = 1/r
		return r;
	}

	/// Slower than sqrt and cbrt !
	template <int32_t N>
	friend bf nth_root(const bf& a_in, bool test = false) { 
		static_assert(N > 1);
		if (a_in.is_zero() || a_in.is_one() || a_in.is_inf() || a_in.is_nan()) { return a_in; }
		using namespace detail;

		bf a = a_in; // we need a copy to adjust range
		auto ep = a_in.exponent();
		ep -= ep % N;
		a.set_exponent(a.exponent() - ep);

		constexpr int32_t n_1 = N - 1;
		const bf nn = bf{ n_1 } / bf{ N };
		const bf an = a / bf{ N };

		//TODO: if N is even then we can not accept negative a !
		const double da = double(a);
		double ax;
		if (da < 0.0) { ax = (-nth_rootd_aprox<N>(-da)); }
		else { ax = nth_rootd_aprox<N>(da);	}

		if (test) printf(" ep %lli  %f  %f \n",ep, da, ax);
		bf x{ax};

		bf new_x;
		bf x_diff;

		limb_t prec = 128; // start precision 128	
		for (int32_t i = 0; i < bf::precision; ++i) {
			/// new_x = nn * x + an / powi_fast<bf, n_1>(x);
			mul(new_x, nn, x, prec);					// new_x = nn * x
			div(x_diff,an,ipow<n_1>(x, prec), prec);	// x_diff = an / x^n_1;
			add(new_x,new_x, x_diff, prec);				// new_x = new_x + x_diff
			sub(x_diff, x, new_x, prec);				// x_diff = x - new_x;

			if(test) printf("i %i, exp %lli <= %lli, len x: %llu\n",i,x_diff.exponent(), -slimb_t(bf::precision-2), x.size());
			if (x_diff.is_zero() || x_diff.exponent() <= -slimb_t(bf::precision-10)) { break; }
			x = new_x;

			// because of quadratic convergence we make precision twice as big for every next iteration
			prec *= 2; if (prec > bf::precision) { /*printf(" prec limit! %lli \n", prec);*/ prec = bf::precision; }
		}
		
		x.set_exponent(x.exponent() + ep/N);
		if (test) printf(" ep %lli \n",x.exponent());
		return x;
	}

}; // class bf

using big_int = cx::bf<BF_PREC_INF>;

template<uint64_t PREC, uint32_t FLAGS>
inline void dump(const char* str, const bf<PREC, FLAGS>& v);

//=--------------------------------------------------------------------------------------------------------------------
template<uint64_t PREC, uint32_t FLAGS>
bf<PREC, FLAGS> ldexp(const bf<PREC, FLAGS>& x, slimb_t exp) {
	auto res = x;
	if (res) { // is not zero ?
		res.m_exponent += exp;
	}
	return res;
}

//=--------------------------------------------------------------------------------------------------------------------
/// arithmetic geometric mean
template<uint64_t PREC, uint32_t FLAGS>
bf<PREC, FLAGS> agm(const bf<PREC, FLAGS>& a, const bf<PREC, FLAGS>& b) 
{
	if (a.is_inf() || a.is_nan()) { return a; } 
	if (b.is_inf() || b.is_nan()) { return b; } 
	using bf = bf<PREC, FLAGS>;

	if (a > b) {
		//TODO: swap(a, b);
	}

	const auto prec = precision(b);
	const auto p = prec + ceil_log2(prec) + 16;
	// printf(" prec: %lli, p: %lli \n",prec,p);

	bf vf, uf;
	bf t, w;

	bf s = (a * b);
	bf u = sqrt(s);
	bf v = (a + b); v /= 2;
	for (int32_t i = 1; i < 64; ++i) {
		t = v - u;
		const auto ep = v.exponent() - (t).exponent();
		//printf(" e: %lli, p: %lli, p/4: %lli \n",ep, p, p/4);
        if (ep > (p-2)) { break; }

		vf = (u + v); vf /= 2;
        if (ep < p/4) {
			uf = (u * v);
			u = sqrt(uf);
		} else {
			s = v - u;
			t = (s * s); t /= 16;
			w = t / vf;
			return vf - w;
		}
		v = vf;
	}
	return v;
}

//=--------------------------------------------------------------------------------------------------------------------
/// log_v2, faster as log() above
//TODO: rename into log() and rename above log() to log_v0()!
template<uint64_t PREC, uint32_t FLAGS>
bf<PREC, FLAGS> log_v2(const bf<PREC, FLAGS>& a, limb_t prec_in = PREC) {
	static_assert(PREC != BF_PREC_INF, "log()  do not work with big_int a.k. bf<BF_PREC_INF> !");
	if (a.is_inf() || a.is_nan()) { return a; } 
	if (a.is_one()) { return {0}; } // log(1) == 0
	if (a.is_zero()) { return {-1}; } //TODO: log(0) == -Inf

	/*constexpr*/const uint64_t PREC2 = prec_in + 64;// PREC + 64;
	using bf = bf<PREC, FLAGS>;

	bf c_pi; // = pi;		 // 3.1415926535897...
	bf_const_pi(&c_pi.data,PREC2, FLAGS);

	bf c_log_2; // = log(2); // 0.6931471805599...
	bf_const_log2(&c_log_2.data,PREC2, FLAGS);

	bf c_one = bf(1);
	bf c_four = bf(4);

	const auto ep_a = a.exponent(); // exponent of a
	const slimb_t prec = precision(a); // wrong !
	//auto p = prec + 2 * ceil_log2(prec) + 12;
	slimb_t p = prec + 16; // test only

	//printf(" ep_a %lli, prec %lli, p %lli \n",ep_a,prec,p);

	bf r, s, t1, t2, diff;

	for (int32_t i = 1; i < 64; ++i) {
		auto m = ((p + 3) / 2) - ep_a;
		// printf("  h %lli\n",p + 3 - ep_a);
		// printf("  m %lli\n",m);
		// printf("  pp %lli\n",(p + 3) / 2);
		s = a; s.set_exponent((p + 3) / 2);
		// println(s);
		// printf("______\n");
#if 1	//TOOD:
		div(t1, c_four, s, PREC2);		// t1 = 4 / s
		t2 = agm(c_one, t1); t2 *= 2;	// t2 = 2 * agm(1, 4 / s)
		div(t2, c_pi, t2, PREC2);		// t2 = pi / t2
		mul(t1, c_log_2, m, PREC2);		// t1 = ln(2) * m
		sub(t1, t2, t1, PREC2);			// t1 = t2 - t1
#else
		t1 = c_four / s;				// t1 = 4 / s
		t2 = agm(c_one, t1); t2 *= 2;	// t2 = 2 * agm(1, 4 / s)
		t2 = c_pi / t2;					// t2 = pi / t2
		t1 = c_log_2 * m;				// t1 = ln(2) * m
		t1 = t2 - t1;					// t1 = t2 - t1
#endif
		const auto cancel = bf_max(0, t2.exponent() - t1.exponent());
		p += cancel + ceil_log2(p);

		diff = t1 - r;
		r = t1;

		if (diff.is_zero()) { /*printf(" zero break\n");*/ break; }
		if (diff.exponent() < -slimb_t(prec - 16)) { /*printf(" prec break! %lli %lli\n",diff.exponent(),-(prec - 16));*/ break; }
	}

	return r;
}


//=--------------------------------------------------------------------------------------------------------------------
// cbrt_v5:  https://github.com/shibatch/sleef/wiki/Divisionless-iterative-approximation-method-of-cube-root
template<uint64_t PREC, uint32_t FLAGS>
bf<PREC, FLAGS> cbrt(const bf<PREC, FLAGS>& a_in, bool test = false)
{
	using bf = bf<PREC, FLAGS>;
	if (a_in.is_zero() || a_in.is_one() || a_in.is_inf() || a_in.is_nan()) { return a_in; }

	bf a = a_in; // we need a copy to adjust range
	/// reduced range to [1, 8)
	auto ep = a_in.exponent() - 1;
	ep -= ep % 3;
	a.set_exponent(a.exponent() - ep); // a = ldexp(a, -e);

	const double da = double(a); // squish big-float into a double
	bf x{ 1.0 / std::cbrt(da) }; // use double 1/cbrt to approximate initial x

	// some temporals
	bf x2, x4; 
	bf new_x;
	bf x_diff;

	limb_t prec = 128; // start precision 

	for (int32_t i = 0; i < 128; i++) { //TODO: remove hard limit of 128!
		/// new_x = (4*x - x^4 * a) / 3

		mul(x2, x, x, prec);   // x2 = x * x;
		new_x = x;
		mul(new_x, 4, prec);   // new_x *= 4;
		mul(x4, x2, x2, prec); // x4 = x2 * x2;
		mul(x2 ,x4, a, prec);  // x_diff = (new_x - x4 * a);
		sub(x_diff, new_x, x2, prec);
		new_x = x_diff;
		div(new_x, 3, prec);   // new_x /= 3;

		sub(x_diff, x, new_x, prec); // x_diff = x - new_x;

		if (test) printf("i %i, exp %lli <= %lli, len x: %llu\n", i, x_diff.exponent(), -slimb_t(bf::precision - 2), x.size()); // test
		if (x_diff.is_zero() || x_diff.exponent() <= -slimb_t(bf::precision - 2)) { break; } // break once we have reached maximum precision
		x = new_x;

		// because of quadratic convergence we make precision twice as big for every next iteration
		prec *= 2; if (prec > bf::precision) { prec = bf::precision; }
	}
	x4 = (a * x) * x; // convert reciprocal-cbrt back to cbrt, multiplication order is very important for final precision !
	x4.set_exponent(x4.exponent() + ep / 3); // return ldexp(...,ep/3)  extend range back 
	return x4;
}

//=--------------------------------------------------------------------------------------------------------------------
/// Newton Raphson (NR) reciprocal.  1/a
/// Note: actually this is slower as simple 1/a !
// x = x*(2-a*x)
template<uint64_t PREC, uint32_t FLAGS>
bf<PREC, FLAGS> rcp(const bf<PREC, FLAGS>& a_in, bool test = false) {
	using bf = bf<PREC, FLAGS>;
	if (a_in.is_zero() || a_in.is_one() || a_in.is_inf() || a_in.is_nan()) { return a_in; }

	bf a = a_in; // we need a copy to adjust range
	auto ep = a_in.exponent(); // -1;
	a.set_exponent(0);
	
	const double da = double(a); // squish big-float into a double
	if (test) printf(" ep %lli  %g \n",ep, da);

	bf x{1.0/da};				 // use double 1/a to approximate initial x
	//if(test) println(x);

	const bf two{2};
	bf new_x;
	bf x_diff;

	limb_t prec = 128; // start precision 
	for (int32_t i = 0; i < 128; i++) {
		/// new_x = (two - a * x) * x; /// x = x*(2-a*x)
		mul(new_x, a, x, prec);			// new_x = a * x
		sub(x_diff, two, new_x, prec);	// x_diff = 2 - new_x
		mul(new_x, x_diff, x, prec);	// new_x = x_diff * x

		sub(x_diff, x, new_x, prec);    // x_diff = x - new_x;
		if(test) printf("i %i, exp %lli <= %lli, len x: %llu\n",i,x_diff.exponent(), -slimb_t(bf::precision-2), x.size());
		if (x_diff.is_zero() || x_diff.exponent() <= -slimb_t(bf::precision - 2)) { break; }
		x = new_x;

		// because of quadratic convergence we make precision twice as big for every next iteration
		prec *= 2; if (prec > bf::precision) { /*printf(" prec limit! %lli \n", prec);*/ prec = bf::precision; }
	}
	x.set_exponent(1 - ep); // return ldexp(...,ep)  extend range back 
	if (test) printf(" ep %lli \n",x.exponent());
	return x;
}


//=--------------------------------------------------------------------------------------------------------------------
// Fibonacci using fast doubling algorithm  https://www.nayuki.io/page/fast-fibonacci-algorithms
// F(2*k+1) = F(k+1)^2 + F(k)^2
template<typename BigT = big_int >
BigT fibonacci(int32_t n) {
	if (n <= 0) return {0};
	BigT a = 0;
	BigT b = 1;
	BigT c, d, e;
	for (int32_t i = 31; i >= 0; i--) {
		d = a * (b*2 - a);
		e = (a*a) + (b*b);
		a = d;
		b = e;
		if ( ((uint32_t(n) >> i) & 1) != 0 ) {
			c = a + b;
			a = b;
			b = c;
		}
	}
	return a;
}


namespace detail {



template<typename BigT>
BigT oddProduct(int32_t m, int32_t len) {
	if (len < 24) {
		auto p = BigT(m);
		for (int32_t k = 2; k <= 2 * (len - 1); k += 2) {
			p *= (m - k);
		}
		return p;
	}
	int32_t hlen = len >> 1;
	return oddProduct<BigT>(m - 2 * hlen, len - hlen) * oddProduct<BigT>(m, hlen);
}

template<typename BigT>
std::array<BigT, 2> oddFactorial(int32_t n) {
	BigT oddFact, sqrOddFact;
	if (n < 69) {
			static const BigT smallOddFactorial[] = {
				{ { 0x8000000000000000 },1 },
				{ { 0x8000000000000000 },1 },
				{ { 0x8000000000000000 },1 },
				{ { 0xc000000000000000 },2 },
				{ { 0xc000000000000000 },2 },
				{ { 0xf000000000000000 },4 },
				{ { 0xb400000000000000 },6 },
				{ { 0x9d80000000000000 },9 },
				{ { 0x9d80000000000000 },9 },
				{ { 0xb130000000000000 },12 },
				{ { 0xdd7c000000000000 },14 },
				{ { 0x9845400000000000 },18 },
				{ { 0xe467e00000000000 },19 },
				{ { 0xb994660000000000 },23 },
				{ { 0xa261d94000000000 },26 },
				{ { 0x983bbbac00000000 },30 },
				{ { 0x983bbbac00000000 },30 },
				{ { 0xa1bf7766c0000000 },34 },
				{ { 0xb5f7665398000000 },37 },
				{ { 0xd815c98344800000 },41 },
				{ { 0x870d9df20ad00000 },44 },
				{ { 0xb141df4dae310000 },48 },
				{ { 0xf3ba930acf836000 },51 },
				{ { 0xaf2e19afc5266d00 },56 },
				{ { 0x83629343d3dcd1c0 },58 },
				{ { 0xcd4a0619fb0907bc },62 },

				{ { 0xa6cc24f51bf75648,0xc000000000000000 },66 },
				{ { 0x8cbc3f2ecf98b0cd,0x6200000000000000 },71 },
				{ { 0xf6496e91eb4b3567,0x6b80000000000000 },73 },
				{ { 0xdf328c343d3c2865,0xb96c000000000000 },78 },
				{ { 0xd13f6370f96865df,0x5dd5400000000000 },82 },
				{ { 0xcab56855719d22b0,0x62e6960000000000 },87 },
				{ { 0xcab56855719d22b0,0x62e6960000000000 },87 },
				{ { 0xd10b13981d2a0bc5,0xe5fdcab000000000 },92 },
				{ { 0xde1bc4d19efcac82,0x445da75b00000000 },96 },
				{ { 0xf2ee5f4545e45cae,0x7ac66f0b88000000 },101 },
				{ { 0x88a61596f7507422,0x250f9e767c800000 },105 },
				{ { 0x9e0008f68df50647,0x7ada0f38fff40000 },110 },
				{ { 0xbba00aa4c892f774,0xe1e2f213aff1c000 },114 },
				{ { 0xe4ab0cf8d4731d96,0x734c9707fe6ea200 },119 },
				{ { 0x8eeae81b84c7f27e,0x080fde64ff052540 },122 },
				{ { 0xb71cf96342202eb1,0x7a5454f166be97ba },127 }, // 41

				{ { 0xf056075246ca3d48,0xf08eaf7cd6da2724,0x2000000000000000 },131}, // 42
				{ { 0xa179cceb478fe12d,0x019fdde7e05a924c,0x4580000000000000 },137}, // 43
				{ { 0xde0779c38265d59d,0xe23bd11ed47c8928,0xdf90000000000000 },140}, // 44
				{ { 0x9c1d419d77af9a33,0x03120f09ad679070,0xbd31400000000000 },146}, // 45
				{ { 0xe06a0e525c0c6da9,0x5469f59de944dfa2,0x0ff6cc0000000000 },150}, // 46
				{ { 0xa4cde2847b992088,0x59fdd05ff74e943b,0x03b93dd000000000 },156}, // 47
				{ { 0xf734d3c6b965b0cc,0x86fcb88ff2f5de58,0x8595dcb800000000 },157}, // 48
				{ { 0xbd44722425f1db5c,0x97597d4e36043e3b,0xc646bcfce0000000 },163}, // 49
				{ { 0x93dd792c3da4f360,0x563de9e51a33509e,0xb2e743a58f000000 },168}, // 50
				{ { 0xeba8f91e823ee3e1,0x8972acc521c1c87c,0xed2093cfdbe80000 },173}, // 51
				{ { 0xbf794a68c9d31927,0x3fad2c602b6d72e5,0x80aa7818e2ac8000 },177}, // 52
				{ { 0x9e90719ec722d0d4,0x80bb68bfa3f6a326,0x0e8d2b749bb6da00 },183}, // 53
				{ { 0x85c9dfddf8056033,0x4c9e2061b25819a8,0x1c471caa636247f0 },188}, // 54

			    { { 0xe5f2f8c582493d58,0x2bafc7a7ea876c18,0xf09a3944dad0eba4,0x8000000000000000 },193}, // 55
			    { { 0xc93499acd20015ad,0x2639ceb2ed367e95,0xd286f21c3f76ce2f,0xf000000000000000 },196}, // 56
			    { { 0xb332d8ddeb08134e,0x360b7c175b4488bd,0x6f802fa12885cfa2,0xb1c0000000000000 },202}, // 57
			    { { 0xa26614891cff517e,0xe0fa68752ab61beb,0xad0c2b2a0cb9442b,0x7116000000000000 },207}, // 58
			    { { 0x95b61aee66bb5f20,0xf766d84c035fe1bd,0x438737cac3bacad8,0x0c40480000000000 },213}, // 59
			    { { 0x8c5ab93f804fa92e,0xe7f06ac74329e3a1,0x6f4ec44e177f1e2a,0x8b7c438000000000 },217}, // 60
			    { { 0x85c67890864bed40,0xb51125c5ec03ecf5,0xde17131a6e6528c0,0x8cf2705600000000 },223}, // 61
			    { { 0x819844cc02198dd6,0xaf689c97bca3cd8e,0x2f265a819af1ff7a,0x888adcd350000000 },228}, // 62
			    { { 0xff23c771a4224f3e,0xa955f44abb627caf,0xecd3822f290c6ef9,0x3cd162c005800000 },233}, // 63
			    { { 0xff23c771a4224f3e,0xa955f44abb627caf,0xecd3822f290c6ef9,0x3cd162c005800000 },233}, // 64
			    { { 0x81902b47b5596c3d,0xd1fda60df3280351,0x5643681bf2d8505a,0x90e2542582cb0000 },240}, // 65
			    { { 0x859caca1f304379f,0xc08d933e62c1436b,0xe0f5835cd26f12dd,0x656966c6aee15800 },245}, // 66
			    { { 0x8be004b98a686a3b,0x3d9436254f625294,0xef8105852c4c47bf,0xc62a5797ff13e820 },251}, // 67
			    { { 0x949e0505230ef0de,0xf16d7987a45877be,0x3e7915dd7f110c3b,0xc28cfd117f0526a2 },255}, // 68
			};

		/*
			// only for  n < 41 !
			static const BigT smallOddFactorial[] = {"0x0000000000000000000000000000001",
			"0x0000000000000000000000000000001", "0x0000000000000000000000000000001",
			"0x0000000000000000000000000000003", "0x0000000000000000000000000000003",
			"0x000000000000000000000000000000f", "0x000000000000000000000000000002d",
			"0x000000000000000000000000000013b", "0x000000000000000000000000000013b",
			"0x0000000000000000000000000000b13", "0x000000000000000000000000000375f",
			"0x0000000000000000000000000026115", "0x000000000000000000000000007233f",
			"0x00000000000000000000000005cca33", "0x0000000000000000000000002898765",
			"0x00000000000000000000000260eeeeb", "0x00000000000000000000000260eeeeb",
			"0x0000000000000000000000286fddd9b", "0x00000000000000000000016beecca73",
			"0x000000000000000000001b02b930689", "0x00000000000000000000870d9df20ad",
			"0x0000000000000000000b141df4dae31", "0x00000000000000000079dd498567c1b",
			"0x00000000000000000af2e19afc5266d", "0x000000000000000020d8a4d0f4f7347",
			"0x000000000000000335281867ec241ef", "0x0000000000000029b3093d46fdd5923",
			"0x0000000000000465e1f9767cc5866b1", "0x0000000000001ec92dd23d6966aced7",
			"0x0000000000037cca30d0f4f0a196e5b", "0x0000000000344fd8dc3e5a1977d7755",
			"0x000000000655ab42ab8ce915831734b", "0x000000000655ab42ab8ce915831734b",
			"0x00000000d10b13981d2a0bc5e5fdcab", "0x0000000de1bc4d19efcac82445da75b",
			"0x000001e5dcbe8a8bc8b95cf58cde171", "0x00001114c2b2deea0e8444a1f3cecf9",
			"0x0002780023da37d4191deb683ce3ffd", "0x002ee802a93224bddd3878bc84ebfc7",
			"0x07255867c6a398ecb39a64b83ff3751", "0x23baba06e131fc9f8203f7993fc1495"};
		}*/

		oddFact = smallOddFactorial[n];
		sqrOddFact = smallOddFactorial[(n / 2)];
	}
	else {
		const auto ofr = oddFactorial<BigT>(n / 2);
		sqrOddFact = ofr[0];
		auto oldOddFact = ofr[1];

		int32_t len = (n - 1) / 4;
		(n % 4) != 2 && (len += 1);
		int32_t high = n - ((n + 1) & 1);
		const auto oddSwing = div(oddProduct<BigT>(high, len), oldOddFact);
		oddFact = (sqrOddFact * sqrOddFact) * oddSwing;
	}
	return std::array<BigT, 2>{oddFact, sqrOddFact};
}

} // namespace detail

//=--------------------------------------------------------------------------------------------------------------------
// Return the factorial of ``n``. Implementation of the swing algorithm using no
// primes. An advanced version based on prime-factorization which is much faster
// is available as the prime-swing factorial. However the claim is that this is
// the fastest algorithm not using prime-factorization. It has the same recursive
// structure as his big brother.
template<typename BigT = big_int >
BigT swing_factorial(int32_t n) {
	if (n <= 1) { return 1; }
	const int32_t sh = n - std::popcount(uint32_t(n));
	auto ofr = detail::oddFactorial<BigT>(n);
	ofr[0] <<= sh;
	return ofr[0];
}

//=--------------------------------------------------------------------------------------------------------------------
/// Not working for negative integers! Ignores sign.
template<uint64_t PREC, uint32_t FLAGS>
constexpr int32_t count_ones(const bf<PREC, FLAGS>& a) {
	const slimb_t n = a.size();
	const limb_t* tab = a.data.tab;
	int32_t cnt = 0;
	for (auto i = n - 1; i >= 0; i--) {
		cnt += std::popcount(tab[i]);
	}
	return cnt;
}

//=--------------------------------------------------------------------------------------------------------------------
/// Number of zeros trailing the binary representation of `x`.
/// trailing_zeros(2) == 1
template<uint64_t PREC, uint32_t FLAGS>
constexpr int32_t trailing_zeros(const bf<PREC, FLAGS>& a) {
	const slimb_t n = a.size();
	const limb_t* tab = a.data.tab;
	if (n > 0) {
		std::countr_zero(tab[0] >> (LIMB_BITS - a.data.expn));
	}
	/*int32_t cnt = 0;
	for (auto i = n - 1; i >= 0; i--) {
		cnt += std::countr_zero(tab[i]);
		printf(" i %i, cnt %i \n ",i,cnt);
	}
	return cnt;*/
	return 0;
}

//=--------------------------------------------------------------------------------------------------------------------
/// println for single bf number.
template<uint64_t PREC, uint32_t FLAGS>
inline void println(const char* str,const bf<PREC, FLAGS>& v, int ndigits = -1, bf_flags_t flags = BF_FTOA_FORMAT_FIXED) {
	if (ndigits < 0) { 
		if (PREC != BF_PREC_INF) { ndigits = (int)ceil(PREC / BITS_PER_DIGIT) - 1; }
		else {
			ndigits = BF_PREC_INF; flags = BF_FTOA_FORMAT_FRAC /*| BF_FTOA_FORMAT_FIXED*/;
		}
	}
	size_t digits_len = 0;
	char* digits = nullptr;
	digits = bf_ftoa(&digits_len, &v.data, 10, ndigits, flags | BF_RNDZ);
	printf("%s%.*s \n",str, (int32_t)digits_len, digits);
	free(digits);
}
//=--------------------------------------------------------------------------------------------------------------------
/// println for single bf number.
template<uint64_t PREC, uint32_t FLAGS>
inline void println(const bf<PREC, FLAGS>& v, int ndigits = -1, bf_flags_t flags = BF_FTOA_FORMAT_FIXED) {
	println("",v,ndigits,flags);
}


//=--------------------------------------------------------------------------------------------------------------------
/// dump single bf number.
template<uint64_t PREC, uint32_t FLAGS>
inline void dump(const char *str, const bf<PREC, FLAGS>& v) {
	printf("%s {{", str);
	const slimb_t n = v.size();
	const limb_t* tab = v.data.tab;
	for (auto i = n - 1; i >= 0; i--) {
		if (i != (n - 1)) printf(",");
		printf("0x" FMT_LIMB, tab[i]);
	}
	printf("},%lli", v.exponent());
	printf("}\n");
}

static_assert(sizeof(bf<>) <= 48);


inline namespace literals {

inline auto operator""_bigf1k(const char* str) { return bf<1024>{str}; }
inline auto operator""_bigf2k(const char* str) { return bf<2048>{str}; }
inline auto operator""_bigf4k(const char* str) { return bf<4096>{str}; }
inline auto operator""_bigf8k(const char* str) { return bf<8192>{str}; }


inline auto operator""_bigi(const char* str) { return bf<BF_PREC_INF>{str}; }
// inline auto operator""_bigi(unsigned long long int n) { return bf<BF_PREC_INF>{n}; }


} // inline namespace literals


}; // namespace


#endif // HUGEFLOAT_RE_HPP_