/*******************************************************************
*
*    Author: Kareem Omar
*    kareem.h.omar@gmail.com
*    https://github.com/komrad36
*
*    Last updated Feb 15, 2021
* 
*    uint128_t replacement !
* 
*******************************************************************/

#pragma once

#include <cstdint>
#include <intrin.h>
#include <iosfwd>

using I8  = int8_t;
using I16 = int16_t;
using I32 = int32_t;
using I64 = int64_t;

using U8  = uint8_t;
using U16 = uint16_t;
using U32 = uint32_t;
using U64 = uint64_t;

#define MAKE_BINARY_OP_HELPERS(op) \
friend auto operator op(const U128& x, U8   y) { return operator op(x, (U128)y); }  \
friend auto operator op(const U128& x, U16  y) { return operator op(x, (U128)y); }  \
friend auto operator op(const U128& x, U32  y) { return operator op(x, (U128)y); }  \
friend auto operator op(const U128& x, U64  y) { return operator op(x, (U128)y); }  \
friend auto operator op(const U128& x, I8   y) { return operator op(x, (U128)y); }  \
friend auto operator op(const U128& x, I16  y) { return operator op(x, (U128)y); }  \
friend auto operator op(const U128& x, I32  y) { return operator op(x, (U128)y); }  \
friend auto operator op(const U128& x, I64  y) { return operator op(x, (U128)y); }  \
friend auto operator op(const U128& x, char y) { return operator op(x, (U128)y); }  \
friend auto operator op(U8   x, const U128& y) { return operator op((U128)x, y); }  \
friend auto operator op(U16  x, const U128& y) { return operator op((U128)x, y); }  \
friend auto operator op(U32  x, const U128& y) { return operator op((U128)x, y); }  \
friend auto operator op(U64  x, const U128& y) { return operator op((U128)x, y); }  \
friend auto operator op(I8   x, const U128& y) { return operator op((U128)x, y); }  \
friend auto operator op(I16  x, const U128& y) { return operator op((U128)x, y); }  \
friend auto operator op(I32  x, const U128& y) { return operator op((U128)x, y); }  \
friend auto operator op(I64  x, const U128& y) { return operator op((U128)x, y); }  \
friend auto operator op(char x, const U128& y) { return operator op((U128)x, y); }

#define MAKE_BINARY_OP_HELPERS_FLOAT(op) \
friend auto operator op(const U128& x, float  y) { return (float)x op y; }   \
friend auto operator op(const U128& x, double y) { return (double)x op y; }  \
friend auto operator op(float  x, const U128& y) { return x op (float)y; }    \
friend auto operator op(double x, const U128& y) { return x op (double)y; }

#define MAKE_BINARY_OP_HELPERS_U64(op) \
friend U128 operator op(const U128& x, U8  n) { return operator op(x, (U64)n); }    \
friend U128 operator op(const U128& x, U16 n) { return operator op(x, (U64)n); }    \
friend U128 operator op(const U128& x, U32 n) { return operator op(x, (U64)n); }    \
friend U128 operator op(const U128& x, I8  n) { return operator op(x, (U64)n); }    \
friend U128 operator op(const U128& x, I16 n) { return operator op(x, (U64)n); }    \
friend U128 operator op(const U128& x, I32 n) { return operator op(x, (U64)n); }    \
friend U128 operator op(const U128& x, I64 n) { return operator op(x, (U64)n); }    \
friend U128 operator op(const U128& x, const U128& n) { return operator op(x, (U64)n); }

// =====================================================================================================================
class U128
{
public:
    friend U128 DivMod(U128 n, U128 d, U128& rem);

    U128() = default;
    U128(U8    x) : m_lo(x), m_hi(0) {}
    U128(U16   x) : m_lo(x), m_hi(0) {}
    U128(U32   x) : m_lo(x), m_hi(0) {}
    U128(U64   x) : m_lo(x), m_hi(0) {}
    U128(I8    x) : m_lo(I64(x)), m_hi(I64(x) >> 63) {}
    U128(I16   x) : m_lo(I64(x)), m_hi(I64(x) >> 63) {}
    U128(I32   x) : m_lo(I64(x)), m_hi(I64(x) >> 63) {}
    U128(I64   x) : m_lo(I64(x)), m_hi(I64(x) >> 63) {}
    U128(U64 hi, U64 lo) : m_lo(lo), m_hi(hi) {}

    // inexact values truncate, as per the Standard [conv.fpint]
    // passing values unrepresentable in the destination format is undefined behavior,
    // as per the Standard, but this implementation saturates
    U128(float x);

    // inexact values truncate, as per the Standard [conv.fpint]
    // passing values unrepresentable in the destination format is undefined behavior,
    // as per the Standard, but this implementation saturates
    U128(double x);

    U128& operator+=(const U128& x)
    {
        static_cast<void>(_addcarry_u64(_addcarry_u64(0, m_lo, x.m_lo, &m_lo), m_hi, x.m_hi, &m_hi));
        return *this;
    }

    friend U128 operator+(const U128& x, const U128& y)
    {
        U128 ret;
        static_cast<void>(_addcarry_u64(_addcarry_u64(0, x.m_lo, y.m_lo, &ret.m_lo), x.m_hi, y.m_hi, &ret.m_hi));
        return ret;
    }

    MAKE_BINARY_OP_HELPERS(+);
    MAKE_BINARY_OP_HELPERS_FLOAT(+);

    U128& operator-=(const U128& x)
    {
        static_cast<void>(_subborrow_u64(_subborrow_u64(0, m_lo, x.m_lo, &m_lo), m_hi, x.m_hi, &m_hi));
        return *this;
    }

    friend U128 operator-(const U128& x, const U128& y)
    {
        U128 ret;
        static_cast<void>(_subborrow_u64(_subborrow_u64(0, x.m_lo, y.m_lo, &ret.m_lo), x.m_hi, y.m_hi, &ret.m_hi));
        return ret;
    }

    MAKE_BINARY_OP_HELPERS(-);
    MAKE_BINARY_OP_HELPERS_FLOAT(-);

    U128& operator*=(const U128& x)
    {
        // ab * cd
        // ==
        // (2^64*a + b) * (2^64*c + d)
        // if a*c == e, a*d == f, b*c == g, b*d == h
        // |ee|ee|  |  |
        // |  |fg|fg|  |
        // |  |  |hh|hh|

        U64 hHi;
        const U64 hLo = _umul128(m_lo, x.m_lo, &hHi);
        m_hi = hHi + m_hi * x.m_lo + m_lo * x.m_hi;
        m_lo = hLo;
        return *this;
    }

    friend U128 operator*(const U128& x, const U128& y)
    {
        U128 ret;
        U64 hHi;
        ret.m_lo = _umul128(x.m_lo, y.m_lo, &hHi);
        ret.m_hi = hHi + y.m_hi * x.m_lo + y.m_lo * x.m_hi;
        return ret;
    }

    MAKE_BINARY_OP_HELPERS(*);
    MAKE_BINARY_OP_HELPERS_FLOAT(*);

    U128& operator/=(const U128& x)
    {
        U128 rem;
        *this = DivMod(*this, x, rem);
        return *this;
    }

    friend U128 operator/(const U128& x, const U128& y)
    {
        U128 rem;
        return DivMod(x, y, rem);
    }

    MAKE_BINARY_OP_HELPERS(/);
    MAKE_BINARY_OP_HELPERS_FLOAT(/);

    U128& operator%=(const U128& x)
    {
        static_cast<void>(DivMod(*this, x, *this));
        return *this;
    }

    friend U128 operator%(const U128& x, const U128& y)
    {
        U128 ret;
        static_cast<void>(DivMod(x, y, ret));
        return ret;
    }

    MAKE_BINARY_OP_HELPERS(%);

    U128& operator&=(const U128& x)
    {
        m_hi &= x.m_hi;
        m_lo &= x.m_lo;
        return *this;
    }

    friend U128 operator&(const U128& x, const U128& y)
    {
        return U128(x.m_hi & y.m_hi, x.m_lo & y.m_lo);
    }

    MAKE_BINARY_OP_HELPERS(&);

    U128& operator|=(const U128& x)
    {
        m_hi |= x.m_hi;
        m_lo |= x.m_lo;
        return *this;
    }

    friend U128 operator|(const U128& x, const U128& y)
    {
        return U128(x.m_hi | y.m_hi, x.m_lo | y.m_lo);
    }

    MAKE_BINARY_OP_HELPERS(|);

    U128& operator^=(const U128& x)
    {
        m_hi ^= x.m_hi;
        m_lo ^= x.m_lo;
        return *this;
    }

    friend U128 operator^(const U128& x, const U128& y)
    {
        return U128(x.m_hi ^ y.m_hi, x.m_lo ^ y.m_lo);
    }

    MAKE_BINARY_OP_HELPERS(^);

    U128& operator>>=(U64 n)
    {
        const U64 lo = __shiftright128(m_lo, m_hi, (U8)n);
        const U64 hi = m_hi >> (n & 63ULL);

        m_lo = n & 64 ? hi : lo;
        m_hi = n & 64 ? 0  : hi;

        return *this;
    }

    friend U128 operator>>(const U128& x, U64 n)
    {
        U128 ret;

        const U64 lo = __shiftright128(x.m_lo, x.m_hi, (U8)n);
        const U64 hi = x.m_hi >> (n & 63ULL);

        ret.m_lo = n & 64 ? hi : lo;
        ret.m_hi = n & 64 ? 0  : hi;

        return ret;
    }

    MAKE_BINARY_OP_HELPERS_U64(>>);

    U128& operator<<=(U64 n)
    {
        const U64 hi = __shiftleft128(m_lo, m_hi, (U8)n);
        const U64 lo = m_lo << (n & 63ULL);

        m_hi = n & 64 ? lo : hi;
        m_lo = n & 64 ? 0 : lo;

        return *this;
    }

    friend U128 operator<<(const U128& x, U64 n)
    {
        U128 ret;

        const U64 hi = __shiftleft128(x.m_lo, x.m_hi, (U8)n);
        const U64 lo = x.m_lo << (n & 63ULL);

        ret.m_hi = n & 64 ? lo : hi;
        ret.m_lo = n & 64 ? 0 : lo;

        return ret;
    }

    MAKE_BINARY_OP_HELPERS_U64(<<);

    friend U128 operator~(const U128& x)
    {
        return U128(~x.m_hi, ~x.m_lo);
    }

    friend U128 operator+(const U128& x)
    {
        return x;
    }

    friend U128 operator-(const U128& x)
    {
        U128 ret;
        static_cast<void>(_subborrow_u64(_subborrow_u64(0, 0, x.m_lo, &ret.m_lo), 0, x.m_hi, &ret.m_hi));
        return ret;
    }

    U128& operator++()
    {
        operator+=(1);
        return *this;
    }

    U128 operator++(int)
    {
        const U128 x = *this;
        operator++();
        return x;
    }

    U128& operator--()
    {
        operator-=(1);
        return *this;
    }

    U128 operator--(int)
    {
        const U128 x = *this;
        operator--();
        return x;
    }

    friend bool operator<(const U128& x, const U128& y)
    {
        U64 unusedLo, unusedHi;
        return _subborrow_u64(_subborrow_u64(0, x.m_lo, y.m_lo, &unusedLo), x.m_hi, y.m_hi, &unusedHi);
    }
    MAKE_BINARY_OP_HELPERS(<);
    MAKE_BINARY_OP_HELPERS_FLOAT(<);

    friend bool operator>(const U128& x, const U128& y) { return y < x; }
    MAKE_BINARY_OP_HELPERS(>);
    MAKE_BINARY_OP_HELPERS_FLOAT(>);

    friend bool operator<=(const U128& x, const U128& y) { return !(x > y); }
    MAKE_BINARY_OP_HELPERS(<=);
    MAKE_BINARY_OP_HELPERS_FLOAT(<=);

    friend bool operator>=(const U128& x, const U128& y) { return !(x < y); }
    MAKE_BINARY_OP_HELPERS(>=);
    MAKE_BINARY_OP_HELPERS_FLOAT(>=);

    friend bool operator==(const U128& x, const U128& y)
    {
        return !((x.m_hi ^ y.m_hi) | (x.m_lo ^ y.m_lo));
    }
    MAKE_BINARY_OP_HELPERS(==);
    MAKE_BINARY_OP_HELPERS_FLOAT(==);

    friend bool operator!=(const U128& x, const U128& y) { return !(x == y); }
    MAKE_BINARY_OP_HELPERS(!=);
    MAKE_BINARY_OP_HELPERS_FLOAT(!=);

    explicit operator bool() const { return m_hi | m_lo; }

    operator U8 () const { return (U8) m_lo; }
    operator U16() const { return (U16)m_lo; }
    operator U32() const { return (U32)m_lo; }
    operator U64() const { return (U64)m_lo; }

    operator I8 () const { return (I8) m_lo; }
    operator I16() const { return (I16)m_lo; }
    operator I32() const { return (I32)m_lo; }
    operator I64() const { return (I64)m_lo; }

    operator char() const { return (char)m_lo; }

    // rounding method is implementation-defined as per the Standard [conv.fpint]
    // this implementation performs IEEE 754-compliant "round half to even" rounding to nearest,
    // regardless of the current FPU rounding mode, which matches the behavior of clang and GCC
    operator float() const;

    // rounding method is implementation-defined as per the Standard [conv.fpint]
    // this implementation performs IEEE 754-compliant "round half to even" rounding to nearest,
    // regardless of the current FPU rounding mode, which matches the behavior of clang and GCC
    operator double() const;

    // caller is responsible for ensuring that buf has space for the U128 AND the null terminator
    // that follows, in the given output base.
    // Common bases and worst-case size requirements:
    // Base  2: 129 bytes (128 + null terminator)
    // Base  8:  44 bytes ( 43 + null terminator)
    // Base 10:  40 bytes ( 39 + null terminator)
    // Base 16:  33 bytes ( 32 + null terminator)
    void ToString(char* buf, U64 base = 10) const;

private:
    U64 m_lo;
    U64 m_hi;
};

#undef MAKE_BINARY_OP_HELPERS
#undef MAKE_BINARY_OP_HELPERS_FLOAT
#undef MAKE_BINARY_OP_HELPERS_U64

// std::ostream& operator<<(std::ostream& os, const U128& x);


// ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

static inline bool FitsHardwareDivL(U64 nHi, U64 nLo, U64 d)
{
	return !(nHi | (d >> 32)) && nLo < (d << 32);
}

static inline U64 HardwareDivL(U64 n, U64 d, U64& rem)
{
	U32 rLo;
	const U32 qLo = _udiv64(n, U32(d), &rLo);
	rem = rLo;
	return qLo;
}

static inline U64 HardwareDivQ(U64 nHi, U64 nLo, U64 d, U64& rem)
{
	nLo = _udiv128(nHi, nLo, d, &nHi);
	rem = nHi;
	return nLo;
}

static inline bool IsPow2(U64 hi, U64 lo)
{
	const U64 T = hi | lo;
	return !((hi & lo) | (T & (T - 1)));
}

static inline U64 CountTrailingZeros(U64 hi, U64 lo)
{
	const U64 nLo = _tzcnt_u64(lo);
	const U64 nHi = 64ULL + _tzcnt_u64(hi);
	return lo ? nLo : nHi;
}

static inline U64 CountLeadingZeros(U64 hi, U64 lo)
{
	const U64 nLo = 64ULL + _lzcnt_u64(lo);
	const U64 nHi = _lzcnt_u64(hi);
	return hi ? nHi : nLo;
}

static inline U128 MaskBitsBelow(U64 hi, U64 lo, U64 n)
{
	return U128(_bzhi_u64(hi, U32(n < 64 ? 0 : n - 64)), _bzhi_u64(lo, U32(n)));
}

inline U128 DivMod(U128 N, U128 D, U128& rem)
{
	if (D > N)
	{
		rem = N;
		return 0;
	}

	U64 nHi = N.m_hi;
	U64 nLo = N.m_lo;
	U64 dHi = D.m_hi;
	U64 dLo = D.m_lo;

	if (IsPow2(dHi, dLo))
	{
		const U64 n = CountTrailingZeros(dHi, dLo);
		rem = MaskBitsBelow(nHi, nLo, n);
		return N >> n;
	}

	if (!dHi)
	{
		if (nHi < dLo)
		{
			U64 remLo;
			U64 Q;
			if (FitsHardwareDivL(nHi, nLo, dLo))
				Q = HardwareDivL(nLo, dLo, remLo);
			else
				Q = HardwareDivQ(nHi, nLo, dLo, remLo);
			rem = remLo;
			return Q;
		}

		U64 remLo;
		const U64 qHi = HardwareDivQ(0, nHi, dLo, remLo);
		const U64 qLo = HardwareDivQ(remLo, nLo, dLo, remLo);
		rem = remLo;
		return U128(qHi, qLo);
	}

	U64 n = _lzcnt_u64(dHi) - _lzcnt_u64(nHi);

	dHi = __shiftleft128(dLo, dHi, U8(n));
	dLo <<= n;

	U64 Q = 0;
	++n;

	do
	{
		U64 tLo, tHi;
		unsigned char carry = _subborrow_u64(_subborrow_u64(0, nLo, dLo, &tLo), nHi, dHi, &tHi);
		nLo = !carry ? tLo : nLo;
		nHi = !carry ? tHi : nHi;
		Q = (Q << 1) + !carry;
		dLo = __shiftright128(dLo, dHi, 1);
		dHi >>= 1;
	} while (--n);

	rem = U128(nHi, nLo);
	return Q;
}

inline U128::U128(float x)
{
	const U32 bits = U32(_mm_cvtsi128_si32(_mm_castps_si128(_mm_set_ss(x))));
	const U32 s = bits >> 31;

	// technically UB but let's be nice
	if (s)
	{
		m_hi = m_lo = 0ULL;
		return;
	}

	const U32 e = (bits >> 23) - 127;
	const U32 m = (bits & ((1U << 23) - 1U)) | (1U << 23);

	// again, technically UB but let's be nice
	if (e >= 128)
	{
		m_hi = m_lo = ~0ULL;
		return;
	}

	if (e >= 23)
		*this = U128(m) << (e - 23);
	else
		*this = m >> (23 - e);
}

inline U128::U128(double x)
{
	const U64 bits = U64(_mm_cvtsi128_si64(_mm_castpd_si128(_mm_set_sd(x))));
	const U64 s = bits >> 63;

	// technically UB but let's be nice
	if (s)
	{
		m_hi = m_lo = 0ULL;
		return;
	}

	const U64 e = (bits >> 52) - 1023;
	const U64 m = (bits & ((1ULL << 52) - 1ULL)) | (1ULL << 52);

	// again, technically UB but let's be nice
	if (e >= 128)
	{
		m_hi = m_lo = ~0ULL;
		return;
	}

	if (e >= 52)
		*this = U128(m) << (e - 52);
	else
		*this = m >> (52 - e);
}

inline U128::operator float() const
{
	if (!*this)
		return 0.0f;

	const U32 numBits = 128U - U32(CountLeadingZeros(m_hi, m_lo));

	U32 bits;

	if (numBits <= 24)
	{
		const U32 m = (U32(m_lo) << (24 - numBits)) & ~(1U << 23);
		const U32 e = numBits + 126;
		bits = (e << 23) | m;
	}
	else
	{
		const U32 s = numBits - 24;
		const U32 m = U32(*this >> s) & ~(1U << 23);
		const U32 G = U32(*this >> (s - 1));
		const U32 R = U32(bool(MaskBitsBelow(m_hi, m_lo, s < 2 ? 0 : s - 2)));
		const U32 e = numBits + 126;
		bits = ((e << 23) | m) + (G & (R | m) & 1U);
	}

	return _mm_cvtss_f32(_mm_castsi128_ps(_mm_cvtsi32_si128((I32)bits)));
}

inline U128::operator double() const
{
	if (!*this)
		return 0.0;

	const U64 numBits = 128ULL - CountLeadingZeros(m_hi, m_lo);

	U64 bits;

	if (numBits <= 53)
	{
		const U64 m = (m_lo << (53 - numBits)) & ~(1ULL << 52);
		const U64 e = numBits + 1022;
		bits = (e << 52) | m;
	}
	else
	{
		const U64 s = numBits - 53;
		const U64 m = U64(*this >> s) & ~(1ULL << 52);
		const U64 G = U64(*this >> (s - 1));
		const U64 R = U64(bool(MaskBitsBelow(m_hi, m_lo, s < 2 ? 0 : s - 2)));
		const U64 e = numBits + 1022;
		bits = ((e << 52) | m) + (G & (R | m) & 1ULL);
	}

	return _mm_cvtsd_f64(_mm_castsi128_pd(_mm_cvtsi64_si128((I64)bits)));
}

inline void U128::ToString(char* buf, U64 base/* = 10*/) const
{
	U64 i = 0;
	if (base >= 2 && base <= 36)
	{
		U128 n = *this;
		U128 r, b = base;
		do
		{
			n = DivMod(n, b, r);
			const char c(r);
			buf[i++] = c + (c >= 10 ? '7' : '0');
		} while (n);

		for (U64 j = 0; j < (i >> 1); ++j)
		{
			const char t = buf[j];
			buf[j] = buf[i - j - 1];
			buf[i - j - 1] = t;
		}
	}
	buf[i] = '\0';
}

/*
inline std::ostream& operator<<(std::ostream& os, const U128& x)
{
	char buf[40];
	x.ToString(buf);
	os << buf;
	return os;
}*/

inline const char* NatVisStr_DebugOnly(const U128& x)
{
	static char buf[40];
	x.ToString(buf);
	return buf;
}

