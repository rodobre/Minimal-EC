#pragma once
#include <vector>
#include <tuple>
#include <cstdint>
#include <algorithm>
#include <ostream>

/// Define the maximum constants allowed by the
/// elliptic curve implementation
/// eq: a y^2 = b x^3 + c x^2 + d x + e
constexpr size_t EllipticCurveConstants = 5u;

/// Math namespace exceptions
enum class MathImplExceptions : uint8_t
{
	MODINV_PARAMS_NOT_COPRIME
};

/// Elliptic curve exceptions
enum class EllipticCurveExceptions : uint8_t
{
	CURVE_IS_SINGULAR,
	POINT_NOT_ON_CURVE
};

/// Binary helpers
namespace BinaryUtils
{
	/// Count the number of bits in an integer
	template <typename T>
	size_t CountBits(const T& num)
	{
		size_t count = 0;
		T n = num;
		T nullT = T(0);

		while (n != nullT)
		{
			++count;
			n >>= 1u;
		}

		return count;
	}

	/// Get bit at position 'i'
	template <typename T>
	T GetBit(const T& num, size_t i)
	{
		T oneT = T(1);
		return static_cast<T>((num & (oneT << i)) >> i);
	}
}

/// Finite field arithmetic namespace
namespace FiniteFieldMathUtils
{
	/// Naive implementation of the greatest
	/// common divisor algorithm
	template <typename T>
	T GCD(const T& a, const T& b)
	{
		T nullT = T(0);

		while (b != nullT)
		{
			T t = b;
			b = a % b;
			a = t;
		}

		return a;
	}

	/// Perform modular addition (account  for overflow)
	template <typename T>
	T ModAdd(const T& a, const T& b, const T& p)
	{
		T lhs = a % p;
		T rhs = b % p;
		T sum = lhs + rhs;

		if (T(0) == rhs) return lhs;

		if (sum >= p || sum < lhs)
			sum -= p;

		return sum;
	}

	/// Perform modular subtraction (account  for overflow)
	template <typename T>
	T ModSub(const T& a, const T& b, const T& p)
	{
		T lhs = a % p;
		T rhs = b % p;

		if (lhs >= rhs)
			return lhs - rhs;
		return p - rhs + lhs;
	}

    /// Perform modular multiplication of a
	/// with respect to p
	template <typename T>
	T ModMul(T a, T b, const T& p)
	{
		T result = T(0);
		T nullT = T(0);
		T oneT = T(1);

		a %= p;

		while (b > nullT)
		{
			if ((b & oneT) != nullT)
				result = ModAdd<T>(result, a, p);

			a = ModAdd<T>(a, a, p);
			b >>= 1;
		}

		return result % p;
	}

	/// Perform modular exponentiation of
	/// a to the power of b mod p
	template <typename T>
	T ModPow(const T& a, const T& b, const T& p)
	{
		T result = T(1);
		T nullT = T(0);
		T oneT = T(1);
		T lhs = a % p;
		T rhs = b;

		while (b != nullT)
		{
			if (b & oneT)
				result = (result * lhs) % p;

			lhs = (lhs * lhs) % p;
			rhs >>= 1;
		}

		return result;
	}

	/// Calculate the greatest common divisor
	/// using the extended Euclidean algorithm
	template <typename T>
	std::tuple<T, T, T> EGCD(const T& a, const T& b)
	{
		T nullT = T(0);

		if (a == nullT)
			return std::tuple<T, T, T>(b, T(0), T(1));
		else
		{
			auto tup = EGCD(b % a, a);
			auto g = std::get<0>(tup);
			auto y = std::get<1>(tup);
			auto x = std::get<2>(tup);

			return std::tuple<T, T, T>(g, x - (b / a) * y, y);
		}
	}

	/// Calculate the Legendre symbol of a
	/// with respect to p
	template <typename T>
	T LegendreSymbol(const T& a, const T& p)
	{
		return ModPow<T>(a, ModSub<T>(p, 1, p) / 2, p);
	}

	/// Returns true if a is a quadratic residue of p
	/// Returns false if a is not a quadratic residue of p
	template <typename T>
	bool QuadraticResidue(const T& a, const T& p)
	{
		return static_cast<bool>(LegendreSymbol<T>(a, p) == T(1));
	}

	/// Computes the multiplicative modular inverse
	/// of a with respect to p (faster due to recursion)
	template <typename T>
	T ModInv(const T& a, const T& p)
	{
		auto result = EGCD<T>(a, p);
		if (std::get<0>(result) == T(1))
			return (p + std::get<1>(result)) % p;
		else
			throw MathImplExceptions::MODINV_PARAMS_NOT_COPRIME;
	}

	/// Computes the multiplicative modular inverse
	/// of a with respect to p (prime)
	template <typename T>
	T ModInvPrime(const T& a, const T& p)
	{
		return ModPow<T>(a, p - T(2), p);
	}

	/// Implementation of the Tonneli-Shanks algorithm
	/// for solving congruences of the form a^2 congr b (mod p)
	/// a.k.a find square root of b modulo p
	template <typename T>
	std::vector<T> TonneliShanks(const T& a, const T& p)
	{
		if (!QuadraticResidue<T>(a, p))
			return {};

		T oneT = T(1);

		T q = p - oneT;
		T s = T(0);

		while (~q & oneT)
		{
			q >>= oneT;
			s += oneT;
		}

		if (s == oneT)
		{
			auto x = ModPow<T>(a, (p + 1) / 4, p);
			return { x, p - x };
		}

		T z = T(0);
		for (T k = oneT; k < p; ++k)
		{
			if (!QuadraticResidue<T>(k, p))
			{
				z = k;
				break;
			}
		}

		T c = ModPow<T>(z, q, p);
		T r = ModPow<T>(a, (q + oneT) / T(2), p);
		T t = ModPow<T>(a, q, p);
		T m = s;

		while (t != oneT)
		{
			T i = oneT;
			T x = (t * t) % p;

			while (x != oneT)
			{
				x = (x * x) % p;
				i += oneT;
			}

			T b = ModPow<T>(c, (oneT << ModSub(ModSub(m, i, p), oneT, p)), p);

			r = ModMul(r, b, p);
			c = ModMul(b, b, p);
			t = ModMul(t, c, p);
			m = i;
		}

		return { r, p - r };
	}

	/// Implementation which returns only one root from
	/// the vector provided by the Tonneli-Shanks algorithm
	template <typename T>
	T ModSqrt(const T& a, const T& p)
	{
		auto congruenceVector = TonneliShanks<T>(a, p);
		if (congruenceVector.empty())
			return T(0);

		return congruenceVector[0];
	}
}

template <class T>
class EllipticCurvePoint;

/// Elliptic curve class
template <class T>
class EllipticCurve
{
private:
	T primeField;
	std::vector<T> curveConstants;

public:
	EllipticCurve()
		:
		primeField(0),
		curveConstants(std::vector<T>(EllipticCurveConstants, 0))
	{}

	EllipticCurve(const T& prime, const T& a, const T& b)
		:
		primeField(prime),
		curveConstants(std::vector<T>({ T(1), T(1), T(0), a, b }))
	{
		using FiniteFieldMathUtils::ModMul;
		using FiniteFieldMathUtils::ModAdd;

		// Elliptic curves cannot be singular
		T lhs_parameter = ModMul<T>(T(4),  ModMul<T>(ModMul<T>(a, a, prime), a, prime), prime);
		T rhs_parameter = ModMul<T>(T(27), ModMul<T>(b, b, prime), prime);

		if(ModAdd<T>(lhs_parameter, rhs_parameter, prime) == T(0))
			throw EllipticCurveExceptions::CURVE_IS_SINGULAR;
	}

	EllipticCurve(const EllipticCurve<T>& other)
		:
		primeField(other.primeField),
		curveConstants()
	{
		curveConstants = other.curveConstants;
	}

	~EllipticCurve() = default;

	/// Evaluate the curve at x
	T evaluatePoint(const T& x)
	{
		const T evaluation = curveConstants[0] * (curveConstants[1] * x * x * x +
			curveConstants[2] * x * x + curveConstants[3] * x + curveConstants[4]);

		return evaluation;
	}

	/// Get the prime over which we have the Galois Field
	const T& getPrime() const
	{
		return primeField;
	}

	/// Get the curve constants
	const std::vector<T>& getConstants() const
	{
		return curveConstants;
	}

	/// Comparsion operator 'is_equal'
	bool operator==(const EllipticCurve<T>& other)
	{
		return ((primeField == other.primeField) && (curveConstants == other.curveConstants));
	}

	/// Comparsion operator 'is_not_equal'
	bool operator!=(const EllipticCurve<T>& other)
	{
		return !((primeField == other.primeField) && (curveConstants == other.curveConstants));
	}

	/// Access the curve constants
	T& operator[](size_t idx) const
	{
		return curveConstants[idx];
	}

	/// Get point at (x, y)
	EllipticCurvePoint<T> operator()(const T& x, const T& y)
	{
		using FiniteFieldMathUtils::ModMul;
		using FiniteFieldMathUtils::ModAdd;

		const T lhs = ModMul<T>(ModMul<T>(curveConstants[0], y, primeField), y, primeField);
		const T rhs = ModAdd<T>(ModAdd<T>(ModAdd<T>(ModMul<T>(ModMul<T>(ModMul<T>(x, x, primeField), x, primeField),
			curveConstants[1], primeField),
			ModMul<T>(ModMul<T>(x, x, primeField), curveConstants[2], primeField), primeField),
			ModMul<T>(x, curveConstants[3], primeField), primeField), curveConstants[4], primeField);

		if (lhs != rhs)
		{
			throw EllipticCurveExceptions::POINT_NOT_ON_CURVE;
			return {};
		}

		return EllipticCurvePoint<T>(x, y, *this);
	}

	/// Pretty print output stream operator
	friend std::ostream& operator<<(std::ostream& out, const EllipticCurve<T>& curve)
	{
		out << "Elliptic curve ";

		const T nullT = T(0), oneT = T(1);

		if (curve.curveConstants[0] == nullT)
			out << "0 = ";
		else if (curve.curveConstants[0] == oneT)
			out << "y ** 2 = ";
		else
			out << curve.curveConstants[0] << " * y ** 2 = ";

		if (curve.curveConstants[1] == nullT)
			;
		else if (curve.curveConstants[1] == oneT)
			out << "x ** 3 + ";
		else
			out << curve.curveConstants[1] << " * x ** 3 + ";

		if (curve.curveConstants[2] == nullT)
			;
		else if (curve.curveConstants[2] == oneT)
			out << "x ** 2 + ";
		else
			out << curve.curveConstants[2] << " * x ** 2 + ";

		if (curve.curveConstants[3] == nullT)
			;
		else if (curve.curveConstants[3] == oneT)
			out << "x + ";
		else
			out << curve.curveConstants[3] << " * x + ";

		if (curve.curveConstants[4] == nullT)
			;
		else
			out << curve.curveConstants[4];

		return out;
	}
};

template <class T>
class EllipticCurvePoint
{
	friend class EllipticCurve<T>;

private:
	T x;
	T y;
	EllipticCurve<T> parentCurve;

	EllipticCurvePoint(const T& _x, const T& _y, const EllipticCurve<T>& _parentCurve)
		:
		x(_x),
		y(_y),
		parentCurve()
	{
		parentCurve = _parentCurve;
	}

public:
	EllipticCurvePoint()
		:
		x(0),
		y(0),
		parentCurve()
	{}

	EllipticCurvePoint(const EllipticCurvePoint<T>& other)
		:
		x(other.x),
		y(other.y),
		parentCurve()
	{
		parentCurve = other.parentCurve;
	}

	~EllipticCurvePoint() = default;

	/// Point doubling
	void Double()
	{
		using FiniteFieldMathUtils::ModMul;
		using FiniteFieldMathUtils::ModInv;
		using FiniteFieldMathUtils::ModSub;
		using FiniteFieldMathUtils::ModAdd;

		T primeField = parentCurve.getPrime();

		T lambda_lhs = ModAdd<T>(ModMul<T>(ModMul<T>(this->x, this->x, primeField), T(3), primeField),
			 parentCurve.getConstants()[3], primeField);
		T lambda_rhs = ModMul<T>(T(2), y, primeField);

		T mod_inv = ModInv<T>(lambda_rhs, primeField);

		T lambda = ModMul<T>(lambda_lhs, mod_inv, primeField);
		T newX = ModSub(ModMul<T>(lambda, lambda, primeField), ModMul(T(2), x, primeField), primeField);
		T newY = ModSub(ModMul<T>(lambda, ModSub<T>(this->x, newX, primeField), primeField), this->y, primeField);

		this->x = newX;
		this->y = newY;
	}

	/// Get inverse of point (mirror against x axis)
	/// ... in short Weierstrass form
	EllipticCurvePoint<T> Inverse()
	{
		using FiniteFieldMathUtils::ModSub;
		return parentCurve(this->x, ModSub<T>(T(0), this->y, parentCurve.getPrime()));
	}

	/// Is the point at infinity ?
	bool AtInfinity()
	{
		T nullT(0);
		return (this->x == nullT) && (this->y == nullT);
	}

	/// Get the order of the point on the curve
	T Order()
	{
		T point_order(1);
		EllipticCurvePoint<T> Q = *this;

		while(!Q.AtInfinity())
		{
			Q += *this;
			++ point_order;
		}

		return point_order;
	}

	/// Discrete logarithm using baby-step giant-step
	/// TODO


	/// Point addition on the elliptic curve
	EllipticCurvePoint<T> operator+(const EllipticCurvePoint<T>& other)
	{
		using FiniteFieldMathUtils::ModMul;
		using FiniteFieldMathUtils::ModInv;
		using FiniteFieldMathUtils::ModSub;

		T primeField = parentCurve.getPrime();

		if (parentCurve != other.parentCurve)
			return {};

		if (*this == other)
		{
			EllipticCurvePoint<T> result = *this;
			result.Double();
			return result;
		}

		// Symmetric points (with respect to x-axis)
		if (this->x == other.x)
			// Return point at infinity
			return EllipticCurvePoint<T>(T(0), T(0), parentCurve);

		T lambda_lhs = ModSub<T>(other.y, this->y, primeField);
		T lambda_rhs = ModSub<T>(other.x, this->x, primeField);

		T mod_inv = ModInv<T>(lambda_rhs, primeField);
		T lambda = ModMul<T>(lambda_lhs, mod_inv, primeField);

		T newX = ModSub<T>(ModSub<T>(ModMul<T>(lambda, lambda, primeField), this->x, primeField), other.x, primeField);
		T newY = ModSub<T>(ModMul<T>(lambda, ModSub<T>(this->x, newX, primeField), primeField), this->y, primeField);

		return EllipticCurvePoint<T>(newX, newY, parentCurve);
	}

	/// Point addition on the elliptic curve
	EllipticCurvePoint<T>& operator+=(const EllipticCurvePoint<T>& other)
	{
		using FiniteFieldMathUtils::ModMul;
		using FiniteFieldMathUtils::ModInv;
		using FiniteFieldMathUtils::ModSub;

		// Get the prime of the finite field
		T primeField = parentCurve.getPrime();

		// Cannot add points on different curves
		if (parentCurve != other.parentCurve)
		{
			throw EllipticCurveExceptions::POINT_NOT_ON_CURVE;
			return *this;
		}

		// Treat doubling case
		if (*this == other)
		{
			this->Double();
			std::cout << "Doubling...\n";
			return *this;
		}

		// Symmetric points (with respect to x-axis)
		if (this->x == other.x)
		{
			// Set point at infinity
			// Potential edge-case
			this->x = 0;
			this->y = 0;

			std::cout << "Returning point at infinity...\n";

			// Return point at infinity
			return *this;
		}

		T lambda_lhs = ModSub<T>(other.y, this->y, primeField);
		T lambda_rhs = ModSub<T>(other.x, this->x, primeField);

		T mod_inv = ModInv<T>(lambda_rhs, primeField);
		T lambda = ModMul<T>(lambda_lhs, mod_inv, primeField);

		T newX = ModSub<T>(ModSub<T>(ModMul<T>(lambda, lambda, primeField), this->x, primeField), other.x, primeField);
		T newY = ModSub<T>(ModMul<T>(lambda, ModSub<T>(this->x, newX, primeField), primeField), this->y, primeField);

		x = newX;
		y = newY;

		return *this;
	}

	/// Use Montgomery ladder to perform scalar multiplication
	template <typename U>
	EllipticCurvePoint<T>& operator*=(const U& scalar)
	{
		if (scalar == U(2))
		{
			this->Double();
			return *this;
		}

		using BinaryUtils::CountBits;
		using BinaryUtils::GetBit;

		EllipticCurvePoint<T> R0 = *this , R1 = *this;
		size_t bit_count = CountBits<U>(scalar);
		R1.Double();

		for (size_t i = bit_count - 1; i > 0; --i)
		{
			if (GetBit<T>(scalar, i - 1) == 0)
			{
				R1 = R0 + R1;
				R0.Double();
			}
			else
			{
				R0 = R0 + R1;
				R1.Double();
			}
		}

		*this = R0;
		return *this;
	}

	/// Access either X or Y (0 or 1)
	T operator[](size_t idx)
	{
		if (idx == 0u)
			return x;
		else if (idx == 1u)
			return y;
		else return {};
	}

	/// Comparsion operator 'is_equal'
	bool operator==(const EllipticCurvePoint<T>& other)
	{
		return ((x == other.x) && (y == other.y) && (parentCurve == other.parentCurve));
	}

	/// Comparsion operator 'is_not_equal'
	bool operator!=(const EllipticCurvePoint<T>& other)
	{
		return !((x == other.x) && (y == other.y) && (parentCurve == other.parentCurve));
	}

	/// Output stream operator
	friend std::ostream& operator<<(std::ostream& out, EllipticCurvePoint<T>& point)
	{
		out << "Point (" << point.x << ", " << point.y << ") on curve " << point.parentCurve;
		return out;
	}
};
