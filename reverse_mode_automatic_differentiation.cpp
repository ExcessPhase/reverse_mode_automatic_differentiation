#include <vector>
#include <iostream>
#include <cmath>
#include <type_traits>
#include <cassert>
#include <optional>


namespace foelsche
{
namespace rmad
{
	/// an environment object keeping values of independent variables
	/// and expecting derivative values
	/// and storing intermediate values
struct environment
{	const std::vector<double> m_sIndependent;
	std::vector<double>&m_rDer;
	typedef std::pair<double, double> doublePair;
	std::vector<std::optional<doublePair> >&m_rValues;
	environment(
		std::vector<double> _s,
		std::vector<double>&_rDer,
		std::vector<std::optional<doublePair> >&_rValues
	)
		:m_sIndependent(std::move(_s)),
		m_rDer(_rDer),
		m_rValues(_rValues)
	{	assert(m_sIndependent.size() == m_rDer.size());
	}
	double getIndependent(const std::size_t _i) const
	{	return m_sIndependent.at(_i);
	}
	void reportDerivative(const std::size_t _i, const double _d) const
	{	m_rDer[_i] += _d;
	}
};
	/// an object with a unique ID within the current thread
struct hasId
{	const std::size_t m_iIndex;
	static thread_local std::size_t s_iNextIndex;
	hasId(void)
		:m_iIndex(s_iNextIndex++)
	{
	}
};
thread_local std::size_t hasId::s_iNextIndex;
template<std::size_t ENUM>
struct independent;
template<std::size_t ENUM>
const independent<ENUM> &X(void);
	/// an independent variable
	/// immutable & unique
template<std::size_t ENUM>
struct independent:hasId
{	double calculate(const environment&_r) const
	{	if (auto &r = _r.m_rValues[m_iIndex])
			return r.value().first;
		else
			return (r = std::make_pair(_r.getIndependent(ENUM), 0.0))->first;
	}
	void backpropagate(const environment&_r, const double _d) const
	{	_r.reportDerivative(ENUM, _d);
	}
	template<std::size_t E>
	friend const independent<E> &X(void);
};
	/// the way to create an independent variable
template<std::size_t ENUM>
const independent<ENUM> &X(void)
{	static const independent<ENUM> s;
	return s;
}
template<typename L, typename R>
struct addition;
template<typename L, typename R>
const addition<L, R> &operator+(const L&_rL, const R&_rR);
	/// the result of adding two values
	/// immutable & unique
template<typename L, typename R>
struct addition:hasId
{	const L &m_sL;
	const R &m_sR;
	addition(const L&_rL, const R&_rR)
		:m_sL(_rL),
		m_sR(_rR)
	{
	}
	double calculate(const environment&_r) const
	{	if (auto &r = _r.m_rValues[m_iIndex])
			return r.value().first;
		else
			return (r = std::make_pair(m_sL.calculate(_r) + m_sR.calculate(_r), 0.0))->first;
	}
	void backpropagate(const environment&_r, const double _d) const
	{	m_sL.backpropagate(_r, _d);
		m_sR.backpropagate(_r, _d);
	}
		/// the way to create an instance
	template<typename L1, typename R1>
	friend const addition<L1, R1> &operator+(const L1&_rL, const R1&_rR)
	{	static const addition<L1, R1> s(_rL, _rR);
		return s;
	}
};
template<typename L, typename R>
struct subtraction;
template<typename L, typename R>
const subtraction<L, R> &operator-(const L&_rL, const R&_rR);
	/// the result of a substraction
	/// immutable & unique
template<typename L, typename R>
struct subtraction:hasId
{	const L &m_sL;
	const R &m_sR;
	subtraction(const L&_rL, const R&_rR)
		:m_sL(_rL),
		m_sR(_rR)
	{
	}
	double calculate(const environment&_r) const
	{	if (auto &r = _r.m_rValues[m_iIndex])
			return r.value().first;
		else
			return (r = std::make_pair(m_sL.calculate(_r) - m_sR.calculate(_r), 0.0))->first;
	}
	void backpropagate(const environment&_r, const double _d) const
	{	m_sL.backpropagate(_r, _d);
		m_sR.backpropagate(_r, -_d);
	}
		/// the way to create an object
	template<typename L1, typename R1>
	friend const subtraction<L1, R1> &operator-(const L1&_rL, const R1&_rR)
	{	static const subtraction<L1, R1> s(_rL, _rR);
		return s;
	}
};
template<typename L, typename R>
struct multiplication;
template<typename L1, typename R1>
const multiplication<L1, R1> &operator*(const L1&_rL, const R1&_rR);
	/// the result of a multiplication
	/// immutable and unique
template<typename L, typename R>
struct multiplication:hasId
{	const L &m_sL;
	const R &m_sR;
	multiplication(const L&_rL, const R&_rR)
		:m_sL(_rL),
		m_sR(_rR)
	{
	}
	double calculate(const environment&_r) const
	{	if (auto &r = _r.m_rValues[m_iIndex])
			return r.value().first;
		else
			return (r = std::make_pair(m_sL.calculate(_r) * m_sR.calculate(_r), 0.0))->first;
	}
	void backpropagate(const environment&_r, const double _d) const
	{	m_sL.backpropagate(_r, _d*m_sR.calculate(_r));
		m_sR.backpropagate(_r, _d*m_sL.calculate(_r));
	}
		/// the way to create an object
	template<typename L1, typename R1>
	friend const multiplication<L1, R1> &operator*(const L1&_rL, const R1&_rR)
	{	static const multiplication<L1, R1> s(_rL, _rR);
		return s;
	}
};
template<typename L, typename R>
struct division;
template<typename L1, typename R1>
const division<L1, R1> &operator/(const L1&_rL, const R1&_rR);
	/// the result of a division
	/// immutable and unique
template<typename L, typename R>
struct division:hasId
{	const L &m_sL;
	const R &m_sR;
	division(const L&_rL, const R&_rR)
		:m_sL(_rL),
		m_sR(_rR)
	{
	}
	double calculate(const environment&_r) const
	{	if (auto &r = _r.m_rValues[m_iIndex])
			return r.value().first;
		else
		{	const auto dInv = 1.0/m_sR.calculate(_r);
			return (r = std::make_pair(m_sL.calculate(_r)*dInv, dInv))->first;
		}
	}
	void backpropagate(const environment&_r, const double _d) const
	{	const auto dInv = _r.m_rValues[m_iIndex]->second;
		const auto d = _d*dInv;
		m_sL.backpropagate(_r, d);
		m_sR.backpropagate(_r, -d*m_sL.calculate(_r)*dInv);
	}
		/// the way to create a division object
	template<typename L1, typename R1>
	friend const division<L1, R1> &operator/(const L1&_rL, const R1&_rR)
	{	static const division<L1, R1> s(_rL, _rR);
		return s;
	}
};
template<environment::doublePair (*P)(const double), typename L>
struct NONLINEAR;
template<environment::doublePair (*P)(const double), typename L>
const NONLINEAR<P, L> &nonlinear(const L&);
	/// a nonlinear function with a single argument
template<environment::doublePair (*P)(const double), typename L>
struct NONLINEAR:hasId
{	const L &m_s;
	NONLINEAR(const L&_r)
		:m_s(_r)
	{
	}
	double calculate(const environment&_r) const
	{	if (auto &r = _r.m_rValues[m_iIndex])
			return r.value().first;
		else
			return (r = P(m_s.calculate(_r)))->first;
	}
	void backpropagate(const environment&_r, const double _d) const
	{	m_s.backpropagate(_r, _d*_r.m_rValues[m_iIndex].value().second);
	}
		/// the way to create an object
	template<environment::doublePair (*P1)(const double), typename L1>
	friend const NONLINEAR<P1, L1> &nonlinear(const L1&_r)
	{	static const NONLINEAR<P1, L1> s(_r);
		return s;
	}
};
	/// helper function to create a nonlinear object
environment::doublePair exp(const double _d)
{	const auto d = std::exp(_d);
	return std::make_pair(d, d);
}
	/// helper function to create a nonlinear object
environment::doublePair sin(const double _d)
{	return std::make_pair(std::sin(_d), std::cos(_d));
}
	/// helper function to create a nonlinear object
environment::doublePair cos(const double _d)
{	return std::make_pair(std::cos(_d), -std::sin(_d));
}
	/// helper function to create a nonlinear object
environment::doublePair sqrt(const double _d)
{	const double d = std::sqrt(_d);
	return std::make_pair(d, 0.5/d);
}
	/// an expression
static const auto &s = nonlinear<sqrt>((X<0>() + X<1>())/(X<0>() - X<1>()));
}
}
	
int main()
{	using namespace foelsche::rmad;
		/// where to keep the derivative vs X0 and X1
	std::vector<double> sDer(2);
		/// where to keep the values
	std::vector<std::optional<environment::doublePair> > sValues(hasId::s_iNextIndex);
		/// first argument are the values for the independent variables
	const environment sEnv({1.2, 1.1}, sDer, sValues);
	std::cout << "value=" << s.calculate(sEnv) << "\n";
	s.backpropagate(sEnv, 1.0);
	for (const auto d : sDer)
		std::cout << "der=" << d << "\n";	
}
