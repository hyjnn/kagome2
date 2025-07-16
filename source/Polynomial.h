#pragma once

#include <iostream>

#include <Eigen/Dense>

#ifdef __GNUG__
#include <quadmath.h>

namespace Eigen {
    template<> struct NumTraits<__float128> : GenericNumTraits<__float128> {
        typedef __float128 Real;
        typedef __float128 NonInteger;
        typedef const __float128& Nested;

        static inline Real epsilon() { return __float128(0); }
        static inline Real dummy_precision() { return __float128(0); }
        static inline int digits10() { return 0; }

        enum {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = Eigen::HugeCost, //As per docs - when in doubt just use HugeCost
            AddCost = Eigen::HugeCost,
            MulCost = Eigen::HugeCost
        };
    };
}

namespace Eigen::numext {
    template<>
    EIGEN_DEVICE_FUNC EIGEN_ALWAYS_INLINE
    typename internal::enable_if<NumTraits<__float128>::IsSigned || NumTraits<__float128>::IsComplex,typename NumTraits<__float128>::Real>::type
    abs(const __float128 &x) {
        EIGEN_USING_STD(abs);
        return fabsq(x);
    }
}

namespace std {
    template <>
    struct formatter<__float128> : formatter<string> {
        auto format(__float128 x, format_context& ctx) const {
            char buffer[64];
            quadmath_snprintf(buffer, sizeof buffer, "%-40.34Qe", x);
            return formatter<string>::format(std::format("{}", buffer), ctx);
        }
    };
}
#endif

//going clockwise, starting in the upper left arrow with (kagome) node positioned horizontally:
//x - 1 1 1 1
//y - 1 0 0 1
//z - 1 1 0 0
class Polynomial {

public:
#ifdef __GNUG__
    using CoeffScalar = __float128;
#endif
#ifdef _MSC_VER
    using CoeffScalar = double;
#endif
    using CoeffMatrix = Eigen::Matrix<CoeffScalar, Eigen::Dynamic, Eigen::Dynamic>;

	CoeffMatrix coefficients;

    Polynomial() = default;
	explicit Polynomial(CoeffMatrix coefList) : coefficients(coefList) {}
	explicit Polynomial(size_t n) : coefficients(CoeffMatrix::Zero(n + 1, n + 1)) {} //This constructor takes the polynomial's degree as its argument

	inline size_t degree() const;

    bool friend operator==(const Polynomial&, const Polynomial&);

	Polynomial& operator+=(const Polynomial&);
	Polynomial& operator-=(const Polynomial&);
	Polynomial& operator*=(const Polynomial&);
	Polynomial& operator*=(double);

	friend Polynomial operator+(Polynomial, const Polynomial&);
	friend Polynomial operator-(Polynomial, const Polynomial&);
	friend Polynomial operator*(const Polynomial&, const Polynomial&);
	friend Polynomial operator*(const Polynomial&, double);
	friend Polynomial operator*(double, const Polynomial&);
	friend std::ostream& operator<<(std::ostream&, const Polynomial&);
};

inline size_t Polynomial::degree() const {
	return coefficients.rows() - 1;
}

namespace Eigen {
    typedef Matrix<Polynomial, Eigen::Dynamic, Eigen::Dynamic> MatrixXPoly;

    template<> struct NumTraits<Polynomial> : GenericNumTraits<Polynomial>
    {
        typedef Polynomial Real;
        typedef Polynomial NonInteger;
        typedef const Polynomial& Nested;

        static inline Real epsilon() { return Polynomial(0); }
        static inline Real dummy_precision() { return Polynomial(0); }
        static inline int digits10() { return 0; }

        enum {
            IsInteger = 0,
            IsSigned = 1,
            IsComplex = 0,
            RequireInitialization = 1,
            ReadCost = Eigen::HugeCost,
            AddCost = Eigen::HugeCost,
            MulCost = Eigen::HugeCost
        };
    };

}