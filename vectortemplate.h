#ifndef VECTORTEMPLATE_H
#define VECTORTEMPLATE_H

#include <functional>
#include <assert.h>
#include <sstream>
#include <cmath>
#include <initializer_list>
#include <array>

/*!
 * Math vector of any size
 */
#include <array_operations.h>
#include "mathutility.h"

namespace math {

/**
 *Math vector implementation
 */
template<class T, size_t N> struct vector_c : std::array<T, N>
{
    using type = vector_c<T, N>;

    /**
     * Sets a vector from a list
     */
    template<class ListIterator>
    struct set_list
    {
        ListIterator& cur_;
        inline set_list(ListIterator& first)
            : cur_(first){}
        inline void operator () (T& x) { x = *(this->cur_++); }
    };

    /**
     * Prints vector element
     */
    struct print_val
    {
        std::ostream& out_;
        inline print_val(std::ostream& out):out_(out){}
        inline void operator () (const T& x) const
        { out_ << x << " "; }
    };

    vector_c() {}

    explicit vector_c(const T& val)
    {
        struct set_val{
            const T& val_;
            inline set_val(const T& val):val_(val){}
            inline void operator()(T& x) const { x = val_; }
        };

        math::array_operations<type, type, N - 1> op;
        op.umap(set_val(val), *this);
    }

    vector_c(const type& v)
    {
        math::array_operations<type, type, N - 1> op;
        op.bmap(set_val<T>(), *this, v);
    }

    vector_c& operator=(const type& v)
    {
        math::array_operations<type, type, N - 1> op;
        op.bmap(set_val<T>(), *this, v);
        return *this;
    }

    vector_c(std::initializer_list<T> list)
    {
        assert(list.size() == N);
        using list_iterator =
            typename std::initializer_list<T>::const_iterator;
        using set_list
            = set_list<list_iterator>;
        list_iterator it = list.begin();
        math::array_operations<type, type, N - 1> op;
        op.umap(set_list(it), *this);
    }
};

/**
 * Vector and number addition
 */
template<class T, size_t N> inline
vector_c<T, N>& operator += (vector_c<T, N>& v, const T& h)
{
    DEF_OPERATION_WITH_VAL(T, +=)
    using type = typename vector_c<T, N>::type;
    math::array_operations<type, type, N - 1> op;
    op.umap(operation(h), v);
    return v;
}

/**
 * Vector vector inplace addition
 */
template<class T, size_t N> inline
vector_c<T, N>& operator += (vector_c<T, N>& vl, const vector_c<T, N>& vr)
{
    using type =  typename vector_c<T, N>::type;
    math::array_operations<type, type, N - 1> op;
    op.bmap(in_place_plus<T>(), vl, vr);
    return vl;
}

/**
 * Vector vector addition
 */
template<class T, size_t N> inline
vector_c<T, N> operator  + (const vector_c<T, N>& vl, const vector_c<T, N>& vr)
{ vector_c<T, N> result(vl); return result+=vr;}

/**
 * vector in place subtraction
 */
template<class T, size_t N> inline
vector_c<T, N>& operator -= (vector_c<T, N>& vl, const vector_c<T, N>& vr)
{
    using type = typename vector_c<T, N>::type;
    math::array_operations<type, type, N - 1> op;
    op.bmap(in_place_sub<T>(), vl, vr);
    return vl;
}

/**
 * vector subtraction
 */
template<class T, size_t N> inline
vector_c<T, N> operator - (const vector_c<T, N>& vl, const vector_c<T, N>& vr)
{ vector_c<T, N> result(vl); return result -= vr; }

/**
 * Vector inplace multiplication by a number
 */
template<class T, size_t N> inline
vector_c<T, N>& operator *= (vector_c<T, N>& vl, const T& h)
{
    DEF_OPERATION_WITH_VAL(T, *=)
    using type = typename vector_c<T, N>::type;
    math::array_operations<type, type, N - 1> op;
    op.umap(operation(h), vl);
    return vl;
}

/**
 * Vector multiplication by a number
 */
template<class T, size_t N> inline
vector_c<T, N> operator * (const vector_c<T, N>& vl, const T& h)
{ vector_c<T, N> result(vl); return result *= h; }
/**
*/
template<class T, size_t N> inline
vector_c<T, N> operator * (const T& h, const vector_c<T, N>& vl)
{ vector_c<T, N> result(vl); return result *= h; }

/**
 * Vector inplace division by a number
 */
template<class T, size_t N> inline
vector_c<T, N>& operator /= (vector_c<T, N>& v, const T& h)
{
    DEF_OPERATION_WITH_VAL(T, /=)
    using type = typename vector_c<T, N>::type;
    math::array_operations<type, type, N - 1> op;
    op.umap(operation(h), v);
    return v;
}
/**
 * Vector division by a number
 */
template<class T, size_t N> inline
vector_c<T, N> operator / (const vector_c<T, N>& v, const T& h)
{ vector_c<T, N> result(v); return result /= h; }

/**
 * Vector dot multiplication
 */
template<class T, size_t N> inline
T operator * (const vector_c<T, N>& vl, const vector_c<T, N>& vr)
{
    using type = typename vector_c<T, N>::type;
    math::array_operations<type, type, N - 1> op;
    T res = 0.0;
    op.bfold(math::in_place_plus<T>(),std::multiplies<T>(),res, vl, vr);
    return res;
}

/**
 * Vector printing
 */
template<class T, size_t N> inline
std::ostream& operator << (std::ostream& out, const vector_c<T, N>& v)
{
    using type = typename vector_c<T, N>::type;
    math::array_operations<type, type, N - 1> op;
    using print_val = typename type::print_val;
    out << "( ";
    op.umap(print_val(out),const_cast<type&>(v));
    out << ")";
    return out;
}

/**
 * Vector's length
 */
template<class T, size_t N>
inline T abs(const vector_c<T, N>& v){ return ::sqrt(v*v); }

/**
 * Sum of all vector elements
 */
template<class T, size_t N>
inline T sum(const vector_c<T, N>& v)
{
    using type = typename vector_c<T, N>::type;
    math::array_operations<type, type, N - 1> op;
    T res = 0.0;
    op.ufold(in_place_plus<T>(), res, v);
    return res;
}

/**
 * Product of all vector elements
 */
template<class T, size_t N>
inline T prod(const vector_c<T, N>& v)
{
    using type = typename vector_c<T, N>::type;
    math::array_operations<type, type, N - 1> op;
    T res = 1.0;
    op.ufold(in_place_mul<T>(), res, v);
    return res;
}

}

#endif // VECTORTEMPLATE_H
