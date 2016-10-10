#ifndef VECTORTEMPLATE_H
#define VECTORTEMPLATE_H

#include <iostream>

#include <assert.h>
#include <sstream>
#include <cmath>
#include <initializer_list>
#include <array>

/*!
 * Math vector of any size
 */
#include "mathutility.h"

namespace math {

using std::size_t;

/**
 *Math vector implementation
 */
template<class T, std::size_t N> struct vector_c;

/**
 * Creates elementwise and folding operations with a vector
 */
template<class T, size_t N>
struct vector_c_op
{
    template<size_t M> using vector = vector_c<T, M>;
    using type = vector_c_op<T, N>;
    using prior_type =  vector_c_op<T, N-1>;

    /**
     * Applies operator op from the left (from the first vector element)
     */
    template<class op, size_t M>
    inline void l_apply (op Op, vector<M>& v) const
    {
        static_assert(M >= N, "Vector index is out of range.");
        prior_type().l_apply(Op, v);
        Op(v[N - 1]);
    }

    template<class op, size_t M>
    inline void apply (op Op, vector<M>& v) const
    {
        static_assert(M >= N, "Vector index is out of range.");
        Op(v[N - 1]);
        prior_type().apply(Op, v);
    }

    template<class op, size_t M>
    inline void apply (op Op, vector<M>& vl, const vector<M>& vr) const
    {
        static_assert(M >= N, "Vector index is out of range.");
        Op(vl[N - 1], vr[N - 1]);
        prior_type().apply(Op, vl, vr);
    }

    /**
     * Applies a one value to an every vector element using binary op
     */
    template<class op, size_t M>
    inline void apply (op Op, vector<M>& vl, const T& val) const
    {
        static_assert(M >= N, "Vector index is out of range.");
        Op(vl[N - 1], val);
        prior_type().apply(Op, vl, val);
    }

    /**
     * Folds all elements
     */
    template<class add, class mul, size_t M> inline
    T fold (add Add, mul Mul, const vector<M>& vl, const vector<M>& vr) const
    {
        static_assert(M >= N, "Vector index is out of range.");
        T z = Mul(vl[N-1], vr[N-1]);
        return Add(z, prior_type().fold(Add, Mul, vl, vr));
    }

    /**
     * Summation trough elements of a vector
     */
    template<class add, size_t M> inline
    T sum (add Add, const vector<M>& v) const
    {
        static_assert(M >= N, "Vector index is out of range.");
        T z = v[N-1];
        return Add(z, prior_type().sum(Add, v));
    }
};

template<class T>
struct vector_c_op<T, 1>
{
    template<size_t N> using vector = vector_c<T, N>;

    template<class op, size_t N>
    inline void l_apply (op Op, vector<N>& v) const
    {
        static_assert(N != 0, "Vector index is out of range.");
        Op(v[0]);
    }

    template<class op, size_t N>
    inline void apply (op Op, vector<N>& v) const
    {
        static_assert(N != 0, "Vector index is out of range.");
        Op(v[0]);
    }

    template<class op, size_t N>
    inline void apply (op Op, vector<N>& vl, const vector<N>& vr) const
    {
        static_assert(N != 0, "Vector index is out of range.");
        Op(vl[0], vr[0]);
    }

    template<class op, size_t N>
    void apply (op Op, vector<N>& vl, const T& val) const
    {
        static_assert(N != 0, "Vector index is out of range.");
        Op(vl[0], val);
    }

    template<class add, class mul, size_t N> inline
    T fold (add /*Add*/, mul Mul, const vector<N>& vl, const vector<N>& vr) const
    {
        static_assert(N != 0, "Vector index is out of range.");
        return Mul(vl[0], vr[0]);
    }

    template<class add, size_t N> inline
    T sum (add /*Add*/, const vector<N>& v)
    {
        static_assert(N != 0, "Vector index is out of range.");
        return v[0];
    }
};

template<class T, std::size_t N> struct vector_c : std::array<T, N>
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
        vector_c_op<T, N>().apply(set_val<T>(), *this, val);
    }

    vector_c(const type& v)
    {
        vector_c_op<T, N>().apply(set_val<T>(), *this, v);
    }

    vector_c& operator = (const type& v)
    {
        vector_c_op<T, N>().apply(set_val<T>(), *this, v);
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
        vector_c_op<T, N>().l_apply(set_list(it), *this);
    }
};

/**
 * Vector and number addition
 */
template<class T, size_t N> inline
vector_c<T, N>& operator += (vector_c<T, N>& v, const T& h)
{
    vector_c_op<T, N>().apply(in_place_plus<T>(), v, h);
    return v;
}

/**
 * Vector vector inplace addition
 */
template<class T, size_t N> inline
vector_c<T, N>& operator += (vector_c<T, N>& vl, const vector_c<T, N>& vr)
{
    vector_c_op<T, N>().apply(in_place_plus<T>(), vl, vr);
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
    vector_c_op<T, N>().apply(in_place_sub<T>(), vl, vr);
    return vl;
}

/**
 * vector subtraction
 */
template<class T, size_t N> inline
vector_c<T, N> operator - (const vector_c<T, N>& vl, const vector_c<T, N>& vr)
{
    vector_c<T, N> result(vl);
    return result -= vr;
}

/**
 * Vector inplace multiplication by a number
 */
template<class T, size_t N> inline
vector_c<T, N>& operator *= (vector_c<T, N>& vl, const T& h)
{
    vector_c_op<T, N>().apply(in_place_mul<T>(), vl, h);
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
    vector_c_op<T, N>().apply(in_place_div<T>(), v, h);
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
#include <iostream>
#include <type_traits>
template<class T, size_t N> inline
T operator * (const vector_c<T, N>& vl, const vector_c<T, N>& vr)
{ return vector_c_op<T, N>().fold(in_place_plus<T>(), mul<T>(), vl, vr); }

/**
 * Vector dot multiplication functor
 */
template<class T, size_t N>
struct vector_dot_mul
{
    using vector = vector_c<T, N>;
    inline T operator () (const vector& v1, const vector& v2) const
    { return v1*v2; }
};

/**
 * Vector printing
 */
template<class T, size_t N> inline
std::ostream& operator << (std::ostream& out, const vector_c<T, N>& v)
{
    using vector = vector_c<T, N>;
    using print_val = typename vector::print_val;
    out << "( ";
    vector_c_op<T, N>().l_apply(print_val(out),const_cast<vector&>(v));
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
{ return vector_c_op<T, N>().sum(in_place_plus<T>(), v); }

/**
 * Product of all vector elements
 */
template<class T, size_t N>
inline T prod(const vector_c<T, N>& v)
{ return vector_c_op<T, N>().sum(in_place_mul<T>(), v); }

}

#endif // VECTORTEMPLATE_H
