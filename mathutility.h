#ifndef MATHUTILITY_H
#define MATHUTILITY_H

/*!
 * Defines some frequently used math functions
 */

namespace math {

/*!
 * Something in square it could be a vector or a number whatsever
 */
template<class sqT, class T> inline sqT sqr(const T& t)
{
    return t*t;
}

/*!
 * Weighting mean using vector::iterator interface
 * accumulator should be initialized with zero value of type T1
 * lengths of v and w are supposed to be equal
 */
template<class T1, class T2, template<class ...> class Vector>
inline void mean
(
        const Vector<T1>& v,
        const Vector<T2>& w,
        T1& accumulator)
{
    typename Vector<T1>::const_iterator it1 = v.begin();
    typename Vector<T2>::const_iterator it2 = w.begin();
    T2 total(0.0); //Total usually has a numeric type hence it can be initialized with 0.0

    for(; it1 != v.end(); ++it1, ++it2)
    {
        total += *it2;
        accumulator += (*it1 * *it2);
    }

    accumulator /= total;
}

}

#endif // MATHUTILITY_H
