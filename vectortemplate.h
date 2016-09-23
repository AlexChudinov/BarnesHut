#ifndef VECTORTEMPLATE_H
#define VECTORTEMPLATE_H

#include <inttypes.h>
#include <sstream>
#include <initializer_list>

/*!
 * Math vector of any size
 */

namespace math {

template <class T, uint64_t N> struct Vector : Vector <T, N-1>
{    
    T x;

    /*!
     *Returns vector components
     */
    template<uint64_t M> T& comp() { return Vector<T, M>::x; }
    template<uint64_t M> const T& comp() const { return Vector<T, M>::x; }
    T& X() { return this->comp<1>(); }
    const T& X() const { return this->comp<1>(); }
    T& Y() { return this->comp<2>(); }
    const T& Y() const { return this->comp<2>(); }
    T& Z() { return this->comp<3>(); }
    const T& Z() const { return this->comp<3>(); }

    //Prints vector
    std::ostream& print(std::ostream& out) const
    {
        static_cast<const Vector<T, N-1>*>(this)->print(out);
        out << ", " <<x;
        return out;
    }

    /*!
     * \brief Vector default constructor fills with zeros
     */
    Vector()
        :
          Vector<T, N-1>(),
          x(0.0)
    {}

    /*!
     * Copy Vector
     */
    Vector(const Vector& v)
        :
          Vector<T, N-1>(static_cast<const Vector<T, N-1>&>(v)),
          x(v.x)
    {}

    /*!
     * \brief Vector initialize list constructor
     * \param list
     */
    Vector(std::initializer_list<T> list)
        :
          Vector<T, N-1>(list),
          x(*(list.begin() + N - 1))
    {}

    /*!
     * Vector in place summation, without looping
     */
    Vector& operator+=(const Vector& v)
    {
        x+=v.x;
        *static_cast<Vector<T, N-1>*>(this)
                += static_cast<const Vector<T, N-1>&>(v);
        return *this;
    }

    /*!
     * \brief operator -= Vector in place subtraction
     * \param v subtractant
     * \return Vector subtraction
     */
    Vector& operator-=(const Vector& v)
    {
        x-=v.x;
        *static_cast<Vector<T, N-1>*>(this)
                -= static_cast<const Vector<T, N-1>&>(v);
        return *this;
    }

    /*!
     * \brief operator *= Vector in place multiplication by a number
     * \param h Number
     * \return Vector muliplied by a number
     */
    Vector& operator*=(const T& h)
    {
        x*=h;
        *static_cast<Vector<T, N-1>*>(this) *= h;
        return *this;
    }

    /*!
     * \brief operator /= Vector in place division by a number
     * \param h Number
     * \return Vector divided by a number
     */
    Vector& operator/=(const T& h)
    {
        x/=h;
        *static_cast<Vector<T, N-1>*>(this) /= h;
        return *this;
    }

    /*!
     * \brief dotMul dot vector product
     * \param v right hand side vector
     * \return dot vector product
     */
    T dotMul(const Vector& v) const
    {
        T result = x*v.x;
        return
                result +=
                static_cast<const Vector<T, N-1>*>(this)->dotMul
                (static_cast<const Vector<T, N-1>&>(v));
    }
};

template <class T> struct Vector<T, 1>
{
    T x;

    /*!
     *Returns vector components
     */

    template<uint64_t M> T& comp() { return Vector<T, M>::x; }
    template<uint64_t M> const T& comp() const { return Vector<T, M>::x; }

    //Prints vector
    std::ostream& print(std::ostream& out) const
    {
        out << x;
        return out;
    }

    /*!
     * \brief Vector default constructor fills with zeros
     */
    Vector(): x(0.0) {}

    /*!
     * Copy Vector
     */
    Vector(const Vector& v) : x(v.x) {}

    /*!
     * \brief Vector initialize list constructor
     * \param list
     */
    Vector(std::initializer_list<T> list)
        :
          x(*list.begin())
    {}

    /*!
     * Vector in place summation, without looping
     */
    Vector& operator+=(const Vector& v)
    {
        x+=v.x;
        return *this;
    }

    /*!
     * \brief operator -= Vector in place subtraction
     * \param v subtractant
     * \return Vector subtraction
     */
    Vector& operator-=(const Vector& v)
    {
        x-=v.x;
        return *this;
    }

    /*!
     * \brief operator *= Vector in place multiplication by a number
     * \param h Number
     * \return Vector muliplied by a number
     */
    Vector& operator*=(const T& h)
    {
        x*=h;
        return *this;
    }

    /*!
     * \brief operator /= Vector in place division by a number
     * \param h Number
     * \return Vector divided by a number
     */
    Vector& operator/=(const T& h)
    {
        x/=h;
        return *this;
    }

    /*!
     * \brief dotMul dot vector product
     * \param v right hand side vector
     * \return dot vector product
     */
    T dotMul(const Vector& v) const
    {
        return x*v.x;
    }
};

template<class T, uint64_t N>
std::ostream& operator<<(std::ostream& out, const Vector<T, N>& v)
{
    out << "("; v.print(out); out << ")";
    return out;
}

/*!
 *Postfix followed by a prefix Vector multiplication by a number
 */
template<class T, uint64_t N>
Vector<T, N> operator*(const Vector<T, N>& v, const T& h)
{
    Vector<T, N> result = v;
    return result *= h;
}

template<class T, uint64_t N>
Vector<T, N> operator*(const T& h, const Vector<T, N>& v)
{
    Vector<T, N> result = v;
    return result *= h;
}

/*!
 *Dot vector multiplication
 */
template<class T, uint64_t N>
T operator*(const Vector<T, N>& v1, const Vector<T, N>& v2)
{
    return v1.dotMul(v2);
}

}

#endif // VECTORTEMPLATE_H
