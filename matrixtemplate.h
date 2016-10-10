#ifndef MATRIXTEMPLATE_H
#define MATRIXTEMPLATE_H

#include <utility>
#include "vectortemplate.h"
#include <vector>

/**
 * Template defines simple math matrix operations
 */

namespace math {

using std::pair;
using std::vector;

template<class T, size_t nRows, size_t nCols> class matrix_c;

template<class T, size_t M, size_t N>
struct matrix_c_op
{
    using type         = matrix_c_op<T, M, N>;
    using prior_column = matrix_c_op<T, M, N-1>;
    using prior_row    = matrix_c_op<T, M-1, N>;
    using prior_diag   = matrix_c_op<T, M-1, N-1>;

    /**
     * Matrix operation with a number
     */
    template<class op, size_t nRows, size_t nCols>
    inline void apply
    (
            op Op,
            matrix_c<T, nRows, nCols>& result,
            const T& Val) const
    {
        static_assert(nRows >= M, "Matrix index is out of range.");
        vector_c_op<T, nRows>().apply(Op, result[M-1], Val);
        prior_row().apply(Op, result, Val);
    }

    /**
     * Folds matrix's rows
     */
    template<class op, size_t nRows, size_t nCols>
    inline void fold
    (
            op Op,
            const matrix_c<T, nRows, nCols>& m,
            const vector_c<T, nCols>& v,
            vector_c<T, nRows>& result) const
    {
        static_assert(nRows >= M, "Matrix element index is out of range.");
        result[M-1] = Op(m[M-1], v);
        prior_row().fold(Op, m, v, result);
    }

    template<class add, class op, size_t nRows, size_t nCols>
    inline void fold
    (
            add Add,
            op Op,
            const vector_c<T, nRows>& v,
            const matrix_c<T, nRows, nCols>& m,
            vector_c<T, nCols>& result) const
    {
        static_assert(nRows >= M, "Matrix element index is out of range.");
        static_assert(nCols >= N, "Matrix element index is out of range.");

        Add(result[N-1], Op(v[M-1], m[M-1][N-1]));
        prior_row().fold(Add, Op, v, m, result);
    }

    template<class add, class op, size_t nRows, size_t nCols, size_t nCols1>
    inline void fold
    (
            add Add,
            op Op,
            const matrix_c<T, nRows, nCols>& m1,
            const matrix_c<T, nCols, nCols1>& m2,
            matrix_c<T, nRows, nCols1>& result) const
    {
        static_assert(nRows >= M, "Matrix element index is out of range.");

        matrix_c_op<T, nCols, nCols1>().fold(Add, Op, m1[M-1], m2, result[M-1]);
        prior_row().fold(Add, Op, m1, m2, result);
    }

    template<size_t K>
    inline void transpose(matrix_c<T, K, K>& m) const
    {
        static_assert(K >= M, "Matrix element index is out of range");
        static_assert(M >= 2, "Matrix are too small");
        std::swap(m[M-2][N-1], m[N-1][M-2]);
        prior_row().transpose(m);
        prior_diag().transpose(m);
    }

    /**
     * Removes vector component projection
     */
    template<size_t nRows, size_t nCols>
    inline void component_remove
    (
        matrix_c<T, nRows, nCols>& m,
        const vector_c<T, nRows>& v) const
    {
        static_assert(nRows >= M,"Matrix index is out of bounds.");
        constexpr size_t Idx = nRows - M;
        m[Idx] = m[Idx] - (m[Idx]*v) * v/(v*v);
        prior_row().component_remove(m, v);
    }

    /**
     * Creates covariance matrix from a vector
     */
    template<size_t nRows>
    inline void cov_matrix_diag
    (
        const vector_c<T, nRows>& v,
        const T& weight,
        matrix_c<T, nRows, nRows>& result) const
    {
        static_assert(nRows >= M && nRows >= N,
                      "Matrix element index is out of range.");
        static_assert(M == N,
                      "This function should goes throug the matrix diagonal");
        result[M-1][N-1] += weight*v[M-1]*v[N-1];
        prior_row().cov_matrix_row(v, weight, result);
        prior_diag().cov_matrix_diag(v, weight, result);
    }

    template<size_t nRows>
    inline void cov_matrix_row
    (
        const vector_c<T, nRows>& v,
        const T& weight,
        matrix_c<T, nRows, nRows>& result) const
    {
        static_assert(nRows >= M && nRows >= N,
                      "Matrix element index is out of range.");
        result[N-1][M-1] = result[M-1][N-1] += weight*v[M-1]*v[N-1];
        prior_row().cov_matrix_row(v, weight, result);
    }
};

template<class T, size_t N>
struct matrix_c_op<T, 1, N>
{
    template<size_t M> using prior_column = matrix_c_op<T, M, N-1>;

    /**
     * Matrix operation with a number
     */
    template<class op, size_t nRows, size_t nCols>
    inline void apply
    (
            op Op,
            matrix_c<T, nRows, nCols>& result,
            const T& Val) const
    { vector_c_op<T, nRows>().apply(Op, result[0], Val); }

    template<class op, size_t nRows, size_t nCols>
    inline void fold
    (
            op Op,
            const matrix_c<T, nRows, nCols>& m,
            const vector_c<T, nCols>& v,
            vector_c<T, nRows>& result) const
    { result[0] = Op(m[0], v); }

    template<class add, class op, size_t nRows, size_t nCols>
    inline void fold
    (
            add Add,
            op Op,
            const vector_c<T, nRows>& v,
            const matrix_c<T, nRows, nCols>& m,
            vector_c<T, nCols>& result) const
    {
        static_assert(nCols >= N, "Matrix element index is out of range.");
        Add(result[N-1], Op(v[0], m[0][N-1]));
        prior_column<nRows>().fold(Add, Op, v, m, result);
    }

    template<class add, class op, size_t nRows, size_t nCols, size_t nCols1>
    inline void fold
    (
            add Add,
            op Op,
            const matrix_c<T, nRows, nCols>& m1,
            const matrix_c<T, nCols, nCols1>& m2,
            matrix_c<T, nRows, nCols1>& result) const
    {
        static_assert(nRows != 0, "Matrix element index is out of range.");
        matrix_c_op<T, nCols, nCols1>().fold(Add, Op, m1[0], m2, result[0]);
    }

    template<size_t K>
    inline void transpose(matrix_c<T, K, K>& /*m*/) const
    { }

    /**
     * Removes vector component projection starting from row nRows0
     */
    template<size_t nRows, size_t nCols>
    inline void component_remove
    (
        matrix_c<T, nRows, nCols>& m,
        const vector_c<T, nRows>& v) const
    {
        constexpr size_t Idx = nRows - 1;
        m[Idx] = m[Idx] - (m[Idx]*v) * v/(v*v);
    }

    /**
     * Creates covariance matrix from a vector
     */
    template<size_t nRows>
    inline void cov_matrix_row
    (
        const vector_c<T, nRows>& v,
        const T& weight,
        matrix_c<T, nRows, nRows>& result) const
    { result[N-1][0] = result[0][N-1] += weight*v[0]*v[N-1]; }
};

template<class T, size_t M>
struct matrix_c_op<T, M, 1>
{
    using prior_row = matrix_c_op<T, M-1, 1>;

    template<class add, class op, size_t nRows, size_t nCols>
    inline void fold
    (
            add Add,
            op Op,
            const vector_c<T, nRows>& v,
            const matrix_c<T, nRows, nCols>& m,
            vector_c<T, nCols>& result) const
    {
        static_assert(nRows >= M, "Matrix element index is out of range.");
        Add(result[0], Op(v[M-1], m[M-1][0]));
        prior_row().fold(Add, Op, v, m, result);
    }
};

template<class T>
struct matrix_c_op<T, 1, 1>
{
    template<class add, class op, size_t nRows, size_t nCols>
    inline void fold
    (
            add Add,
            op Op,
            const vector_c<T, nRows>& v,
            const matrix_c<T, nRows, nCols>& m,
            vector_c<T, nCols>& result) const
    { Add(result[0], Op(v[0], m[0][0])); }

    template<size_t K>
    inline void transpose(matrix_c<T, K, K>& /*m*/) const
    {}

    /**
     * Creates covariance matrix from a vector
     */
    template<size_t nRows>
    inline void cov_matrix_diag
    (
        const vector_c<T, nRows>& v,
        const T& weight,
        matrix_c<T, nRows, nRows>& result) const
    { result[0][0] += weight*v[0]*v[0]; }
};

template <class T, size_t nRows, size_t nCols = nRows>
struct matrix_c : vector_c<vector_c<T, nCols>, nRows>
{
    using base_type = vector_c<vector_c<T, nCols>, nRows>;

    matrix_c(){}

    explicit matrix_c(T const& val)
        :
          base_type(vector_c<T, nCols>(val))
    {}
};

/**
 * Matrix multiplication by a vector
 */
template<class T, size_t M, size_t N>
vector_c<T, M> operator * (const matrix_c<T, M, N>& m, const vector_c<T, N>& v)
{
    vector_c<T, M> result;
    matrix_c_op<T, M, N>().fold(vector_dot_mul<T, N>(), m, v, result);
    return result;
}

template<class T, size_t M, size_t N>
vector_c<T, N> operator * (const vector_c<T, M>& v, const matrix_c<T, M, N>& m)
{
    vector_c<T, N> result(0.0);
    matrix_c_op<T, M, N>().fold(in_place_plus<T>(), mul<T>(), v, m, result);
    return result;
}

/**
 * Matrix multiplication
 */
template<class T, size_t M, size_t K, size_t N>
matrix_c<T, M, N> operator * (const matrix_c<T, M, K>& m1, const matrix_c<T, K, N>& m2)
{
    matrix_c<T, M, N> result(0.0);
    matrix_c_op<T, M, N>().fold(in_place_plus<T>(), mul<T>(), m1, m2, result);
    return result;
}

/**
 * Matrix division by a number
 */
template<class T, size_t M, size_t N>
matrix_c<T, M, N>& operator /= (matrix_c<T, M, N>& m, const T& h)
{
    matrix_c_op<T, M, N>().apply(in_place_div<T>(), m, h);
    return m;
}
template<class T, size_t M, size_t N>
matrix_c<T, M, N> operator / (const matrix_c<T, M, N>& m, const T& h)
{
    matrix_c<T, M, N> result(m);
    return result /= h;
}

/**
 * Matrix multiplication by a number
 */
template<class T, size_t M, size_t N>
matrix_c<T, M, N>& operator *= (matrix_c<T, M, N>& m, const T& h)
{
    matrix_c_op<T, M, N>().apply(in_place_mul<T>(), m, h);
    return m;
}
template<class T, size_t M, size_t N>
matrix_c<T, M, N> operator * (const matrix_c<T, M, N>& m, const T& h)
{
    matrix_c<T, M, N> result(m);
    return result *= h;
}
template<class T, size_t M, size_t N>
matrix_c<T, M, N> operator / (const T& h, const matrix_c<T, M, N>& m)
{
    matrix_c<T, M, N> result(m);
    return result *= h;
}

/**
 * Matrix transposition
 */
template<class T, size_t M>
matrix_c<T, M, M> transpose(const matrix_c<T, M, M>& m)
{
    matrix_c<T, M, M> result(m);
    matrix_c_op<T, M, M>().transpose(result);
    return result;
}

/**
 * Calculates covariance matrix
*/
template<class T, std::size_t N>
matrix_c<T, N, N> cov
(
        const vector<vector_c<T, N>>& vectors,
        const vector<T>& w
)
{
    assert(vectors.size() == w.size());

    matrix_c<T, N, N> covMatrix(0.0);

    typename vector<T>::const_iterator it = w.begin();

    T total = 0.0;

    for (const vector_c<T, N>& v:vectors)
    {
        total += *it;
        matrix_c_op<T, N, N>().cov_matrix_diag(v, *(it++), covMatrix);
    }

    return (covMatrix/total) *
            static_cast<double>(w.size())
            /static_cast<double>(w.size()-1);
}

/**
 * Calculates first principal component
 */
template<class T, std::size_t N>
vector_c<T, N> pc1
(
        const vector<vector_c<T, N>>& vectors,
        const vector<T>& w,
        T relTol = 1.0e-10,
        size_t maxItter = 1000
)
{
    typedef vector_c<T, N> vector_type;
    typedef matrix_c<T, N, N> matrix_type;

    //Make first approximation
    vector_type eigen_vector(0.0);
    typename vector<T>::const_iterator it = w.begin();
    for(const vector_type& v : vectors)
    {
        if (eigen_vector*v < 0.0)
            eigen_vector -= (v*v)*v * *(it++);
        else
            eigen_vector += (v*v)*v * *(it++);
    }

    if(abs(eigen_vector) == 0) return eigen_vector;
    eigen_vector /= abs(eigen_vector);

    matrix_type covMatrix = cov(vectors,w);
    eigen_vector = covMatrix * eigen_vector;
    T disp0, disp1 = abs(eigen_vector);
    size_t iter = 0;
    do
    {
        disp0 = disp1;
        eigen_vector /= disp1;
        eigen_vector = covMatrix * eigen_vector;
        disp1 = abs(eigen_vector);
    } while(std::fabs(disp1-disp0)/disp0 > relTol
            && (++iter) != maxItter);

    return eigen_vector;
}

}

#endif // MATRIXTEMPLATE_H
