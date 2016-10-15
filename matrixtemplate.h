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

template<class T, size_t m, size_t n> class matrix_c;
template<class T, size_t m, size_t n> class proxy_matrix_row;
template<class T, size_t m, size_t n> class proxy_matrix_col;

template <class T, size_t m, size_t n>
struct matrix_c : vector_c<vector_c<T, n>, m>
{
    using base_type = vector_c<vector_c<T, n>, m>;
    using type = matrix_c<T, m, n>;

    matrix_c(){}

    explicit matrix_c(T const& val)
        :
          base_type(vector_c<T, n>(val))
    {}

    matrix_c(std::initializer_list<vector_c<T,n>> list)
        : base_type(list)
    {}

    inline proxy_matrix_col<T, m, n> column(size_t idx_col);

    inline proxy_matrix_row<T, m, n> row(size_t idx_row);
};

template<class T, size_t m, size_t n>
struct proxy_matrix_col
{
    using matrix = matrix_c<T, m, n>;
    matrix& A_;
    size_t idx_col_;

    inline proxy_matrix_col(matrix& A, size_t idx_col)
        :
          A_(A), idx_col_(idx_col){}

    inline T& operator[](size_t idx_row){ return A_[idx_row][idx_col_]; }
    inline const T& operator[](size_t idx_row) const
    { return A_[idx_row][idx_col_]; }
};

template<class T, size_t m, size_t n>
proxy_matrix_col<T, m, n> matrix_c<T, m, n>::column(size_t idx_col)
{
    return proxy_matrix_col<T, m, n>(*this, idx_col);
}

template<class T, size_t m, size_t n>
struct proxy_matrix_row
{
    using matrix = matrix_c<T, m, n>;
    matrix& A_;
    size_t idx_row_;

    inline proxy_matrix_row(matrix& A, size_t idx_row)
        :
          A_(A), idx_row_(idx_row){}

    inline T& operator[](size_t idx_col){ return A_[idx_row_][idx_col]; }
    inline const T& operator[](size_t idx_col) const
    { return A_[idx_row_][idx_col]; }
};

template<class T, size_t m, size_t n>
proxy_matrix_row<T, m, n> matrix_c<T, m, n>::row(size_t idx_row)
{
    return proxy_matrix_row<T, m, n>(*this, idx_row);
}

/**
 * Matrix multiplication by a vector
 */
template<class T, size_t m, size_t n>
vector_c<T, m> operator * (const matrix_c<T, m, n>& A, const vector_c<T, n>& x)
{
    using matrix = matrix_c<T, m, n>;
    using vector_row = vector_c<T, n>;
    using vector_col = vector_c<T, m>;

    vector_col res;

    struct fold_matrix_row
    {
        const matrix& A_;
        const vector_row& x_;
        size_t idx_row_;

        inline fold_matrix_row(const matrix& A, const vector_row& x)
            :
              A_(A), x_(x), idx_row_(0) {}

        inline void operator()(T& y)
        {
            y = A_[idx_row_++]*x_;
        }
    } fold_row(A, x);

    math::array_operations<vector_col, vector_col, m - 1> op;
    op.umap(std::ref(fold_row), res);

    return res;
}

template<class T, size_t m, size_t n>
vector_c<T, N> operator * (const vector_c<T, M>& v, const matrix_c<T, M, N>& m)
{
    vector_c<T, N> result(0.0);
    matrix_c_op<T, M, N>().fold(in_place_plus<T>(), mul<T>(), v, m, result);
    return result;
}

///**
// * Matrix multiplication
// */
//template<class T, size_t M, size_t K, size_t N>
//matrix_c<T, M, N> operator * (const matrix_c<T, M, K>& m1, const matrix_c<T, K, N>& m2)
//{
//    matrix_c<T, M, N> result(0.0);
//    matrix_c_op<T, M, N>().fold(in_place_plus<T>(), mul<T>(), m1, m2, result);
//    return result;
//}

///**
// * Matrix division by a number
// */
//template<class T, size_t M, size_t N>
//matrix_c<T, M, N>& operator /= (matrix_c<T, M, N>& m, const T& h)
//{
//    matrix_c_op<T, M, N>().apply(in_place_div<T>(), m, h);
//    return m;
//}
//template<class T, size_t M, size_t N>
//matrix_c<T, M, N> operator / (const matrix_c<T, M, N>& m, const T& h)
//{
//    matrix_c<T, M, N> result(m);
//    return result /= h;
//}

///**
// * Matrix multiplication by a number
// */
//template<class T, size_t M, size_t N>
//matrix_c<T, M, N>& operator *= (matrix_c<T, M, N>& m, const T& h)
//{
//    matrix_c_op<T, M, N>().apply(in_place_mul<T>(), m, h);
//    return m;
//}
//template<class T, size_t M, size_t N>
//matrix_c<T, M, N> operator * (const matrix_c<T, M, N>& m, const T& h)
//{
//    matrix_c<T, M, N> result(m);
//    return result *= h;
//}
//template<class T, size_t M, size_t N>
//matrix_c<T, M, N> operator * (const T& h, const matrix_c<T, M, N>& m)
//{
//    matrix_c<T, M, N> result(m);
//    return result *= h;
//}

///**
// * Matrix transposition
// */
//template<class T, size_t M>
//matrix_c<T, M, M> transpose(const matrix_c<T, M, M>& m)
//{
//    matrix_c<T, M, M> result(m);
//    matrix_c_op<T, M, M>().transpose(result);
//    return result;
//}

///**
// * Calculates covariance matrix
//*/
//template<class T, std::size_t N>
//matrix_c<T, N, N> cov
//(
//        const vector<vector_c<T, N>>& vectors,
//        const vector<T>& w
//)
//{
//    assert(vectors.size() == w.size());

//    matrix_c<T, N, N> covMatrix(0.0);

//    typename vector<T>::const_iterator it = w.begin();

//    T total = 0.0;

//    for (const vector_c<T, N>& v:vectors)
//    {
//        total += *it;
//        matrix_c_op<T, N, N>().cov_matrix_diag(v, *(it++), covMatrix);
//    }

//    return (covMatrix/total) *
//            static_cast<double>(w.size())
//            /static_cast<double>(w.size()-1);
//}

///**
// * Calculates first principal component
// */
//template<class T, std::size_t N>
//vector_c<T, N> pc1
//(
//        const vector<vector_c<T, N>>& vectors,
//        const vector<T>& w,
//        T relTol = 1.0e-10,
//        size_t maxItter = 1000
//)
//{
//    typedef vector_c<T, N> vector_type;
//    typedef matrix_c<T, N, N> matrix_type;

//    //Make first approximation
//    vector_type eigen_vector(0.0);
//    typename vector<T>::const_iterator it = w.begin();
//    for(const vector_type& v : vectors)
//    {
//        if (eigen_vector*v < 0.0)
//            eigen_vector -= (v*v)*v * *(it++);
//        else
//            eigen_vector += (v*v)*v * *(it++);
//    }

//    if(abs(eigen_vector) == 0) return eigen_vector;
//    eigen_vector /= abs(eigen_vector);

//    matrix_type covMatrix = cov(vectors,w);
//    eigen_vector = covMatrix * eigen_vector;
//    T disp0, disp1 = abs(eigen_vector);
//    size_t iter = 0;
//    do
//    {
//        disp0 = disp1;
//        eigen_vector /= disp1;
//        eigen_vector = covMatrix * eigen_vector;
//        disp1 = abs(eigen_vector);
//    } while(std::fabs(disp1-disp0)/disp0 > relTol
//            && (++iter) != maxItter);

//    return eigen_vector;
//}

///**
// * Matrix determinant, test function
// */
//template<class T, size_t N>
//inline T det(const matrix_c<T, N, N>& m)
//{
//    matrix_c<T, N, N> _tm(m);
//    T result = 1.0;
//    if(matrix_c_op<T, N, N>().tri(_tm) == N)
//        matrix_c_op<T, N, N>().fold_diag
//                (in_place_mul<T>(),_tm,result);
//    else result = 0.0;
//    return result;
//}

///**
// * Solves a linear equation system Ax=b,
// * returns a rank of matrix A. If it is smaller
// * than a size of A then x is unidentified
// */
//template<class T, size_t N>
//inline vector_c<T, N> solve
//(
//        const matrix_c<T, N, N>& A,
//        const vector_c<T, N>& b)
//{
//    vector_c<T, N> x = b;
//    matrix_c<T, N, N> _tA(A);
//    matrix_c_op<T, N, N>().tri(_tA, x);
//    matrix_c_op<T, N, N>().backward_gauss_sweep(_tA, x);
//    return x;
//}

}

#endif // MATRIXTEMPLATE_H
