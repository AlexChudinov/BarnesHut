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

template<class T, size_t m, size_t n> struct matrix_c;
template<class T, size_t m, size_t n> struct proxy_matrix_row;
template<class T, size_t m, size_t n> struct const_proxy_matrix_row;
template<class T, size_t m, size_t n> struct proxy_matrix_col;
template<class T, size_t m, size_t n> struct const_proxy_matrix_col;
template<class T, size_t M, size_t N> struct proxy_matrix_diag;
template<class T, size_t M, size_t N> struct const_proxy_matrix_diag;

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

    inline proxy_matrix_diag<T, m, n> diag();

    inline const_proxy_matrix_col<T, m, n> column(size_t idx_col) const;

    inline const_proxy_matrix_row<T, m, n> row(size_t idx_row) const;

    inline const_proxy_matrix_diag<T, m, n> diag() const;
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
};

template<class T, size_t m, size_t n>
struct const_proxy_matrix_col
{
    using matrix = matrix_c<T, m, n>;
    const matrix& A_;
    size_t idx_col_;

    inline const_proxy_matrix_col(const matrix& A, size_t idx_col)
        :
          A_(A), idx_col_(idx_col){}

    inline const_proxy_matrix_col(const proxy_matrix_col<T, m, n>& col)
        :
          A_(col.A_), idx_col_(col.idx_col_){}

    inline const T& operator[](size_t idx_row) const
    { return A_[idx_row][idx_col_]; }
};

template<class T, size_t m, size_t n> inline
proxy_matrix_col<T, m, n> matrix_c<T, m, n>::column(size_t idx_col)
{
    return proxy_matrix_col<T, m, n>(*this, idx_col);
}

template<class T, size_t m, size_t n> inline
const_proxy_matrix_col<T, m, n> matrix_c<T, m, n>::column(size_t idx_col) const
{
    return const_proxy_matrix_col<T, m, n>(*this, idx_col);
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
};

template<class T, size_t m, size_t n>
struct const_proxy_matrix_row
{
    using matrix = matrix_c<T, m, n>;
    const matrix& A_;
    size_t idx_row_;

    inline const_proxy_matrix_row(const matrix& A, size_t idx_row)
        :
          A_(A), idx_row_(idx_row)
    {}

    inline const_proxy_matrix_row(const proxy_matrix_row<T, m, n>& row)
        :
          A_(row.A_), idx_row_(row.idx_row_)
    {}

    inline const T& operator[](size_t idx_col) const
    { return A_[idx_row_][idx_col]; }
};

template<class T, size_t m, size_t n> inline
proxy_matrix_row<T, m, n> matrix_c<T, m, n>::row(size_t idx_row)
{
    return proxy_matrix_row<T, m, n>(*this, idx_row);
}

template<class T, size_t m, size_t n> inline
const_proxy_matrix_row<T, m, n> matrix_c<T, m, n>::row(size_t idx_row) const
{
    return const_proxy_matrix_row<T, m, n>(*this, idx_row);
}

template<class T, size_t M, size_t N>
struct proxy_matrix_diag
{
    using matrix = matrix_c<T, M, N>;
    matrix& A_;
    inline proxy_matrix_diag(matrix& A) : A_(A){}
    inline T& operator[](size_t diag_idx){ return A_[diag_idx][diag_idx]; }
};

template<class T, size_t M, size_t N>
struct const_proxy_matrix_diag
{
    using matrix = matrix_c<T, M, N>;
    const matrix& A_;
    inline const_proxy_matrix_diag(const matrix& A) : A_(A){}
    inline const T& operator[](size_t diag_idx) const { return A_[diag_idx][diag_idx]; }
};

template<class T, size_t M, size_t N> inline
proxy_matrix_diag<T, M, N> matrix_c<T, M, N>::diag()
{
    return proxy_matrix_diag<T, M, N>(*this);
}

template<class T, size_t M, size_t N> inline
const_proxy_matrix_diag<T, M, N> matrix_c<T, M, N>::diag() const
{
    return const_proxy_matrix_diag<T, M, N>(*this);
}

/**
 * Folds a matrix row with a matrix column
 */
template<class T, size_t m, size_t k, size_t n>
T operator *
(
        const_proxy_matrix_row<T, m, k>& row,
        const_proxy_matrix_col<T, k, n>& col)
{
    using row_vector = const_proxy_matrix_row<T, m, k>;
    using col_vector = const_proxy_matrix_col<T, k, n>;
    T res(0.0);
    math::array_operations<row_vector, col_vector, k-1> op;
    op.bfold(math::in_place_plus<T>(), std::multiplies<T>(), res, row, col);
    return res;
}
/**
 * Folds a matrix row with a matrix column (rvalue variant)
 */
template<class T, size_t m, size_t k, size_t n>
T operator *
(
        const_proxy_matrix_row<T, m, k>&& row,
        const_proxy_matrix_col<T, k, n>&& col)
{
    using row_vector = const_proxy_matrix_row<T, m, k>;
    using col_vector = const_proxy_matrix_col<T, k, n>;
    T res(0.0);
    math::array_operations<row_vector, col_vector, k-1> op;
    op.bfold(math::in_place_plus<T>(), std::multiplies<T>(), res, row, col);
    return res;
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

    vector_col res(0.0);

    struct fold_matrix_row
    {
        const matrix& A_;
        const vector_row& x_;
        vector_col& y_;

        inline fold_matrix_row
                (
                    const matrix& A,
                    const vector_row& x,
                    vector_col& y)
            :
              A_(A), x_(x), y_(y) {}

        inline void operator()(size_t idx_row)
        {
            y_[idx_row] = A_[idx_row]*x_;
        }
    };

    math::For<0, m, true>().Do(fold_matrix_row(A, x, res));

    return res;
}

/**
 * Matrix multiplication
 */
template<class T, size_t m, size_t k, size_t n>
matrix_c<T, m, n> operator *
(
        const matrix_c<T, m, k>& m1,
        const matrix_c<T, k, n>& m2)
{
    using lmatrix_type = matrix_c<T, m, k>;
    using rmatrix_type = matrix_c<T, k, n>;
    using result_type  = matrix_c<T, m, n>;
    using result_row   = proxy_matrix_row<T, m, n>;

    result_type res(0.0);

    struct one_result_row
    {
        const lmatrix_type& m1_;
        const rmatrix_type& m2_;
        result_type& res_;
        inline one_result_row
                (
                    const lmatrix_type& m1,
                    const rmatrix_type& m2,
                    result_type& res) : m1_(m1), m2_(m2), res_(res)
        {}

        inline void operator () (size_t row_idx)
        {
            struct one_result_elem
            {
                const lmatrix_type& m1_;
                const rmatrix_type& m2_;
                result_row res_;

                inline one_result_elem
                        (
                            const lmatrix_type& m1,
                            const rmatrix_type& m2,
                            result_row res)
                    :m1_(m1), m2_(m2), res_(res)
                {}

                inline void operator()(size_t col_idx)
                {
                    res_[col_idx] = m1_.row(res_.idx_row_) * m2_.column(col_idx);
                }
            };

            math::For<0, n, true>()
                    .Do(one_result_elem(m1_, m2_, res_.row(row_idx)));
        }
    };

    math::For<0, m, true>().Do(one_result_row(m1, m2, res));

    return res;
}

/**
 * Matrix division by a number
 */
template<class T, size_t m, size_t n>
matrix_c<T, m, n>& operator /= (matrix_c<T, m, n>& M, const T& h)
{
    using row_type = vector_c<T, n>;
    using matrix = matrix_c<T, m, n>;
    DEF_OPERATION_WITH_VAL_1(T, row_type, /=);
    math::array_operations<matrix, matrix, m-1> op;
    op.umap(operation(h), M);
    return M;
}
template<class T, size_t m, size_t n>
matrix_c<T, m, n> operator / (const matrix_c<T, m, n>& M, const T& h)
{
    matrix_c<T, m, n> result(M);
    return result /= h;
}

/**
 * Matrix multiplication by a number
 */
template<class T, size_t m, size_t n>
matrix_c<T, m, n>& operator *= (matrix_c<T, m, n>& M, const T& h)
{
    using row_type = vector_c<T, n>;
    using matrix = matrix_c<T, m, n>;
    DEF_OPERATION_WITH_VAL_1(T, row_type, *=);
    math::array_operations<matrix, matrix, m-1> op;
    op.umap(operation(h), M);
    return M;
}
template<class T, size_t m, size_t n>
matrix_c<T, m, n> operator * (const matrix_c<T, m, n>& M, const T& h)
{
    matrix_c<T, m, n> result(M);
    return result *= h;
}
template<class T, size_t m, size_t n>
matrix_c<T, m, n> operator * (const T& h, const matrix_c<T, m, n>& M)
{
    matrix_c<T, m, n> result(M);
    return result *= h;
}

/**
 * Matrix transposition
 */
template<class T, size_t m, size_t n>
matrix_c<T, n, m> transpose(const matrix_c<T, m, n>& M)
{
    using matrix = matrix_c<T, m, n>;
    using matrix_result = matrix_c<T, n, m>;
    using row_vector = const_proxy_matrix_row<T, m, n>;
    using col_vector = proxy_matrix_col<T, n, m>;

    matrix_result result;

    struct set_elements
    {
        matrix_result& res_;
        const matrix& M_;
        inline set_elements(matrix_result& res, const matrix& M)
            : res_(res), M_(M)
        {}
        inline void operator()(size_t res_row_idx)
        {
            math::array_operations<col_vector, row_vector, n-1> op;
            col_vector col(res_, res_row_idx);
            op.bmap(set_val<T>(), col, M_.row(res_row_idx));
        }
    };

    math::For<0, m, true>().Do(set_elements(result, M));

    return result;
}

/**
* Calculates covariance matrix
*/
template<class T, std::size_t N>
matrix_c<T, N, N> cov
(
        const vector<vector_c<T, N>>& vectors,
        const vector<T>& w,
        const vector_c<T, N>& v_mean = 0.0
)
{
    using matrix   = matrix_c<T, N, N>;
    using vector_c = vector_c<T, N>;

    assert(vectors.size() == w.size());

    matrix covMatrix(0.0);
    size_t diag_idx = 0;

    struct set_lower_matrix_triangle
    {
        matrix& m_;
        const vector_c& v_;
        size_t& diag_idx_;
        inline set_lower_matrix_triangle(matrix& m, const vector_c& v, size_t& diag_idx)
            : m_(m), v_(v), diag_idx_(diag_idx) {}
        inline void operator()(size_t row_idx, T w)
        {
            m_[row_idx][diag_idx_] += w * v_[row_idx] * v_[diag_idx_];
        }
    };

    struct set_upper_matrix_triangle
    {
        matrix& m_;
        const vector_c& v_;
        const size_t& diag_idx_;
        inline set_upper_matrix_triangle(matrix& m, const vector_c& v, const size_t& diag_idx)
            : m_(m), v_(v), diag_idx_(diag_idx) {}
        inline void operator()(size_t row_idx)
        {
            m_[diag_idx_][row_idx] = m_[row_idx][diag_idx_];
        }
    };

    struct set_diag_idx
    {
        size_t& diag_idx_;
        inline set_diag_idx(size_t& diag_idx) : diag_idx_(diag_idx) {}
        inline void operator()(size_t diag_idx){ diag_idx_ = diag_idx; }
    };

    typename vector<T>::const_iterator it = w.begin();

    T total = 0.0;

    for (const vector_c& v:vectors)
    {
        total += *it;
        vector_c vv = v - v_mean;
        math::For<0, N, true>()
                .Do_for_triangle(
                    set_diag_idx(diag_idx),
                    set_lower_matrix_triangle(covMatrix, vv, diag_idx),*it);
        math::For<0, N, true>()
                .Do_for_triangle(
                    set_diag_idx(diag_idx),
                    set_upper_matrix_triangle(covMatrix, vv, diag_idx));
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
    using vector_c = vector_c<T, N>;
    using matrix = matrix_c<T, N, N>;

    //Calculate cloud center
    vector_c v0(0.0);
    math::mean(vectors, w, v0);

    //Make first approximation
    vector_c eigen_vector(0.0);
    typename vector<T>::const_iterator it = w.begin();
    for(const vector_c& v : vectors)
    {
        vector_c vv = v - v0;
        if (eigen_vector*vv < 0.0)
            eigen_vector -= (vv*vv)*vv * *(it++);
        else
            eigen_vector += (vv*vv)*vv * *(it++);
    }

    if(abs(eigen_vector) == 0) return eigen_vector;
    eigen_vector /= abs(eigen_vector);

    matrix covMatrix = cov(vectors, w, v0);
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

/**
 * Matrix determinant
 */
template<class T, size_t N>
inline T det(const matrix_c<T, N, N>& m)
{
    using matrix = matrix_c<T, N, N>;

    matrix tri(m);
    size_t row_num;
    T coef;

    struct set_column
    {
        matrix& tri_;
        const size_t& row_num_;
        inline do_triangle
                (
                    matrix& tri,
                    const size_t& row_num)
            : tri_(tri), row_num_(row_num){}
        inline void operator()(size_t col_num, const T& coef)
        {
            tri_[row_num_][col_num] -= (coef*tri_[row_num_-1][col_num]);
        }
    };

    struct set_row
    {
        matrix& tri_;
        size_t row_num_;
        T& coef_;
        inline set_row(matrix& tri, size_t row_num, T& coef)
            : tri_(tri), row_num_(row_num), coef_(coef) {}
        inline void operator()(size_t row_num1)
        {
            coef_ = tri_[row_num1][row_num_]/tri_[row_num_][row_num_];
            tri_[row_num1][row_num_] = 0.0;
        }
    };

    struct diag_prod
    {
        const matrix& m_;
        inline diag_prod(const matrix& m) : m_(m) {}
        inline void operator()(size_t diag_idx, T& res) const
        {
            res *= m_[diag_idx][diag_idx];
        }
    };

    For<1, N, true>()
            .Do_for_triangle(
                set_row(tri,row_num,coef),
                do_triangle(tri,row_num), std::cref(coef));

    T result = 1.0;

    For<0, N, true>().Do(diag_prod(tri),std::ref(result));

    return result;
}

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
