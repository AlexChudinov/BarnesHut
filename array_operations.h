#ifndef ARRAY_OPERATIONS_H
#define ARRAY_OPERATIONS_H

#include <cstddef>

namespace math {
    /**
     * Defines operations over array using recursive templates
     */
    template
    <class array_type_a,
     class array_type_b,
     size_t N>
    struct array_operations
    {
        using prior = array_operations<array_type_a, array_type_b, N - 1>;

        template<class op> inline
        void umap(op Op, array_type_a& a) const
        {
            prior().umap(Op, a);
            Op(a[N]);
        }

        template<class op> inline
        void bmap(op Op, array_type_a& a, const array_type_b& b) const
        {
            prior().bmap(Op, a, b);
            Op(a[N], b[N]);
        }

        template<class add, class T> inline
        void ufold(add Add, T& a, const array_type_a& b) const
        {
            prior().ufold(Add, a, b);
            Add(a, b[N]);
        }

        template<class add, class mul, class T> inline
        void bfold(
                add Add,
                mul Mul,
                T& a,
                const array_type_a& b,
                const array_type_b& c) const
        {
            prior().bfold(Add, Mul, a, b, c);
            Add(a, Mul(b[N], c[N]));
        }
    };

    template
    <class array_type_a,
      class array_type_b>
    struct array_operations<array_type_a, array_type_b, 0>
    {
        template<class op> inline
        void umap(op Op, array_type_a& a) const { Op(a[0]); }

        template<class op> inline
        void bmap(op Op, array_type_a& a, const array_type_b& b) const
        { Op(a[0], b[0]); }

        template<class add, class T> inline
        void ufold(add Add, T& a, const array_type_a& b) const
        { Add(a, b[0]); }

        template<class add, class mul, class T> inline
        void bfold(
                add Add,
                mul Mul,
                T& a,
                const array_type_a& b,
                const array_type_b& c) const
        {
            Add(a, Mul(b[0], c[0]));
        }
    };

    /**
     * Templates implementing recursive for-loop
     */
    template<size_t begin, size_t end, bool stop> struct For;

    template<size_t begin, size_t end>
    struct For<begin, end, true>
    {
        using next = For<begin + 1, end, (begin + 1) < end>;

        template<class operation, typename ... T>
        inline void Do(operation Op, T ... args) const
        {
            Op(begin, args...);
            next().Do(Op, args...);
        }

        template<class op1, class op2, typename ... T>
        inline void Do_for_triangle(op1 Op1, op2 Op2, T ... args) const
        {
            Op1(begin);
            For<begin, end, begin < end>().Do(Op2, args ...);
            next().Do_for_triangle(Op1, Op2, args ...);
        }
    };

    template<size_t begin, size_t end>
    struct For<begin, end, false>
    {
        template<class operation, typename ... T>
        inline void Do(operation /*Op*/, T ...) const{}
        template<class op1, class op2, typename ... T>
        inline void Do_for_triangle(op1 /*Op1*/, op2 /*Op2*/, T ...) const{}
    };

//One value applicator to a hole array
#define DEF_OPERATION_WITH_VAL(value_type, op) \
    struct operation{ \
        const value_type& val_;\
        inline operation(const value_type& val) : val_(val) {} \
        inline void operator()(value_type& x) const { x op val_; } \
    };
#define DEF_OPERATION_WITH_VAL_1(type1, type2, op)\
    struct operation{ \
        const type1& val_;\
        inline operation(const type1& val) : val_(val) {} \
        inline void operator()(type2& x) const { x op val_; } \
    };
}

#endif // ARRAY_OPERATIONS_H
