#ifndef MATRIXTEMPLATE_H
#define MATRIXTEMPLATE_H

#include "vectortemplate.h"

/**
 * Template defines simple math matrix operations
 */

namespace math {
    template <class T, uint64_t nRows, uint64_t nCols>
    struct Matrix : Vector<Vector<T, nCols>, nRows>
    {
        typedef Vector<Vector<T, nCols>, nRows> base_type;

        /**
         * Acces to a matrix components
         */
        template<uint64_t M, uint64_t N> T& comp()
        { return static_cast<base_type*>(this)->comp<M>().comp<N>(); }
        template<uint64_t M, uint64_t N> const T& comp() const
        {
            return static_cast<const base_type*>(this)->
                        comp<M>().comp<N>();
        }


    };
}

#endif // MATRIXTEMPLATE_H
