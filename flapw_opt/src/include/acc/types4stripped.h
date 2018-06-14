#ifndef _TEST_TYPES
#define _TEST_TYPES

#include "cx_double.h"

class myTMat {
    public:
        int aa_size;
        int ab_size;
        int bb_size;
        const cx_double *aa;
        const cx_double *ab;
        const cx_double *bb;

        myTMat( ) {};
        void set( int _aa_size, int _ab_size, int _bb_size,
                const cx_double *src_aa, const cx_double *src_ab, const cx_double *src_bb )
        {
            aa_size = _aa_size;
            ab_size = _ab_size;
            bb_size = _bb_size;
            aa = src_aa;
            ab = src_ab;
            bb = src_bb;
        }

        ~myTMat( )
        {
        }
};

#endif // _TEST_TYPES
