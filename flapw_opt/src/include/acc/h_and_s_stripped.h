#ifndef _H_AND_S_STRIPPED
#define _H_AND_S_STRIPPED

#include <string>
#include "cx_double.h"
#include "types4stripped.h"

/*double create_HS2_stripped( int L, int a, int G,*/
/*cx_double *A, cx_double *B,*/
/*cx_double *H, cx_double *S,*/
/*myTMat *T,*/
/*int *l_maxs,*/
/*int max_lmax,*/
/*double *u_norms, */
/*std::string method );*/
void create_HS2_stripped( int L, int a, int G,
                   cx_double *A, cx_double *B,
                   cx_double *H, cx_double *S,
                   myTMat *T,
                   int *l_maxs,
                   int max_lmax,
                   double *u_norms, 
                   std::string method, int tile_size );
void create_HS2_stripped_v2( int L, int a, int G,
        cx_double *A, cx_double *B,
        /*cx_double *A_ptr, cx_double *A_pinned,*/
        /*cx_double *B_ptr, cx_double *B_pinned,*/
                   cx_double *H, cx_double *S,
                   myTMat *T,
                   int *l_maxs,
                   int max_lmax,
                   double *u_norms, 
                   std::string method,
                   cx_double *W );

#endif // _H_AND_S_STRIPPED
