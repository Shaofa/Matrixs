#ifndef __MATRIX_H__
#define __MATRIX_H__

#define MAX_MAT_SIZE 289
#define ROW 0
#define COL 1
#define ELEMENT(M,col,i,j)	*((M)+(col)*(i)+(j))

#include "typedefs.h"

#ifdef __cplusplus
extern "C"
{
#endif

    int16 m_add(Uint16 row, Uint16 col, float32*left, float32*right, float32*out);

    int16 m_sub(Uint16 row, Uint16 col, float32*left, float32*right, float32*out);

    int16 m_mul_scale(Uint16 row, Uint16 col, float32*M, float32 a);

    int16 m_mul(float32*left, Uint16 row_left, Uint16 col_left, float32*right, Uint16 row_right, Uint16 col_right, float32*out);
    
    float32 m_det(float32* M, Uint16 n);

    int16 m_joint(float32* M, Uint16 n, float32*ans);

    int16 m_inv(float32 *M, Uint16 n, float32* M_inv);

    int16 m_trans(float32*M, Uint16 row, Uint16 col, float32 *M_tran);

    int16 m_repmat(Uint16 row, Uint16 col, float32* M, Uint16 num, Uint16 ROW_or_COL, float32* out);
	
    int16 m_copy(const float32* src, float32* dst, const Uint16 row, const Uint16 col);
    
    int16 m_chol(Uint16 row, Uint16 col, float32* M);

    int16 m_qr(const float32* mat, const Uint16 row, const Uint16 col, float32* out);

    int16 m_cholupdate(float32* mat, const float32* vec, const Uint16 row, const Uint16 col, const char sign);

    void m_print(const float32* mat, const Uint16 m, const Uint16 n, const char *prompt);

#ifdef __cplusplus
}
#endif

#endif
