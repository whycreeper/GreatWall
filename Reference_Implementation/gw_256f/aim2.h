// SPDX-License-Identifier: MIT

#ifndef AIM2_H
#define AIM2_H

#include "field.h"
#include "params.h"
#include <stdint.h>

static const GF constants[AIM2_NUM_INPUT_Matrix] =
{
  {0x24a19947b3916cf7,0xba7c9045f12c7f99,0xb8e1afed6a267e96,0x2ffd72dbd01adfb7},
  {0x0d95748f728eb658,0xa458fea3f4933d7e,0x636920d871574e69,0x0801f2e2858efc16}
};

#define GF_exp_invmer_e_3 AIMER_NAMESPACE(GF_exp_invmer_e_3)
void GF_exp_invmer_e_3(GF out, const GF in);

#define generate_matrices_L_and_U AIMER_NAMESPACE(generate_matrices_L_and_U)
void generate_matrices_L_and_U(
        GF matrix_L[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD],
        GF matrix_U[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD],
        const uint8_t iv[AIM2_IV_SIZE]);

#define generate_matrix_LU AIMER_NAMESPACE(generate_matrix_LU)
void generate_matrix_LU(GF matrix_A[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD],
                        const uint8_t iv[AIM2_IV_SIZE]);

#define aim2_2_sbox_outputs AIMER_NAMESPACE(aim2_2_sbox_outputs)
void aim2_2_sbox_outputs(GF sbox_outputs[AIM2_NUM_INPUT_Matrix], const GF pt, const GF matrix_A[AIM2_NUM_BITS_FIELD]);

#define aim2 AIMER_NAMESPACE(aim2)
void aim2(uint8_t ct[AIM2_NUM_BYTES_FIELD],
          const uint8_t pt[AIM2_NUM_BYTES_FIELD],
          const uint8_t iv[AIM2_IV_SIZE]);

#endif // AIM2_H
