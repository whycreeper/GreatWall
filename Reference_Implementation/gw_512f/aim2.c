// SPDX-License-Identifier: MIT

#include "aim2.h"
#include "field.h"
#include "hash.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>

// inverse Mersenne S-box with e = 3
// (2 ^ 3 - 1) ^ (-1) mod (2 ^ 512 - 1)
// = 0x49249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249249
// =   100 100 100 100 ... 100 1
// 170 times '100' and '1'
void GF_exp_invmer_e_3(GF out, const GF in)
{
  GF t1 = {0,};
  GF table_4 = {0,}, table_4_16 = {0,}, table_4_17_2 = {0,}, table_4_17_8 = {0,};

  // table_4 = in ^ 4
  GF_sqr_s(table_4, in);
  GF_sqr_s(table_4, table_4);

  // table_4_16
  GF_sqr_s(t1, table_4);
  for(int i = 1 ; i < 3; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4);

  GF_sqr_s(t1, table_4_16);
  for(int i = 1 ; i < 6; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4_16);

  GF_sqr_s(t1, table_4_16);
  for(int i = 1 ; i < 12; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4_16);

  GF_sqr_s(t1, table_4_16);
  for(int i = 1 ; i < 24; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4_16);

  // table_4_17
  GF_sqr_s(t1, table_4_16);
  for(int i = 1 ; i < 3; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4);

  // table_4_17_2
  GF_sqr_s(t1, table_4_16);
  for(int i = 1 ; i < 51; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_17_2, t1, table_4_16);

  // table_4_17_8
  GF_sqr_s(t1, table_4_17_2);
  for(int i = 1 ; i < 102; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_17_8, t1, table_4_17_2);

  GF_sqr_s(t1, table_4_17_8);
  for(int i = 1 ; i < 204; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_17_8, t1, table_4_17_8);

  // table_4_17_10
  GF_sqr_s(t1, table_4_17_8);
  for(int i = 1 ; i < 102; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(t1, t1, table_4_17_2);

  GF_sqr_s(t1, t1);
  GF_mul_s(out, t1, in);
}


void generate_matrices_L_and_U(
        GF matrix_L[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD],
        GF matrix_U[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD],
        const uint8_t iv[AIM2_IV_SIZE])
{
  uint8_t buf[AIM2_NUM_BYTES_FIELD];
  uint64_t ormask, lmask, umask;
  hash_instance ctx;
  GF temp = {0,};

  // initialize hash
  hash_init(&ctx);
  hash_update(&ctx, iv, AIM2_IV_SIZE);
  hash_final(&ctx);

  for (size_t num = 0; num < AIM2_NUM_INPUT_Matrix; num++)
  {
    for (size_t row = 0; row < AIM2_NUM_BITS_FIELD; row++)
    {
      hash_squeeze(&ctx, buf, AIM2_NUM_BYTES_FIELD);
      GF_from_bytes(temp, buf);

      ormask = ((uint64_t)1) << (row % 64);
      lmask = ((uint64_t)-1) << (row % 64);
      umask = ~lmask;

      size_t inter = row / 64;
      size_t col_word;
      for (col_word = 0; col_word < inter; col_word++)
      {
        // L is zero, U is full
        matrix_L[num][row][col_word] = 0;
        matrix_U[num][row][col_word] = temp[col_word];
      }
      matrix_L[num][row][inter] = (temp[inter] & lmask) | ormask;
      matrix_U[num][row][inter] = (temp[inter] & umask) | ormask;
      for (col_word = inter + 1; col_word < AIM2_NUM_WORDS_FIELD; col_word++)
      {
        // L is full, U is zero
        matrix_L[num][row][col_word] = temp[col_word];
        matrix_U[num][row][col_word] = 0;
      }
    }
  }

  hash_ctx_release(&ctx);
}

void generate_matrix_LU(GF matrix_A[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD],
                        const uint8_t iv[AIM2_IV_SIZE])
{
  GF matrix_L[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD];
  GF matrix_U[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD];

  generate_matrices_L_and_U(matrix_L, matrix_U, iv);

  for (size_t num = 0; num < AIM2_NUM_INPUT_Matrix; num++)
  {
    for (size_t i = 0; i < AIM2_NUM_BITS_FIELD; i++)
    {
      GF_transposed_matmul(matrix_A[num][i], matrix_U[num][i],
                           (const GF *)matrix_L[num]);
    }
  }
}

void aim2(uint8_t ct[AIM2_NUM_BYTES_FIELD],
          const uint8_t pt[AIM2_NUM_BYTES_FIELD],
          const uint8_t iv[AIM2_IV_SIZE])
{
  GF matrix_L[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD];
  GF matrix_U[AIM2_NUM_INPUT_Matrix][AIM2_NUM_BITS_FIELD];

  GF state;
  GF pt_GF = {0,}, ct_GF = {0,};
  GF_from_bytes(pt_GF, pt);

  generate_matrices_L_and_U(matrix_L, matrix_U, iv);

  GF_exp_invmer_e_3(state, pt_GF);
  GF_transposed_matmul(state, state, (const GF *)matrix_U[0]);
  GF_transposed_matmul(state, state, (const GF *)matrix_L[0]);
  GF_add(state, state, constants[0]);
  GF_add(state, state, pt_GF);

  GF_exp_invmer_e_3(state, state);
  GF_transposed_matmul(state, state, (const GF *)matrix_U[1]);
  GF_transposed_matmul(state, state, (const GF *)matrix_L[1]);
  GF_add(state, state, constants[1]);
  GF_add(state, state, pt_GF);

  GF_exp_invmer_e_3(state, state);
  GF_add(state, state, pt_GF);

  GF_to_bytes(ct, state);
}

void aim2_2_sbox_outputs(GF sbox_outputs[AIM2_NUM_INPUT_Matrix], const GF pt, const GF matrix_A[AIM2_NUM_BITS_FIELD])
{
  GF_exp_invmer_e_3(sbox_outputs[0], pt);
  
  GF_transposed_matmul(sbox_outputs[1], sbox_outputs[0], (const GF *)matrix_A);
  GF_add(sbox_outputs[1], sbox_outputs[1], constants[0]);
  GF_add(sbox_outputs[1], sbox_outputs[1], pt);

  GF_exp_invmer_e_3(sbox_outputs[1], sbox_outputs[1]);
}
