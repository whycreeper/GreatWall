// SPDX-License-Identifier: MIT

#include "aim2.h"
#include "field.h"
#include "hash.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>
void GF_print(const GF x)
{
    printf("GF = 0x%016llx%016llx\n",
           (unsigned long long)x[1],   // high 64 bits
           (unsigned long long)x[0]);  // low 64 bits
}
// inverse Mersenne S-box with e = 5
// (2 ^ 5 - 1) ^ (-1) mod (2 ^ 192 - 1)
// = 0x5294a5294a5294a5294a5294a5294a5294a5294a5294a529
// =   01010 01010 01010 ... 01010 01
// 38 times '01010' and '01'
void GF_exp_invmer_e_3(GF out, const GF in)
{
  size_t i;
  GF t1 = {0,};
  GF table_a = {0,}, table_a_2 = {0,}, table_a_4 = {0,}, table_a_32 = {0,};

  // table_a = in ^ a
  GF_sqr_s(table_a, in);
  GF_sqr_s(table_a, table_a);
  GF_mul_s(table_a, table_a, in);
  GF_sqr_s(table_a, table_a);

  // table_a_2
  GF_sqr_s(t1, table_a);
  for(int i = 1 ; i < 5; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_a_2, t1, table_a);

  // table_a_4
  GF_sqr_s(t1, table_a_2);
  for(int i = 1 ; i < 10; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_a_4, t1, table_a_2);
  
  // table_a_32
  GF_sqr_s(t1, table_a_4);
  for(int i = 1 ; i < 20; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_a_32, t1, table_a_4);

  GF_sqr_s(t1, table_a_32);
  for(int i = 1 ; i < 40; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_a_32, t1, table_a_32);

  GF_sqr_s(t1, table_a_32);
  for(int i = 1 ; i < 80; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_a_32, t1, table_a_32);

  // table_a_36
  GF_sqr_s(t1, table_a_32);
  for(int i = 1 ; i < 20; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_a_32, t1, table_a_4);

  // table_a_38
  GF_sqr_s(t1, table_a_32);
  for(int i = 1 ; i < 10; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_a_32, t1, table_a_2);

  GF_sqr_s(table_a_32, table_a_32);
  GF_sqr_s(table_a_32, table_a_32);
  GF_mul_s(out, table_a_32, in);
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
