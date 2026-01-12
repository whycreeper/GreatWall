// SPDX-License-Identifier: MIT

#include "pylon.h"
#include "field.h"
#include "hash.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>

// inverse Mersenne S-box with e = 3
// (2 ^ 3 - 1) ^ (-1) mod (2 ^ 128 - 1)
// = 0x49249249249249249249249249249249
// =   100 100 100 100 ... 100 1
// 42 times '100' and '1'
void GF_exp_invmer_e_3(GF out, const GF in)
{
  size_t i;
  GF t1 = {0,};
  GF table_4 = {0,}, table_4_4 = {0,}, table_4_16 = {0,};

  // table_4 = in ^ 4
  GF_sqr_s(table_4, in);
  GF_sqr_s(table_4, table_4);

  // table_4_4
  GF_sqr_s(t1, table_4);
  for(int i = 1 ; i < 3; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_4, t1, table_4);

  GF_sqr_s(t1, table_4_4);
  for(int i = 1 ; i < 6; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_4, t1, table_4_4);

  // table_4_16
  GF_sqr_s(t1, table_4_4);
  for(int i = 1 ; i < 12; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4_4);

  GF_sqr_s(t1, table_4_16);
  for(int i = 1 ; i < 24; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4_16);

  // table_4_20
  GF_sqr_s(t1, table_4_16);
  for(int i = 1 ; i < 12; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4_4);
  
  // table_4_21
  GF_sqr_s(t1, table_4_16);
  for(int i = 1 ; i < 3; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4);

  // table_4_21
  GF_sqr_s(t1, table_4_16);
  for(int i = 1 ; i < 63; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_4_16, t1, table_4_16);

  GF_sqr_s(t1, table_4_16);
  GF_mul_s(out, t1, in);
}



void pylon(uint8_t ct[PYLON_NUM_BYTES_FIELD],
          const uint8_t pt[PYLON_NUM_BYTES_FIELD])
{
  GF state;
  GF pt_GF = {0,}, ct_GF = {0,};
  GF_from_bytes(pt_GF, pt);

  GF_exp_invmer_e_3(state, pt_GF);
  GF_transposed_matmul(state, state, (const GF *)m1_transposed);
  GF_add(state, state, constants[0]);
  GF_add(state, state, pt_GF);

  GF_exp_invmer_e_3(state, state);
  GF_transposed_matmul(state, state, (const GF *)m2_transposed);
  GF_add(state, state, constants[1]);
  GF_add(state, state, pt_GF);

  GF_exp_invmer_e_3(state, state);
  GF_add(state, state, constants[2]);
  GF_add(state, state, pt_GF);

  GF_to_bytes(ct, state);
}

void pylon_2_sbox_outputs(GF sbox_outputs[PYLON_NUM_INPUT_Matrix + 1], const GF pt)
{
  GF_exp_invmer_e_3(sbox_outputs[0], pt);
  
  GF_transposed_matmul(sbox_outputs[1], sbox_outputs[0], (const GF *)m1_transposed);
  GF_add(sbox_outputs[2], sbox_outputs[1], constants[0]);
  GF_add(sbox_outputs[1], sbox_outputs[2], pt);

  GF_exp_invmer_e_3(sbox_outputs[1], sbox_outputs[1]);
}
