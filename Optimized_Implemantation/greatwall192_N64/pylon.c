// SPDX-License-Identifier: MIT

#include "pylon.h"
#include "field.h"
#include "hash.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>

// inverse Mersenne S-box  with e = 5
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
  GF_sqr(table_a, in);
  GF_sqr(table_a, table_a);
  GF_mul(table_a, table_a, in);
  GF_sqr(table_a, table_a);

  // table_a_2
  GF_sqr(t1, table_a);
  for(int i = 1 ; i < 5; i++){
    GF_sqr(t1, t1);
  }
  GF_mul(table_a_2, t1, table_a);

  // table_a_4
  GF_sqr(t1, table_a_2);
  for(int i = 1 ; i < 10; i++){
    GF_sqr(t1, t1);
  }
  GF_mul(table_a_4, t1, table_a_2);
  
  // table_a_32
  GF_sqr(t1, table_a_4);
  for(int i = 1 ; i < 20; i++){
    GF_sqr(t1, t1);
  }
  GF_mul(table_a_32, t1, table_a_4);

  GF_sqr(t1, table_a_32);
  for(int i = 1 ; i < 40; i++){
    GF_sqr(t1, t1);
  }
  GF_mul(table_a_32, t1, table_a_32);

  GF_sqr(t1, table_a_32);
  for(int i = 1 ; i < 80; i++){
    GF_sqr(t1, t1);
  }
  GF_mul(table_a_32, t1, table_a_32);

  // table_a_36
  GF_sqr(t1, table_a_32);
  for(int i = 1 ; i < 20; i++){
    GF_sqr(t1, t1);
  }
  GF_mul(table_a_32, t1, table_a_4);

  // table_a_38
  GF_sqr(t1, table_a_32);
  for(int i = 1 ; i < 10; i++){
    GF_sqr(t1, t1);
  }
  GF_mul(table_a_32, t1, table_a_2);

  GF_sqr(table_a_32, table_a_32);
  GF_sqr(table_a_32, table_a_32);
  GF_mul(out, table_a_32, in);
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
  GF_add(state, state, pt_GF);
  GF_add(state, state, constants[2]);

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
