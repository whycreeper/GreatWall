// SPDX-License-Identifier: MIT

#include "pylon.h"
#include "field.h"
#include "hash.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>

// inverse Mersenne S-box with e = 5
// (2 ^ 5 - 1) ^ (-1) mod (2 ^ 384 - 1)
// = 0x108421084210842108421084210842108421084210842108421084210842108421084210842108421084210842108421
// =   10842 10842 10842 ... 10842 1
// 19 times '10842' and '1'
// 76 times '10000' and '1'
void GF_print(const GF x)
{
    uint8_t buf[PYLON_NUM_BYTES_FIELD];
    GF_to_bytes(buf, x);

    printf("GF bytes = ");
    for (int i = 0; i < PYLON_NUM_BYTES_FIELD; i++)
        printf("%02x", buf[i]);
    printf("\n");
}
void GF_exp_invmer_e_3(GF out, const GF in)
{
  GF t1 = {0,};
  GF table_16 = {0,}, table_16_4 = {0,}, table_16_8 = {0,}, table_16_64 = {0,};

  // table_16 = in ^ 16
  GF_sqr_s(table_16, in);
  GF_sqr_s(table_16, table_16);
  GF_sqr_s(table_16, table_16);
  GF_sqr_s(table_16, table_16);
  // table_16_4
  GF_sqr_s(t1, table_16);
  for(int i = 1 ; i < 5; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_16_4, t1, table_16);

  GF_sqr_s(t1, table_16_4);
  for(int i = 1 ; i < 10; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_16_4, t1, table_16_4);

  // table_16_8
  GF_sqr_s(t1, table_16_4);
  for(int i = 1 ; i < 20; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_16_8, t1, table_16_4);

  // table_16_64
  GF_sqr_s(t1, table_16_8);
  for(int i = 1 ; i < 40; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_16_64, t1, table_16_8);
  
  GF_sqr_s(t1, table_16_64);
  for(int i = 1 ; i < 80; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_16_64, t1, table_16_64);

  GF_sqr_s(t1, table_16_64);
  for(int i = 1 ; i < 160; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(table_16_64, t1, table_16_64);

  // table_16_76
  GF_sqr_s(t1, table_16_64);
  for(int i = 1 ; i < 40; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(t1, t1, table_16_8);

  GF_sqr_s(t1, t1);
  for(int i = 1 ; i < 20; i++){
    GF_sqr_s(t1, t1);
  }
  GF_mul_s(t1, t1, table_16_4);

  GF_sqr_s(t1, t1);
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
