// SPDX-License-Identifier: MIT

#include "api.h"
#include "pylon.h"
#include "field.h"
#include "hash.h"
#include "params.h"
#include "sign.h"
#include "tree.h"
#include "common/crypto_declassify.h"
#include "common/rng.h"
#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
void GF_print(const GF x)
{
    uint8_t buf[PYLON_NUM_BYTES_FIELD];
    GF_to_bytes(buf, x);

    printf("GF bytes = ");
    for (int i = 0; i < PYLON_NUM_BYTES_FIELD; i++)
        printf("%02x", buf[i]);
    printf("\n");
}
void commit_and_expand_tape_x4(tape_t *tapes, uint8_t *commits,
                               const hash_instance_x4 *ctx_precom,
                               const uint8_t *seeds, size_t rep, size_t party)
{
  hash_instance_x4 ctx;
  uint8_t bufs[4][GREATWALL_SEED_SIZE + 2];
  const uint8_t *in_ptrs[4] = {bufs[0], bufs[1], bufs[2], bufs[3]};
  uint8_t *out_ptrs[4];

  memcpy(&ctx, ctx_precom, sizeof(hash_instance_x4));

  bufs[0][0] = (uint8_t)(rep);
  bufs[1][0] = (uint8_t)(rep);
  bufs[2][0] = (uint8_t)(rep);
  bufs[3][0] = (uint8_t)(rep);

  bufs[0][1] = (uint8_t)(party + 0);
  bufs[1][1] = (uint8_t)(party + 1);
  bufs[2][1] = (uint8_t)(party + 2);
  bufs[3][1] = (uint8_t)(party + 3);

  memcpy(&bufs[0][2], seeds + 0 * GREATWALL_SEED_SIZE, GREATWALL_SEED_SIZE);
  memcpy(&bufs[1][2], seeds + 1 * GREATWALL_SEED_SIZE, GREATWALL_SEED_SIZE);
  memcpy(&bufs[2][2], seeds + 2 * GREATWALL_SEED_SIZE, GREATWALL_SEED_SIZE);
  memcpy(&bufs[3][2], seeds + 3 * GREATWALL_SEED_SIZE, GREATWALL_SEED_SIZE);

  hash_update_x4(&ctx, in_ptrs, GREATWALL_SEED_SIZE + 2);
  hash_final_x4(&ctx);

  out_ptrs[0] = commits + 0 * GREATWALL_COMMIT_SIZE;
  out_ptrs[1] = commits + 1 * GREATWALL_COMMIT_SIZE;
  out_ptrs[2] = commits + 2 * GREATWALL_COMMIT_SIZE;
  out_ptrs[3] = commits + 3 * GREATWALL_COMMIT_SIZE;

  hash_squeeze_x4(&ctx, out_ptrs, GREATWALL_COMMIT_SIZE);

  out_ptrs[0] = (uint8_t *)(tapes);
  out_ptrs[1] = (uint8_t *)(tapes + 1);
  out_ptrs[2] = (uint8_t *)(tapes + 2);
  out_ptrs[3] = (uint8_t *)(tapes + 3);

  hash_squeeze_x4(&ctx, out_ptrs, sizeof(tape_t));
}

void pylon_mpc_N(mult_chk_N_t *mult_chk, const GF ct_GF)
{
  GF_transposed_matmul_add_N(mult_chk->Delta_share, mult_chk->x_shares[0], m1_transposed);
  GF_add(mult_chk->Delta_share[0], mult_chk->Delta_share[0], constants[0]);

  GF_sqr_N(mult_chk->z_shares[0], (const GF *)mult_chk->x_shares[0]);
  for (size_t i = 1; i < 3; i++)
  {
    GF_sqr_N(mult_chk->z_shares[0], (const GF *)mult_chk->z_shares[0]); 
  }

  GF_sqr_N(mult_chk->z_shares[1], (const GF *)mult_chk->x_shares[1]);
  for (size_t i = 1; i < 3; i++)
  {
    GF_sqr_N(mult_chk->z_shares[1], mult_chk->z_shares[1]); 
  }

  GF ct_c2;
  GF_add(ct_c2, ct_GF, constants[2]);
  for(size_t j = 0 ; j < GREATWALL_N; j++){
    GF_copy(mult_chk->z_shares[2][j], mult_chk->pt_share[j]);
  }
  GF_add(mult_chk->z_shares[2][0], mult_chk->z_shares[2][0],ct_c2);
  for (size_t i = 1; i <= 3; i++)
  {
    GF_sqr_N(mult_chk->z_shares[2], mult_chk->z_shares[2]);
  }

  for(size_t j = 0 ; j < GREATWALL_N; j++){
    GF_copy(mult_chk->x_shares[2][j],mult_chk->pt_share[j]);
  }
  GF_transposed_matmul_add_N(mult_chk->x_shares[2], mult_chk->x_shares[1], m2_transposed);
  GF_add(mult_chk->x_shares[2][0], mult_chk->x_shares[2][0], constants[1]);

  GF_mul_add_N(mult_chk->z_shares[2], mult_chk->x_shares[2], ct_c2);
}

// committing to the seeds and the execution views of the parties
void run_phase_1(signature_t *sign,
                 uint8_t commits[GREATWALL_T][GREATWALL_N][GREATWALL_COMMIT_SIZE],
                 uint8_t nodes[GREATWALL_T][2 * GREATWALL_N - 1][GREATWALL_SEED_SIZE],
                 mult_chk_N_t mult_chk[GREATWALL_T],
                 GF alpha_v_shares[GREATWALL_T][3][GREATWALL_N],
                 const uint8_t *sk, const uint8_t *m, size_t mlen)
{
  GF pt_GF = {0,}, ct_GF = {0,};
  GF_from_bytes(pt_GF, sk);
  GF_from_bytes(ct_GF, sk + PYLON_NUM_BYTES_FIELD);

  // message pre-hashing
  hash_instance ctx;
  hash_init_prefix(&ctx, HASH_PREFIX_0);
  hash_update(&ctx, sk + PYLON_NUM_BYTES_FIELD, PYLON_NUM_BYTES_FIELD);
  hash_update(&ctx, m, mlen);
  hash_final(&ctx);

  uint8_t mu[GREATWALL_COMMIT_SIZE];
  hash_squeeze(&ctx, mu, GREATWALL_COMMIT_SIZE);

  GF sbox_outputs[GREATWALL_L];//sbox_outputs[2] not sbox output
  pylon_2_sbox_outputs(sbox_outputs, pt_GF);
  // GF_print(sbox_outputs[0]);
  // GF_print(sbox_outputs[1]);
  // GF_print(sbox_outputs[2]);
  // GF tmpp,tmpp2;
  // GF_copy(tmpp, sbox_outputs[1]);
  // for(int i=1;i<7;i++)GF_mul(tmpp,tmpp,sbox_outputs[1]);
  // GF_add(tmpp, tmpp, constants[0]);
  // GF_add(tmpp, tmpp, pt_GF);
  // GF_print(tmpp);printf("should equal to\n");
  // GF_print(ct_GF);
  // GF_add(tmpp,ct_GF, constants[2]);
  // GF_add(tmpp,tmpp, pt_GF);
  // GF_copy(tmpp2, tmpp);
  // for(int i=1;i<7;i++)GF_mul(tmpp2,tmpp2,tmpp);
  // GF_add(tmpp2,tmpp2,pt_GF);
  // GF_add(tmpp2,tmpp2,constants[1]);
  // GF_print(tmpp2);printf("should equal to\n");
  // GF_transposed_matmul(tmpp2, sbox_outputs[1], m2_transposed);
  // GF_print(tmpp2);
  // GF_add(tmpp,tmpp,pt_GF);
  // GF_add(tmpp,tmpp,constants[1]);
  // GF_exp_invmer_e_3(tmpp,tmpp);
  // GF_add(tmpp,tmpp,pt_GF);
  // GF_add(tmpp,tmpp,constants[2]);
  // GF_print(tmpp);


  // generate per-signature randomness
  uint8_t random[SECURITY_BYTES];
  randombytes(random, SECURITY_BYTES);

  // generate salt
  hash_init_prefix(&ctx, HASH_PREFIX_3);
  hash_update(&ctx, sk, PYLON_NUM_BYTES_FIELD);
  hash_update(&ctx, mu, GREATWALL_COMMIT_SIZE);
  hash_update(&ctx, random, SECURITY_BYTES);
  hash_final(&ctx);
  hash_squeeze(&ctx, sign->salt, GREATWALL_SALT_SIZE);

  // generate root seeds and expand seed trees
  for (size_t rep = 0; rep < GREATWALL_T; rep++)
  {
    hash_squeeze(&ctx, nodes[rep][0], GREATWALL_SEED_SIZE);
  }
  expand_trees(nodes, sign->salt);

  // hash_instance for h_1
  hash_init_prefix(&ctx, HASH_PREFIX_1);
  hash_update(&ctx, mu, GREATWALL_COMMIT_SIZE);
  hash_update(&ctx, sign->salt, GREATWALL_SALT_SIZE);

  hash_instance_x4 ctx_precom;
  hash_init_prefix_x4(&ctx_precom, HASH_PREFIX_5);
  hash_update_x4_1(&ctx_precom, sign->salt, GREATWALL_SALT_SIZE);
  
  memset(mult_chk, 0, sizeof(mult_chk_N_t) * GREATWALL_T);

  for (size_t rep = 0; rep < GREATWALL_T; rep++)
  {
    // initialize adjustment values
    tape_t tapes[4];
    tape_t delta;
    memset(&delta, 0, sizeof(tape_t));


    for (size_t party = 0; party < GREATWALL_N; party += 4)
    {
      commit_and_expand_tape_x4(tapes, commits[rep][party], &ctx_precom,
                             nodes[rep][party + GREATWALL_N - 1], rep, party);
      hash_update(&ctx, commits[rep][party], 4 * GREATWALL_COMMIT_SIZE);
      
      for (size_t j = 0; j < 4 && party + j <= GREATWALL_N - 1; j++)
      {
        GF_add(delta.pt_share, delta.pt_share, tapes[j].pt_share);
        GF_add(delta.t_shares[0], delta.t_shares[0], tapes[j].t_shares[0]);
        GF_add(delta.t_shares[1], delta.t_shares[1], tapes[j].t_shares[1]);
        GF_add(delta.a_share, delta.a_share, tapes[j].a_share);
        GF_add(delta.a2_share, delta.a2_share, tapes[j].a2_share);
        GF_add(delta.c_share, delta.c_share, tapes[j].c_share);

        if (party + j == GREATWALL_N - 1)
        {
          GF_add(delta.pt_share, delta.pt_share, pt_GF);// true pt
          GF_add(delta.t_shares[0], delta.t_shares[0], sbox_outputs[0]);// true t0
          GF_add(delta.t_shares[1], delta.t_shares[1], sbox_outputs[1]);// true t1
          
          GF_mul_add(delta.c_share, pt_GF, delta.a_share);
          GF_mul_add(delta.c_share, sbox_outputs[2], delta.a2_share);// true c

          GF_to_bytes(sign->proofs[rep].delta_pt_bytes, delta.pt_share);// delta_pt into sign
          GF_to_bytes(sign->proofs[rep].delta_ts_bytes[0], delta.t_shares[0]);// delta_t0 into sign
          GF_to_bytes(sign->proofs[rep].delta_ts_bytes[1], delta.t_shares[1]);// delta_t1 into sign
          GF_to_bytes(sign->proofs[rep].delta_c_bytes, delta.c_share);// delta_c into sign

          GF_add(tapes[j].pt_share, delta.pt_share, tapes[j].pt_share);//adjust
          GF_add(tapes[j].t_shares[0], delta.t_shares[0], tapes[j].t_shares[0]);
          GF_add(tapes[j].t_shares[1], delta.t_shares[1], tapes[j].t_shares[1]);
          GF_add(tapes[j].c_share, delta.c_share, tapes[j].c_share);
        }
        
        GF_copy(mult_chk[rep].pt_share[party + j], tapes[j].pt_share);
        GF_copy(mult_chk[rep].x_shares[0][party + j], tapes[j].t_shares[0]);
        GF_copy(mult_chk[rep].x_shares[1][party + j], tapes[j].t_shares[1]);

        GF_add(alpha_v_shares[rep][0][party + j], tapes[j].a_share, tapes[j].a2_share);
        GF_copy(alpha_v_shares[rep][1][party + j], tapes[j].a2_share);
        GF_copy(alpha_v_shares[rep][2][party + j], tapes[j].c_share);
      }
    }
    // for (size_t party = 0; party < GREATWALL_N; party++){
      
      // GF_print(alpha_v_shares[rep][0][party]);
      // GF_print(alpha_v_shares[rep][1][party]);
      // GF_print(alpha_v_shares[rep][2][party]);
    // }

    pylon_mpc_N(&mult_chk[rep], ct_GF);

    // GF tmp2;
    //GF_set0(tmp2);
    // for (size_t party = 0; party < GREATWALL_N; party++){
    //   GF_add(tmp2,tmp2,mult_chk[rep].pt_share[party]);
    // }
    // GF_print(tmp2);printf("=\n");
    // GF_print(pt_GF);
    // GF_set0(tmp2);
    // for (size_t party = 0; party < GREATWALL_N; party++){
    //   GF_add(tmp2,tmp2,mult_chk[rep].x_shares[0][party]);
    // }
    // GF_print(tmp2);printf("=\n");
    // GF_print(sbox_outputs[0]);
    // GF_set0(tmp2);
    // for (size_t party = 0; party < GREATWALL_N; party++){
    //   GF_add(tmp2,tmp2,mult_chk[rep].x_shares[1][party]);
    // }
    // GF_print(tmp2);printf("=\n");
    // GF_print(sbox_outputs[1]);
    // GF_set0(tmp2);
    // for (size_t party = 0; party < GREATWALL_N; party++){
    //   GF_add(tmp2,tmp2,mult_chk[rep].x_shares[2][party]);
    // }
    // GF_print(tmp2);printf("=\n");
    // GF_transposed_matmul(tmp2, sbox_outputs[1], m2_transposed);
    // GF_add(tmp2,tmp2,constants[1]);
    // GF_add(tmp2,tmp2,pt_GF);
    // GF_print(tmp2);
    // GF_set0(tmp2);
    // for (size_t party = 0; party < GREATWALL_N; party++){
    //   GF_add(tmp2,tmp2,mult_chk[rep].z_shares[0][party]);
    // }
    // GF_print(tmp2);printf("=\n");
    // GF_copy(tmp2, sbox_outputs[0]);
    // for (int i=1;i<8;i++){
    //   GF_mul(tmp2,tmp2,sbox_outputs[0]);
    // }
    // GF_print(tmp2);
    // GF_set0(tmp2);
    // for (size_t party = 0; party < GREATWALL_N; party++){
    //   GF_add(tmp2,tmp2,mult_chk[rep].z_shares[1][party]);
    // }
    // GF_print(tmp2);printf("=\n");
    // GF_copy(tmp2, sbox_outputs[1]);
    // for (int i=1;i<8;i++){
    //   GF_mul(tmp2,tmp2,sbox_outputs[1]);
    // }
    // GF_print(tmp2);
    // GF_set0(tmp2);
    // for (size_t party = 0; party < GREATWALL_N; party++){
    //   GF_add(tmp2,tmp2,mult_chk[rep].z_shares[2][party]);
    // }
    // GF_print(tmp2);printf("=\n");
    // GF_copy(tmp2, pt_GF);GF_add(tmp2,tmp2,ct_GF);GF_add(tmp2,tmp2,constants[2]);
    // GF_sqr(tmp2,tmp2);GF_sqr(tmp2,tmp2);GF_sqr(tmp2,tmp2);
    // GF tmp_h;
    // GF_add(tmp_h, ct_GF, constants[2]);
    // GF tmp_s;
    // GF_set0(tmp_s);
    // for (size_t party = 0; party < GREATWALL_N; party++){
    //   GF_add(tmp_s,tmp_s,mult_chk[rep].x_shares[2][party]);
    // }
    // GF_mul_add(tmp2, tmp_h, tmp_s);
    // GF_print(tmp2);
    // GF_mul(tmp2, tmp_s, pt_GF);
    // GF_print(tmp2);



    // NOTE: depend on the order of values in proof_t
    hash_update(&ctx, sign->proofs[rep].delta_pt_bytes,
                PYLON_NUM_BYTES_FIELD * (GREATWALL_L + 1));
  }

  // commit to salt, (all commitments of parties' seeds,
  // delta_pt, delta_t, delta_c) for all repetitions
  hash_final(&ctx);
  hash_squeeze(&ctx, sign->h_1, GREATWALL_COMMIT_SIZE);

}

void run_phase_2_and_3(signature_t *sign,
                       GF alpha_v_shares[GREATWALL_T][3][GREATWALL_N],
                       const mult_chk_N_t mult_chk[GREATWALL_T])
{
  hash_instance ctx_e;
  hash_init(&ctx_e);
  hash_update(&ctx_e, sign->h_1, GREATWALL_COMMIT_SIZE);
  hash_final(&ctx_e);

  hash_instance ctx;
  hash_init_prefix(&ctx, HASH_PREFIX_2);
  hash_update(&ctx, sign->h_1, GREATWALL_COMMIT_SIZE);
  hash_update(&ctx, sign->salt, GREATWALL_SALT_SIZE);


  for (size_t rep = 0; rep < GREATWALL_T; rep++)
  {
    GF alpha = {0,};
    GF alpha2 = {0,};
    GF epsilons[GREATWALL_L];
    hash_squeeze(&ctx_e, (uint8_t *)epsilons, sizeof(epsilons));

    GF alpha_v_shares_hi[3][GREATWALL_N];
    memset(alpha_v_shares_hi, 0, sizeof(alpha_v_shares_hi));

    POLY_mul_add_N(alpha_v_shares[rep][0], alpha_v_shares_hi[0],
                   mult_chk[rep].x_shares[0], epsilons[0]);
    POLY_mul_add_N(alpha_v_shares[rep][2], alpha_v_shares_hi[2],
                   mult_chk[rep].z_shares[0], epsilons[0]);

    POLY_mul_add_N(alpha_v_shares[rep][1], alpha_v_shares_hi[1],
                   mult_chk[rep].x_shares[1], epsilons[1]);
    POLY_mul_add_N(alpha_v_shares[rep][2], alpha_v_shares_hi[2],
                   mult_chk[rep].z_shares[1], epsilons[1]);

    POLY_mul_add_N(alpha_v_shares[rep][0], alpha_v_shares_hi[0],
                   mult_chk[rep].x_shares[2], epsilons[2]);
    POLY_mul_add_N(alpha_v_shares[rep][2], alpha_v_shares_hi[2],
                   mult_chk[rep].z_shares[2], epsilons[2]);

    POLY_red_N(alpha_v_shares[rep][0], (const GF *)alpha_v_shares_hi[0]);
    POLY_red_N(alpha_v_shares[rep][1], (const GF *)alpha_v_shares_hi[1]);
    POLY_red_N(alpha_v_shares[rep][2], (const GF *)alpha_v_shares_hi[2]);
    
    for (size_t party = 0; party < GREATWALL_N; party++)
    {
      // GF_print(mult_chk[rep].pt_share[party]);
      // GF_print(mult_chk[rep].Delta_share[party]);
      // GF_print(alpha_v_shares[rep][2][party]);
        GF_add(alpha_v_shares[rep][0][party], 
        alpha_v_shares[rep][0][party], alpha_v_shares[rep][1][party]);
      GF_add(alpha, alpha, alpha_v_shares[rep][0][party]);
      GF_add(alpha2, alpha2, alpha_v_shares[rep][1][party]);
    }
    
    
    // GF_print(alpha);;
    // GF_print(alpha2);


    GF_mul_add_N(alpha_v_shares[rep][2], mult_chk[rep].pt_share, alpha);
    GF_mul_add_N(alpha_v_shares[rep][2], mult_chk[rep].Delta_share, alpha2);

    // GF v_sum={0,};
    // for (size_t party = 0; party < GREATWALL_N; party++){
    //   GF_add(v_sum,v_sum,alpha_v_shares[rep][2][party]);
    // }
    // printf("%lld:",rep);GF_print(v_sum);

    hash_update(&ctx, (const uint8_t *)alpha_v_shares[rep],
                PYLON_NUM_BYTES_FIELD * 3 * GREATWALL_N);
  
  }
  hash_final(&ctx);
  hash_squeeze(&ctx, sign->h_2, GREATWALL_COMMIT_SIZE);

}

////////////////////////////////////////////////////////////////////////////////
int crypto_sign_keypair(uint8_t *pk, uint8_t *sk)
{
  if (!pk || !sk)
  {
    return -1;
  }

  randombytes(sk, PYLON_NUM_BYTES_FIELD);

  pylon(pk, sk);
  memcpy(sk + PYLON_NUM_BYTES_FIELD, pk, PYLON_NUM_BYTES_FIELD);

  return 0;
}

int crypto_sign_signature(uint8_t *sig, size_t *siglen,
        const uint8_t *m, size_t mlen,
        const uint8_t *sk)
{
  hash_instance ctx;
  signature_t *sign = (signature_t *)sig;

  //////////////////////////////////////////////////////////////////////////
  // Phase 1: Committing to the seeds and the execution views of parties. //
  //////////////////////////////////////////////////////////////////////////

  // nodes for seed trees
  uint8_t nodes[GREATWALL_T][2 * GREATWALL_N - 1][GREATWALL_SEED_SIZE];

  // commitments for seeds
  uint8_t commits[GREATWALL_T][GREATWALL_N][GREATWALL_COMMIT_SIZE];

  // multiplication check inputs
  mult_chk_N_t mult_chk[GREATWALL_T];

  // multiplication check outputs
  GF alpha_v_shares[GREATWALL_T][3][GREATWALL_N];

  // commitments for phase 1
  run_phase_1(sign, commits, 
    (uint8_t (*)[2 * GREATWALL_N - 1][GREATWALL_SEED_SIZE])nodes, mult_chk,
     alpha_v_shares, sk, m, mlen);

  /////////////////////////////////////////////////////////////////
  // Phase 2, 3: Challenging and committing to the simulation of //
  //             the multiplication checking protocol.           //
  /////////////////////////////////////////////////////////////////

  // compute the commitment of phase 3
  run_phase_2_and_3(sign, alpha_v_shares, mult_chk);

  //////////////////////////////////////////////////////
  // Phase 4: Challenging views of the MPC protocols. //
  //////////////////////////////////////////////////////

  hash_init(&ctx);
  hash_update(&ctx, sign->h_2, GREATWALL_COMMIT_SIZE);
  hash_final(&ctx);

  uint8_t indices[GREATWALL_T]; // GREATWALL_N <= 256
  hash_squeeze(&ctx, indices, GREATWALL_T);
  for (size_t rep = 0; rep < GREATWALL_T; rep++)
  {
    indices[rep] &= (1 << GREATWALL_LOGN) - 1;
  }

  //////////////////////////////////////////////////////
  // Phase 5: Opening the views of the MPC protocols. //
  //////////////////////////////////////////////////////

  crypto_declassify(indices, sizeof(indices));
  for (size_t rep = 0; rep < GREATWALL_T; rep++)
  {
    size_t i_bar = indices[rep];
    reveal_all_but(sign->proofs[rep].reveal_path,
                   (const uint8_t (*)[GREATWALL_SEED_SIZE])nodes[rep], i_bar);
    memcpy(sign->proofs[rep].missing_commitment, commits[rep][i_bar],
           GREATWALL_COMMIT_SIZE);
    GF_to_bytes(sign->proofs[rep].missing_alpha_share_bytes,
                alpha_v_shares[rep][0][i_bar]);
    GF_to_bytes(sign->proofs[rep].missing_alpha2_share_bytes,
                alpha_v_shares[rep][1][i_bar]);
  }
  *siglen = CRYPTO_BYTES;

  return 0;
}

int crypto_sign(uint8_t *sm, size_t *smlen,
        const uint8_t *m, size_t mlen,
        const uint8_t *sk)
{
  crypto_sign_signature(sm + mlen, smlen, m, mlen, sk);

  memcpy(sm, m, mlen);
  *smlen += mlen;

  return 0;
}

int crypto_sign_verify(const uint8_t *sig, size_t siglen,
        const uint8_t *m, size_t mlen,
        const uint8_t *pk)
{

  if (siglen != CRYPTO_BYTES)
  {
    return -1;
  }

  const signature_t *sign = (const signature_t *)sig;

  GF ct_GF = {0,};
  GF_from_bytes(ct_GF, pk);
  // derive the binary matrix and the vector from the initial vector

  hash_instance ctx_e, ctx_h1, ctx_h2;

  // indices = Expand(h_2)
  hash_init(&ctx_e);
  hash_update(&ctx_e, sign->h_2, GREATWALL_COMMIT_SIZE);
  hash_final(&ctx_e);

  uint8_t indices[GREATWALL_T]; // GREATWALL_N <= 256
  hash_squeeze(&ctx_e, indices, GREATWALL_T);
  for (size_t rep = 0; rep < GREATWALL_T; rep++)
  {
    indices[rep] &= (1 << GREATWALL_LOGN) - 1;
  }

  // epsilons = Expand(h_1)
  hash_init(&ctx_e);
  hash_update(&ctx_e, sign->h_1, GREATWALL_COMMIT_SIZE);
  hash_final(&ctx_e);

  // message pre-hashing
  uint8_t mu[GREATWALL_COMMIT_SIZE];
  hash_init_prefix(&ctx_h1, HASH_PREFIX_0);
  hash_update(&ctx_h1, pk, PYLON_NUM_BYTES_FIELD);
  hash_update(&ctx_h1, m, mlen);
  hash_final(&ctx_h1);
  hash_squeeze(&ctx_h1, mu, GREATWALL_COMMIT_SIZE);

  // ready for computing h_1' and h_2'
  hash_init_prefix(&ctx_h1, HASH_PREFIX_1);
  hash_update(&ctx_h1, mu, GREATWALL_COMMIT_SIZE);
  hash_update(&ctx_h1, sign->salt, GREATWALL_SALT_SIZE);

  hash_init_prefix(&ctx_h2, HASH_PREFIX_2);
  hash_update(&ctx_h2, sign->h_1, GREATWALL_COMMIT_SIZE);
  hash_update(&ctx_h2, sign->salt, GREATWALL_SALT_SIZE);

  hash_instance_x4 ctx_precom;
  hash_init_prefix_x4(&ctx_precom, HASH_PREFIX_5);
  hash_update_x4_1(&ctx_precom, sign->salt, GREATWALL_SALT_SIZE);

  for (size_t rep = 0; rep < GREATWALL_T; rep++)
  {
    size_t i_bar = indices[rep];
    uint8_t nodes[2 * GREATWALL_N - 2][GREATWALL_SEED_SIZE];
    memset(nodes, 0, sizeof(nodes));

    reconstruct_tree(nodes, sign->salt, sign->proofs[rep].reveal_path,
                     rep, i_bar);

    mult_chk_N_t mult_chk;
    memset(&mult_chk, 0, sizeof(mult_chk_N_t));

    GF epsilons[GREATWALL_L];
    hash_squeeze(&ctx_e, (uint8_t *)epsilons, sizeof(epsilons));

    GF alpha = {0,}, alpha2 = {0,};
    GF alpha_v_shares[3][GREATWALL_N], alpha_v_shares_hi[3][GREATWALL_N];
    memset(alpha_v_shares_hi, 0, sizeof(alpha_v_shares_hi));
    GF_set0(alpha_v_shares[2][i_bar]);

    
    for (size_t party = 0; party < GREATWALL_N; party += 4)
    {
      tape_t tapes[4];
      uint8_t commits[4][GREATWALL_COMMIT_SIZE];
      commit_and_expand_tape_x4(tapes, commits[0], &ctx_precom,
                                nodes[party + GREATWALL_N - 2], rep, party);

      if (party / 4 == i_bar / 4)
      {
        memcpy(commits[i_bar % 4], sign->proofs[rep].missing_commitment,
               GREATWALL_COMMIT_SIZE);
      }
      hash_update(&ctx_h1, commits[0], 4 * GREATWALL_COMMIT_SIZE);

      if (party == GREATWALL_N - 4)
      {
        GF temp = {0,};

        GF_from_bytes(temp, sign->proofs[rep].delta_pt_bytes);
        GF_add(tapes[3].pt_share, tapes[3].pt_share, temp);

        GF_from_bytes(temp, sign->proofs[rep].delta_ts_bytes[0]);
        GF_add(tapes[3].t_shares[0], tapes[3].t_shares[0], temp);
        GF_from_bytes(temp, sign->proofs[rep].delta_ts_bytes[1]);
        GF_add(tapes[3].t_shares[1], tapes[3].t_shares[1], temp);

        GF_from_bytes(temp, sign->proofs[rep].delta_c_bytes);
        GF_add(tapes[3].c_share, tapes[3].c_share, temp);
      }

      for (size_t j = 0; j < 4 && party + j <= GREATWALL_N - 1; j++)
      {
        GF_copy(mult_chk.pt_share[party + j], tapes[j].pt_share);
        GF_copy(mult_chk.x_shares[0][party + j], tapes[j].t_shares[0]);
        GF_copy(mult_chk.x_shares[1][party + j], tapes[j].t_shares[1]);
        GF_add(alpha_v_shares[0][party + j], tapes[j].a_share, tapes[j].a2_share);
        GF_copy(alpha_v_shares[1][party + j], tapes[j].a2_share);
        GF_copy(alpha_v_shares[2][party + j], tapes[j].c_share);
      }
    }
    for (size_t party = 0; party < GREATWALL_N; party++){
      
      // GF_print(alpha_v_shares[0][party]);
      // GF_print(alpha_v_shares[1][party]);
      // GF_print(alpha_v_shares[2][party]);
    }

    pylon_mpc_N(&mult_chk, ct_GF);
    hash_update(&ctx_h1, sign->proofs[rep].delta_pt_bytes,
                PYLON_NUM_BYTES_FIELD * (GREATWALL_L + 1));

    POLY_mul_add_N(alpha_v_shares[0], alpha_v_shares_hi[0],
                   (const GF *)mult_chk.x_shares[0], epsilons[0]);
    POLY_mul_add_N(alpha_v_shares[2], alpha_v_shares_hi[2],
                   (const GF *)mult_chk.z_shares[0], epsilons[0]);

    POLY_mul_add_N(alpha_v_shares[1], alpha_v_shares_hi[1],
                   (const GF *)mult_chk.x_shares[1], epsilons[1]);
    POLY_mul_add_N(alpha_v_shares[2], alpha_v_shares_hi[2],
                   (const GF *)mult_chk.z_shares[1], epsilons[1]);

    POLY_mul_add_N(alpha_v_shares[0], alpha_v_shares_hi[0],
                   (const GF *)mult_chk.x_shares[2], epsilons[2]);
    POLY_mul_add_N(alpha_v_shares[2], alpha_v_shares_hi[2],
                   (const GF *)mult_chk.z_shares[2], epsilons[2]);

    POLY_red_N(alpha_v_shares[0], (const GF *)alpha_v_shares_hi[0]);
    POLY_red_N(alpha_v_shares[1], (const GF *)alpha_v_shares_hi[1]);
    POLY_red_N(alpha_v_shares[2], (const GF *)alpha_v_shares_hi[2]);
    GF_from_bytes(alpha_v_shares[0][i_bar],
                  sign->proofs[rep].missing_alpha_share_bytes);
    GF_add(alpha_v_shares[0][i_bar],alpha_v_shares[0][i_bar],
                  sign->proofs[rep].missing_alpha2_share_bytes);
    GF_from_bytes(alpha_v_shares[1][i_bar],
                  sign->proofs[rep].missing_alpha2_share_bytes);
    GF_set0(alpha);GF_set0(alpha2);
    for (size_t party = 0; party < GREATWALL_N; party++)
    {
      // GF_print(mult_chk.pt_share[party]);
      // GF_print(mult_chk.Delta_share[party]);
      // GF_print(alpha_v_shares[2][party]);
      GF_add(alpha_v_shares[0][party], 
        alpha_v_shares[0][party], alpha_v_shares[1][party]);
      GF_add(alpha, alpha, alpha_v_shares[0][party]);
      GF_add(alpha2, alpha2, alpha_v_shares[1][party]);
    }
    
    // alpha is opened, so we can finish calculating v_share
    // GF_print(alpha);
    // GF_print(alpha2);
    GF_mul_add_N(alpha_v_shares[2], (const GF *)mult_chk.pt_share, alpha);
    GF_mul_add_N(alpha_v_shares[2], (const GF *)mult_chk.Delta_share, alpha2);

    GF_set0(alpha);
    for (size_t party = 0; party < GREATWALL_N; party++)
    {
      GF_add(alpha, alpha, alpha_v_shares[2][party]);
    }
    GF_add(alpha_v_shares[2][i_bar], alpha, alpha_v_shares[2][i_bar]);
    // GF_print(alpha_v_shares[2][10]);
    hash_update(&ctx_h2, (const uint8_t *)alpha_v_shares,
                PYLON_NUM_BYTES_FIELD * 3 * GREATWALL_N);
  }

  uint8_t h_1_prime[GREATWALL_COMMIT_SIZE];
  hash_final(&ctx_h1);
  hash_squeeze(&ctx_h1, h_1_prime, GREATWALL_COMMIT_SIZE);

  uint8_t h_2_prime[GREATWALL_COMMIT_SIZE];
  hash_final(&ctx_h2);
  hash_squeeze(&ctx_h2, h_2_prime, GREATWALL_COMMIT_SIZE);

  if (memcmp(h_1_prime, sign->h_1, GREATWALL_COMMIT_SIZE) != 0 ||
      memcmp(h_2_prime, sign->h_2, GREATWALL_COMMIT_SIZE) != 0)
  {
    return -1;
  }

  return 0;
}

int crypto_sign_open(uint8_t *m, size_t *mlen,
        const uint8_t *sm, size_t smlen,
        const uint8_t *pk)
{
  if (smlen < CRYPTO_BYTES)
  {
    return -1;
  }

  const size_t message_len = smlen - CRYPTO_BYTES;
  const uint8_t *message = sm;
  const uint8_t *signature = sm + message_len;
  if (crypto_sign_verify(signature, CRYPTO_BYTES, message, message_len, pk))
  {
    return -1;
  }

  memmove(m, message, message_len);
  *mlen = message_len;

  return 0;
}
