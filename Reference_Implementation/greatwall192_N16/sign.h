// SPDX-License-Identifier: MIT

#ifndef SIGN_H
#define SIGN_H

#include "field.h"
#include "hash.h"
#include "params.h"
#include <stddef.h>
#include <stdint.h>

typedef struct tape_t
{
  GF pt_share;
  GF t_shares[GREATWALL_L - 1];
  GF a_share;
  GF a2_share;
  GF c_share;
} tape_t;

typedef struct mult_chk_t
{
  GF pt_share;
  GF Delta_share;
  GF x_shares[GREATWALL_L];
  GF z_shares[GREATWALL_L];
} mult_chk_t;

typedef struct proof_t
{
  uint8_t reveal_path[GREATWALL_LOGN][GREATWALL_SEED_SIZE];
  uint8_t missing_commitment[GREATWALL_COMMIT_SIZE];
  uint8_t delta_pt_bytes[PYLON_NUM_BYTES_FIELD];
  uint8_t delta_c_bytes[PYLON_NUM_BYTES_FIELD];
  uint8_t delta_ts_bytes[GREATWALL_L - 1][PYLON_NUM_BYTES_FIELD];
  uint8_t missing_alpha_share_bytes[PYLON_NUM_BYTES_FIELD];
  uint8_t missing_alpha2_share_bytes[PYLON_NUM_BYTES_FIELD];
} proof_t;

typedef struct signature_t
{
  uint8_t salt[GREATWALL_SALT_SIZE];
  uint8_t h_1[GREATWALL_COMMIT_SIZE];
  uint8_t h_2[GREATWALL_COMMIT_SIZE];
  proof_t proofs[GREATWALL_T];
} signature_t;

#define pylon_mpc GREATWALL_NAMESPACE(pylon_mpc)
void pylon_mpc(mult_chk_t *mult_chk,
              const GF ct_GF, int party);

#define commit_and_expand_tape GREATWALL_NAMESPACE(commit_and_expand_tape)
void commit_and_expand_tape(tape_t *tape, uint8_t *commit,
                            const hash_instance *ctx_precom,
                            const uint8_t seed[GREATWALL_SEED_SIZE],
                            size_t rep, size_t party);

#define run_phase_1 GREATWALL_NAMESPACE(run_phase_1)
void run_phase_1(signature_t *sign,
                 uint8_t commits[GREATWALL_T][GREATWALL_N][GREATWALL_COMMIT_SIZE],
                 uint8_t nodes[GREATWALL_T][2 * GREATWALL_N - 1][GREATWALL_SEED_SIZE],
                 mult_chk_t mult_chk[GREATWALL_T][GREATWALL_N],
                 GF alpha_v_shares[GREATWALL_T][3][GREATWALL_N],
                 const uint8_t *sk, const uint8_t *m, size_t mlen);

#define run_phase_2_and_3 GREATWALL_NAMESPACE(run_phase_2_and_3)
void run_phase_2_and_3(signature_t *sign,
                       GF alpha_v_shares[GREATWALL_T][3][GREATWALL_N],
                       const mult_chk_t mult_chk[GREATWALL_T][GREATWALL_N]);

#endif // SIGN_H
