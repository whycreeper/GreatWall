// SPDX-License-Identifier: MIT

#ifndef PARAMS_H
#define PARAMS_H

#define GREATWALL_NAMESPACE(s)          GREATWALL128_N256_avx2_##s

#define SECURITY_BITS               128                  // security parameter
#define SECURITY_BYTES              (SECURITY_BITS / 8)  // byte size of security parameter

#define PYLON_NUM_BITS_FIELD         SECURITY_BITS        // number of bits in field element
#define PYLON_NUM_BYTES_FIELD        SECURITY_BYTES       // number of bytes in field element
#define PYLON_NUM_WORDS_FIELD        (SECURITY_BITS / 64) // number of 64-bit words in element
#define PYLON_NUM_BITS_WORD          64                   // number of bits in word

#define PYLON_NUM_INPUT_SBOX         3                    // number of PYLON input S-boxes
#define PYLON_NUM_INPUT_Matrix       2                    // number of PYLON input Matrix

#define GREATWALL_SALT_SIZE             SECURITY_BYTES       // byte size of salt
#define GREATWALL_SEED_SIZE             SECURITY_BYTES       // byte size of seed
#define GREATWALL_COMMIT_SIZE           (SECURITY_BYTES * 2) // byte size of commitment

#define GREATWALL_L                     PYLON_NUM_INPUT_SBOX
#define GREATWALL_T                     17                   // number of parallel repetitions (Tau)
#define GREATWALL_N                     256                   // number of MPC parties (N)
#define GREATWALL_LOGN                  8                    // log_2(N)

#endif // PARAMS_H
