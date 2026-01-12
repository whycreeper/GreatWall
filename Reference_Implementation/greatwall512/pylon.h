// SPDX-License-Identifier: MIT

#ifndef PYLON_H
#define PYLON_H

#include "field.h"
#include "params.h"
#include <stdint.h>

static const GF constants[PYLON_NUM_INPUT_SBOX] =
{
  {0xb6ad000cebff7fe0,0xc2656b9df9f54a81,0xa33fa3010a33d18b,0x2deb8fad1f2471ef,0x939bab528ea20bc9,0xaf62ca6e898dbb7c,0xbfc6792535598d4f,0x36d6e6dc08a843fb},
  {0x2b5d61431f59f74f,0xd5002f5b05b73994,0x3c506361627f6bf6,0x79067f00511d021f,0x0454c531acb924cb,0xb88a986d6e24664f,0x41cd34eed42ee4cf,0x2691a69b30d73c0e},
  {0x8c5ac7f49a32b469,0x111089d07b03556d,0x9b505aadc52506c0,0xc6fb51e2a7fe2581,0xf074ee179197ff67,0x00dbc7f712ba4cc8,0xfa048f4ca013bce3,0xe6d46b12397a0fec}
};
#define GF_exp_invmer_e_3 GREATWALL_NAMESPACE(GF_exp_invmer_e_3)
void GF_exp_invmer_e_3(GF out, const GF in);


#define pylon_2_sbox_outputs GREATWALL_NAMESPACE(pylon_2_sbox_outputs)
void pylon_2_sbox_outputs(GF sbox_outputs[PYLON_NUM_INPUT_Matrix + 1], const GF pt);

#define pylon GREATWALL_NAMESPACE(pylon)
void pylon(uint8_t ct[PYLON_NUM_BYTES_FIELD],
          const uint8_t pt[PYLON_NUM_BYTES_FIELD]);

#endif // PYLON_H
