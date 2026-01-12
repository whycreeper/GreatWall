// SPDX-License-Identifier: MIT

#ifndef PYLON_H
#define PYLON_H

#include "field.h"
#include "params.h"
#include <stdint.h>

static const GF constants[PYLON_NUM_INPUT_SBOX] =
{
  {0x787d77617b34735d,0xc6fd9c99eb3903ec,0xbe6581828743823a,0x075cbd797975e486,0xa7ec1a700e769b34,0x260d91e3c7a2d419},
  {0x10dbc0329b73651f,0xf167444bda04fde9,0xe39d39f2961b59e8,0x1626e60e1db93941,0xcc5a8de0be6aae71,0xce9718650984b890},
  {0xaebfee879304a4b2,0x211795397f79c994,0x81fd60f7059c6480,0x1dd3e27d1ba2644c,0xe29cf0afba62f3c5,0x05a8a37ab135ed7e}
};

#define GF_exp_invmer_e_3 GREATWALL_NAMESPACE(GF_exp_invmer_e_3)
void GF_exp_invmer_e_3(GF out, const GF in);


#define pylon_2_sbox_outputs GREATWALL_NAMESPACE(pylon_2_sbox_outputs)
void pylon_2_sbox_outputs(GF sbox_outputs[PYLON_NUM_INPUT_Matrix + 1], const GF pt);

#define pylon GREATWALL_NAMESPACE(pylon)
void pylon(uint8_t ct[PYLON_NUM_BYTES_FIELD],
          const uint8_t pt[PYLON_NUM_BYTES_FIELD]);

#endif // PYLON_H
