// SPDX-License-Identifier: MIT

#ifndef FIELD_H
#define FIELD_H

#include "params.h"
#include <stddef.h>
#include <stdint.h>

typedef uint64_t GF[3];

#define GF_set0 GREATWALL_NAMESPACE(GF_set0)
void GF_set0(GF a);
#define GF_copy GREATWALL_NAMESPACE(GF_copy)
void GF_copy(GF out, const GF in);
#define GF_to_bytes GREATWALL_NAMESPACE(GF_to_bytes)
void GF_to_bytes(uint8_t *out, const GF in);
#define GF_from_bytes GREATWALL_NAMESPACE(GF_from_bytes)
void GF_from_bytes(GF out, const uint8_t *in);

#define GF_add GREATWALL_NAMESPACE(GF_add)
void GF_add(GF c, const GF a, const GF b);
#define GF_mul GREATWALL_NAMESPACE(GF_mul)
void GF_mul(GF c, const GF a, const GF b);
#define GF_mul_N GREATWALL_NAMESPACE(GF_mul_N)
void GF_mul_N(GF c[GREATWALL_N], const GF a[GREATWALL_N], const GF b);
#define GF_mul_add GREATWALL_NAMESPACE(GF_mul_add)
void GF_mul_add(GF c, const GF a, const GF b);
#define GF_mul_add_N GREATWALL_NAMESPACE(GF_mul_add_N)
void GF_mul_add_N(GF c[GREATWALL_N], const GF a[GREATWALL_N], const GF b);
#define GF_sqr GREATWALL_NAMESPACE(GF_sqr)
void GF_sqr(GF c, const GF a);
#define GF_sqr_N GREATWALL_NAMESPACE(GF_sqr_N)
void GF_sqr_N(GF c[GREATWALL_N], const GF a[GREATWALL_N]);

#define GF_transposed_matmul GREATWALL_NAMESPACE(GF_transposed_matmul)
void GF_transposed_matmul(GF c, const GF a, const GF b[PYLON_NUM_BITS_FIELD]);
#define GF_transposed_matmul_add_N GREATWALL_NAMESPACE(GF_transposed_matmul_add_N)
void GF_transposed_matmul_add_N(GF c[GREATWALL_N], const GF a[GREATWALL_N],
                                const GF b[PYLON_NUM_BITS_FIELD]);

#define POLY_mul_add_N GREATWALL_NAMESPACE(POLY_mul_add_N)
void POLY_mul_add_N(GF lo[GREATWALL_N], GF hi[GREATWALL_N],
                    const GF a[GREATWALL_N], const GF b);
#define POLY_red_N GREATWALL_NAMESPACE(POLY_red_N)
void POLY_red_N(GF lo[GREATWALL_N], const GF hi[GREATWALL_N]);

#endif // FIELD_H
