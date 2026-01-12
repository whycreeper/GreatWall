// SPDX-License-Identifier: MIT

#ifndef TREE_H
#define TREE_H

#include "params.h"
#include <stddef.h>
#include <stdint.h>

#define expand_trees GREATWALL_NAMESPACE(expand_trees)
void expand_trees(uint8_t nodes[GREATWALL_T][2 * GREATWALL_N - 1][GREATWALL_SEED_SIZE],
                  const uint8_t salt[GREATWALL_SALT_SIZE]);

#define reveal_all_but GREATWALL_NAMESPACE(reveal_all_but)
void reveal_all_but(uint8_t reveal_path[GREATWALL_LOGN][GREATWALL_SEED_SIZE],
                    const uint8_t nodes[2 * GREATWALL_N - 1][GREATWALL_SEED_SIZE],
                    size_t cover_index);

#define reconstruct_tree GREATWALL_NAMESPACE(reconstruct_tree)
void reconstruct_tree(uint8_t nodes[2 * GREATWALL_N - 2][GREATWALL_SEED_SIZE],
                      const uint8_t salt[GREATWALL_SALT_SIZE],
                      const uint8_t reveal_path[GREATWALL_LOGN][GREATWALL_SEED_SIZE],
                      size_t rep_index,
                      size_t cover_index);

#endif // TREE_H
