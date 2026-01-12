// SPDX-License-Identifier: MIT

#include <stdio.h>
#include <string.h>
#include "../field.h"
#include "../pylon.h"

int main()
{
  int succ = 1, i;

  // Inverted order compared to Sage!
  uint8_t pt[16] = {
    0xc0, 0x7a, 0x49, 0x16, 0xda, 0x53, 0x6d, 0x85,
    0xdd, 0x8a, 0x31, 0x3d, 0x89, 0xa4, 0x90, 0x5b
  };

const uint8_t ct_expected[16] = {
0x54, 0x63, 0x80, 0xf6,
0x21, 0xb6, 0x1e, 0x31,
0x26, 0x9b, 0xfe, 0xc3,
0xcd, 0xa9, 0x3b, 0x16
};

  uint8_t ct[16] = {0,};

  pylon(ct, pt);

  printf("PLAINTEXT                 : ");
  for (i = (int)sizeof(pt) - 1; i >= 0; i--)
  {
    printf("%02x", pt[i]);
  }
  printf("\n");

  printf("CIPHERTEXT                : ");
  for (i = (int)sizeof(ct) - 1; i >= 0; i--)
  {
    printf("%02x", ct[i]);
  }
  printf("\n");

  if (memcmp(ct_expected, ct, sizeof(GF)))
  {
    succ = 0;
    printf("[-] ct != ct_expected\n");
  }

  if (succ)
    printf("Passed!\n");
  else
    printf("Error!\n");
  return 0;
}
