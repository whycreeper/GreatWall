// SPDX-License-Identifier: MIT

#include <stdio.h>
#include <string.h>
#include "../field.h"
#include "../pylon.h"

int main()
{
  int succ = 1, i;

  // Inverted order compared to Sage!
  uint8_t pt[32] = {
    0x17, 0xa3, 0xab, 0xc5, 0x1f, 0xaa, 0x51, 0x1a,
    0xbb, 0x46, 0x1b, 0x66, 0x6a, 0x63, 0x94, 0x19,
    0x3c, 0x0d, 0x86, 0x30, 0x49, 0x56, 0x3f, 0xc6,
    0x9b, 0x14, 0xc0, 0xec, 0xd8, 0xec, 0xad, 0x1e
  };


const uint8_t ct_expected[32] = {
0x04, 0xf6, 0xf9, 0x1f, 0x38, 0x38, 0x3a, 0xf6,
0xda, 0xe4, 0xba, 0xa7, 0x7a, 0xa9, 0x11, 0x12,
0x5e, 0x75, 0xa0, 0x07, 0x8d, 0x44, 0xb2, 0x8e,
0xd7, 0x07, 0x1b, 0x29, 0xf1, 0xc5, 0x10, 0x6a
};


  uint8_t ct[32] = {0,};

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

  // For  extracting test data -------------------------------------------------
  // printf("CIPHERTEXT(little endian) : ");
  // for (i = 0; i < (int)sizeof(ct) - 1; i++)
  // {
  //   printf("0x%02x, ", ct[i]);
  // }
  // printf("0x%02x\n", ct[i]);
  // For  extracting test data -------------------------------------------------

  return 0;
}
