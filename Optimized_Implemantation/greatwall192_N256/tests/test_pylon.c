// SPDX-License-Identifier: MIT

#include <stdio.h>
#include <string.h>
#include "../field.h"
#include "../pylon.h"

int main()
{
  int succ = 1, i;

  // Inverted order compared to Sage!
  uint8_t pt[24] ={
  0x37, 0x3b, 0x19, 0xcf, 0xc8, 0xc7, 0xff, 0xf4,
  0x79, 0xb9, 0xd9, 0xcf, 0x35, 0x7b, 0x5e, 0xce,
  0xe2, 0x66, 0x7c, 0x95, 0x97, 0x26, 0x74, 0x5b
};

  const uint8_t ct_expected[24] = {
  0x99, 0xfd, 0xd3, 0x33, 0xcb, 0x06, 0x36, 0x48,
  0x14, 0xe5, 0xa4, 0xce, 0xb2, 0xd2, 0xfa, 0x30,
  0x53, 0x47, 0xad, 0x2f, 0x6a, 0x2d, 0xf4, 0xa8
};
  uint8_t ct[24] = {0,};

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
