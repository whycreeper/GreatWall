// SPDX-License-Identifier: MIT

#include <stddef.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "../api.h"
#include "../common/cpucycles.h"
#include "../common/speed_print.h"
#include "../common/rng.h"

#define NTESTS 100

uint64_t t[NTESTS];
double t_ms[NTESTS];

static void print_results_double(const char *s, double *t, size_t n)
{
  size_t i;
  double min = t[0], max = t[0], sum = 0.0;

  for(i = 0; i < n; i++)
  {
    if(t[i] < min) min = t[i];
    if(t[i] > max) max = t[i];
    sum += t[i];
  }

  printf("%s\n", s);
  printf("median :  %.6f ms\n", t[n/2]);
  printf("average:  %.6f ms\n", sum / n);
}
#define CPU_MHZ 2419.201
int main()
{
  uint64_t freq = (uint64_t)(CPU_MHZ * 1e6); 
  int i,kkk;
  int ret = 0;

  size_t mlen = 32;
  size_t smlen = 0;
  uint8_t message[mlen];
  uint8_t m[mlen];
  uint8_t pk[CRYPTO_PUBLICKEYBYTES];
  uint8_t sk[CRYPTO_SECRETKEYBYTES];
  uint8_t sm[mlen + CRYPTO_BYTES];

  uint64_t time_begin, time_end;

  randombytes((unsigned char*)message, 32);

  for(i = 0; i < NTESTS; i++)
  {
    time_begin = cpucycles();
    ret = crypto_sign_keypair(pk, sk);
    time_end   = cpucycles();
    t[i] = time_end - time_begin;
    t_ms[i] = (double)t[i] * 1000.0 / (double)freq;
  }
  print_results("greatwall_keypair: ", t, NTESTS);//print_results_double("greatwall_keypair (ms)    : ", t_ms, NTESTS);

  for(i = 0; i < NTESTS; i++)
  {
    time_begin = cpucycles();
    ret = crypto_sign(sm, &smlen, message, mlen, sk);
    time_end   = cpucycles();
    t[i] = time_end - time_begin;
    t_ms[i] = (double)t[i] * 1000.0 / (double)freq;
  }
  print_results("greatwall_sign   : ", t, NTESTS);//print_results_double("greatwall_sign (ms)    : ", t_ms, NTESTS);

  for(i = 0; i < NTESTS; i++)
  {
    time_begin = cpucycles();
    ret = crypto_sign_open(m, &mlen, sm, smlen, pk);
    time_end   = cpucycles();
    t[i] = time_end - time_begin;
    t_ms[i] = (double)t[i] * 1000.0 / (double)freq;
  }
  print_results("greatwall_verify : ", t, NTESTS);//print_results_double("greatwall_verify (ms)    : ", t_ms, NTESTS);

  return ret;
}
