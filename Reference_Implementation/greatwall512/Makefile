# SPDX-License-Identifier: MIT

CC = gcc
CFLAGS += -I. -O3 -Wall -Wextra -march=native -fomit-frame-pointer -w
NISTFLAGS = -Wno-sign-compare -Wno-unused-but-set-variable -Wno-unused-parameter -Wno-unused-result

EXECUTABLE_TESTPYLON  = tests/test_pylon
EXECUTABLE_TESTSIGN  = tests/test_sign
EXECUTABLE_TESTSPEED = tests/test_speed
EXECUTABLE_KAT       = PQCgenKAT_sign

.PHONY: all

all: $(EXECUTABLE_TESTPYLON) $(EXECUTABLE_TESTSIGN) $(EXECUTABLE_TESTSPEED) $(EXECUTABLE_KAT)

.c.o:
	$(CC) -c $(CFLAGS) $< -o $@

$(EXECUTABLE_TESTPYLON): $(EXECUTABLE_TESTPYLON).c pylon.c field.c hash.c common/fips202.c
	$(CC) $(CFLAGS) $^ -o $@

$(EXECUTABLE_TESTSIGN): $(EXECUTABLE_TESTSIGN).c pylon.c field.c hash.c sign.c tree.c common/aes.c common/fips202.c common/rng.c
	$(CC) $(CFLAGS) $^ -o $@

$(EXECUTABLE_TESTSPEED): $(EXECUTABLE_TESTSPEED).c pylon.c field.c hash.c sign.c tree.c common/aes.c common/speed_print.c common/cpucycles.c common/fips202.c common/rng.c
	$(CC) $(CFLAGS) $^ -o $@

$(EXECUTABLE_KAT): common/$(EXECUTABLE_KAT).c pylon.c field.c hash.c sign.c tree.c common/aes.c common/fips202.c common/rng.c
	$(CC) $(CFLAGS) $^ -o $@

clean:
	rm -f $(wildcard *.o) $(EXECUTABLE_TESTPYLON) $(EXECUTABLE_TESTSIGN) $(EXECUTABLE_TESTSPEED) $(EXECUTABLE_KAT) PQCsignKAT_*.rsp PQCsignKAT_*.req
