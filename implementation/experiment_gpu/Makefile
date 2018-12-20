CC=g++
CFLAGS=
CDBG_FLAGS=-g
INCLUDE=include
OPT=-O3
linear_gkr_opt_parallel: src/linear_gkr/main_fast_para_track.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -std=c++11 src/linear_gkr/verifier_fast_para_track.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover_fast_para_track.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_fast_para_track.cpp -o bin/main $(OPT) -lgmp -lgmpxx -lpthread
