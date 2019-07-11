CC=clang++
CFLAGS=-Wall -march=native
CDBG_FLAGS=-g
INCLUDE=-I include -I include/ate-pairing/include
LIB=include/ate-pairing/lib
OPT=-O3

#main program (Linear GKR + Polynomial delegation + ZK)
linear_gkr_zk: src/linear_gkr/main_zk.cpp
	$(CC) $(CFLAGS) $(OPT) $(INCLUDE) -L $(LIB) -std=c++11 src/VPD/inputvpd.cpp src/VPD/vpdR.cpp src/VPD/vpd_test.cpp src/linear_gkr/zk_verifier.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/zk_prover.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_zk.cpp -o bin/main_zk  -lgmp -lgmpxx -lzm
#GKR protocol without zero knowledge and polynomial delegation (Pure Linear GKR)
linear_gkr_opt: src/linear_gkr/main_fast_track.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -std=c++11 src/linear_gkr/verifier_fast_track.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover_fast_track.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_fast_track.cpp -o bin/main $(OPT) -lgmp -lgmpxx

#Debug purpose
linear_gkr: src/linear_gkr/main_fast_track.cpp
	$(CC) $(CFLAGS) $(CDBG_FLAGS) $(INCLUDE) -std=c++11 src/linear_gkr/verifier_fast_track.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover_fast_track.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_fast_track.cpp -o bin/main -lgmp -lgmpxx

#Debug purpose
linear_gkr_dbg: src/linear_gkr/main_fast_track.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -g -std=c++11 src/linear_gkr/verifier_fast_track.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover_fast_track.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_fast_track.cpp -o bin/main -lgmp -lgmpxx

#Debug purpose
linear_gkr_slow_track: src/linear_gkr/main_slow_track.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -std=c++11 src/linear_gkr/verifier.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_slow_track.cpp -o bin/main_slow -O2 -lgmp -lgmpxx

#clean
clean:
	rm bin/main

#Debug purpose
linear_gkr_zk_dbg: src/linear_gkr/main_zk.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -L $(LIB) -std=c++11 src/VPD/inputvpd.cpp src/VPD/vpdR.cpp src/VPD/vpd_test.cpp src/linear_gkr/zk_verifier.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/zk_prover.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_zk.cpp -g -o bin/main_zk -lgmp -lgmpxx -lzm


#Debug purpose
vpd_test:
	$(CC) $(CFLAGS) $(INCLUDE) -L $(LIB) -std=c++11 src/VPD/inputvpd.cpp -g -o inputvpd $(OPT) -lzm -lgmp -lgmpxx

#Debug purpose
clogc: src/traditional_gkr/ClogC.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -std=c++11 src/traditional_gkr/verifier_traditional.cpp src/linear_gkr/prime_field.cpp src/traditional_gkr/prover_clogc.cpp src/linear_gkr/polynomial.cpp src/traditional_gkr/ClogC.cpp -g -o bin/main_clogc $(OPT) -lgmp -lgmpxx

#Debug purpose
bclogc: src/traditional_gkr/bclogc.cpp
	$(CC) $(CFLAGS) $(INCLUDE) -std=c++11 src/traditional_gkr/verifier_traditional.cpp src/linear_gkr/prime_field.cpp src/traditional_gkr/prover_clogc.cpp src/linear_gkr/polynomial.cpp src/traditional_gkr/bclogc.cpp -g -o bin/main_bclogc $(OPT) -lgmp -lgmpxx

#Debug purpose
bc_clogc:
	$(CC) $(CFLAGS) $(INCLUDE) -std=c++11 src/traditional_gkr/verifier_bc_clogc.cpp src/linear_gkr/prime_field.cpp src/traditional_gkr/prover_bc_clogc.cpp src/linear_gkr/polynomial.cpp src/traditional_gkr/bc_clogc.cpp -g -o bin/main_bc_clogc -lgmp -lgmpxx



#linear_gkr_zk: 
#	$(CC) $(CFLAGS) $(CDBG_FLAGS) -I $(INCLUDE) -std=c++11 -pg src/linear_gkr/zk_verifier.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/zk_prover.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_zk.cpp -o bin/main -lgmp -lgmpxx
	
#linear_gkr_opt_parallel: src/linear_gkr/main_fast_para_track.cpp
#	$(CC) $(CFLAGS) -I $(INCLUDE) -std=c++11 src/linear_gkr/verifier_fast_para_track.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover_fast_para_track.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_fast_para_track.cpp -o bin/main $(OPT) -lgmp -lgmpxx -lpthread
