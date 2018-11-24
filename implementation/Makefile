CC=g++
CFLAGS=
CDBG_FLAGS=-g
INCLUDE=include
OPT=-O2

linear_gkr: src/linear_gkr/main_fast_track.cpp
	$(CC) $(CFLAGS) $(CDBG_FLAGS) -I $(INCLUDE) -std=c++11 src/linear_gkr/verifier_fast_track.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover_fast_track.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_fast_track.cpp -o bin/main -lgmp -lgmpxx
linear_gkr_opt: src/linear_gkr/main_fast_track.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -std=c++11 src/linear_gkr/verifier_fast_track.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover_fast_track.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_fast_track.cpp -o bin/main $(OPT) -lgmp -lgmpxx
linear_gkr_slow_track: src/linear_gkr/main_slow_track.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -std=c++11 src/linear_gkr/verifier.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_slow_track.cpp -o bin/main_slow -O2 -lgmp -lgmpxx
clean:
	rm bin/main
linear_gkr_prof:
	$(CC) $(CFLAGS) $(CDBG_FLAGS) -I $(INCLUDE) -std=c++11 -pg src/linear_gkr/verifier_fast_track.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover_fast_track.cpp src/linear_gkr/polynomial.cpp src/linear_gkr/main_fast_track.cpp -o bin/main -lgmp -lgmpxx