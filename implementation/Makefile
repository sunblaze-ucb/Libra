CC=g++
CFLAGS=
CDBG_FLAGS=-g
INCLUDE=include
OPT=-O2

linear_gkr: src/linear_gkr/main.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) -std=c++11 src/linear_gkr/verifier.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover.cpp src/linear_gkr/main.cpp -o bin/main -lgmp -lgmpxx
linear_gkr_opt: src/linear_gkr/main.cpp
	$(CC) $(CFLAGS) $(CDBG_FLAGS) -I $(INCLUDE) -std=c++11 src/linear_gkr/verifier.cpp src/linear_gkr/prime_field.cpp src/linear_gkr/prover.cpp src/linear_gkr/main.cpp -o bin/main $(OPT) -lgmp -lgmpxx
clean:
	rm bin/*