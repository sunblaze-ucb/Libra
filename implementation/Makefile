CC=g++
CFLAGS=
CDBG_FLAGS=-g
INCLUDE=include
OPT=-O2

linear_gkr: src/linear_gkr/main.cpp
	$(CC) $(CFLAGS) -I $(INCLUDE) src/linear_gkr/main.cpp -o bin/main
linear_gkr_opt: src/linear_gkr/main.cpp
	$(CC) $(CFLAGS) $(CDBG_FLAGS) -I $(INCLUDE) src/linear_gkr/main.cpp -o bin/main $(OPT)
clean:
	rm bin/*