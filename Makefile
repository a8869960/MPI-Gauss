FLAGS = -isystem /opt/toolchain/devel3/impi-4.1.3.049/intel64/include/ -O3 -mfpmath=sse -fstack-protector-all -g -W -Wall -Wextra -Wunused -Wcast-align -Werror -pedantic -pedantic-errors -Wfloat-equal -Wpointer-arith -Wformat-security -Wmissing-format-attribute -Wformat=1 -Wwrite-strings -Wcast-align -Wno-long-long -Woverloaded-virtual -Wnon-virtual-dtor -Wcast-qual -Wno-suggest-attribute=format
CC = mpicxx

all: main.o matrix_input.o matrix_output.o mpi_functions.o mpi_gauss.o matrix_operations.o for_gauss.o inverse_matrix.o results.o
	$(CC) main.o matrix_input.o matrix_output.o mpi_functions.o mpi_gauss.o matrix_operations.o for_gauss.o inverse_matrix.o results.o

main.o: main.cpp functions.h
	$(CC) $(FLAGS) -c main.cpp

matrix_input.o: matrix_input.cpp functions.h
	$(CC) $(FLAGS) -c matrix_input.cpp

matrix_output.o: matrix_output.cpp functions.h
	$(CC) $(FLAGS) -c matrix_output.cpp

mpi_functions.o: mpi_functions.cpp functions.h
	$(CC) $(FLAGS) -c mpi_functions.cpp

mpi_gauss.o: mpi_gauss.cpp functions.h
	$(CC) $(FLAGS) -c mpi_gauss.cpp

matrix_operations.o: matrix_operations.cpp functions.h
	$(CC) $(FLAGS) -c matrix_operations.cpp

for_gauss.o: for_gauss.cpp functions.h
	$(CC) $(FLAGS) -c for_gauss.cpp

inverse_matrix.o: inverse_matrix.cpp functions.h
	$(CC) $(FLAGS) -c inverse_matrix.cpp

results.o: results.cpp functions.h
	$(CC) $(FLAGS) -c results.cpp

clean:
	rm -f *.out *.o *.gch
