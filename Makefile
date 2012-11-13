MPICC=mpicc

all: matrix_power matrix_power

matrix_power: matrix_power.c
			$(MPICC) -lm -o matrix_power matrix_power.c

clean:
	    rm -rf matrix_power
