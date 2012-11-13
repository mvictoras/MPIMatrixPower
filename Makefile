CC=mpicc

all: matrix_power matrix_power_serial

matrix_power: matrix_power.c
			$(CC) -lm -o matrix_power matrix_power.c

matrix_power_serial: matrix_power_serial.c
			$(CC) -lm -o matrix_power_serial matrix_power_serial.c

clean:
	    rm -rf matrix_power matrix_power_serial
