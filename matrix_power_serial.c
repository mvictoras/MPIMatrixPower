/**
* Copyright (c) <2012>, <Victor Mateevitsi> www.vmateevitsi.com
* All rights reserved.
* 
* Redistribution and use in source and binary forms, with or without
* modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the following disclaimer in the
*    documentation and/or other materials provided with the distribution.
* 3. All advertising materials mentioning features or use of this software
*    must display the following acknowledgement:
*    This product includes software developed by the <organization>.
* 4. Neither the name of the <organization> nor the
*    names of its contributors may be used to endorse or promote products
*    derived from this software without specific prior written permission.

* THIS SOFTWARE IS PROVIDED BY <COPYRIGHT HOLDER> ''AS IS'' AND ANY
* EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
* WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
* DISCLAIMED. IN NO EVENT SHALL <COPYRIGHT HOLDER> BE LIABLE FOR ANY
* DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
* (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
*	ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
*	(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
*	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*/

#include <getopt.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

int n, *local_coords;
float one_prob, minus_one_prob;
char s_local_coords[255];
int computerStats = 0;

// Timing
double gen_time, proc_time, comm_time, total_time;
int hasReceivedAkk = -1;
typedef enum { 
  serial, 
  mesh
} TYPE;
TYPE type;

int parse_arguments(int argc, char **argv);
float *generate_array(float one_prob, float minus_one_prob);

float *generate_serial();

void serial_deter();

void one_d_partitioning(float *A);
float *matrix_mult_serial(float *A);
float *block_mat_mult(float *A, int q);

void process_row_and_column(float *A, float *left_row, float *top_row, int startingRow, int k, int endingRow, int startingColumn, int endingColumn, int numRows, int numColumns, int local_k);

void printMatrix(float *A, int nElem);

int main(int argc, char **argv) {
  
	double t_start, t_end;
  int *ptr_gen_array, elem_per_node;
  float *local_array;  
  int i, num_procs, local_rank, name_len;
  gen_time = 0.0; proc_time = 0.0; comm_time = 0.0; total_time = 0.0;
 
  // Parse the arguments
 printf("YES\n"); 
  if( parse_arguments(argc, argv) ) return 1;
  /*
  ///t_start = MPI_Wtime();
  if( type == serial ) {
    local_array = generate_serial(&elem_per_node);
    float *C = matrix_mult_serial(local_array);
    one_d_partitioning(C);
    free(C);
  } 
  
  ///t_end = MPI_Wtime();
  total_time = t_end - t_start;

  if( computerStats ) {
    printf("%d\tg\t%s\t%d\t%f\n", n, s_local_coords, num_procs, gen_time);
    printf("%d\tp\t%s\t%d\t%f\n", n, s_local_coords, num_procs, proc_time);
    printf("%d\tc\t%s\t%d\t%f\n", n, s_local_coords, num_procs, comm_time);
    printf("%d\tt\t%s\t%d\t%f\n", n, s_local_coords, num_procs, total_time);
  }

  free(local_array); */
	return 0;
}

int parse_arguments(int argc, char **argv) {
	int c;
  int option_index = 0;
  static struct option long_options[] =
  {
    {"input_method1", required_argument,  0, 'q'},
    {"input_method2", required_argument,  0, 'w'},
    {0, 0, 0, 0}
  };

  char *result = NULL;
  char delims[] = ",";
	while( (c = getopt_long (argc, argv, "n:t:q:w:c", long_options, &option_index)) != -1 ) {
		switch(c) {
      case 'q':
        printf("optarg:%s\n", optarg);
        //result = strtok(optarg, delims);
        //one_prob = atof(result);
        //result = strtok(NULL, delims);
        //minus_one_prob = atof(result);
        break;
      case 'w':
        break;
			case 'n':
				n = atoi(optarg);
				break;
      case 'c':
        computerStats = 1;
        break;
      case 't':
        if( strcmp(optarg, "serial" ) == 0 ) type = serial;
        else if( strcmp(optarg, "mesh" ) == 0 ) type = mesh;
        else {
          fprintf( stderr, "Option -%c %s in incorrect. Allowed values are: serial, mesh\n", optopt, optarg);
          return 1;
        }
        break;
      case '?':
				if( optopt == 'n' )
					fprintf (stderr, "Option -%c requires an argument.\n", optopt);
				else if (isprint (optopt))
					fprintf (stderr, "Unknown option `-%c'.\n", optopt);
				else
					fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
				return 1;
			default:
				fprintf(stderr, "Usage: %s -n <number of numbers> \n", argv[0]);
				fprintf(stderr, "\tExample: %s -n 1000\n", argv[0]);
				return 1;
		}
	}
	return 0;
}

float *generate_array(float one_prob, float minus_one_prob) {
	unsigned int iseed = (unsigned int)time(NULL);
  int i,j; 
  float *gen_array;
  double start, end, dt;
 
  srand (iseed);
	gen_array = (float *)malloc(sizeof(float) * (n * n));
  
  ///start = MPI_Wtime();
  float a = 0;
  int c = 0;
  int d = 0;
	for( i = 0; i < n; i++ ) {
    for( j = 0; j < n; j++ ) {
      double number = (double)rand() / RAND_MAX;
      //printf("%f\n", number);
      if( number <= one_prob ) { a = 1; c++; }
      else if( number - one_prob <= minus_one_prob ) { a = -1; d++; }
      else a = 0;
		  gen_array[i * n + j] = a;
    }
  }
  ///end = MPI_Wtime();
  dt = end - start;
  gen_time = dt;
  int c_n = n * n * one_prob;
  int d_n = n * n * minus_one_prob;

  return gen_array;
}


float *generate_serial(char *proc_name) {
  float *local_array;
  double start, end, dt;
  int i, j;

  local_array = generate_array(one_prob, minus_one_prob);
 
  return local_array;
}


float *matrix_mult_serial(float *A) {
  return block_mat_mult(A, n);
}

float *block_mat_mult(float *A, int q) {
  int i, j, k;
  float *C;
  C = (float*)malloc(sizeof(float) * q * q);
  for( i = 0; i < q; i++ ) {
    for( j = 0; j < q; j++ ) { 
      bzero(C + i * q, q);
      for( k = 0; k < q; k++ ) {
        C[i * q + j] = C[i * q + j] + A[i * q + k] * A[k * q + j];
      }
    }
  }
  return C;
}

void one_d_partitioning(float *A) {
  int k, i, j;
  long double determinant;
  double start, end, dt;

  ///start = MPI_Wtime();
  for( k = 0;  k < n; k++ ) {
    for( j = k + 1; j < n; j++ ) {
      A[k * n + j] = A[k * n + j] / A[k * n + k];
    }
    for( i = k + 1; i < n; i++ ) {
      for( j = k + 1; j < n; j++ ) {
        A[i * n + j] -= A[i * n + k] * A[k * n + j];
      }
      A[i * n + k] = 0;
    }
  }
  ///end = MPI_Wtime();
  dt = end - start;
  proc_time += dt;

  ///start = MPI_Wtime();
  
  determinant = 1.0f;
   
  for( i = 0; i < n; i++ ) {
    printf("%f " , A[i * n + i]);
    determinant = determinant * A[i * n + i];
  }
  printf("\nDet is: %Lf\n", determinant);
  ///end = MPI_Wtime();
  dt = end - start;
  proc_time += dt;
}

void process_row_and_column(float *A, float *left_row, float *top_row, int k,
                            int startingRow, int endingRow, int startingColumn, int endingColumn, 
                            int numRows, int numColumns, int local_k) {
  int i, j;
  int index_row = 0;
  int index_column = 0;
  int starting_x = 0;
  int starting_y = 0;
  if( k >= startingColumn && k < endingColumn ) {
    starting_y = local_k + 1;
  }
  if( k >= startingRow && k < endingRow ) {
    starting_x = local_k + 1;
  }

  char str[255];
  for( i = starting_x; i < numRows; i++ ) {
    index_column = 0;
    for( j = starting_y; j < numColumns; j++ ) {
      A[i * numRows + j] -= left_row[index_row] * top_row[index_column];    
      index_column++;;
    }
    index_row++;
  }
}

void printMatrix(float *A, int nElem) {
  int i, j;
  for( i = 0; i < nElem; i++ ) {
    for( j = 0; j < nElem; j++ ) {
      printf("%f ", A[i * nElem + j]);
    }
    printf("\n");
  }
  printf("\n");
}

