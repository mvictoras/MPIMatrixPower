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
#include <mpi.h>
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
float *generate_array(int num_procs, char *proc_name, int local_rank, float one_prob, float minus_one_prob);
void create_2dmesh_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs);
void create_ring_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs);

float *generate_serial(MPI_Comm *comm_new, int local_rank, int num_procs, char *proc_name, int *elem_per_node);
float *generate_mesh(MPI_Comm *comm_new, int local_rank, int num_procs, char *proc_name, int *elem_per_node);

void serial_deter();

void one_d_partitioning(MPI_Comm *comm_new, float *A, int local_rank, int num_procs);
void two_d_partitioning(MPI_Comm *comm_new, float *A, int local_rank, int num_procs);
float *matrix_mult_serial(float *A);
float *block_mat_mult(float *A, int q);

void process_row_and_column(float *A, float *left_row, float *top_row, int startingRow, int k, int endingRow, int startingColumn, int endingColumn, int numRows, int numColumns, int local_k);

void printMatrix(float *A, int nElem);
void send_to(MPI_Comm *comm, int direction, float *A, int size, int row, int column, int n);
void receive_from_left(MPI_Comm *comm, int direction, float *A, int size, int row, int column, int n, int k);

int main(int argc, char **argv) {
	double t_start, t_end;
  int *ptr_gen_array, elem_per_node;
  float *local_array;  
  int i, num_procs, local_rank, name_len;
  	
	char proc_name[MPI_MAX_PROCESSOR_NAME];
	MPI_Comm comm_new;
  gen_time = 0.0; proc_time = 0.0; comm_time = 0.0; total_time = 0.0;
 
  // Parse the arguments
  if( parse_arguments(argc, argv) ) return 1;
  // Initialize MPI
  MPI_Init(&argc, &argv); 
	MPI_Get_processor_name(proc_name, &name_len);
  
  // Initially create the topology
  if( type == serial ) {
    create_ring_topology(&comm_new, &local_rank, &num_procs);
  } else {
    create_2dmesh_topology(&comm_new, &local_rank, &num_procs);
  }

  t_start = MPI_Wtime();
  if( type == serial ) {
    local_array = generate_serial(&comm_new, local_rank, num_procs, proc_name, &elem_per_node);
    float *C = matrix_mult_serial(local_array);
    one_d_partitioning(&comm_new, C, local_rank, num_procs);
    free(C);
  } else {
    local_array = generate_mesh(&comm_new, local_rank, num_procs, proc_name, &elem_per_node);
  }

  /*
  if( type == ring ) {
    one_d_partitioning(&comm_new, local_array, local_rank, num_procs);
  } else {
    two_d_partitioning(&comm_new, local_array, local_rank, num_procs);
  }
*/
  t_end = MPI_Wtime();
  total_time = t_end - t_start;

  if( computerStats ) {
    printf("%d\tg\t%s\t%d\t%f\n", n, s_local_coords, num_procs, gen_time);
    printf("%d\tp\t%s\t%d\t%f\n", n, s_local_coords, num_procs, proc_time);
    printf("%d\tc\t%s\t%d\t%f\n", n, s_local_coords, num_procs, comm_time);
    printf("%d\tt\t%s\t%d\t%f\n", n, s_local_coords, num_procs, total_time);
  }

  free(local_array);
	//MPI_Comm_free(&comm_new);
	MPI_Finalize(); // Exit MPI
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
  char delims[] = "m";
	while( (c = getopt_long (argc, argv, "n:t:q:w:c", long_options, &option_index)) != -1 ) {
		switch(c) {
      case 'q':
        //printf("optarg:%s\n", optarg);
        result = strtok(optarg, delims);
        one_prob = atof(result);
        result = strtok(NULL, delims);
        minus_one_prob = atof(result);
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


void create_ring_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs) {
	int dims[1], periods[1];
 
  MPI_Comm_size(MPI_COMM_WORLD, num_procs);

  dims[0] = *num_procs;
  periods[0] = 1;
  local_coords = (int *) malloc(sizeof(int) * 1);
  // Create the topology
	MPI_Cart_create(MPI_COMM_WORLD, 1, dims, periods, 0, comm_new);
	MPI_Comm_rank(*comm_new, local_rank);
	MPI_Comm_size(*comm_new, num_procs);
	MPI_Cart_coords(*comm_new, *local_rank, 1, local_coords);
  sprintf(s_local_coords, "[%d]", local_coords[0]);
}

void create_2dmesh_topology(MPI_Comm *comm_new, int *local_rank, int *num_procs) {
  int *dims, i, *periods, nodes_per_dim;
  
  MPI_Comm_size(MPI_COMM_WORLD, num_procs);
	MPI_Comm_rank(MPI_COMM_WORLD, local_rank);
 
  int dimension = 2;
  nodes_per_dim = (int) sqrt( (double) *num_procs );
  local_coords = (int *) malloc(sizeof(int) * dimension);
  dims = (int *) malloc(sizeof(int) * dimension);
  periods = (int *) malloc(sizeof(int) * dimension);
  for( i = 0; i < dimension; i++ ) {
    dims[i] = nodes_per_dim;
    periods[i] = 0;
  }

  MPI_Cart_create(MPI_COMM_WORLD, dimension, dims, periods, 0, comm_new);
	MPI_Comm_size(*comm_new, num_procs);
  MPI_Cart_coords(*comm_new, *local_rank, dimension, local_coords);
  sprintf(s_local_coords, "[%d][%d]", local_coords[0], local_coords[1]);
}

float *generate_array(int num_procs, char *proc_name, int local_rank, 
                      float one_prob, float minus_one_prob) {
	unsigned int iseed = (unsigned int)time(NULL);
  int i,j; 
  float *gen_array;
  double start, end, dt;
 
  srand (iseed);
	gen_array = (float *)malloc(sizeof(float) * (n * n));
  
  start = MPI_Wtime();
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
  end = MPI_Wtime();
  dt = end - start;
  gen_time = dt;
  int c_n = n * n * one_prob;
  int d_n = n * n * minus_one_prob;

  if( !computerStats )
    printf("(%s(%d/%d)%s: %d random numbers generated in %1.8fs one:%d/%d minus one:%d/%d\n", proc_name, local_rank, num_procs, s_local_coords, n * n, dt, c, c_n, d, d_n);
  
  return gen_array;
}


float *generate_serial(MPI_Comm *comm_new, int local_rank, int num_procs, 
                 char *proc_name, int *elem_per_node) {
  float *local_array;
  double start, end, dt;
  int i, j;

  if( local_rank == 0 ) {
    local_array = generate_array(num_procs, proc_name, local_rank, one_prob, minus_one_prob);
  } else local_array = (float *) malloc(sizeof(float) * n * n);
 
  if( !computerStats )
    printf("(%s(%d/%d)%s: It took %1.8fs to receive the sub-array\n", proc_name, local_rank, num_procs, s_local_coords, dt);

  return local_array;
}


float *generate_mesh(MPI_Comm *comm_new, int local_rank, int num_procs, 
                 char *proc_name, int *elem_per_node) {
  float *local_array;
  float *tmp_array;
  double start, end, dt;
  int i, j;
  int nProc = n / sqrt(num_procs);
  MPI_Status status;

  local_array = (float *) malloc(sizeof(float) * nProc * nProc);
  if( local_rank == 0 ) {
    tmp_array = generate_array(num_procs, proc_name, local_rank, one_prob, minus_one_prob);
  }
    /*
    int i, j, k, l;
    for( i = 0; i < sqrt(num_procs); i++ ) {
      for( j = 0; j < sqrt(num_procs); j++ ) {
        if( i == 0 && j == 0 ) {
          int index = 0;
          for( k = 0; k < nProc; k++ ) {
            for( l = 0; l < nProc; l++ ) {
              local_array[index++] = tmp_array[k * n + l];
            }
          }
        }
        else {
          float *buff_to_send = (float*)malloc(sizeof(float) * nProc * nProc);
          int startingRow = i * nProc;
          int startingColumn = j * nProc;
          int index = 0;
          for( k = startingRow; k < startingRow + nProc; k++ ) {
            for( l = startingColumn; l < startingColumn + nProc; l++ ) {
              buff_to_send[index++] = tmp_array[k * n + l];
            }
          }
          MPI_Send(buff_to_send, nProc * nProc, MPI_FLOAT, j * sqrt(num_procs) + i, 0, *comm_new);
          free (buff_to_send);
        }
      }
    }
    free(tmp_array);
  } else {
    MPI_Recv(local_array, nProc * nProc, MPI_FLOAT, 0, 0, *comm_new, &status);
  }
*/
  if( !computerStats )
    printf("(%s(%d/%d)%s: It took %1.8fs to receive the sub-array\n", proc_name, local_rank, num_procs, s_local_coords, dt);

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

void one_d_partitioning(MPI_Comm *comm_new, float *A, int local_rank, int num_procs) {
  MPI_Status status;
  int k, i, j, startingRow, endingRow, numRows;
  long double determinant;
  double start, end, dt;

  numRows = n / num_procs;
  startingRow = local_rank * numRows;
  endingRow = startingRow + numRows;

  start = MPI_Wtime();
  if( local_rank != 0 ) {
    MPI_Recv(A, n * n, MPI_FLOAT, local_rank - 1, 0, *comm_new, &status);  
  }
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;
 
  start = MPI_Wtime();
  for( k = startingRow;  k < endingRow; k++ ) {
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
  end = MPI_Wtime();
  dt = end - start;
  proc_time += dt;

  start = MPI_Wtime();
  if( local_rank != num_procs - 1 )
    MPI_Send(A, n * n, MPI_FLOAT, local_rank + 1, 0, *comm_new);
  
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;

  determinant = 1.0f;
  if( !computerStats && local_rank == num_procs - 1) {
   
    for( i = 0; i < n; i++ ) {
      printf("%f " , A[i * n + i]);
      determinant = determinant * A[i * n + i];
    }
    printf("\nDet is: %Lf\n", determinant);
  }
  end = MPI_Wtime();
  dt = end - start;
  proc_time += dt;
}

void two_d_partitioning(MPI_Comm *comm_new, float *A, int local_rank, int num_procs) {
  MPI_Status status;
  int k, i, j, startingRow, endingRow, numRows, startingColumn, endingColumn, numColumns;
  int n_startingRow, n_startingColumn, n_local_coords[2];
  //long double determinant;
  double start, end, dt;
  int p = (int) sqrt(num_procs);
  int dis, left_rank, right_rank, up_rank, down_rank;
  MPI_Request req;
  
  numRows = n / p;
  numColumns = numRows;

  startingRow = local_coords[1] * numRows;
  endingRow = startingRow + numRows;

  startingColumn = local_coords[0] * numRows;
  endingColumn = startingColumn + numColumns;

  
  start = MPI_Wtime();
  
  for( k = 0; k < n; k++ ) {
    float Akk[1];
    int local_k = k % numRows;
    
    // Send A(k,k) to the right
    start = MPI_Wtime();
    if( k >= startingColumn && k < endingColumn && k >= startingRow && k < endingRow ) {
      send_to(comm_new, 0, A, 1, local_k, local_k, numRows);
      Akk[0] = A[local_k * numRows + local_k];
    } else if( k < startingColumn && k >= startingRow && k < endingRow ) {
      receive_from_left(comm_new, 0, Akk, 1, 0, 0, numRows, k);
    }
    end = MPI_Wtime();
    dt = end - start;
    comm_time += dt;

    // Now calculate the row
    start = MPI_Wtime();
    if( k >= startingColumn && k < endingColumn && k >= startingRow && k < endingRow ) {
      for( j = local_k + 1; j < numColumns; j++ ) {
        A[local_k * numRows + j] /= Akk[0];
      }
    } else if( k >= startingRow && k < endingRow && k < startingColumn ) {
      for( j = 0; j < numColumns; j++ ) {
        A[local_k * numRows + j] /= Akk[0];
      }
    }
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;

    // Now calculate the box
    int m, bOutside = 1; 
    float top_row[numRows]; 

    start = MPI_Wtime();
    // k is West of this Partition
    if( k >= startingRow && k < endingRow & k < startingColumn ) {
      send_to(comm_new, 1, A, numColumns, local_k, 0, numRows);
      
      for( m = 0; m < numColumns; m++ ) {
        top_row[m] = A[local_k * numRows + m];
      }
      bOutside = -1;
    } 
    // k is in this BOX
    else if( k >= startingRow && k < endingRow && k >= startingColumn && k < endingColumn ) {
      int size = numColumns - (local_k + 1);
      if( size != 0 ) {
        send_to(comm_new, 1, A, size, local_k, local_k + 1, numRows);

        for( m = 0; m < size; m++ ) {
          top_row[m] = A[local_k * numRows + local_k + 1 + m];
        }
        bOutside = -1;
      }
    } // k is NW of this box 
    else if( k < startingRow && k < startingColumn ) {
      int sender_row = k / numRows;
      int sender_column = k / numColumns;
      int sender_rank = local_coords[0] * sqrt(num_procs) + sender_row;
      
      MPI_Recv(top_row, numColumns, MPI_FLOAT, sender_rank, 0, *comm_new, &status);
      
      bOutside = -1;
    } 
    // k is N of this box
    else if( k < startingRow && k >= startingColumn && k < endingColumn ) {
      int sender_row = k / numRows;
      int sender_column = k / numColumns;
      int sender_rank = sender_column * sqrt(num_procs) + sender_row;
      int size = numColumns - (local_k + 1);
      
      if( size != 0 ) { 
        //top_row = (float *)malloc(sizeof(float) * numberToReceive);
        //printf("%d Waiting to receive from:%d\n", local_rank, sender_rank);
        MPI_Recv(top_row, size, MPI_FLOAT, sender_rank, 0, *comm_new, &status);
        
        bOutside = -1;
      } 
    }
    float left_row[numRows];
    
    // k is N of this Box
    if( k >= startingColumn && k < endingColumn & k < startingRow ) {
      for(m = 0; m < numRows; m++ ) {
        left_row[m] = A[m * numColumns + local_k];
      }
      send_to(comm_new, 0, left_row, numRows, 0, 0, 0);
      
      bOutside = -1;
    } 
    // k is IN this box 
    else if( k >= startingRow && k < endingRow && k >= startingColumn && k < endingColumn ) {
      //int local_k = k % numRows;
      int size = numColumns - (local_k + 1);
      if( size != 0 ) {
        for(m = 0; m < size; m++ ) {
          left_row[m] = A[(local_k + 1) * numColumns + local_k];
        }
        send_to(comm_new, 0, left_row, size, 0, 0, 0);
        
        bOutside = -1;
      }
    } 
    // k is SW from this box
    else if( k < startingRow && k < startingColumn ) {
      int sender_row = k / numRows;
      int sender_column = k / numColumns;
      int sender_rank = sender_column * sqrt(num_procs) + local_coords[1];
      
      MPI_Recv(left_row, numColumns, MPI_FLOAT, sender_rank, 0, *comm_new, &status);
      
      bOutside = -1;
    } 
    // k is W of this box
    else if( k < startingColumn && k >= startingRow && k < endingRow ) {
      int sender_row = k / numRows;
      int sender_column = k / numColumns;
      int sender_rank = sender_column * sqrt(num_procs) + local_coords[1];
      int local_k = k % numRows;
      int numberToReceive = numColumns - (local_k + 1);
      
      if( numberToReceive != 0 ) { 
        MPI_Recv(left_row, numberToReceive, MPI_FLOAT, sender_rank, 0, *comm_new, &status);
        bOutside = -1;
      }
    }
    end = MPI_Wtime();
    dt = end - start;
    comm_time += dt;
   

    // Now process the box
    if( bOutside < 0 ) {
      start = MPI_Wtime();
      process_row_and_column(A, left_row, top_row, k, startingRow, endingRow, startingColumn, endingColumn, numRows, numColumns, local_k);
      end = MPI_Wtime();
      dt = end - start;
      proc_time += dt;
    }
  } // end for

  float determinant[1];
  float result[1];
  determinant[0] = 1;
  if( local_coords[0] == local_coords[1] ) {
    start = MPI_Wtime();
    for(i = 0; i < numRows; i++ ) {
      determinant[0] *= A[i * numRows + i];
    }
    end = MPI_Wtime();
    dt = end - start;
    proc_time += dt;
  }
  
  start = MPI_Wtime();
  MPI_Reduce(determinant, result, 1, MPI_FLOAT, MPI_PROD, 0, *comm_new);
  end = MPI_Wtime();
  dt = end - start;
  comm_time += dt;
  if( !computerStats && local_rank == 0 ) {
    printf("Determinant is %f\n", result[0]);
  }
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

/*
 * int direction (0 horizontal, 1 vertical)
 * int distance  
 */
void send_to(MPI_Comm *comm, int direction, float *A, int size, int row, int column, int n) {
  int prev_rank, next_rank;
  int distance = 1;
  MPI_Cart_shift(*comm, direction, distance, &prev_rank, &next_rank);
  while(next_rank >= 0) {
    MPI_Send(A + row * n + column, size, MPI_FLOAT, next_rank, 0, *comm);
    MPI_Cart_shift(*comm, direction, ++distance, &prev_rank, &next_rank);
  }
} 
    
void receive_from_left(MPI_Comm *comm, int direction, float *A, int size, int row, int column, int n, int k) {
  int prev_rank, next_rank, coords[2], startingRow, startingColumn;
  int distance = 1;
  MPI_Status status;

  MPI_Cart_shift(*comm, direction, distance, &prev_rank, &next_rank);
  while(prev_rank >= 0) {
    MPI_Cart_coords(*comm, prev_rank, 2, coords);
    
    startingRow     = coords[1] * n; 
    startingColumn  = coords[0] * n;
    
    if( k >= startingColumn && k < startingColumn + n ) {
      MPI_Recv(A, size, MPI_FLOAT, prev_rank, 0, *comm, &status);
    }
    MPI_Cart_shift(*comm, 0, ++distance, &prev_rank, &next_rank);
  }
}
