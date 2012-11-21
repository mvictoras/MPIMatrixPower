Description
============

For this project I implemented a program that calculates the |Xk| in parallel using various parallel formulations. The pdf of the assignment can be found here: File:MPIMatrixPower.pdf
[edit] Source code

You can download the source code of the application here: https://github.com/mvictoras/MPIMatrixPower
git clone https://github.com/mvictoras/MPIMatrixPower
[edit] Compilation

I have included a Makefile to compile the project.
make
[edit] Running

The project can be run on ARGO, that is a University of Illinois at Chicago cluster available for students and staff. If you are on ARGO (or any cluster that processes are submitted using TORQUE), you need to run:
./submit.sh 
-n <# of numbers> 
-p <# of nodes> 
-l <# of processors per node> 
-t <formulation: serial, strassen, cannon, dns>
-q <Input Method1: 2 numbers: the probability of −1, and the probability of +1 with m as a delimiter>
-w <Input Method2: Sequence of 4 numbers from the range {-1, 0, 1} with m as a delimiter>
-k <#power of matrix>
[-a] <use associativity>
If your cluster allows you to run directly programs (without process queues, then you run:
mpirun -np <num_of_nodes> matrix_power
-n <# of numbers> 
-p <# of nodes> 
-t <formulation: serial, strassen, cannon, dns>
-q <Input Method1: 2 numbers: the probability of −1, and the probability of +1 with m as a delimiter>
-w <Input Method2: Sequence of 4 numbers from the range {-1, 0, 1} with m as a delimiter>
-k <#power of matrix>
[-a] <use associativity> 
[-c]
All the arguments are mandatory except from the -c and -a options that are optional. By default the program prints the results into a human readable format, -c outputs in a tab separated format, so that the data can be easily read and processed by other processes. -a enables the associativity.
