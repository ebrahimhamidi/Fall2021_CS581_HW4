/* 
   Sample Solution to the game of life program.
   Author: Purushotham Bangalore
   Date: Feb 17, 2010

   Use -DDEBUG1 for output at the start and end.
   Use -DDEBUG2 for output at each iteration.

   To compile: gcc -Wall -O -o life life.c
   To run: ./life <problem size> <max iterations> <nthread> <output dir>
           ./life 1000 1000 1 .                   (on your local system)
           ./life 5000 5000 1 /scratch/$USER/     (on DMC at ASC)
*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h> //new

#define DIES   0
#define ALIVE  1

/* function to measure time taken */
double gettime(void) {
  struct timeval tval;

  gettimeofday(&tval, NULL);

  return( (double)tval.tv_sec + (double)tval.tv_usec/1000000.0 );
}

/* allocate row-major two-dimensional array */
int **allocarray(int P, int Q) {
  int i, *p, **a;

  p = (int *)malloc(P*Q*sizeof(int));
  a = (int **)malloc(P*sizeof(int*));
  for (i = 0; i < P; i++)
    a[i] = &p[i*Q]; 

  return a;
}

/* free allocated memory */
void freearray(int **a) {
  free(&a[0][0]);
  free(a);
}

/* print arrays in 2D format */
void printarray(int **a, int N, int k) {
  int i, j;
  printf("Life after %d iterations:\n", k) ;
  for (i = 1; i < N+1; i++) {
    for (j = 1; j< N+1; j++)
      printf("%d ", a[i][j]);
    printf("\n");
  }
  printf("\n");
}

/* write array to a file (including ghost cells) */
void writefile(int **a, int N, FILE *fptr) {
  int i, j;
  for (i = 0; i < N; i++) {
    for (j = 0; j< N+2; j++)
      fprintf(fptr, "%d ", a[i][j]);
    fprintf(fptr, "\n");
  }
}

/* update each cell based on old values */
int compute(int **life, int **temp, int new_mycount, int N) {
  int i, j, value, flag=0;

  for (i = 1; i < new_mycount+1; i++) {
    for (j = 1; j < N+1; j++) {
      /* find out the value of the current cell */
      value = life[i-1][j-1] + life[i-1][j] + life[i-1][j+1]
	+ life[i][j-1] + life[i][j+1]
	+ life[i+1][j-1] + life[i+1][j] + life[i+1][j+1] ;
      
      /* check if the cell dies or life is born */
      if (life[i][j]) { // cell was alive in the earlier iteration
	if (value < 2 || value > 3) {
	  temp[i][j] = DIES ;
	  flag++; // value changed 
	}
	else // value must be 2 or 3, so no need to check explicitly
	  temp[i][j] = ALIVE ; // no change
      } 
      else { // cell was dead in the earlier iteration
	if (value == 3) {
	  temp[i][j] = ALIVE;
	  flag++; // value changed 
	}
	else
	  temp[i][j] = DIES; // no change
      }
    }
  }

  return flag;
}


int main(int argc, char **argv) {
  int N, NTIMES,  **life=NULL, **local_life=NULL, **temp=NULL, **ptr, NTHREADS;
  int i, j, k, flag=1, myflag=1;
  int myN, rank, size, remain,  mycount, *counts=NULL, *displs=NULL, new_mycount; //new  
  //int rank, size;
  int *x=NULL, *y; //new
  double t1, t2;
  char filename[BUFSIZ];
  FILE *fptr;


  
  MPI_Status status;
  MPI_Init(&argc, &argv);//new
  MPI_Comm_size(MPI_COMM_WORLD, &size);//new
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);//new
  
  //  if (argc != 4) {
  //  printf("Usage: %s <size> <max. iterations> <nthreads> <output dir>\n", argv[0]);
  //  exit(-1);
 //}  

//  N = atoi(argv[1]);//new


//  NTIMES = atoi(argv[2]);
 // sprintf(filename,"%s/output.%d.%d.%d", argv[3], N, NTIMES, size);
  
  
  
    if (argc != 5) {
    printf("Usage: %s <size> <max. iterations> <nthreads> <output dir>\n", argv[0]);
    exit(-1);
  }

  N = atoi(argv[1]);
  NTIMES = atoi(argv[2]);
  NTHREADS = atoi(argv[3]);
  sprintf(filename,"%s/output.%d.%d.%d", argv[4], N, NTIMES, NTHREADS);

    myN = (N) / size;
  remain = (N) % size;
  
  	if (rank == 0) {
  
  if ((fptr = fopen(filename, "w")) == NULL) {
     printf("Error opening file %s for writing\n", argv[3]);
     perror("fopen");
     exit(-1);
  }
  
 // if ((fptr = fopen(filename, "w")) == NULL) {
   //  printf("Error opening file %s for writing\n", argv[2]);
    // perror("fopen");
    // exit(-1);
 // }

	  /* Allocate memory for both arrays */
	  life = allocarray(N,N+2);	  

	  /* Initialize the boundaries of the life matrix */
	  for (i = 0; i < N; i++) {		
		life[i][0] = life[i][N+1] = DIES ;
	  }

	  /* Initialize the life matrix */
	  for (i = 0; i < N; i++) {
		srand48(54321|i);
		for (j = 1; j< N+1; j++)
		  if (drand48() < 0.5) 
		life[i][j] = ALIVE ;
		  else
		life[i][j] = DIES ;
	  }	  	  
	  
		x = &life[0][0];  
	  
	  	counts = malloc(sizeof(int)*size);
	    displs = malloc(sizeof(int)*size);
	    
	    for (i = 0; i < size; i++)
			counts[i] = myN + ((i < remain)?1:0);
			
	    displs[0] = 0;
	    for (i = 1; i < size; i++)
			displs[i] = displs[i-1] + counts[i-1];
			
	    for (i = 1; i < size; i++)
			displs[i] = displs[i] * (N+2);
			
	    for (i = 0; i < size; i++)
			counts[i] = counts[i] *(N+2);			
	}
	mycount = (myN + ((rank < remain)?1:0)) *(N+2);
    
#ifdef DEBUG
	if (rank == 0) {
	   for (i = 0; i < size; i++)
		   printf("counts[%d] = %d, displs[%d] = %d\n", i, counts[i], i, displs[i]);
	} else {
	  printf("rank = %d, mycount = %d\n", rank, mycount);
	}
#endif
    local_life = allocarray((mycount/(N+2))+2,N+2);
	y = &local_life[1][0];

	MPI_Scatterv(x, counts, displs, MPI_INT, y, mycount, MPI_INT, 0, MPI_COMM_WORLD);
	

#ifdef DEBUG1
  /* Display the initialized life matrix */
  printarray(life, N, 0);
#endif

  temp = allocarray((mycount/(N+2))+2,N+2);
  t1 = MPI_Wtime(); 
  /* Play the game of life for given number of iterations */
  for (k = 0; k < NTIMES; k++) {

	y = &local_life[0][0]; 
    int up_nbr, down_nbr;

    up_nbr = rank + 1;
    if (up_nbr >= size) up_nbr = MPI_PROC_NULL;
    down_nbr = rank - 1;
    if (down_nbr < 0) down_nbr = MPI_PROC_NULL;

	MPI_Sendrecv(&y[(N+2)], (N+2), MPI_INT, down_nbr, 0, 
		      &y[0], (N+2), MPI_INT, down_nbr, 0, 
		      MPI_COMM_WORLD, &status );

	MPI_Sendrecv(&y[(mycount/(N+2))*(N+2)], (N+2), MPI_INT, up_nbr, 0,
		      &y[((mycount/(N+2))+1)*(N+2)], (N+2), MPI_INT, up_nbr, 0, 
		      MPI_COMM_WORLD, &status );			  

    myflag = 0;
	//temp = allocarray((mycount/(N+2))+2,N+2);
	new_mycount = (mycount/(N+2));			
	myflag = compute(local_life, temp, new_mycount, N);

    MPI_Allreduce(&myflag, &flag, 1, MPI_INT, MPI_SUM,
         MPI_COMM_WORLD);
   
	  if (flag != 0) {
		ptr = local_life;		
		local_life = temp;
		temp = ptr;

#ifdef DEBUG2
    /* Print no. of cells alive after the current iteration */
    printf("No. of cells whose value changed in iteration %d = %d\n",k+1,flag) ;

    /* Display the life matrix */
    printarray(life, N, k+1);
#endif
   //MPI_Barrier(MPI_COMM_WORLD);
	  		flag = 0;
	  }else
		break;
  }
  t2 = MPI_Wtime(); 
  
  MPI_Gatherv(&y[N+2], mycount, MPI_INT, x, counts, displs, MPI_INT, 0, MPI_COMM_WORLD);	
		
#ifdef DEBUG1
  /* Display the life matrix after k iterations */
  printarray(life, N, k);
#endif
  if (rank == 0) {
	  
  printf("Time taken %f seconds for %d iterations\n", t2 - t1, k);

  /* Write the final array to output file */
  printf("Writing output to file: %s\n", filename);
  printf("Program terminates normally\n");
  writefile(life, N, fptr);
  fclose(fptr);
  freearray(life);
  free(counts);
  free(displs);
 // free(x);
   }
 // freearray(local_life);
  freearray(temp);
  //freearray(ptr);
 // free(y);
 
  MPI_Finalize();//new
  return 0;
}
