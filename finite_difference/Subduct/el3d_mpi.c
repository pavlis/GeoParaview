/*###################################################################
	Parallel poroelastice wave propagation in 3D Version 1.0
	April, 2005
  
	Author:	Dong-Hoon Sheen
		Ph.D. student at Seoul National University
	       	Research associate at Indiana University
	
	Reference: 
	Parallel implementation of a velocity-stress
	    staggered-grid finite-difference method 
	    for poroelastic wave propagation
	    Sheen, D.-H., Tuncay, K., Baag, C.-E. and Ortoleva, P. J.
###################################################################*/

#include "el3d_mpi.h"

int main(int argc,char *argv[])
{
   int n,n1,n2,n3;
   int numProcs;
   double t1,t2;

   MPI_Init(&argc,&argv);
   
   MPI_Comm_rank(MPI_COMM_WORLD, &Myrank);
   MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

   elmdls.rank[3]   = Myrank;
   elmdls.nProcs[3] = numProcs;

   n1 = (pow(numProcs,1./3.)+0.1);
   n2 = n1;
   n3 = n1;
   find_opt(&n1,&n2,&n3,numProcs);
   if( Myrank == 0 ) 
   {
	if( n1 == 0 )
	    sherror("\n# Check the number of processors : %d\n",numProcs);

	printf( "\n# The number of processors in this run is %d.\n", numProcs );

	if( argc > 1 ) Read_setup(argv[1]);
	else Read_setup("input_pewave3d.txt");
   }

   ErrorCheck("Error from Setup Processers in main\n");
   comm_distribute_setup(Myrank);
   ErrorCheck("Error from Distribute Processers in main\n");

   Initialize(n1,n2,n3);

   Make_PML_layer();
   ErrorCheck("Error from Make PML in main\n");

   Set_source();
   Set_receiver();
   ErrorCheck("Initialize in main\n");

   Print_Configuration();
   ErrorCheck("Error from Print config in main\n");

   Alloc_Mem_Func();
   ErrorCheck("Error from Memory allocation in main\n");

   Init_model_parameter();
   ErrorCheck("Init para in main\n");

   Write_mdl_parameter();

   t1 = MPI_Wtime();

   Initialize_fp();
   Update();

   t2 = (MPI_Wtime()-t1);
   printf("%d:Time:%f\n",Myrank,t2);

   comm_write_seismogram();

   MPI_Finalize();
   return 0;
}

int prime(int number)
{
    register int i;
    for(i=number/2;i>=2;i--)
    {
            if(number%i == 0 ) number = i;
    }
    return number;
}

void find_opt(int *m1,int *m2,int *m3,int number)
{
    register int i,j,k,n;
    if(*m1 < 2)
    {
        *m1 = 0;
        *m2 = 0;
        *m3 = 0;
        return;
    }
    n = *m1 * prime( *m1 );
    i=*m1;
    for(j=i;j<=n;j++)
    {
        for(k=j;k<=n;k++)
        {
            if( i*j*k == number )
            {
                *m1 = i;
                *m2 = j;
                *m3 = k;
                return;
            }
            if( i*j*k > number ) break;
        }
    }
    *m1 -= 1;
    *m2 -= 1;
    *m3 -= 1;
    find_opt(m1,m2,m3,number);
    return;
}

