#include <stdio.h>
#include <stdlib.h>
#include<time.h>

#include <mpi.h>
#define N 8

int* createForest(int n,int m);
void initializeForestRandom(int* f);
void print_Forest(int* f,int n, int m);
void scatter_Forest(int* f,int n,int m);
int* compute(int* f,int rank, int world_size);


MPI_Comm comm_world;
int main(int argc, char** argv) {
    int world_size;
    int rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    /**
     * 
     * Inserire il controllo del numero dei Processori
     * 
    */
    

    MPI_Status status;
    int msg[N];
    int *forest;
    int *sendcount;
    int *starting; 
    int t;
    MPI_Datatype row;
    MPI_Type_contiguous( N , MPI_INT , &row);
    MPI_Type_commit( &row);
    sendcount=malloc(sizeof(int)*world_size);
    starting=malloc(sizeof(int)*world_size);
    int recvcount;
    if (rank==0)
    {
    
    forest=createForest(N,N);
    /*
    Inizialzzazione Foresta
    */
    srand(time(0));
    for (size_t i = 0; i < N*N; i++)
    {
        forest[i]= rand()%3;
    }

    print_Forest(forest,N,N);
    /*
    Stampa foresta
    */
    for (size_t i = 0; i < N*N; i++)
    {
        if (i%N==0)
        {
         printf("\n");   
        }
        switch (forest[i])
        {
        
        case 0:
            printf(" \U0001F7E9 ");
            break;
        case 1:
            printf(" \U0001F333 ");
            break;
        case 2:
            printf(" \U0001F525 ");
            break;

        default:
            break;
        }
        
    }
    printf("\n");
    }

    scatter_Forest(forest,N,N);
/*Sending Forest method*/
    int execess_row=N%world_size;
    int row_for_each_process=N/world_size;
    int recv_count;
    for (size_t i = 0; i < world_size; i++)
    {
        sendcount[i]=row_for_each_process;
        if (execess_row>0)
        {
            execess_row-=1;
            sendcount[i]+=1;
        }
        if (i==0)
        {
            starting[i]=0;
        }
        else{
            starting[i]=starting[i-1]+sendcount[i-1];
        }     
        if (rank==i)
        {
            recv_count=sendcount[i];
        }
        if (rank==0)
        {
            MPI_Send(&forest[starting[i]*N] , sendcount[i] , row , i , 10 , MPI_COMM_WORLD);
        }
        
    }


    MPI_Recv( &msg , recv_count , row, 0 , 10 , MPI_COMM_WORLD , &status);
    
    msg=compute(msg,rank,world_size);
/*Compute*/

    for (size_t i = 0; i < N-1 ; i++)
    {
        int isempty=-1;
        int number_of_element=recv_count*N;
        for (size_t j = 0; j < number_of_element; j++)
        {
    
            if (msg[j]==2)
            {

                int nrow=j/N;
                int position=j%N;
                if ( (position!=0) && (msg[j-1]==1) )
                {
                    msg[j-1]=2;
                }
                if ( (position!=N-1) && (msg[j+1]==1) )
                {
                    msg[j+1]=2;
                }
                if ( (nrow!=0) && (msg[j-N]==1) )
                {
                    msg[j-N]=2;
                }
                if ( (nrow!=recv_count-1) && (msg[j+N]==1) )
                {
                    msg[j+N]=2;
                }
                if ( (rank!=0) && (nrow==0) )
                {
                    isempty=1;
                    MPI_Send(&position,1,MPI_INT,rank-1,10,MPI_COMM_WORLD);
                }
                if ( (rank!=world_size-1) && (nrow==recv_count-1) )
                {
                    isempty=1;
                    MPI_Send(&position,1,MPI_INT,rank+1,10,MPI_COMM_WORLD);
                }                      
            }
            
        }
        if (isempty==-1)
        {
            if (rank!=0)
                {
                    MPI_Send(&isempty,1,MPI_INT,rank-1,10,MPI_COMM_WORLD);
                }
                if (rank!=world_size-1)
                {
                    MPI_Send(&isempty,1,MPI_INT,rank+1,10,MPI_COMM_WORLD);
                }  
        }
        if (rank!=0)
        {
            MPI_Recv(&t, 1, MPI_INT, rank-1 , 10 , MPI_COMM_WORLD ,&status);
            if (t!=-1)
            {
                if (msg[t]==1)
                {
                    msg[t]=2;
                }
            }
            
        }
        if (rank!=world_size-1)
        {
            MPI_Recv(&t, 1, MPI_INT, rank+1 , 10 , MPI_COMM_WORLD ,&status);
            if (t!=-1)
            {
                int pn=N*recv_count-1;
                pn+=t;
                if (msg[pn]==1)
                {
                    msg[pn]=2;
                }
            }
        }
    }
        
        
    /**
     * 
     * RESULTS
     * 
    */
    int *final;
    if (rank==0)
    {
        final=malloc(sizeof(int)*N*N);
    }
    MPI_Gatherv( &msg, recv_count , row , final , sendcount , starting , row , 0 , MPI_COMM_WORLD);
    if (rank==0)
    {

    print_Forest(forest,N,N)
    /**
     * Stampa foresta
    */
    for (size_t i = 0; i < N*N; i++)
    {
        if (i%N==0)
        {
         printf("\n");   
        }
        switch (final[i])
        {
        
        case 0:
            printf(" \U0001F7E9 ");
            break;
        case 1:
            printf(" \U0001F333 ");
            break;
        case 2:
            printf(" \U0001F525 ");
            break;

        default:
            break;
        }
        
    }
    printf("\n");
    
    }

    MPI_Finalize();
}