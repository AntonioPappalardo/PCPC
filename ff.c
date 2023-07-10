#include <stdio.h>
#include <stdlib.h>
#include<time.h>

#include <mpi.h>
#define N 8

typedef struct{
    int* matrix; 
    int rows;
    int columns;
}Forest;

Forest* createForest(int n,int m);
void initializeForestRandom(Forest* f);
void print_Forest(Forest* f);
void print_Forest_emoji(Forest* f);
Forest* scatter_Forest(int n,int m);
Forest* compute(Forest* f,int rank, int world_size);

int* sendcount;
int* starting; 
int recv_count;

MPI_Comm comm_world;
MPI_Status status;
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
    int* msg;
    Forest* forest;
    int t;
    MPI_Datatype row;
    MPI_Type_contiguous( N , MPI_INT , &row);
    MPI_Type_commit( &row);
    sendcount=malloc(sizeof(int)*world_size);
    starting=malloc(sizeof(int)*world_size);
    
    forest= scatter_Forest(N,N);
    Forest* ff_compiled;
    ff_compiled=forest;
    ff_compiled=compute(ff_compiled,rank,world_size); 
    /**
     * 
     * RESULTS
     * 
    */
    Forest* final_forest;
    int* pointer_matrix;
    if (rank==0)
    {
        final_forest=createForest(N,N);
        pointer_matrix=&final_forest->matrix[0];
    }
    printf("\n");
    MPI_Gatherv( &(ff_compiled->matrix[0]), sendcount[rank] , row , pointer_matrix , sendcount , starting , row , 0 , MPI_COMM_WORLD);
    if (rank==0)
    {
        print_Forest_emoji(final_forest);
    }
    MPI_Finalize();
}

/**
 * Function for create a forest with n row and m column
*/
Forest* createForest(int n,int m){
    Forest* f=(Forest*)malloc(sizeof(Forest));
    f->matrix=(int*)malloc(sizeof(int)*n*m);
    f->columns=m;
    f->rows=n;
    return f;
}

/**
 * Function for Initialize a Forest
 * where:
 * 0 is an empty cell
 * 1 is a cell with tree
 * 2 is a cell with fire
*/
void initializeForestRandom(Forest* f){
    srand(time(0));
    for (size_t i = 0; i < f->rows*f->columns; i++)
    {
        f->matrix[i]= rand()%3;
    }
}

/**
 * Function for print the forest with emoji
*/
void print_Forest_emoji(Forest *f){
    for (size_t i = 0; i < f->rows; i++)
    {
        for (size_t j = 0; j < f->columns; j++)
        {
            int position=(i*f->rows)+j;
            switch (f->matrix[position])
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
}
void print_Forest(Forest *f){
    for (size_t i = 0; i < f->rows; i++)
    {
        for (size_t j = 0; j < f->columns; j++)
        {
            int position=(i*f->rows)+j;
            printf(" %d ",f->matrix[position]);
        }
        printf("\n");
    }
}
Forest* scatter_Forest(int n, int m)
{
    MPI_Datatype row;
    MPI_Type_contiguous( N , MPI_INT , &row);
    MPI_Type_commit( &row);
    int world_size;
    int rank;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    int execess_row=n%world_size;
    int row_for_each_process=n/world_size;
    Forest* f;
    
    if (rank==0)
    {
        f=createForest(n,m);
        initializeForestRandom(f);
        print_Forest_emoji(f);
    }
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
        if (rank==i && rank!=0)
        {
            recv_count=sendcount[i];
        }     
    }
    if (rank!=0)
    {
        f=createForest(recv_count,m);

    }
    MPI_Scatterv( f->matrix, sendcount, starting, row, f->matrix, recv_count, row, 0, MPI_COMM_WORLD);
    
    if(rank==0)
    {
        f->rows=sendcount[0];
        f->matrix=realloc(f->matrix, sizeof(int)*sendcount[0]*m);
    }
    return f;
}
Forest* compute(Forest* f,int rank, int world_size)
{
    int number_of_element=f->rows*f->columns;
    int t;
    Forest* toreturn=createForest(f->rows,f->columns);
    
    
    for (int i = 0; i < f->rows; i++)
    {
        int isempty=-1;
        for (int j = 0; j < f->columns; j++)
        {
            int ap=j+(i*f->columns);
            if ( toreturn->matrix[ap]!=2)
            {
                toreturn->matrix[ap]=f->matrix[ap];
            }
            
            if (f->matrix[ap]==2)
            {
                if ( (j!=0) && (f->matrix[ap-1]==1) )
                {
                    toreturn->matrix[ap-1]=2;
                }
                if ( (j!=f->columns-1) && (f->matrix[ap+1]==1) )
                {
                    toreturn->matrix[ap+1]=2;
                }
                if ( (i!=0) && (f->matrix[ap-f->columns]==1) )
                {
                    toreturn->matrix[ap-f->columns]=2;
                }
                if ( (i!=f->rows-1) && (f->matrix[ap+f->columns]==1) )
                {
                    toreturn->matrix[ap+f->columns]=2;
                }
                if ( (rank!=0) && (i==0) )
                {
                    isempty=1;
                    MPI_Send(&j,1,MPI_INT,rank-1,10,MPI_COMM_WORLD);
                }
                if ( (rank!=world_size-1) && (i==f->rows-1) )
                {
                    isempty=1;
                    MPI_Send(&j,1,MPI_INT,rank+1,10,MPI_COMM_WORLD);
                }                      
            }
            else
            {
                if (i==f->rows-1)
                {
                    if (rank!=world_size-1)
                    {
                        MPI_Send(&isempty,1,MPI_INT,rank+1,10,MPI_COMM_WORLD);
                    }  
                }
                if(i==0)
                {
                    if (rank!=0)
                    {
                        MPI_Send(&isempty,1,MPI_INT,rank-1,10,MPI_COMM_WORLD);
                    }    
                }
            }

            if (i==0 && rank!=0)
            {
                
                MPI_Recv(&t, 1, MPI_INT, rank-1, 10, MPI_COMM_WORLD, &status);
                if (t!=-1)
                {
                    if (f->matrix[t]==1)
                    {
                        toreturn->matrix[t]=2;
                    }
                }
            }
            if (i==f->rows-1 && rank!=world_size-1)
            {
                MPI_Recv(&t, 1, MPI_INT, rank+1, 10, MPI_COMM_WORLD, &status);
                if (t!=-1)
                {
                    if (f->matrix[t+(f->columns*i)]==1)
                    {
                        toreturn->matrix[t+(f->columns*i)]=2;
                    }
                }
            }
        }
    }
    return toreturn;
}