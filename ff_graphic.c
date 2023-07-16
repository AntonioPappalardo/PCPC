#include <stdio.h>
#include <stdlib.h>
#include<time.h>

#include <mpi.h>
#define N 8
#define pp 10
#define pf 3
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

MPI_Group used_group;
MPI_Group world_group;
MPI_Comm used_comm;
MPI_Status status;
int main(int argc, char** argv) {
    int world_size;
    int rank;
    MPI_Init(NULL, NULL);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_group(MPI_COMM_WORLD, &world_group);

    int in_rows=atoi(argv[1]);
    int in_columns=atoi(argv[2]);
    int generations=atoi(argv[3]);

    //if the number is less than the number of processes, the program uses only n_rows processes.
    int * ranks;
    if(in_rows<world_size){
        int ranks[in_rows];
        for(int i=0;i<in_rows;i++){
            ranks[i]=i;
        }
    
        MPI_Group_incl(world_group,in_rows,ranks,&used_group);    
        MPI_Comm_create(MPI_COMM_WORLD,used_group,&used_comm);        
        
    }else{
        used_comm=MPI_COMM_WORLD;
    }
    if(used_comm==MPI_COMM_NULL){
            MPI_Finalize();
            exit(0);
    }
    MPI_Comm_size(used_comm,&world_size);
    MPI_Comm_rank(used_comm,&rank);

    int* msg;
    Forest* forest;
    int t;
    MPI_Datatype row;
    MPI_Type_contiguous( in_columns , MPI_INT , &row);
    MPI_Type_commit( &row);
    sendcount=malloc(sizeof(int)*world_size);
    starting=malloc(sizeof(int)*world_size);
    
    forest= scatter_Forest(in_rows,in_columns);
    Forest* ff_compiled;
    ff_compiled=forest;
    for (size_t i = 0; i < generations; i++)
    {
        ff_compiled=compute(ff_compiled,rank,world_size); 
    }
    
    Forest* final_forest;
    int* pointer_matrix;
    if (rank==0)
    {
        final_forest=createForest(in_rows,in_columns);
        pointer_matrix=&final_forest->matrix[0];
    }
    printf("\n");
    MPI_Gatherv( &(ff_compiled->matrix[0]), sendcount[rank] , row , pointer_matrix , sendcount , starting , row , 0 , used_comm);
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
                case 3:
                    printf(" \U0001F4A8 ");
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
    MPI_Type_contiguous( m , MPI_INT , &row);
    MPI_Type_commit( &row);
    int world_size;
    int rank;
    MPI_Comm_size(used_comm, &world_size);
    MPI_Comm_rank(used_comm, &rank);
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
        if (rank==i)
        {
            recv_count=sendcount[i];
        }     
    }
    if (rank!=0)
    {
        f=createForest(recv_count,m);

    }
    MPI_Scatterv( f->matrix, sendcount, starting, row, f->matrix, recv_count, row, 0, used_comm);
    
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
    for (int i = 0; i < f->rows*f->columns; i++)
    {
        toreturn->matrix[i]=f->matrix[i];
    }
    
    
    for (int i = 0; i < f->rows; i++)
    {
        for (int j = 0; j < f->columns; j++)
        {
            int notburning=-1;
            int ap=j+(i*f->columns);
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
                    notburning=1;
                    MPI_Send(&j,1,MPI_INT,rank-1,10,used_comm);
                }
                if ( (rank!=world_size-1) && (i==f->rows-1) )
                {
                    notburning=1;
                    MPI_Send(&j,1,MPI_INT,rank+1,10,used_comm);
                }                   
            }
            else
            {
                srand(time(0));
                if (f->matrix[ap]==1 && (rand()%100<pf))
                {
                    f->matrix[ap]=2;
                }
                if (f->matrix[ap]==0 && (rand()%100<pp))
                {
                    f->matrix[ap]=1;
                }
                if (i==f->rows-1)
                {
                    if (rank!=world_size-1)
                    {
                        MPI_Send(&notburning,1,MPI_INT,rank+1,10,used_comm);
                    }
                }
                if(i==0)
                {
                    if (rank!=0)
                    {
                        MPI_Send(&notburning,1,MPI_INT,rank-1,10,used_comm);
                    }
                }
            }

            if (i==0 && rank!=0)
            {
                
                MPI_Recv(&t, 1, MPI_INT, rank-1, 10, used_comm, &status);
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
                MPI_Recv(&t, 1, MPI_INT, rank+1, 10, used_comm, &status);

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