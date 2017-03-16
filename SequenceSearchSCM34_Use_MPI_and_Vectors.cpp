#include <iostream>
#include <omp.h>
#include <fstream>
#include "mpi.h"
#include <time.h>
#include <ctime>
#include <vector>

#pragma GCC optimize ("unroll-loops")
inline std::vector<int> SequenceSearch(int ThNum, int n, int a1, int a2, int a3, int a4, int a5, int a6, int a7, int a8, int a9, int a10, int rank, int numprocs, std::vector<int> LIST){

    bool GOOD;  //If good then we found a "sequence"
    int DotSum=0;  //dot product sum
    int Seq[n]={0};  //initialize the sequence to all 0
    int Seq_Write[n]={0};
    Seq[1]=1;  //the first entry is 0, the second is 1
    int COUNT=0;

    int M[2*n]={0};  //Set up the matrix
    int C[n-1]={0};  //set up the "correlation". Need one less entry than sequence length

    Seq[2]=a1;
    Seq[3]=a2;
    Seq[4]=a3;
    Seq[5]=a4;
    Seq[6]=a5;
    Seq[7]=a6;
    Seq[8]=a7;
    Seq[9]=a8;
    Seq[10]=a9;
    Seq[11]=a10;


    for (int a11=-1;a11<2;a11=a11+2){
        Seq[12]=a11;
        for (int a12=-1;a12<2;a12=a12+2){
            Seq[13]=a12;
            for (int a13=-1;a13<2;a13=a13+2){
                Seq[14]=a13;
                for (int a14=-1;a14<2;a14=a14+2){
                    Seq[15]=a14;
                    for (int a15=-1;a15<2;a15=a15+2){
                        Seq[16]=a15;
                        for (int a16=-1;a16<2;a16=a16+2){
                             Seq[17]=a16;
                            for (int a17=-1;a17<2;a17=a17+2){
                                Seq[18]=a17;
                                for (int a18=-1;a18<2;a18=a18+2){
                                     Seq[19]=a18;
                                    for (int a19=-1;a19<2;a19=a19+2){
                                        Seq[20]=a19;
                                        for (int a20=-1;a20<2;a20=a20+2){
                                            Seq[21]=a20;
                                            for (int a21=-1;a21<2;a21=a21+2){
                                                Seq[22]=a21;
                                                for (int a22=-1;a22<2;a22=a22+2){
                                                    Seq[23]=a22;
                                                    for (int a23=-1;a23<2;a23=a23+2){
                                                        Seq[24]=a23;
                                                        for (int a24=-1;a24<2;a24=a24+2){
                                                            Seq[25]=a24;
                                                            for (int a25=-1;a25<2;a25=a25+2){
                                                                Seq[26]=a25;
                                                                for (int a26=-1;a26<2;a26=a26+2){
                                                                    Seq[27]=a26;
                                                                    for (int a27=-1;a27<2;a27=a27+2){
                                                                        Seq[28]=a27;
                                                                        for (int a28=-1;a28<2;a28=a28+2){
                                                                            Seq[29]=a28;
                                                                                for (int a29=-1;a29<2;a29=a29+2){
                                                                                    Seq[30]=a29;
                                                                                    for (int a30=-1;a30<2;a30=a30+2){
                                                                                        Seq[31]=a30;
                                                                                            for (int a31=-1;a31<2;a31=a31+2){
                                                                                                Seq[32]=a31;
                                                                                                for (int a32=-1;a32<2;a32=a32+2){
                                                                                                    Seq[33]=a32;
                //Place First Row in the matrix
                for (int i=0;i<2*n;i++)
                {
                  M[i]=Seq[i%n];
                }

                GOOD=1;  //assume good
                for (int i=1;i<n;i++)
                {
                    C[i-1]=0;
                    for (int j=0;j<n;j++)
                    {
                        C[i-1] += M[j]*M[j+i];
                    }
                    if (abs(C[i-1])>2){
                        GOOD=0;
                        break;
                    }
                }

                if(GOOD){   //only if the dot products are less than 2
                    COUNT++;

                    #pragma omp critical  //avoid race conditions
                    {//Write to vector
                        for (int i=0;i<n;i++){
                           LIST.push_back(Seq[i]);
                        }
                    }
                }

}}}}}}}}}}}}}}}}}}}}}}  //End Sequence Loop

    return LIST;
}


int main(int argc, char* argv[])
{
	clock_t t1,t2;
	t1=clock();

	std::time_t startt, endt;
	long delta = 0;
	startt = std::time(NULL);

	int numprocs, rank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int iam =0, np = 1;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);

    const int n = 34;
    const int CORES = 4;

    int COUNTS[1024] = {0};
    std::vector<int> LIST;
    LIST.reserve(1);
    LIST[0]=5;

    #pragma omp parallel num_threads(CORES) //Run on 16 threads
    {
        std::vector<int> LIST2;
        LIST2.reserve(1);
        LIST2[0]=5;

		int A[1024][10]={0};
		int count=0;
		for (int A1=-1;A1<2;A1=A1+2){
		for (int A2=-1;A2<2;A2=A2+2){
		for (int A3=-1;A3<2;A3=A3+2){
		for (int A4=-1;A4<2;A4=A4+2){
		for (int A5=-1;A5<2;A5=A5+2){
		for (int A6=-1;A6<2;A6=A6+2){
		for (int A7=-1;A7<2;A7=A7+2){
		for (int A8=-1;A8<2;A8=A8+2){
		for (int A9=-1;A9<2;A9=A9+2){
		for (int A10=-1;A10<2;A10=A10+2){
			A[count][0]=A1;
			A[count][1]=A2;
			A[count][2]=A3;
			A[count][3]=A4;
			A[count][4]=A5;
			A[count][5]=A6;
			A[count][6]=A7;
			A[count][7]=A8;
			A[count][8]=A9;
			A[count][9]=A10;
			count++;
		}}}}}}}}}}  //Finish Matrix A to Parallelize on

        //get the current thread # and send the respective starting bits
        int thNum = omp_get_thread_num();

        for (int k=thNum+4*rank;k<1024;k=k+24)
            {
                LIST2 = SequenceSearch(k, n, A[k][0], A[k][1], A[k][2],  A[k][3], A[k][4], A[k][5], A[k][6], A[k][7], A[k][8], A[k][9], rank, numprocs, LIST2);
            }

        #pragma omp critical
        LIST.insert(LIST.end(),LIST2.begin(),LIST2.end());
    } //EXIT THREAD OPENMP

    int ListSize[numprocs];
    int temp;

    if (rank==0){
        ListSize[0]=LIST.size();
        for (int z=1; z<numprocs; z++){
            MPI_Recv(&temp,1,MPI_INT,z,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            ListSize[z]=temp;
        }
    } else {
        int temp2;
        temp2 = LIST.size();
        MPI_Send(&temp2,1,MPI_INT,0,1,MPI_COMM_WORLD);
    }

    if (rank==0){
        for (int i=0;i<numprocs;i++){
            std::cout << ListSize[i]/n << "\n";
        }
    }

    std::vector<int> TempList;
    std::vector<int> TempList2;

    if (rank==0){
        //Print the current list
        std::cout << "\nThe Good Vectors are\n";
        for (int i=0;i<ListSize[0];i++)
        {
            if (i%n==0){
                std::cout<<"\n";
            }
            std::cout << LIST[i];
        }
        //Recv the others one at a time and print
        for (int z=1; z<numprocs; z++){
            TempList2.resize(ListSize[z]);
            MPI_Recv(&TempList2.front(),ListSize[z],MPI_INT,z,1,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
            std::cout<<"\n";
            for (int i=0; i<ListSize[z]; i++){
                if (i>0 && i%n==0){
                    std::cout<<"\n";
                }
                std::cout << TempList2[i];
            }
        }
    } else {  //Send the LIST
        MPI_Send(&LIST.front(),LIST.size(),MPI_INT,0,1,MPI_COMM_WORLD);
    }

	endt = std::time(NULL);

	t2=clock();
	float diff ((float)t2-(float)t1);
	delta = endt-startt;
    if (rank==0){
        std::cout << '\n' << delta << '\n';
    }

	 return 0;
}

