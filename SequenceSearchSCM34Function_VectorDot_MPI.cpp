#include <iostream>
#include <omp.h>
#include <fstream>
#include "mpi.h"

#pragma GCC optimize ("unroll-loops")
inline int SequenceSearch(int ThNum, int n, int a1, int a2, int a3, int a4, int a5){
    bool GOOD;  //If good then we found a "sequence"
    int DotSum=0;  //dot product sum
    int Seq[n]={0};  //initialize the sequence to all 0
    Seq[1]=1;  //the first entry is 0, the second is 1
    int COUNT=0;

    int M[2*n]={0};  //Set up the matrix
    int C[n-1]={0};  //set up the "correlation". Need one less entry than sequence length

    //Open a file based on the core # as to avoid parallel "race conditions"
    std::ofstream DataFile;
    DataFile.open("SCM" +  std::to_string(ThNum) + ".txt",std::fstream::app);

                                        Seq[2]=a1;
                                        Seq[3]=a2;
                                        Seq[4]=a3;
                                        Seq[5]=a4;
                                        Seq[6]=a5;


    #pragma omp critical  //avoid race conditions
    {  //write to command line
        std::cout << "\n" << ThNum+1 << "/32   Entering Search With " << Seq[0] << Seq[1] << Seq[2] << Seq[3] << Seq[4] << Seq[5] << Seq[6] << "\n";
    }


                //for (int a4=-1;a4<2;a4=a4+2){
                   // for (int a5=-1;a5<2;a5=a5+2){
                        for (int a6=-1;a6<2;a6=a6+2){
                                Seq[7]=a6;
                            for (int a7=-1;a7<2;a7=a7+2){
                                Seq[8]=a7;
                                for (int a8=-1;a8<2;a8=a8+2){
                                    Seq[9]=a8;
                                    for (int a9=-1;a9<2;a9=a9+2){
                                        Seq[10]=a9;
                                            for (int a10=-1;a10<2;a10=a10+2){
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
                                                //disp seq
                                                #pragma omp critical  //avoid race conditions
                                                {  //write to command line
                                                    std::cout << "\n";
                                                    for (int i=0;i<n;i++){
                                                        std::cout << Seq[i] <<" " ;
                                                    }
                                                     std::cout << "\n";
                                                    for (int i=0;i<n-1;i++){
                                                        std::cout << C[i] <<" " ;
                                                    }
                                                    //std::cout << "\n";
                                                }

                                                //write to file
                                                for (int i=0;i<n;i++){
                                                    DataFile << Seq[i] <<" " ;
                                                }
                                                DataFile << "\n";


                                            }
                                                                                                                                            }
                                                                                                                                        }
                                                                                                                                }
                                                                                                                            }
                                                                                                                    }
                                                                                                                }
                                                                                                            }
                                                                                                        }
                                                                                                    }
                                                                                                }
                                                                                            }
                                                                                        }
                                                                                    }
                                                                                }
                                                                            }
                                                                        }
                                                                    }
                                                                }
                                                            }
                                                        }
                                                    }
                                                }
                                            }
                                    }
                                }
                            }
                        }
                    //}
                //}
                DataFile.close();
    return COUNT;
}





int main(int argc, char* argv[])
{

	int numprocs, rank, namelen;
	char processor_name[MPI_MAX_PROCESSOR_NAME];
	int iam =0, np = 1;
	
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Get_processor_name(processor_name, &namelen);
	
//while(1){	
const int n = 34;
const int CORES = 4;
//Direct 8 threads
int COUNTS[32] = {0};
    #pragma omp parallel num_threads(CORES) //Run on 16 threads
    {

    //set up the 16 threads to start at different 4-bit positions
    float A[32][5]=
    {
    {1,1,1,1,1},
    {1,1,1,1,-1},
    {1,1,1,-1,1},
    {1,1,1,-1,-1},
    {1,1,-1,1,1},
    {1,1,-1,1,-1},
    {1,1,-1,-1,1},
    {1,1,-1,-1,-1},
    {1,-1,1,1,1},
    {1,-1,1,1,-1},
    {1,-1,1,-1,1},
    {1,-1,1,-1,-1},
    {1,-1,-1,1,1},
    {1,-1,-1,1,-1},
    {1,-1,-1,-1,1},
    {1,-1,-1,-1,-1},
    {-1,1,1,1,1},
    {-1,1,1,1,-1},
    {-1,1,1,-1,1},
    {-1,1,1,-1,-1},
    {-1,1,-1,1,1},
    {-1,1,-1,1,-1},
    {-1,1,-1,-1,1},
    {-1,1,-1,-1,-1},
    {-1,-1,1,1,1},
    {-1,-1,1,1,-1},
    {-1,-1,1,-1,1},
    {-1,-1,1,-1,-1},
    {-1,-1,-1,1,1},
    {-1,-1,-1,1,-1},
    {-1,-1,-1,-1,1},
    {-1,-1,-1,-1,-1},
    };

    //get the current thread # and send the respective starting bits
    int thNum = omp_get_thread_num();

    //for (int k=0;k<32/CORES;k++)
    //for (int k=rank*8+thNum;k<=((rank*8)+7);k=k+4)
    for (int k=rank*6+thNum;k<((rank*6)+6);k=k+4)
        {
            //COUNTS[thNum+k*CORES] = SequenceSearch(thNum+k*CORES, n, A[thNum+k*CORES][0], A[thNum+k*CORES][1], A[thNum+k*CORES][2],  A[thNum+k*CORES][3], A[thNum+k*CORES][4]);
            if (k<32){
            COUNTS[k] = SequenceSearch(k, n, A[k][0], A[k][1], A[k][2],  A[k][3], A[k][4]);
			}
        }


    }

    int TOTAL = 0;
    for (int i=0;i<32;i++){
        TOTAL += COUNTS[i];
    }
    std::cout<< "\n" << TOTAL;

    //For loop threading.
/*
    int A[32][5]=
    {
    {1,1,1,1,1},
    {1,1,1,1,-1},
    {1,1,1,-1,1},
    {1,1,1,-1,-1},
    {1,1,-1,1,1},
    {1,1,-1,1,-1},
    {1,1,-1,-1,1},
    {1,1,-1,-1,-1},
    {1,-1,1,1,1},
    {1,-1,1,1,-1},
    {1,-1,1,-1,1},
    {1,-1,1,-1,-1},
    {1,-1,-1,1,1},
    {1,-1,-1,1,-1},
    {1,-1,-1,-1,1},
    {1,-1,-1,-1,-1},
    {-1,1,1,1,1},
    {-1,1,1,1,-1},
    {-1,1,1,-1,1},
    {-1,1,1,-1,-1},
    {-1,1,-1,1,1},
    {-1,1,-1,1,-1},
    {-1,1,-1,-1,1},
    {-1,1,-1,-1,-1},
    {-1,-1,1,1,1},
    {-1,-1,1,1,-1},
    {-1,-1,1,-1,1},
    {-1,-1,1,-1,-1},
    {-1,-1,-1,1,1},
    {-1,-1,-1,1,-1},
    {-1,-1,-1,-1,1},
    {-1,-1,-1,-1,-1},
    };

    #pragma omp parallel for num_threads(8) schedule(dynamic,1)
    for (int StartPoint=0;StartPoint<n;StartPoint++){
        SequenceSearch(StartPoint, n, A[StartPoint][0], A[StartPoint][1], A[StartPoint][2], A[StartPoint][3], A[StartPoint][4]);
    }
*/
//}
 return 0;

}

