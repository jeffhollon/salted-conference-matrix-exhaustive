#include <iostream>
#include <omp.h>
#include <fstream>

void SequenceSearch(int ThNum, int n, int a1, int a2, int a3){
    bool GOOD;  //If good then we found a "sequence"
    int DotSum;  //dot product sum
    int Seq[n]={0};  //initialize the sequence to all 0
    Seq[1]=1;  //the first entry is 0, the second is 1

    int M[n][n]={0};  //Set up the matrix
    int C[n-1]={0};  //set up the "correlation". Need one less entry than sequence length

    //Open a file based on the core # as to avoid parallel "race conditions"
    std::ofstream DataFile;
    DataFile.open("SCM" +  std::to_string(ThNum) + ".txt",std::fstream::app);

    /*
    I know this is ugly, but this is the simplest form for searching.
    There are other/better ways,but this was easiest for me to code in C++.
    This area loops through all possible sequences of entries {-1,+1}
    */
                for (int a4=-1;a4<2;a4=a4+2){
                    for (int a5=-1;a5<2;a5=a5+2){
                        for (int a6=-1;a6<2;a6=a6+2){
                            for (int a7=-1;a7<2;a7=a7+2){
                                for (int a8=-1;a8<2;a8=a8+2){
                                    for (int a9=-1;a9<2;a9=a9+2){
                                            for (int a10=-1;a10<2;a10=a10+2){
                                                for (int a11=-1;a11<2;a11=a11+2){
                                                    for (int a12=-1;a12<2;a12=a12+2){
                                                        for (int a13=-1;a13<2;a13=a13+2){
                                                            for (int a14=-1;a14<2;a14=a14+2){
                                                                for (int a15=-1;a15<2;a15=a15+2){
                                                                    for (int a16=-1;a16<2;a16=a16+2){
                                                                        for (int a17=-1;a17<2;a17=a17+2){
                                                                            for (int a18=-1;a18<2;a18=a18+2){
                                                                                for (int a19=-1;a19<2;a19=a19+2){
                                                                                    for (int a20=-1;a20<2;a20=a20+2){
                                                                                        for (int a21=-1;a21<2;a21=a21+2){
                                                                                            for (int a22=-1;a22<2;a22=a22+2){
                                                                                                for (int a23=-1;a23<2;a23=a23+2){
                                                                                                    for (int a24=-1;a24<2;a24=a24+2){
                                                                                                        for (int a25=-1;a25<2;a25=a25+2){
                                                                                                            for (int a26=-1;a26<2;a26=a26+2){
                                                                                                                for (int a27=-1;a27<2;a27=a27+2){
                                                                                                                    for (int a28=-1;a28<2;a28=a28+2){

                                        //Build the Sequence
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
                                        Seq[12]=a11;
                                        Seq[13]=a12;
                                        Seq[14]=a13;
                                        Seq[15]=a14;
                                        Seq[16]=a15;
                                        Seq[17]=a16;
                                        Seq[18]=a17;
                                        Seq[19]=a18;
                                        Seq[20]=a19;
                                        Seq[21]=a20;
                                        Seq[22]=a21;
                                        Seq[23]=a22;
                                        Seq[24]=a23;
                                        Seq[25]=a24;
                                        Seq[26]=a25;
                                        Seq[27]=a26;
                                        Seq[28]=a27;
                                        Seq[29]=a28;


                                        //Place First Row in the matrix
                                          for (int i=1;i<n;i++)
                                          {
                                              M[0][i]=Seq[i];
                                          }

                                        //Place remaining entries as cyclic shifts. For "sequences"
                                        //this is equivalent to right cyclic shifts of all entries.
                                          for (int i=1;i<n;i++)
                                          {
                                              for (int j=0;j<n;j++)
                                              {
                                                      if (j>0){
                                                        M[i][j]=M[i-1][j-1];
                                                      }

                                                      if (j==0){
                                                        M[i][j]=M[i-1][n-1];
                                                      }
                                              }
                                          }

                                          //Dot product of first row with all other rows
                                            GOOD=1;  //assume good
                                            for (int i=1;i<n;i++)
                                            {
                                                DotSum = 0;

                                                for (int j=0;j<n;j++)
                                                {
                                                    DotSum = DotSum + M[0][j]*M[i][j];
                                                }
                                                if (abs(DotSum)>2){
                                                    GOOD=0;  //if any of the dot products is larger than 2 then stop searching
                                                    break;
                                                }
                                                C[i]=DotSum;
                                            }


                                            if(GOOD){   //only if the dot products are less than 2
                                                //disp seq
                                                #pragma omp critical  //avoid race conditions
                                                {  //write to command line
                                                std::cout << "\n";
                                                for (int i=0;i<n;i++){
                                                    std::cout << Seq[i] <<" " ;
                                                }
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
                DataFile.close();
}





int main(int argc, char* argv[])
{


    #pragma omp parallel num_threads(8) //Run on 8 threads
    {

    //set up the 8 threads to start at different 3-bit positions
    int A[8][4]=
    {
    {1,1,1},
    {1,1,-1},
    {1,-1,1},
    {1,-1,-1},
    {-1,1,1},
    {-1,1,-1},
    {-1,-1,1},
    {-1,-1,-1},
    };

    //get the current thread # and send the respective starting bits
    int thNum = omp_get_thread_num();
    SequenceSearch(thNum, 30, A[thNum][0], A[thNum][1], A[thNum][2]);


    }

 return 0;

}

