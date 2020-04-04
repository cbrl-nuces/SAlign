#include "Alignment.h"
#include <iostream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>

using namespace std;

exception er;

int main(int argc, char* argv[])
{
    double lambda = 0.2; //a is alpha that controlls the factor of edgeweights
    double alpha=0.1;
    int degree = 10; //controlls the step of making the skeleton
    string n;
    double beta = 0.7;

    char* name1; //name of the first network
    char* name2; //name of second network
    char* blastFile;
    char* structuralFile;
    try
    {
        if(argc < 3) {
            cout << "There should be two files as input!" <<endl;
            return -1;
        }
        else //input arguments
        {
            int i = 1; //counter for input parameters
            name1 = argv[ i++ ]; //first network
            name2 = argv[ i++ ]; //second network
	    cout<<name1<<" :: "<<name2<<endl;

			while (i<argc) //check all the input parameters
            {
                if ( ( strlen(argv[i]) == 2 ) && ( argv[i][0]=='-' ) && ( i + 1 < argc) ) //wether or not the parameter has started with '-' and has a value
                {
                    i++; //to read the value of input parameter
                    if( argv[i-1][1]=='l' )
                    {
                        lambda = atof(argv[i]);
                        if( lambda <0 || lambda > 1) //the value of a should be between zero and one
                        {
                            cout << "Error : value of a must be between 0 and 1" << endl;
                            return -1;
                        }
                    }
                    else if( argv[i-1][1]=='d')
                    {
                        degree = atoi(argv[i]);
                        if( degree > 100)
                        {
                            cout << "Error : value of t must be between 0 and 100" << endl;
                            return -1;
                        }
                    }
                    else if( argv[i-1][1]=='n')
                    {
                        n = argv[i];
                    }
                    else if( argv[i-1][1]=='t')
                    {
                        beta = atof(argv[i]);
                        if( beta > 1 || beta < 0)
                        {
                            cout << "Error : value of beta must be between 0 and 1" << endl;
                            return -1;
                        }
                    }
                    else if(argv[i-1][1]=='b')
                    {
                        blastFile=argv[i];
                    }
                    else if(argv[i-1][1]=='s')
                    {
                        structuralFile=argv[i];
                    }
                    else if(argv[i-1][1]=='a')
                    {
                        alpha = atof(argv[i]);
                    }

                    i++;// to reach the next input parameter if there is
                }
                else
                {
                    cout <<  strlen(argv[i])<<i+1 << argc << endl;
                    cout << "Error in argument : " << argv[i] << endl;
                    return -1;
                }
            }
		} //end else

        cout << "\n=============== SAlign ==================\n";
        cout << "alignment num: " << n << endl;
        //making the name for output file
        string strm;
        strm = n + '-' + name1 + '_' + name2;
        cout << "output file name: " << strm << endl;
        cout << "==========================================\n" << endl;


        //construct the networks
        Network network1(name1);
        Network network2(name2);
        bool reverse = false; //it means size of first input network is bigger than second one

        if(network1.size > network2.size)
            reverse = true;
        //making the skeletons of the networks

        network1.makeSkeleton(degree);
        network2.makeSkeleton(degree);



        //align two networks with each other
        Alignment alignment( network1, network2);
        if(alpha!=1) {
            alignment.readblast(blastFile);
            alignment.readTmalignScore(structuralFile);
        }

        alignment.align(lambda, alpha, beta, name1, name2);

        alignment.outputEvaluation(strm, lambda, alpha, beta);
        alignment.outputAlignment(strm);

    }
	catch(exception &e)
	{
		cout << "Error in arguments or input files!" << endl;
		e.what();
		return -1;
	}
    return 0;
}
