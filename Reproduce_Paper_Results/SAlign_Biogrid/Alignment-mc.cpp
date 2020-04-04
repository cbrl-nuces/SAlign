// monte carlo chnages ==> 306 to 416 (code in block "changes started" to "changes end" in rewritten)
// file path is mentioned in 641 677 (where output alignment files can be stored)

#include "Alignment.h"
#include <fstream>
#include <cstdlib>
#include <iomanip>
#include <cmath>
#include <ctime> 

#include <chrono>
#include <cstdint>
#include <iostream>

using namespace std;
//constructor
//finds the smaller network and the maximum degree of the input networks
//Inputs are two files of networks net1 and net2
Alignment::Alignment( Network net1, Network net2)
{
    //compare networks to find the biggest one
    if( net1.size > net2.size )
	{
		reverse = true;
		network1 = net2;
		network2 = net1;
	}
	else
	{
		reverse = false;
		network1 = net1;
		network2 = net2;
	}

	//maximum degree of the network
    if(network1.maxDeg > network2.maxDeg)
		maxDeg = network1.maxDeg;
	else
		maxDeg = network2.maxDeg;

    blast = new float*[network1.size];
    structralScore = new float*[network1.size];
    bioScore = new float*[network1.size];

    for (int c=0; c<network1.size; c++) {
        blast[c]=new float[network2.size];
        structralScore[c] = new float[network2.size];
        bioScore[c] = new float[network2.size];
    }
    for (int c1=0; c1<network1.size; c1++) {
        for (int c2=0; c2<network2.size; c2++) {
            blast[c1][c2]=0;
            structralScore[c1][c2] = 0;
            bioScore[c1][c2] = 0;
        }
    }
}

void Alignment::readblast(string blastFile) {
  const string BLAST_FILES_PATH = "blast/";

    float ** temp = new float*[network1.size];
    for (int c=0; c<network1.size; c++) {
        temp[c]=new float[network2.size];
    }
    for (int c1=0; c1<network1.size; c1++) {
        for (int c2=0; c2<network2.size; c2++) {
            temp[c1][c2]=0;
        }
    }

    float max = 0 ;
    //blast values
    string blast_file_path = BLAST_FILES_PATH + blastFile;
    ifstream inputFile;
    string token1,token2,line;
    float token3;
    inputFile.open(blast_file_path.c_str());
    while (getline(inputFile, line)) {
        istringstream tokenizer(line);
        getline(tokenizer, token1, '\t');
        getline(tokenizer, token2, '\t');
        tokenizer >> token3;
        if(max<token3) max = token3;
        temp[network1.mapName[token1]][network2.mapName[token2]]=token3;
    }

    //normalize between zero and 1
    for (int c1=0; c1<network1.size; c1++)
        for (int c2=0; c2<network2.size; c2++)
            blast[c1][c2] = temp[c1][c2]/max;
}

void Alignment::readTmalignScore(string structuralFile) {
  const string TM_ALIGN_PATH = "tmalign/";

    float ** temp = new float*[network1.size];
    for (int c=0; c<network1.size; c++) {
        temp[c]=new float[network2.size];
    }
    for (int c1=0; c1<network1.size; c1++) {
        for (int c2=0; c2<network2.size; c2++) {
            temp[c1][c2]=0;
        }
    }

    float max = 0 ;
    //blast values
    string scoring_file_path = TM_ALIGN_PATH + structuralFile;
    ifstream inputFile;
    string token1,token2,line;
    float token3;
    inputFile.open(scoring_file_path.c_str());
    while (getline(inputFile, line)) {
        istringstream tokenizer(line);
        getline(tokenizer, token1, '\t');
        getline(tokenizer, token2, '\t');
        tokenizer >> token3;
        if(max<token3) max = token3;
        temp[network1.mapName[token1]][network2.mapName[token2]]=token3;
    }

    //normalize between zero and 1
    for (int c1=0; c1<network1.size; c1++){
        for (int c2=0; c2<network2.size; c2++){
            structralScore[c1][c2] = temp[c1][c2];
          }
        }
}

uint64_t timeSinceEpochMillisec() {
  using namespace std::chrono;
  return duration_cast<microseconds>(system_clock::now().time_since_epoch()).count();
}

//produce a mapping between nodes of two network with respect to input parameter a.
//Input parameter a acontrols the factor edgeWeight in assigning the scores to the nodes. a should be between 0 and 1.
void Alignment::align(double lambda, double alpha, double beta, string a, string b)
{
    bool flag;  //check wether or not all the nodes of the smaller network are aligned?
   // cout << "Beta: " << beta << endl;
    //temporary
    float temp;
    float a1,a11;
    float a2,a22;
    float MINSCORE = -100000;
    int coeff;

    if(network2.numOfEdge>network1.numOfEdge) {
        coeff = network2.numOfEdge/network1.numOfEdge;
    }
    else {
        coeff = network1.numOfEdge/network2.numOfEdge;

    }
    int maxNode; // node with max score
    bool *alignNodes1 = new bool[network1.size]; //aligned nodes of the smaller network
    bool *alignNodes2 = new bool[network2.size]; //aligned nodes of the bigger network
    alignment = new int[network1.size]; //alignment array
    float *nodeScore1 = new float[network1.size]; //scores of nodes of smaller network
    float *nodeScore2 = new float[network2.size]; //scores of nodes of bigger network
    double **alignScore = new double*[network1.size]; //this matrix contains the score of each matching pair
    double **similarityScore = new double*[network1.size];
    int *best = new int[network1.size]; //array of best align scores
    float ss;
    //initial values
    for(int c1=0; c1<network1.size; c1++)
        alignScore[c1]=new double[network2.size];
    for(int c1=0; c1<network1.size; c1++)
        similarityScore[c1]=new double[network2.size];
    for(int c1=0; c1<network1.size; c1++)
        alignNodes1[c1]=false;
    for(int c1=0; c1<network2.size; c1++)
        alignNodes2[c1]=false;
    for(int c1=0; c1<network1.size; c1++)
        alignment[c1]=-1;
    for(int c1=0; c1<network1.size; c1++)
        best[c1]=0;

    ofstream NS;
    //initialize nodeScore fro both networks
    for(int c1=0; c1< network1.size; c1++)
        nodeScore1[c1]=(1-lambda)*network1.nodeWeight[c1];
    for(int c1=0; c1< network2.size; c1++)
        nodeScore2[c1]=(1-lambda)*network2.nodeWeight[c1];

    //find max score
    //finding the nodescore
    for (int c1=0; c1<network1.size; c1++){
        for (int c2=0; c2<network1.size; c2++)
            nodeScore1[c1]+= lambda*network1.edgeWeight[c1][c2];
    }
    for (int c1=0; c1<network2.size; c1++){
        for (int c2=0; c2<network2.size; c2++)
            nodeScore2[c1] += lambda*network2.edgeWeight[c1][c2];
     }

    //======first network
    float max = -10000;
    for (int c1=0; c1<network1.size; c1++) {
        if (max < nodeScore1[c1]) {
            max = nodeScore1[c1];
        }
    }

    //====== second network
    for (int c1=0; c1<network2.size; c1++) {
        if (max < nodeScore2[c1]) {
            max = nodeScore2[c1];
        }
    }

    //normalize with respect to max
    for (int c1=0; c1<network1.size; c1++) {
        nodeScore1[c1] = nodeScore1[c1]/max;
    }

    for (int c1=0; c1<network2.size; c1++) {
        nodeScore2[c1] = nodeScore2[c1]/max;
    }
    //END of normalization

    //finding the alignscore
    for(int c1=0; c1<network1.size; c1++){
        for(int c2=0; c2<network2.size; c2++){
          alignScore[c1][c2] = (nodeScore1[c1]>nodeScore2[c2])? nodeScore2[c2]:nodeScore1[c1];
          alignScore[c1][c2] = alpha * (alignScore[c1][c2]);

	  //cout<<alignScore[c1][c2]<<endl;

          if ((network1.getName(c1) == "P01823" && network2.getName(c2) == "A2NI60") || (network1.getName(c1) == "P01642" && network2.getName(c2) == "Q6GMX6") || (network1.getName(c1) == "Q5R1F5" && network2.getName(c2) == "A0JD37")) {
            cout << network1.getName(c1) << " " << network2.getName(c2) << " " << alignScore[c1][c2]<< " " << structralScore[c1][c2] << " " << bioScore[c1][c2] << endl;
          }
        }
      }

    for(int c1=0; c1<network1.size; c1++){
        for(int c2=0; c2<network2.size; c2++){
          bioScore[c1][c2] += beta*blast[c1][c2]; //adding similarity
          bioScore[c1][c2] += (1-beta)*structralScore[c1][c2]; //adding structure
          bioScore[c1][c2] = (1-alpha)*bioScore[c1][c2];
          alignScore[c1][c2] += bioScore[c1][c2];
        }
      }
      // cout << "Network1:  " << network1.name << network1.size << endl;
      // cout << "Network2:  " << network2.name << network2.size << endl;
      // exit(0);
    int counter = 0;
    flag=true;
    int temp1,temp2;
    int best_of = 5; // can be 5 or 10
    int array [best_of]; // contains the index of top 5 or 10 Alignment scores
    float arr_2 [best_of]; // contains top 5 or 10 Alignment scores
    int counteralign=0;
    int selected = 0; // selected index
    float number = 0.0; // fpr random number storage
    int select_zero = 0; // how many time zeroth index is selected (best one)
    int var = 0; // best of 10 selection
    int index = 0; // index of best selection
    cout<<"\n Alignment Start...\n\n";
    while(flag){ //there is some unaligned nodes in determined iteration
        //find the maximum value of each row of alignscore and save it in array "best"
	for(int c1=0; c1<network1.size; c1++)
        {
            if(!alignNodes1[c1]){
                temp=MINSCORE;
                for(int c2=0; c2<=network2.size; c2++){
                    if(temp<alignScore[c1][c2] && !alignNodes2[c2]){
                        if(alignScore[c1][c2]==temp) {
                            temp1 = (network1.deg[c1]>network2.deg[c2]) ? network2.deg[c2]/network1.deg[c1]:network1.deg[c1]/network2.deg[c2];
                            temp2 = (network1.deg[c1]>network2.deg[best[c1]]) ? network2.deg[best[c1]]/network1.deg[c1]:network1.deg[c1]/network2.deg[best[c1]];
                            if(temp1 > temp2) {
                                best[c1]=c2;
                                temp = alignScore[c1][c2];
                            }
                        }
                        else {
                            best[c1]=c2;
                            temp = alignScore[c1][c2];
                        }
                    }
                }
            }
        }

        //doing the alignment
        //find the maximum value of array "best" that means the best score in matrix "alignScore"
        temp=MINSCORE;
        flag=false;

	for(int c1=0; c1<network1.size; c1++){
            if(temp<alignScore[c1][best[c1]] && !alignNodes1[c1] && !alignNodes2[best[c1]]){ //=
                flag=true; //still there is node that is not aligned
                if(alignScore[c1][best[c1]]==temp) {

                    if(network1.deg[c1] > network1.deg[maxNode]) {
                         maxNode = c1;
                        temp = alignScore[c1][best[c1]];
                    }
                }
                else {
                    temp = alignScore[c1][best[c1]];
                    maxNode = c1;
                }
            }
         }

//   Top n elements code  
	 std::fill_n(array, best_of, -1);
         for(int c2=0; c2<=network2.size; c2++)
	 {
             if(!alignNodes2[c2])
             {
		 for (int m=0; m<best_of; m++)
		 {
		      if(alignScore[maxNode][c2] >= alignScore[maxNode][array[m]])
		      { 
		          for (int k=best_of-1; k>m; k--)
			  {
			      array[k] = array[k-1];
			  }
			  array[m] = c2;
			  break;
			}
		   } 
	       }
	 }


	std::fill_n(arr_2, best_of, 0.0);
	float sum_of_best = 0.0;
	float best_of_best = 0;

	// Array filling and getting best of top 5 or top 10
	for (int i=0; i<best_of; i++)
	{
	    if (alignScore[maxNode][array[i]] >= 0.0)
	    {
	        arr_2[i]=alignScore[maxNode][array[i]];
		if (arr_2[i] > best_of_best)
		{
		    best_of_best = arr_2[i];
		}
	    }
	}


	//Monte-carlo 
	for (int i=0; i<best_of; i++)
	{
	    arr_2[i] = exp((arr_2[i]-best_of_best)/0.10); // exp(best) - exp(best-each_element)
	    sum_of_best += arr_2[i];
	}

	//normalization between 0 to max to get the sum = 1 (probablities of selection)
	for (int i=0; i<best_of; i++)
	{
	    arr_2[i]=arr_2[i]/sum_of_best;
	}

	// summation (cummulative sum) e.f 0.1, 0.3, 0.6 ==> 0.1, 0.4, 1.0 
	for (int i=1; i<best_of; i++)
	{
	    arr_2[i]+=arr_2[i-1];
	}

	int selection_count [best_of] = {0}; // this array stores, how many time a particular index is selected	
	int check=1;	
	for (int j=0; j<10; j++)
	{
		number = (timeSinceEpochMillisec() % 101)/100.0; // number is not index here
		check=1;
		for (int i=1;i<best_of;i++)
		{	
		    if (number <= arr_2[i] and number > arr_2[i-1])
		    {
			check=0;
			selection_count[i] += 1;
			break;
		    }
		}

		if (check == 1)
		{
		    selection_count[0] += 1;
		}
	}

	// Getting the index that is selected more time
	var = 0;
	index = 0;
	for (int i=0; i<best_of; i++)
	{
	    if (selection_count[i] > var)
	    {
		var = selection_count[i];
		index = i;
	    }
	}
	selected = array[index];
	if (index == 0){select_zero += 1;}


//////////////////////////////////////////  MC End  ////////////////////////////////////////////////////////

	best[maxNode] = selected;
        if(flag){ //there is some node in first network that are not still aligned

        if (alignScore[maxNode][best[maxNode]] > 0.0) {
          alignment[maxNode]=best[maxNode]; //align two nodes;

            alignNodes1[maxNode]=true;
            alignNodes2[best[maxNode]]=true;
            //align degree one neighbors together
            for(int j=0; j<network1.deg[maxNode]; j++){
                for(int k=0; k<network2.deg[best[maxNode]]; k++)
                    if( !alignNodes1[network1.neighbor[maxNode][j]] && !alignNodes2[network2.neighbor[best[maxNode]][k]])
                    {
                        if(network1.deg[network1.neighbor[maxNode][j]]==1 && network2.deg[network2.neighbor[best[maxNode]][k]]==1)
                        {
                            alignment[network1.neighbor[maxNode][j]] = network2.neighbor[best[maxNode]][k];

                            alignNodes1[network1.neighbor[maxNode][j]] = true;
                            alignNodes2[network2.neighbor[best[maxNode]][k]] = true;
                        }
                    }
                  }


            //update the align scores
            for(int c1=0; c1 <network1.deg[maxNode]; c1++)
                for(int c2=0; c2<network2.deg[best[maxNode]]; c2++)
                    alignScore[ network1.neighbor[maxNode][c1]][network2.neighbor[best[maxNode]][c2]]=(alignScore[ network1.neighbor[maxNode][c1]][network2.neighbor[best[maxNode]][c2]]+(coeff/max));
        }
      }
        counter = counter + 1;
        if ( counter % 100 == 0)
	{
            cout << counter << endl;
	    //exit(0);
	}
    }//end flag
    //cout<<"\n Zero Selection : "<<select_zero<<endl;

    //memory leak
    delete [] alignNodes1;
    delete [] alignNodes2;
    delete [] nodeScore1;
    delete [] nodeScore2;
    delete [] structralScore;
    delete [] bioScore;
    delete [] best;

    for(int j=0; j<network1.size; j++)
    {
        delete [] alignScore[j];
    }
    delete [] alignScore;


    evaluate(); //calculate the measurment evaluations
}

//calculate the evaluation measurments EC (Edge Correctness), IC (Interaction Correctness), NC (Node Correctness), CCCV and CCCE (largest Common Connected subraph with recpect to Vertices and Edges)
void Alignment::evaluate(void)
{
	CCCV = getCCCV(); //calculate CCCV
	CCCE = getCCCE(); //calculate CCCE
	EC = getEC();     //calculate Edge Correctness
    S3 = getS3();      //calculate S3
}

//calculate CCCV
//return the number of vertices of largest common connected subgraph of the alignment
int Alignment::getCCCV(void)
{
    int *subGraph;
    int compNum = 1; //number of connected components
	int *q = new int[network1.size]; //nodes that are already processed
	comp = new int[network1.size]; //dtermines the connected component each node belongs to.
    for(int i=0; i<network1.size; i++)
	{
		comp[i] = network1.size;
		q[i] = i;
	}

	int last = 0;

	//for each node of the network
    for(int i=0; i<network1.size; i++)
	{
		if(comp[i]==network1.size)
		{
			q[0] = i;
			comp[i] = compNum;
			compNum++;
			last = 1;
            		//finds all connected nodes tho the node i that is not alredy in a connected component
			for(int k=0; k<last; k++)
				for(int j=0; j<network1.deg[q[k]]; j++)
                   			//the node is not already processed
					if( comp[q[k]] < comp[network1.neighbor[q[k]][j]])
					{
                        if (alignment[q[k]] != -1)
                            for( int l=0; l < network2.deg[alignment[q[k]]]; l++ )
                                if(network2.neighbor[alignment[q[k]]][l] == alignment[network1.neighbor[q[k]][j]])
                                {
                                    comp[network1.neighbor[q[k]][j]] = comp[q[k]];
                                    q[last] = network1.neighbor[q[k]][j];
                                    last++;
                                }
					}
		}
	}

	subGraph = new int[compNum-1]; //array of connected components
	for(int i=0; i<compNum-1; i++)
		subGraph[i] = 0;
	for(int i=0; i<network1.size; i++)
		subGraph[comp[i]-1]++; //number of nodes in a definit connected component

    //find the component with maximum nodes
    maxComp = 0;
	for(int i=0; i<compNum-1; i++)
	{
		if(subGraph[maxComp] < subGraph[i])
			maxComp = i;
	}

    int temp = subGraph[maxComp];

    //memory leak
    delete [] subGraph;
    delete [] q;

	return temp;
}

//calculate the evaluation measurment CCCE
//return the number of edges of largest common connected subgraph of the alignment
int Alignment::getCCCE(void)
{
	int edgeComp = 0;
    ofstream CC;
    //for each node of first network
	for(int i=0; i<network1.size; i++)
	{
        //for each neighbor of node i
		for(int j=0; j<network1.deg[i]; j++)
            //for each neighbor l of a node in second network that is aligned with node i
			if (alignment[i] != -1)
                for( int l=0; l < network2.deg[alignment[i]]; l++ )
                    if(network2.neighbor[ alignment[i] ][l] == alignment[network1.neighbor[i][j]])
                        if( comp[i]-1 == maxComp){
                            edgeComp++;
                        }
	}

	return ( edgeComp / 2 );
}

//calculate the evaluation measurment EC
//returns the percent of edges that are mapped correctly in alignment
float Alignment::getEC(void)
{
	int totalScore=0;

	//for each node i in first network
    for(int i=0; i<network1.size; i++)
	{
        //for each neighbor j of node i
        for(int j=0; j<network1.deg[i]; j++)
			//for each neighbor l of a node in second network that is aligned with node i
            if (alignment[i] != -1)
                for( int l=0; l < network2.deg[alignment[i]]; l++ ) {
                    if(network2.neighbor[ alignment[i] ][l] == alignment[ network1.neighbor[i][j] ]) {
                        totalScore++;
                    }
                }
	}

	//minimum number of edges of two networks
    int minEdge = ( network1.numOfEdge > network2.numOfEdge)? network2.numOfEdge : network1.numOfEdge;
    //calculate EC(edge correctness)
	return ( (float) totalScore ) / ( 2 * minEdge );
}

float Alignment::getS3(void)
{
	int totalScore=0;
    int* alignnodes = new int[network1.size];
    int num_edge_net2=0;

	//for each node i in first network
    for(int i=0; i<network1.size; i++)
	{
        alignnodes[i]=alignment[i];
        //for each neighbor j of node i
        for(int j=0; j<network1.deg[i]; j++)
			//for each neighbor l of a node in second network that is aligned with node i
            if (alignment[i] != -1)
                for( int l=0; l < network2.deg[alignment[i]]; l++ ) {
                    if(network2.neighbor[ alignment[i] ][l] == alignment[ network1.neighbor[i][j] ]) {
                        totalScore++;
                    }
                }
	}
    totalScore=totalScore/2;

    for(int i=0; i<network1.size; i++)
        if (alignment[i] != -1)
            for(int j=0; j<network2.deg[alignnodes[i]]; j++)
                for(int l=0; l<network1.size; l++)
                    if(network2.neighbor[alignnodes[i]][j]==alignnodes[l])
                        num_edge_net2++;
    num_edge_net2=num_edge_net2/2;
	//minimum number of edges of two networks
    int minEdge = ( network1.numOfEdge > network2.numOfEdge)? network2.numOfEdge : network1.numOfEdge;
    //calculate EC(edge correctness)
	return ( (float) totalScore ) / ( minEdge + float(num_edge_net2) - totalScore );
}

//print the evaluation measurments in a file with input parameter name
//Input parameter name determines the file that result are to be written in.
void Alignment::outputEvaluation(string name, double lambda, double alpha, double beta)
{
	string outFile = "alignments/mc_alignments/" + name; // path of alignment file
    //add a definit suffix to the file
	outFile.append(".eval");
	ofstream outputFile( outFile.c_str());

    //print in console

	outputFile << "++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++" << endl;
	outputFile << "*** CONNECTED COMPONENTS SIZE : " << endl;
	outputFile << "Nodes = " << CCCV << endl;
	outputFile << "Edges = " << CCCE << endl;
	outputFile << "===============================================================" << endl;
  outputFile << "Topological Score : " << alpha << endl;
  outputFile << "Biological Score : " << 1-alpha << endl;
  outputFile << "Structure Weight : " << 1-beta << endl;
  outputFile << "Sequence Weight : " << beta << endl;

  if(reverse)
	{
		outputFile << "G1:  Nodes : " << network2.size << "  - Edges : " << network2.numOfEdge << endl;
		outputFile << "G2:  Nodes : " << network1.size << "  - Edges : " << network1.numOfEdge << endl;
	}
	else
	{
		outputFile << "G1:  Nodes : " << network1.size << "  - Edges : " << network1.numOfEdge << endl;
		outputFile << "G2:  Nodes : " << network2.size << "  - Edges : " << network2.numOfEdge << endl;
	}

	outputFile << "EC : " << EC << endl;
	outputFile << "S3 : " << S3 << endl;
}

//print the alignment(mapping) in a file with input parameter name
//Input parameter name determines the file that mapping is to be written in.
void Alignment::outputAlignment(string name)
{
	string alignFile = "alignments/mc_alignments/" + name; // path of alignment file

	alignFile.append(".alignment");


	ofstream alignmentFile( alignFile.c_str());
	if(reverse)
        {
	  for(int i=0; i<network1.size; i++)
          {
            if (network1.getName( alignment[ i ] ).empty())
            {
             continue;
            }
            alignmentFile << network1.getName( alignment[ i ] ) << ' ' << network2.getName( i )<< endl;
          }
	}
        else 
        {
	  for(int i=0; i<network1.size; i++){
            if (network2.getName( alignment[ i ] ).empty()) {
              continue;
            }

            alignmentFile << network1.getName( i ) << ' ' << network2.getName( alignment[ i ] )<< endl;
          }
        }
        cout<<"\nDONE WITH ALIGNMENT. ALIGNMENT FILE IS SAVED IN ALIGNMENT FOLDER\n"<<endl;
}
//instructor
Alignment::Alignment(void)
{
}
//destructor
Alignment::~Alignment(void)
{
}
