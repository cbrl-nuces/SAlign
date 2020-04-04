
////////////////////////////// HOW TO RUN ///////////////////////////////////

##Execution:

1-a) Cammand 1-a (for SAlign):
	g++ SAlign.cpp Network.cpp Alignment.cpp -o SAlign

1-b) Cammand 1-b (for SAlign_mc):
	g++ SAlign.cpp Network.cpp Alignment-mc.cpp -o SAlign

2) Cammand 2:
	./SAlign net1.txt net2.txt -l 0.1 -a 0.1 -d 10 -b all_files.bitscore -s net1_net2.tmscore -t 0.7 -n 1

General format of cammand 2:
./SAlign Network_1 Network_2 -l 0.1 -a alpha_value -d 10 -b sequence_similarity_file -s strcuture_similarity_file -t beta_value -n alignment_number


--------------------------------------------------------------------------------------------------------------------------------------------
NOTE: change the alignment number for every run of Monte-Carlo based alignment (we used average of 10 runs for Monte-Carlo based alignments)
--------------------------------------------------------------------------------------------------------------------------------------------


## PARAMETERS:

all_files.bitscore --> pre-computed sequence similarity scores
net1_net2.tmscore --> pre-computed structure similarity scores

l 0.1 ==> node score weightage taken from hubalign
a 0.1 ==> alpha = 0.1 --> 10% topological score and 90% biological score
t 0.7 ==> beta = 0.7 --> 70% sequence and 30% structure
d = 10 ==> degree threshold is 10 used in topologiocal score calculation
n = 1 ==> alignment number is 1


////////////////////////////// FUNCTIONALITY AND FLOW //////////////////////////////////

Run these cammands to calculate alignment of two networks.

First cammand is used for the compilation of all files ==> SAlign will be the resultant compiled file of first cammand which is used in the second cammad (./SAlign ...) 
the second cammad is used to run the alignment file. Alignments are generated in the same folder,,, last digit of second cammand is the alignment number

After alignment, go to the semnatic folder and run the files gosemsem.r using cammand "r filename.r"
Then run average_semantic_similarity.py (Command: python3 filename.py or python filename.py) file to calculate average AFS and number of aligned nodes.


---------------------------------------------------------------------------------------------------
##Cautions:
See the file paths of gosemsim-greedy.r and average_semantic_similarity.py to avoid any path error.
---------------------------------------------------------------------------------------------------


### CONTACT-US

for any confusion/query please email us at "hammad.naveed@nu.edu.pk" or "umair.ayub@nu.edu.pk"
