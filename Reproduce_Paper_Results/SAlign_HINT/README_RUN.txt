
////////////////////////////// HOW TO RUN ///////////////////////////////////

##Execution:

1-a) Cammand 1 (for SAlign):
	g++ SAlign.cpp Network.cpp Alignment.cpp -o SAlign

1-b) Cammand 1 (for SAlign_mc):
	g++ SAlign.cpp Network.cpp Alignment-mc.cpp -o SAlign

2) Cammand 2:
	./SAlign MusMusculus_htb_hq.txt HomoSapiens_htb_hq.txt -l 0.1 -a 0.1 -d 10 -b all_files.bitscore -s MusMusculus_htb_hq_HomoSapiens_htb_hq.tmscore -t 0.7 -n 1

General format of cammand 2:
./SAlign Network_1 Network_2 -l 0.1 -a alpha_value -d 10 -b sequence_similarity_file -s strcuture_similarity_file -t beta_value -n alignment_number

NOTE: change the alignment number for every run of Monte-Carlo based alignment (we used average of 10 runs for Monte-Carlo based alignments)

##PARAMETERS:
all_files.bitscore --> pre-computed sequence similarity scores -> please see "blast" folder
MusMusculus_htb_hq_HomoSapiens_htb_hq.tmscore --> pre-computed structure similarity scores -> please see "tmalign" folder 
l 0.1 ==> node score weightage taken from hubalign
a 0.1 ==> alpha = 0.1 --> 10% topological score and 90% biological score
t 0.7 ==> beta = 0.7 --> 70% sequence and 30% structure
d = 10 ==> degree threshold is 10 used in topologiocal score calculation
n = 1 ==> alignment number is 1

##IN CASE OF PATH ERRORS:
Set the paths of the files. You need to update the follwing files:

1) Network.cpp:
      Path of the files "MusMusculus_htb_hq.txt" "HomoSapiens_htb_hq.txt" --> HINT_dataset/ (update line # 57 of the code)

2) Alignment.cpp:
      Path of the file all_files.bitscore --> blast (update line # 57 of the code)
      Path of the file XXXX.tmscore --> tmalign (update line # 91 of the code)


////////////////////////////// FUNCTIONALITY AND FLOW //////////////////////////////////

Run these cammands to calculate alignment of two pairs...
Result will be written in the alignments folder (see the alignment file that you are using to see where it writes the alignment files. last functions of every file is used to write the alignmnent files)

First cammand is used for compilation of all files ==> SAlign will be the resultant complied file of first cammand which is used in the second cammad (./SAlign ...) 
the second cammad is used to run the alignment file. Alignments are generated in alignments folder,,, last digit of second cammand is the alignment number

After alignment, go to the semnatic folder and run the files gosemsem-greedy.r using cammand "r filename.r"
Then run average_semantic_similarity.py (Command: python3 filename.py or python filename.py) file to calculate average AFS and number of aligned nodes.

##Cautions:
See the file paths of gosemsim-greedy.r and average_semantic_similarity.py to avoid any path error.


////////// CONTACT-US /////////

for any confusion/query please email us at "hammad.naveed@nu.edu.pk" or "umair.ayub@nu.edu.pk"
