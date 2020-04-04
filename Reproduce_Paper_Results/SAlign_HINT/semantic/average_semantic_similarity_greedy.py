for i in range(0,2):
    sum_one = 0.0
    avg_one = 0.0
    node_one = 0.0
    
    if i==0:
        f = open("semantic_files/greedy_alignments/MusMusculus_htb_hq_SaccharomycesCerevisiaeS288C_htb_hq_mf.txt")
    if i==1:
        f = open("semantic_files/greedy_alignments/MusMusculus_htb_hq_SaccharomycesCerevisiaeS288C_htb_hq_bp.txt")
    mf_lines = f.readlines()
    f.close()

    for j in range(len(mf_lines)):

        l = mf_lines[j].strip("\n")
        l = l.split("\t")
        if l[-1] != 'NA':
            sum_one+=float(l[-1])
            node_one+=1

    avg_one = sum_one/float(node_one)
    if i==0:
        print ("MF Results:",avg_one, node_one)
    if i==1:
        print ("BP Results:",avg_one, node_one)

