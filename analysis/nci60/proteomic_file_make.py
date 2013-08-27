#!/usr/bin/python

# We generate correlations between recon2 genes' 
# mRNA expression profiles and proteomic profiles, 
# output the consolidated data for analysis, and 
# output protein only and protein + scaled mRNA data
# as input to FALCON.

# INPUT:
# Recon 2 gene list (argv[1])
# CORE to Proteomic NCI-60 cell line name map (argv[2]; NCI60_labels.csv)
# Entrez ID to IPI ID map file (argv[3]; Entrez_and_IPI_unique.csv)
# Proteomic expression file (argv[4]; Expanded_IPI_Listing.csv)
# mRNA expression file (argv[5]; Gholami_Table_S8_mRNA.csv)

# Output:
# RUNDIR/nci60prot/<individual_exp_files.csv>
# RUNDIR/nci60prot_mRNA/<individual_exp_files.csv>

# 
# Example run in directory
# cd ~/FBA/models/Analysis/CancerExpression/NCI60/
# ~/FBA/FALCON/analysis/nci60/proteomic_file_make.py ~/FBA/models/rec2.genes NCI60_labels.csv Entrez_and_IPI_unique.csv Expanded_IPI_Listing.csv Gholami_Table_S8_mRNA.csv


import sys
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy
import scipy.stats

if len(sys.argv) != 6:
    sys.exit('ERROR: Usage %s model_genes_file' % sys.argv[0])


modgenesList = []
modgenes = {}
modgenesFI = open(sys.argv[1],'r')
lineno = -1
while True:
    line = modgenesFI.readline()
    line = line.strip()
    if len(line) == 0:
        break
    lineno = lineno + 1
    modgenesList = modgenesList + [line]
    modgenes[line] = lineno
modgenesFI.close()
numgene = len(modgenesList)


n60_COREtoPROT = {}
n60_PROTtoCORE = {}
nci60labelFI = open(sys.argv[2],'r')
header = nci60labelFI.readline()
while True:
    line = nci60labelFI.readline()
    line = line.strip()
    if len(line) == 0:
        break
    colvals = line.split("\t")
    colvals[0] = colvals[0].strip()
    colvals[1] = colvals[1].strip()
    n60_COREtoPROT[colvals[0]] = colvals[1]
    n60_PROTtoCORE[colvals[1]] = colvals[0]
nci60labelFI.close()

EntrezToIPI = {}
IPIToEntrez = {}
entrezIPIFI = open(sys.argv[3],'r')    
# no header
while True:
    line = entrezIPIFI.readline()
    line = line.strip()
    if len(line) == 0:
        break
    colvals = line.split("\t")
    colvals[0] = colvals[0].strip()
    colvals[1] = colvals[1].strip()
    EntrezToIPI[colvals[0]] = colvals[1]
    IPIToEntrez[colvals[1]] = colvals[0]
entrezIPIFI.close()


#In CPC file with best matching proteins, we have
#column 28 (index 27) should be the first cell line
#column 1 (index 0) has the IPI id of the best match
#
#In the extended ipi listing, 
#column 2 (index 1) should be the first cell line
#
cell_line_start = 1

#Read in Proteomic data
#Use a dict of dict of values:
#      [cell line][gene]->expression

PROTD = {}
PROTDcols = {}
IPIdataFI = open(sys.argv[4],'r')    
header = IPIdataFI.readline()
#first get column indices for different cell lines
col_names = header.split(",")
num_cols = len(col_names)
for i in range(cell_line_start,num_cols):
    col = col_names[i]
    col = col.strip()
    col = col.replace(" ","_")
    col = col.replace("CPC_","",1)
    if n60_PROTtoCORE.has_key(col):
        PROTDcols[i] = col
        PROTD[col] = {}
    else:
        raise Exception("cell line label unmatched in Protein data.")
#read in the data
while True:
    line = IPIdataFI.readline()
    line = line.strip()
    if len(line) == 0:
        break
    colvals = line.split(",")
    IPI = colvals[0].strip()
    for i in range(cell_line_start, num_cols):
        cv = colvals[i].strip()
        if cv == "NA" or cv == "":
            cv = "nan"
        cvf = float(cv)
        PROTD[PROTDcols[i]][IPI] = cvf    
IPIdataFI.close()

#Now handle the mRNA data similarly
MRNAD = {}
MRNADcols = {}
mrnaFI = open(sys.argv[5],'r')
#column 9 (index 8) should be the first cell line
#column 2 (index 1) has the IPI id, but sometimes it may not exist
#column 8 (index 7) whether the probe maps to proteomics
#Notes: take the average of expression for probes from different 
#arrays (after exponentiation)?
header = mrnaFI.readline()
header2 = mrnaFI.readline()
#first get column indices for different cell lines
col_names = header.split("\t")
for i in range(8,len(col_names)):
    col = col_names[i]
    col = col.strip()
    if n60_PROTtoCORE.has_key(col):
        MRNADcols[i] = col
        MRNAD[col] = {}
    else:
        raise Exception("cell line label unmatched in mRNA data.")
#read in the data
while True:
    line = mrnaFI.readline()
    line = line.strip()
    if len(line) == 0:
        break
    colvals = line.split("\t")
    IPI = colvals[1].strip()
    uniqSel = colvals[6].strip()
    map2prot = colvals[7].strip()
    if (
            IPIToEntrez.has_key(IPI) 
            and 
            (
                uniqSel == "TRUE" 
                or 
                map2prot == "TRUE"
                )
            ):
            for i in range(8, len(colvals)):
                    cv = colvals[i].strip()
                    if cv == "NA":
                            cv = "nan"
                    cvf = pow(2,float(cv))
                    if MRNAD[MRNADcols[i]].has_key(IPI):
                            MRNAD[MRNADcols[i]][IPI].append(cvf)
                    else:
                            MRNAD[MRNADcols[i]][IPI] = [cvf]    
mrnaFI.close()    
for cl in MRNAD.keys():
    for g in MRNAD[cl].keys():
        MRNAD[cl][g] = sum(MRNAD[cl][g])/len(MRNAD[cl][g])
        
#Construct correlation and mRNA+Prot matrix 
# y: prot such that: prot AND model AND mRNA
# x: mRNA such that: prot AND model AND mRNA
# y_all: prot such that: prot AND mRNA
# x_all: mRNA such that: prot AND mRNA
y = []
x = []
y_all = []
x_all = []    
    
MPmissed = 0    
for cl in PROTD.keys():
    if MRNAD.has_key(cl):
        for g in PROTD[cl].keys():
            if MRNAD[cl].has_key(g):
                if ((MRNAD[cl][g] == MRNAD[cl][g])
                    and (PROTD[cl][g] == PROTD[cl][g])):
                    x_all.append(MRNAD[cl][g])
                    y_all.append(PROTD[cl][g])
                    if IPIToEntrez.has_key(g):
                        ENTREZ = IPIToEntrez[g]
                        if modgenes.has_key(ENTREZ):
                            x.append(MRNAD[cl][g])
                            y.append(PROTD[cl][g])
                    else:
                        MPmissed += 1
    else:
        raise Exception("mRNA <-> protein cell line " 
                        + cl + " label mismatch!")
print("Found " + str(len(x)) + 
      " (mRNA, Prot) data points for model genes.")
print("m-P mismatches: " + str(MPmissed))
x = np.array(x)
y = np.array(y)
A = np.vstack([x, np.ones(len(x))]).T
m, b = np.linalg.lstsq(A, y)[0]
R = np.corrcoef(x,y)
print("y = " + str(m) + "x + " + str(b) + 
      "with Pearson's r = " + str(R[0][1]))
print("Spearman's rho: " + str(scipy.stats.stats.spearmanr(x, y)[0]))
print("Kendall's tau: " + str(scipy.stats.stats.kendalltau(x, y)[0]))
fig = plt.figure()
## left panel
ax1 = fig.add_subplot(121)
ax1.set_title('Metabolic Genes')
ax1.scatter(x,y,color='blue',s=5,edgecolor='none')
ax1.set_aspect(1./ax1.get_data_ratio()) # make axes square
ax1.set_xlabel('mRNA intensity')
ax1.set_ylabel('protein copy number')
# Now for all (including non-model) genes
x_all = np.array(x_all)
y_all = np.array(y_all)
A = np.vstack([x_all, np.ones(len(x_all))]).T
m, b = np.linalg.lstsq(A, y_all)[0]
R = np.corrcoef(x_all,y_all)
print("y_all = " + str(m) + "x_all + " + str(b) + 
      "with Pearson's r = " + str(R[0][1]))
#print("Spearman's rho: " + str(scipy.stats.stats.spearmanr(x_all, y_all)[0]))
#print("Kendall's tau: " + str(scipy.stats.stats.kendalltau(x_all, y_all)[0]))
## right panel
ax2 = fig.add_subplot(122)
ax2.set_title('All Genes')
ax2.scatter(x_all,y_all,color='blue',s=5,edgecolor='none')
ax2.set_aspect(1./ax2.get_data_ratio()) # make axes square
ax2.set_xlabel('mRNA intensity')
#ax2.set_ylabel('protein copy number')
fig.savefig('mRNA_protein_correlation.png')
# 3: prot + scaled(mRNA) such that: prot AND model, else mRNA AND model



#Need to redo this for the new data:
exit(0)


rnaseqFI = open('nci60-chiron-expression-data_genes.tsv','r')
colnumheader = rnaseqFI.readline().strip()
header = rnaseqFI.readline().strip()
col_names = header.split("\t")
rnaseqFI.close()

for col in range(5,63):
    tissue = col_names[col].strip()
    tissue = tissue.replace("(","_")
    tissue = tissue.replace(")","_")
    tissue = tissue.replace(" ","_")    
    tissue = tissue.replace("/","_")
    tissue = tissue.replace("-","_")
    
    print(tissue)
    print(col)
    # Stupidly open file repeatedly
    rnaseqFI = open('nci60-chiron-expression-data_genes.tsv','r')
    outFI = open('NCI60exp/' + tissue + '.csv','w')
    outFI.write("gene\tmean\tvar\n")
    colnumheader = rnaseqFI.readline().strip()
    header = rnaseqFI.readline().strip()
    while True:
        line = rnaseqFI.readline()
        if len(line) < 2:
            break
        ll = line.split('\t')
        exp = ll[col].strip()
        entrez = ll[2].strip()
        if modgenes.has_key(entrez):
            outlist = [entrez, exp, "1"] 
            outFI.write("\t".join(outlist)+"\n")
    outFI.close()
    rnaseqFI.close()
