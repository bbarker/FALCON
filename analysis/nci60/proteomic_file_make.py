#!/usr/bin/python

# Copyright Brandon Barker
# August 2013

# We generate correlations between recon2 genes' 
# mRNA expression profiles and proteomic profiles, 
# output the consolidated data for analysis, and 
# output protein only and protein + mRNA data 
# as input to FALCON.

# INPUT: (Please edit if necessary)
# Recon 2 gene list 
rec2genes = '/home/brandon/FBA/models/rec2.genes'
# CORE to Proteomic NCI-60 cell line name map (NCI60_labels.csv)
nci60labels = '/home/brandon/FBA/models/Analysis/CancerExpression/NCI60/NCI60_labels.csv'
# Entrez ID to IPI ID map file (Entrez_and_IPI_unique.csv)
geneIDdb = '/home/brandon/FBA/models/Analysis/CancerExpression/NCI60/Entrez_and_IPI_unique.csv'
# Proteomic expression file (protLFQ.csv or Expanded_IPI_Listing.csv)
protEXP = '/home/brandon/FBA/models/Analysis/CancerExpression/NCI60/protLFQ.csv'
# protEXP = '/home/brandon/FBA/models/Analysis/CancerExpression/NCI60/prot_iBAQ.csv'
# Deep Proteomic expression file ( )
deepProtEXP = '/home/brandon/FBA/models/Analysis/CancerExpression/NCI60/DeepProtLFQ.csv'
# mRNA expression file (Gholami_Table_S8_mRNA.csv)
mrnaEXP = '/home/brandon/FBA/models/Analysis/CancerExpression/NCI60/Gholami_Table_S8_mRNA.csv'


# OUTPUT:
#
# scatter between proteomic and deep proteomic data:
# RUNDIR/prot_DeepProt_correlation.*
#
# scatter between merged proteomic and mRNA data:
# RUNDIR/mRNA_protein_correlation.*
#
# 3d surface of pearson correlations used for noise analysis
# RUNDIR/Intensity_corr_mesh_*
#
# 2d plot of pearson correlation for noise analysis
# RUNDIR/Intensity_corr_line.*
#
# Histogram of protein and microarray absolute expression
# RUNDIR/expression_dists.*
#
# Histogram of protein and microarray intensities subsequently
# mapped to zero.
# RUNDIR/expression_zero_dists.*
#
# mRNA files for input to FALCON:
# RUNDIR/nci60mRNA/<individual_exp_files.csv>
#
# merged proteomic files for input to FALCON:
# RUNDIR/nci60prot/<individual_exp_files.csv>
#
# merged proteomic and mRNA files for input to FALCON:
# RUNDIR/nci60prot_mRNA/<individual_exp_files.csv>
#
# table of all 3 data types from all cell lines in the model:
modDatSave = 'model_expression.csv'
#
#
# mRNA files for input to FALCON:
# RUNDIR/nci60mRNA_thresh/threshval/<individual_exp_files.csv>
#
# merged proteomic files for input to FALCON:
# RUNDIR/nci60prot_thresh/threshval/<individual_exp_files.csv>
#
# merged proteomic and mRNA files for input to FALCON:
# RUNDIR/nci60prot_mRNA_thresh/threshval/<individual_exp_files.csv>


import multiprocessing

# 
# Example run in directory
# cd ~/FBA/models/Analysis/CancerExpression/NCI60/
# ~/FBA/FALCON/analysis/nci60/proteomic_file_make.py

import os
import sys
import copy
import numpy as np
#import pylab as P
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import scipy
import scipy.stats
from sklearn import linear_model
import pickle
#import Bio.Cluster
#import rpy2

NaN = float("nan")

def product(*args):
    if not args:
        return iter(((),)) # yield tuple()
    return (items + (item,)
            for items in product(*args[:-1]) for item in args[-1])

def veccorr(x, y):
    R = np.corrcoef(x,y)
    return R[0,1]

def plotScatterCorr(ax, x, y, fig_title, x_title, y_title, txtpos, intercept=True):
    (m,b) = (NaN,NaN)
    if intercept:
        clf = linear_model.LinearRegression()
        A = np.vstack([x, np.ones(len(x))]).T
        # m, b = np.linalg.lstsq(A, y)[0]
        clf.fit(A,y)
        m = clf.coef_[0]
        b = clf.intercept_
    else:
        clf = linear_model.LinearRegression(fit_intercept=False)
        A = np.vstack([x, np.zeros(len(x))]).T
        # m, b = np.linalg.lstsq(A, y)[0]
        clf.fit(A,y)
        m = clf.coef_[0]
        b = 0                    
    R = np.corrcoef(x,y)
    xmin = min(x)
    xmax = max(x)
    ymin = min(y)
    ymax = max(y)
    eps_x = (xmax-xmin)/100
    eps_y = (ymax-ymin)/100
    ax.set_title(fig_title, fontsize=10)    
    ax.axis([xmin-eps_x, xmax+eps_x, ymin-eps_y, ymax+eps_y])
    ax.scatter(x,y,color='blue',s=0.5,edgecolor='none')
    ax.set_aspect(1./ax.get_data_ratio()) # make axes square
    if len(x_title) > 0:
        ax.set_xlabel(x_title)
    if len(y_title) > 0:    
        ax.set_ylabel(y_title)
    y_lin = m*x+b
    ax.plot(x,y_lin,color='red')
    ax.text(txtpos[0], txtpos[1], "y = " + str(m)[0:4] + 
            "x + " + str(b)[0:4], fontsize=8)
    ax.text(txtpos[2], txtpos[3], "r = " + str(R[0][1])[0:4], fontsize=8)
    return (m, b, R[0][1])
    
def getProt_mRNA_pairs(MRNAD, PROTD, IPIToEntrez):
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
                        and (PROTD[cl][g] == PROTD[cl][g])
                        and PROTD[cl][g] > 0.0):
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
    # print("Found " + str(len(x)) + 
    #      " (mRNA, Prot) data points for model genes.")
    # print("m-P mismatches: " + str(MPmissed))
    return(x, y, x_all, y_all)

# Apparently a closure won't work with multiprocessing.Pool
# def makesublistcorr(x, y):        
#     def sublistcorr(tup):
#         xc = x[tup[0]:end+1-tup[1]]     
#         yc = y[tup[0]:end+1-tup[1]]
#         return veccorr(xc, yc)
#     return sublistcorr

# Use a class instead
class sublistcorr:
    def __init__(self, x, y, end):
        self.x = x
        self.y = y
        self.end = end
    def __call__(self, ij):
        end = self.end
        xc = self.x[ij[0]:end+1-ij[1]]     
        yc = self.y[ij[0]:end+1-ij[1]]
        return veccorr(xc, yc)
    
def lineCorrPlot(x, y, m, ival, fig, axpos):
    # Need to normalize x and y
    xmean = np.mean(x)
    ymean = np.mean(y)
    xymid = np.mean([xmean, ymean])
    x = xymid/xmean * x
    y = xymid/ymean * y
    irange = range(0, m+1, ival)
    ilen = len(irange)
    z = np.zeros([ilen, 1])
    xfull = np.zeros([ilen,1])
    # Need to sort (x, y)
    # min may be biased against larger values.
    # It will be interesting to check the analysis
    # with independently sorted indices, which
    # also removes any scaling issues.
    (x,y) = zip(*sorted(zip(x,y), key=np.min)) 
    end = len(x)-1
    print([x[0], y[0], x[end], y[end]]) 
    i_by_j = product(irange, [0])
    len_i_by_j = ilen 
    pool = multiprocessing.Pool(30)
    Rvals = pool.map(sublistcorr(x, y, end), i_by_j)
    pool.close()
    i_by_j = product(irange, [0])
    for ij_idx in range(0,len_i_by_j):
        ij = i_by_j.next()
        z[int(ij[0]/ival)] = Rvals[int(ij[0]/ival)]
        xfull[int(ij[0]/ival)] = ij[0]
    ax = fig.add_subplot(axpos)
    ax.plot(xfull, z, lw=2)
    ax.set_xlabel('removed from bottom')
    ax.set_ylabel("Pearson's r")    
    return (ax, Rvals)    

def centerMeshCorrPlot(x, y, m, n, ival):
    # Need to normalize x and y
    xmean = np.mean(x)
    ymean = np.mean(y)
    xymid = np.mean([xmean, ymean])
    x = xymid/xmean * x
    y = xymid/ymean * y
    irange = range(0, m+1, ival)
    jrange = range(0, n+1, ival)
    ilen = len(irange)
    jlen = len(jrange)
    z = np.zeros([ilen, jlen])
    xfull = np.zeros([ilen, jlen])
    yfull = np.zeros([ilen, jlen])
    # Need to sort (x, y)
    (x,y) = zip(*sorted(zip(x,y), key=np.min)) 
    end = len(x)-1
    print([x[0], y[0], x[end], y[end]]) 
    i_by_j = product(irange, jrange)
    len_i_by_j = ilen * jlen
    pool = multiprocessing.Pool(30)
    Rvals = pool.map(sublistcorr(x, y, end), i_by_j)
    print(type(Rvals))
    pool.close()
    i_by_j = product(irange, jrange)
    for ij_idx in range(0,len_i_by_j):
        ij = i_by_j.next()
        z[int(ij[0]/ival)][int(ij[1]/ival)] = Rvals[ij_idx]
        xfull[int(ij[0]/ival)][int(ij[1]/ival)] = ij[0]
        yfull[int(ij[0]/ival)][int(ij[1]/ival)] = ij[1]        
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.plot_surface(xfull, yfull, z, cmap=cm.coolwarm, linewidth=0.2)
    ax.set_xlabel('removed from bottom')
    ax.set_ylabel('removed from top')
    ax.set_zlabel("Pearson's r")    
    return (fig, ax)    

modgenesList = []
modgenes = {}
modgenesFI = open(rec2genes,'r')
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
nci60labelFI = open(nci60labels,'r')
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
entrezIPIFI = open(geneIDdb,'r')    
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
IPIdataFI = open(protEXP,'r')    
header = IPIdataFI.readline()
#first get column indices for different cell lines
col_names = header.split(",")
num_cols = len(col_names)
for i in range(cell_line_start,num_cols):
    col = col_names[i]
    col = col.strip()
    col = col.replace(" ","_")
    #col = col.replace("CPC_","",1)
    #col = col.replace("iBAQ_","",1)
    col = col.replace("LFQ_","",1)    
    if n60_PROTtoCORE.has_key(col):
        PROTDcols[i] = col
        PROTD[col] = {}
    else:
        raise Exception("cell line label " + col + " unmatched in Protein data.")
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


DPROTD = {}
DPROTDcols = {}
IPIdataFI = open(deepProtEXP,'r')    
header = IPIdataFI.readline()
#first get column indices for different cell lines
col_names = header.split(",")
num_cols = len(col_names)
for i in range(cell_line_start,num_cols):
    col = col_names[i]
    col = col.strip()
    col = col.replace(" ","_")
    #col = col.replace("CPC_","",1)
    #col = col.replace("iBAQ_","",1)
    col = col.replace("LFQ_","",1)    
    if n60_PROTtoCORE.has_key(col):
        DPROTDcols[i] = col
        DPROTD[col] = {}
    else:
        raise Exception("cell line label " + col + " unmatched in Protein data.")
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
        DPROTD[DPROTDcols[i]][IPI] = cvf    
IPIdataFI.close()

#Now handle the mRNA data similarly
MRNAD = {}
MRNADcols = {}
mrnaFI = open(mrnaEXP,'r')
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
            # These options appear to be equivalent here:
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



# Construct correlation between proteomic and deep Proteomic data.
y_notDeep = []
x_Deep = []
y_ModnotDeep = []
x_ModDeep = []
y_notDeep0 = []
x_Deep0 = []

    
for cl in DPROTD.keys():
    if PROTD.has_key(cl):
        for g in DPROTD[cl].keys():
            if PROTD[cl].has_key(g):
                if ((DPROTD[cl][g] == DPROTD[cl][g])
                    and (PROTD[cl][g] == PROTD[cl][g])
                    and PROTD[cl][g] > 0.0
                    and DPROTD[cl][g] > 0.0
                    ):
                    x_Deep.append(DPROTD[cl][g])
                    y_notDeep.append(PROTD[cl][g])
                    if IPIToEntrez.has_key(g):
                        ENTREZ = IPIToEntrez[g]
                        if modgenes.has_key(ENTREZ):
                            x_ModDeep.append(DPROTD[cl][g])
                            y_ModnotDeep.append(PROTD[cl][g])

                if ((DPROTD[cl][g] == DPROTD[cl][g])
                    and (PROTD[cl][g] == PROTD[cl][g])
                    ):
                    y_notDeep0.append(PROTD[cl][g])
                    x_Deep0.append(DPROTD[cl][g])

                    
print("Found " + str(len(y_notDeep)) + 
      " (Prot, Deep-Prot) data points.")

fig = plt.figure()

y_notDeep0 = np.array(y_notDeep0)
x_Deep0 = np.array(x_Deep0)
ax1 = fig.add_subplot(131)
(m, b, r) = plotScatterCorr(ax1, x_Deep0, y_notDeep0, 'With Zeros', 'deep proteomic intensity',
                            'proteomic intensity', [4.5, 1.7, 4.5, 1.2])
print("y_notDeep0 = " + str(m) + "x_Deep0 + " + str(b) + 
      " with Pearson's r = " + str(r))
print("Spearman's rho: " + str(scipy.stats.stats.spearmanr(x_Deep0, y_notDeep0)[0]))
print("Kendall's tau: " + str(scipy.stats.stats.kendalltau(x_Deep0, y_notDeep0)[0]))
print("")

y_notDeep = np.array(y_notDeep)
x_Deep = np.array(x_Deep)
ax2 = fig.add_subplot(132)
(m, b, r) = plotScatterCorr(ax2, x_Deep, y_notDeep, 'No Zeros', 'deep proteomic intensity',
                            '', [4, 8.5, 4, 8.2])
print("y_notDeep = " + str(m) + "x_Deep + " + str(b) + 
      " with Pearson's r = " + str(r))
print("Spearman's rho: " + str(scipy.stats.stats.spearmanr(x_Deep, y_notDeep)[0]))
print("Kendall's tau: " + str(scipy.stats.stats.kendalltau(x_Deep, y_notDeep)[0]))
print("")

y_ModnotDeep = np.array(y_ModnotDeep)
x_ModDeep = np.array(x_ModDeep)
ax3 = fig.add_subplot(133)
(m, b, r) = plotScatterCorr(ax3, x_ModDeep, y_ModnotDeep, 'Model Genes; No Zeros', 'deep proteomic intensity',
                            '', [4.3, 8.4, 4.3, 8.1])
print("y_ModnotDeep = " + str(m) + "x_ModDeep + " + str(b) + 
      " with Pearson's r = " + str(r))
print("Spearman's rho: " + str(scipy.stats.stats.spearmanr(x_ModDeep, y_ModnotDeep)[0]))
print("Kendall's tau: " + str(scipy.stats.stats.kendalltau(x_ModDeep, y_ModnotDeep)[0]))
print("")
(mp, bp) = (m, b)

fig.tight_layout()
fig.savefig('prot_DeepProt_correlation.png', bbox_inches='tight',
            dpi=300)

print("Total number of Deep vs notDeep pairs, excluding zeros: " +
      str(len(x_Deep)))


# This is compute intensive, so it should be commented out usually
if False:
    (fig, ax) = centerMeshCorrPlot(x_Deep, y_notDeep, 20000, 300, 1)
    ax.view_init(azim=72)
    fig.tight_layout()
    fig.savefig('Intensity_corr_mesh_72', bbox_inches='tight', dpi=300)
    ax.view_init(azim=162)
    fig.tight_layout()
    fig.savefig('Intensity_corr_mesh_162', bbox_inches='tight', dpi=300)

    if not os.path.isdir('test_corr_fig'):                                   
        os.mkdir('test_corr_fig')                                            
    for i in range(0,360, 2):                                                
        ax.view_init(azim=i)                                                 
        fig.tight_layout()                                                   
        fig.savefig('test_corr_fig/Intensity_corr_mesh_' + str(i) + '.png',  
                    bbox_inches='tight', dpi=100)                            
    
fig = plt.figure()
xfmt = mpl.ticker.ScalarFormatter()
xfmt.set_powerlimits((-3, 2))
xfmt.set_scientific(True)
(ax1, Rbott) = lineCorrPlot(x_Deep, y_notDeep, len(x_Deep)-3, 1, fig, 121)
ax1.xaxis.set_major_formatter(xfmt)
(ax2, Rbottzoom) = lineCorrPlot(x_Deep, y_notDeep, 7000, 1, fig, 122)
ax2.xaxis.set_major_formatter(xfmt)
ax1.set_title("all pairs")
ax2.set_title("zoomed to bottom 7000 pairs")

# Now based on this analysis, we can define a new
# zero point for the protein intensities. 
(x_Deep,y_notDeep) = zip(*sorted(zip(x_Deep,y_notDeep), key=np.min))
PIntensZeroVal = 10000
for i in range(1500, 4000):
    # scaling shouldn't matter in this low range;
    # but later, we scale to the non-deep values
    if  y_notDeep[i] < x_Deep[i]:
        PIntensZeroVal = y_notDeep[i]
        print("Found intensity zero value " + str(PIntensZeroVal) + 
              " at index " + str(i))
        break
ax2.annotate('LFQ threshold for zeros: ' + str(PIntensZeroVal)[0:4], 
             xy=(i, Rbott[i]), xytext=(i-400, Rbott[i]+0.02), 
             arrowprops=dict(facecolor='black', shrink=0.05),)
fig.tight_layout()
fig.savefig('Intensity_corr_line.png', bbox_inches='tight', dpi=300)

    

#output = open('Intensity_corr_mesh.pickle', 'wb')
#pickle.dump(fig,output)
#output.close()


# Make copies of main expression dictionaries
# for later use in zero-threshold testing.
OMRNAD = copy.deepcopy(MRNAD)
OPROTBoth = copy.deepcopy(PROTD)
for cl in DPROTD.keys():
    for g in DPROTD[cl].keys():
        if DPROTD[cl][g] > 0:
            OPROTBoth[cl][g] = mp*DPROTD[cl][g] + bp
# At this point, the proteins have been normalized,
# but not the microarray data.



        
# Combine proteomic data and set low intensity
# values to zero.
# Also, convert existing 0s in proteomic data to NaNs
# based on the observation observed above that
# these 0s are most certainly "missing data"
#
# We also want
# to get a ballpark zero-value for microarray intensities
# based on this, but since there may be considerable variation
# the best we can think of is to assume the median value
# for those microarray intensities matching protein intensities
# below the threshold. 
mrnaZeroCandidates = []
protZeroCandidates = []    
for cl in PROTD.keys():
    for g in PROTD[cl].keys():
        if PROTD[cl][g] == 0:
            PROTD[cl][g] = NaN
        elif PROTD[cl][g] < PIntensZeroVal:
            protZeroCandidates.append(PROTD[cl][g])
            PROTD[cl][g] = 0
            if MRNAD[cl].has_key(g):
                if MRNAD[cl][g] == MRNAD[cl][g]:
                    mrnaZeroCandidates.append(MRNAD[cl][g])
        else:
            PROTD[cl][g] = PROTD[cl][g] - PIntensZeroVal    
for cl in DPROTD.keys():
    for g in DPROTD[cl].keys():
        if DPROTD[cl][g] == 0:
            DPROTD[cl][g] = NaN        
        DPROTD[cl][g] = mp*DPROTD[cl][g] + bp
        if DPROTD[cl][g] < PIntensZeroVal:
            protZeroCandidates.append(DPROTD[cl][g])
            DPROTD[cl][g] = 0
            if MRNAD[cl].has_key(g):
                if MRNAD[cl][g] == MRNAD[cl][g]:
                    if PROTD[cl].has_key(g):
                        if PROTD[cl][g] != PROTD[cl][g]:
                            mrnaZeroCandidates.append(MRNAD[cl][g])
                    else:
                        mrnaZeroCandidates.append(MRNAD[cl][g])                    
        else:
            DPROTD[cl][g] = DPROTD[cl][g] - PIntensZeroVal
mrnaZeroCandidates = np.array(mrnaZeroCandidates)
protZeroCandidates = np.array(protZeroCandidates)
MIntensMed = np.median(mrnaZeroCandidates)
MIntensMin = np.min(mrnaZeroCandidates)
MIntensMax = np.max(mrnaZeroCandidates)
MIntens1q = np.percentile(mrnaZeroCandidates,25)
MIntens3q = np.percentile(mrnaZeroCandidates,75)
# Since the r^2 value is approximately 0.25, 
# probably the most we can assume is that the first
# quartile is zero or very small, but this is not
# a safe analysis and more work is needed.
MIntensZeroVal = MIntens1q
print("Min, Quartiles, Max from mRNA zero-cutoff values: " + 
      str([MIntensMin, MIntens1q, MIntensMed, MIntens3q, 
           MIntensMax]))
for cl in MRNAD.keys():
    for g in MRNAD[cl].keys():
        if MRNAD[cl][g] < MIntensZeroVal:
            MRNAD[cl][g] = 0
        else:
            MRNAD[cl][g] = MRNAD[cl][g] - MIntensZeroVal    
            
            
PROTBoth = copy.deepcopy(PROTD)
for cl in DPROTD.keys():
    for g in DPROTD[cl].keys():
        if DPROTD[cl][g] > 0:
            PROTBoth[cl][g] = DPROTD[cl][g]

# Construct correlation and mRNA+Prot matrix 
# y: prot such that: prot AND model AND mRNA
# x: mRNA such that: prot AND model AND mRNA
# y_all: prot such that: prot AND mRNA
# x_all: mRNA such that: prot AND mRNA
# b versions include deep proteomic data as well
# d versions are deep only

(x, y, x_all, y_all) = getProt_mRNA_pairs(MRNAD, PROTD, IPIToEntrez)
(x_b, y_b, x_ball, y_ball) = getProt_mRNA_pairs(MRNAD, PROTBoth, IPIToEntrez)
(x_d, y_d, x_dall, y_dall) = getProt_mRNA_pairs(MRNAD, DPROTD, IPIToEntrez)

fig = plt.figure()

# Analayze metabolic gene correlation
x = np.array(x)
y = np.array(y)
ax1 = fig.add_subplot(221)
(m, b, r) = plotScatterCorr(ax1, x, y, 'Metabolic Genes', 'mRNA intensity',
                            'protein intensity', [0.7, 3.7, 0.7, 3.45], intercept=False)
print("y = " + str(m) + "x + " + str(b) + 
      " with Pearson's r = " + str(r))
print("Spearman's rho: " + str(scipy.stats.stats.spearmanr(x, y)[0]))
print("Kendall's tau: " + str(scipy.stats.stats.kendalltau(x, y)[0]))
#print("Kendall's tau: " + str(rpy.r.cor(x, y, method="kendall")))
#print("Kendall's tau: " + str(1 - Bio.Cluster.distancematrix((x,y), dist="k")[1][0]))
print("")

# Analyze all (including non-model) gene correlation
x_all = np.array(x_all)
y_all = np.array(y_all)
ax2 = fig.add_subplot(222)
(m, b, r) = plotScatterCorr(ax2, x_all, y_all, 'All Genes', 'mRNA intensity',
                            'protein intensity', [0.7, 3.8, 0.7, 3.5], intercept=False)
print("y_all = " + str(m) + "x_all + " + str(b) + 
      " with Pearson's r = " + str(r))
print("Spearman's rho: " + str(scipy.stats.stats.spearmanr(x_all, y_all)[0]))
print("")
#print("Kendall's tau: " + str(1 - Bio.Cluster.distancematrix((x_all,y_all), dist="k")[1][0]))


# Analayze metabolic gene correlation
x_d = np.array(x_d)
y_d = np.array(y_d)
ax3 = fig.add_subplot(223)
(m, b, r) = plotScatterCorr(ax3, x_d, y_d, 'Metabolic Genes', 'mRNA intensity',
                            'deep protein intensity', [0.7, 2.75, 0.7, 2.55], intercept=False)
print("y_d = " + str(m) + "x_d + " + str(b) + 
      " with Pearson's r = " + str(r))
print("Spearman's rho: " + str(scipy.stats.stats.spearmanr(x_d, y_d)[0]))
print("Kendall's tau: " + str(scipy.stats.stats.kendalltau(x_d, y_d)[0]))
print("")

# Analyze all (including non-model) gene correlation
x_b = np.array(x_b)
y_b = np.array(y_b)
ax4 = fig.add_subplot(224)
(m, b, r) = plotScatterCorr(ax4, x_b, y_b, 'Metabolic Genes', 'mRNA intensity',
                            'protein intensity (both)', [0.7, 3.7, 0.7, 3.4], intercept=False)
print("y_b = " + str(m) + "x_b + " + str(b) + 
      " with Pearson's r = " + str(r))
print("Spearman's rho: " + str(scipy.stats.stats.spearmanr(x_b, y_b)[0]))
print("")
(m_MtoP, b_MtoP) = (m, b)

fig.tight_layout()
fig.savefig('mRNA_protein_correlation.png', bbox_inches='tight',
            dpi=300)


# Combine proteomic and mRNA data
PROTMRNA = copy.deepcopy(PROTBoth)
for cl in MRNAD.keys():
    for g in MRNAD[cl].keys():        
        if PROTMRNA[cl].has_key(g):
            if PROTMRNA[cl][g] != PROTMRNA[cl][g]:
                if MRNAD[cl][g] == MRNAD[cl][g]:
                    PROTMRNA[cl][g] = m_MtoP*MRNAD[cl][g] + \
                    b_MtoP

        else:
            PROTMRNA[cl][g] = m_MtoP*MRNAD[cl][g] + \
            b_MtoP

# Plot histograms for different datatypes 
def flattenExpression(DATA):
    x = []
    x_mod = []
    for cl in DATA.keys():
        for g in DATA[cl].keys():
            if DATA[cl][g] == DATA[cl][g]:
                x.append(DATA[cl][g])
                if IPIToEntrez.has_key(g):
                    ENTREZ = IPIToEntrez[g]
                    if modgenes.has_key(ENTREZ):
                        x_mod.append(DATA[cl][g])
    return (np.array(x), np.array(x_mod))                        
                        
(xM, xM_mod) = flattenExpression(MRNAD)
xM = m_MtoP*xM + b_MtoP
xM_mod = m_MtoP*xM_mod + b_MtoP
(xP, xP_mod) = flattenExpression(PROTBoth)
(xPM, xPM_mod) = flattenExpression(PROTMRNA)
fig = plt.figure()
minx = np.min(np.concatenate((xM,xP,xPM)))
maxx = np.max(np.concatenate((xM,xP,xPM)))
minxmod = np.min(np.concatenate((xM_mod,xP_mod,xPM_mod)))
maxxmod = np.min(np.concatenate((xM_mod,xP_mod,xPM_mod)))
bins = np.linspace(minx, maxx, 50)
binsmod = np.linspace(minxmod, maxxmod, 50)
ax1 = fig.add_subplot(121)
ax1.hist(xM, bins, alpha=0.5, normed=1, label='mRNA')
ax1.hist(xP, bins, alpha=0.5, normed=1, label='Protein')
ax1.set_xlabel('shifted intensity')
ax1.set_ylabel('normalized density')
ax1.set_title('all genes')
ax1.legend()
ax2 = fig.add_subplot(122)
ax2.hist(xM_mod, bins, alpha=0.5, normed=1, label='mRNA')
ax2.hist(xP_mod, bins, alpha=0.5, normed=1, label='Protein')
ax2.set_xlabel('shifted intensity')
ax2.set_title('model genes')
ax2.legend()
fig.tight_layout()
fig.savefig('expression_dists.png', bbox_inches='tight', dpi=300)

fig = plt.figure()
mrnaZeroCandidates = m_MtoP*mrnaZeroCandidates + b_MtoP
minx = np.min(np.concatenate((mrnaZeroCandidates,protZeroCandidates)))
maxx = np.max(np.concatenate((mrnaZeroCandidates,protZeroCandidates)))
bins = np.linspace(minx, maxx, 50)
ax1 = fig.add_subplot(111)
ax1.hist(mrnaZeroCandidates, bins, alpha=0.5, normed=1, label='mRNA')
ax1.hist(protZeroCandidates, bins, alpha=0.5, normed=1, label='Protein')
ax1.set_xlabel('shifted intensity')
ax1.set_title("All zero-shifted proteins' intensities and corresponding" + 
              " mRNA intensities")
ax1.legend()
fig.tight_layout()
fig.savefig('expression_zero_dists.png', bbox_inches='tight', dpi=300)


# Begin the process of summing isoform data as input to FALCON.
modMcount = 0
modPcount = 0
modMPcount = 0
modelDataFI = open(modDatSave,'w')
header = "Entrez Gene\tIPI\tCell Line (Proteomics Label)\t" + \
  "Cell Line (CoRe Label)\tmRNA Data\tProtein Data\tProtein + mRNA\n"
modelDataFI.write(header)
if not os.path.isdir('nci60mRNA'):
    if os.path.exists('nci60mRNA'):
        raise Exception("Specified output directory is a file!")
    os.mkdir('nci60mRNA')
if not os.path.isdir('nci60prot'):
    if os.path.exists('nci60prot'):
        raise Exception("Specified output directory is a file!")
    os.mkdir('nci60prot')
if not os.path.isdir('nci60prot_mRNA'):
    if os.path.exists('nci60prot_mRNA'):
        raise Exception("Specified output directory is a file!")
    os.mkdir('nci60prot_mRNA')

DBGgenes = open('DBGgenes.txt','w')
IPIEntvals = {}
testCL = PROTMRNA.keys()[0]
nCL = len(PROTMRNA.keys())
ModMRNAD = {}
ModPROTBoth = {}
ModPROTMRNA = {}
print("Number of cell lines: " + str(nCL))
for clPROT in n60_PROTtoCORE.keys():
    ModMRNAD[clPROT] = {}
    ModPROTBoth[clPROT] = {}
    ModPROTMRNA[clPROT] = {}
    clCORE = n60_PROTtoCORE[clPROT]
    tissue = clCORE
    tissue = tissue.replace("(","_")
    tissue = tissue.replace(")","_")
    tissue = tissue.replace(" ","_")    
    tissue = tissue.replace("/","_")
    tissue = tissue.replace("-","_")    
    for g in PROTMRNA[clPROT]:
        ENTREZ = "UNKNOWN ENTREZ ID"
        if IPIToEntrez.has_key(g):
            ENTREZ = IPIToEntrez[g]  
            if modgenes.has_key(ENTREZ): 
                mpEXP = PROTMRNA[clPROT][g]
                pEXP = NaN
                mEXP = NaN
                if mpEXP == mpEXP:
                    if ModPROTMRNA[clPROT].has_key(ENTREZ):
                        ModPROTMRNA[clPROT][ENTREZ].append(mpEXP)
                    else:
                        ModPROTMRNA[clPROT][ENTREZ] = [mpEXP]                    
                    DBGgenes.write(ENTREZ + "\n")
                    if clPROT == testCL:
                        if IPIEntvals.has_key(ENTREZ):
                            IPIEntvals[ENTREZ].append([g,mpEXP])
                        else:
                            IPIEntvals[ENTREZ] = [[g,mpEXP]]    
                if PROTBoth[clPROT].has_key(g):
                    pEXP = PROTBoth[clPROT][g]
                    if PROTBoth[clPROT][g] == PROTBoth[clPROT][g]:
                        if ModPROTBoth[clPROT].has_key(ENTREZ):
                            ModPROTBoth[clPROT][ENTREZ].append(pEXP)
                        else:
                            ModPROTBoth[clPROT][ENTREZ] = [pEXP]
                if MRNAD[clPROT].has_key(g):    
                    mEXP = m_MtoP*MRNAD[clPROT][g] + b_MtoP
                    if MRNAD[clPROT][g] == MRNAD[clPROT][g]:
                        if ModMRNAD[clPROT].has_key(ENTREZ):
                            ModMRNAD[clPROT][ENTREZ].append(mEXP)
                        else:
                            ModMRNAD[clPROT][ENTREZ] = [mEXP]
DBGgenes.close()


###################################################
# It would appear that the IPI to Entrez mapping may be 
# many to one based on the output to DBGgenes.
#
# Let's look at the values for different IPIs for the same
# ENTREZ id.
DBGIPIEnt = open('DBGIPIEnt.txt','w')
for ENTREZ in IPIEntvals.keys():
    DBGIPIEnt.write(ENTREZ + "\t" + str(IPIEntvals[ENTREZ]) + "\n")
DBGIPIEnt.close()
# Based on inspection of this result, the typical case seems
# to be taht different isoforms map to the same ENTREZ id,
# so the most likely action will be to sum the values when passed
# to the model.
###################################################

for clPROT in n60_PROTtoCORE.keys():
    clCORE = n60_PROTtoCORE[clPROT]
    tissue = clCORE
    tissue = tissue.replace("(","_")
    tissue = tissue.replace(")","_")
    tissue = tissue.replace(" ","_")    
    tissue = tissue.replace("/","_")
    tissue = tissue.replace("-","_")    
    MoutFI = open('nci60mRNA/' + tissue + '.csv','w')
    PoutFI = open('nci60prot/' + tissue + '.csv','w')
    MPoutFI = open('nci60prot_mRNA/' + tissue + '.csv','w')
    MoutFI.write("gene\tmean\tvar\n")
    PoutFI.write("gene\tmean\tvar\n")
    MPoutFI.write("gene\tmean\tvar\n")    
    for g in ModPROTMRNA[clPROT]:
        mpEXP = np.array(ModPROTMRNA[clPROT][g])
        mpEXP = mpEXP[np.logical_not(np.isnan(mpEXP))]
        pEXP = NaN
        mEXP = NaN
        if len(mpEXP) > 0:
            mpEXP = sum(mpEXP)
            modMPcount = modMPcount + 1
        else:
            mpEXP = NaN    
        if ModPROTBoth[clPROT].has_key(g):
            pEXP = np.array(ModPROTBoth[clPROT][g])
            pEXP = pEXP[np.logical_not(np.isnan(pEXP))]
            if len(pEXP) > 0:
                pEXP = sum(pEXP)
                modPcount = modPcount + 1
            else:
                pEXP = NaN    
        if ModMRNAD[clPROT].has_key(g):    
            mEXP = np.array(ModMRNAD[clPROT][g])
            mEXP = mEXP[np.logical_not(np.isnan(mEXP))]
            if len(mEXP) > 0:
                mEXP = sum(mEXP)
                modMcount = modMcount + 1
            else:
                mEXP = NaN
            # This was already done above:    
            # mEXP = m_MtoP*mEXP + b_MtoP
        mpEXP = str(mpEXP)
        mEXP = str(mEXP)
        pEXP = str(pEXP)                
        outlist = [g, EntrezToIPI[g], clPROT, clCORE,
                   mEXP, pEXP, mpEXP]
        modelDataFI.write("\t".join(outlist)+"\n")
        outlist = [g, mEXP, "1"]
        MoutFI.write("\t".join(outlist)+"\n")
        outlist = [g, pEXP, "1"]
        PoutFI.write("\t".join(outlist)+"\n")
        outlist = [g, mpEXP, "1"]
        MPoutFI.write("\t".join(outlist)+"\n")
    MoutFI.close()
    PoutFI.close()
    MPoutFI.close()   
modelDataFI.close()

modMcount = modMcount/float(nCL)
modPcount = modPcount/float(nCL)
modMPcount = modMPcount/float(nCL)

print("Average # of model genes with data (mRNA, Protein, Both): (" + \
      str(modMcount) + ", " + str(modPcount) + ", " + str(modMPcount) + ").") 




# Now we write a global function that does much of the above anlaysis
# but no figures, for the purpose of systematically testing different
# zero thresholds.       
# Percentiles and zero values can be calculated in the calling loop.
(Mflat, Mflat_mod) = flattenExpression(OMRNAD)
Mflat = sorted(Mflat)
def zeroAdjustExpression(OMICRD, OPROTD, pz, mz):
    p_zero = y_notDeep[0]
    m_zero = Mflat[0]
    if pz >= 0:
        p_zero = y_notDeep[pz]
    if mz >= 0:    
        m_zero = Mflat[mz]
    LMICRD = copy.deepcopy(OMICRD)
    LPROTD = copy.deepcopy(OPROTD)
    # Zero adjust expression values
    for cl in LPROTD.keys():
        for g in LPROTD[cl].keys():
            if LPROTD[cl][g] == 0:
                LPROTD[cl][g] = NaN
            elif LPROTD[cl][g] < p_zero:
                LPROTD[cl][g] = 0
            else:
                LPROTD[cl][g] = LPROTD[cl][g] - p_zero    
    for cl in LMICRD.keys():
        for g in LMICRD[cl].keys():
            if LMICRD[cl][g] < m_zero:
                LMICRD[cl][g] = 0
            else:
                LMICRD[cl][g] = LMICRD[cl][g] - m_zero

    # Recalculate line of best fit for micr -> prot
    (x_b, y_b, x_ball, y_ball) = getProt_mRNA_pairs(LMICRD, LPROTD, IPIToEntrez)
    clf = linear_model.LinearRegression(fit_intercept=False)
    A = np.vstack([x_b, np.zeros(len(x_b))]).T
    clf.fit(A,y_b)
    m = clf.coef_[0]
                
    # Combine proteomic and mRNA data
    LPROTMRNA = copy.deepcopy(LPROTD)
    for cl in LMICRD.keys():
        for g in LMICRD[cl].keys():        
            if LPROTMRNA[cl].has_key(g):
                if LPROTMRNA[cl][g] != LPROTMRNA[cl][g]:
                    if LMICRD[cl][g] == LMICRD[cl][g]:
                        LPROTMRNA[cl][g] = m*LMICRD[cl][g]

            else:
                LPROTMRNA[cl][g] = m*LMICRD[cl][g]

    # Begin the process of summing isoform data as input to FALCON.
    thresh_label = str(m_zero)[0:5] + "_" + str(p_zero)[0:5] + "_" + str(m)[0:5]    
    if pz == -1:
        if not os.path.isdir('nci60mRNA_thresh'):
            if os.path.exists('nci60mRNA_thresh'):
                raise Exception("Specified output directory is a file!")
            os.mkdir('nci60mRNA_thresh')
        if not os.path.isdir('nci60mRNA_thresh/' + thresh_label):
            if os.path.exists('nci60mRNA_thresh/' + thresh_label):
                raise Exception("Specified output directory is a file!")
            os.mkdir('nci60mRNA_thresh/' + thresh_label)
    if mz == -1:
        if not os.path.isdir('nci60prot_thresh'):
            if os.path.exists('nci60prot_thresh'):
                raise Exception("Specified output directory is a file!")
            os.mkdir('nci60prot_thresh')
        if not os.path.isdir('nci60prot_thresh/' + thresh_label):
            if os.path.exists('nci60prot_thresh/' + thresh_label):
                raise Exception("Specified output directory is a file!")
            os.mkdir('nci60prot_thresh/' + thresh_label)
    if not os.path.isdir('nci60prot_mRNA_thresh'):
        if os.path.exists('nci60prot_mRNA_thresh'):
            raise Exception("Specified output directory is a file!")
        os.mkdir('nci60prot_mRNA_thresh')
    if not os.path.isdir('nci60prot_mRNA_thresh/' + thresh_label):
        if os.path.exists('nci60prot_mRNA_thresh/' + thresh_label):
            raise Exception("Specified output directory is a file!")
        os.mkdir('nci60prot_mRNA_thresh/' + thresh_label)

        
    IPIEntvals = {}
    testCL = LPROTMRNA.keys()[0]
    nCL = len(LPROTMRNA.keys())
    ModLMICRD = {}
    ModLPROTD = {}
    ModLPROTMRNA = {}
    for clPROT in n60_PROTtoCORE.keys():
        ModLMICRD[clPROT] = {}
        ModLPROTD[clPROT] = {}
        ModLPROTMRNA[clPROT] = {}
        clCORE = n60_PROTtoCORE[clPROT]
        tissue = clCORE
        tissue = tissue.replace("(","_")
        tissue = tissue.replace(")","_")
        tissue = tissue.replace(" ","_")    
        tissue = tissue.replace("/","_")
        tissue = tissue.replace("-","_")    
        for g in LPROTMRNA[clPROT]:
            ENTREZ = "UNKNOWN ENTREZ ID"
            if IPIToEntrez.has_key(g):
                ENTREZ = IPIToEntrez[g]  
                if modgenes.has_key(ENTREZ): 
                    mpEXP = LPROTMRNA[clPROT][g]
                    pEXP = NaN
                    mEXP = NaN
                    if mpEXP == mpEXP:
                        if ModLPROTMRNA[clPROT].has_key(ENTREZ):
                            ModLPROTMRNA[clPROT][ENTREZ].append(mpEXP)
                        else:
                            ModLPROTMRNA[clPROT][ENTREZ] = [mpEXP]
                        if clPROT == testCL:
                            if IPIEntvals.has_key(ENTREZ):
                                IPIEntvals[ENTREZ].append([g,mpEXP])
                            else:
                                IPIEntvals[ENTREZ] = [[g,mpEXP]]
                    if mz == -1:
                        if LPROTD[clPROT].has_key(g):
                            pEXP = LPROTD[clPROT][g]
                            if LPROTD[clPROT][g] == LPROTD[clPROT][g]:
                                if ModLPROTD[clPROT].has_key(ENTREZ):
                                    ModLPROTD[clPROT][ENTREZ].append(pEXP)
                                else:
                                    ModLPROTD[clPROT][ENTREZ] = [pEXP]
                    if pz == -1:
                        if LMICRD[clPROT].has_key(g):    
                            mEXP = m*LMICRD[clPROT][g]
                            if LMICRD[clPROT][g] == LMICRD[clPROT][g]:
                                if ModLMICRD[clPROT].has_key(ENTREZ):
                                    ModLMICRD[clPROT][ENTREZ].append(mEXP)
                                else:
                                    ModLMICRD[clPROT][ENTREZ] = [mEXP]

    for clPROT in n60_PROTtoCORE.keys():
        clCORE = n60_PROTtoCORE[clPROT]
        tissue = clCORE
        tissue = tissue.replace("(","_")
        tissue = tissue.replace(")","_")
        tissue = tissue.replace(" ","_")    
        tissue = tissue.replace("/","_")
        tissue = tissue.replace("-","_")
        MoutFI = 0
        PoutFI = 0
        if pz == -1:
            MoutFI = open('nci60mRNA_thresh/' + thresh_label + '/'
                          + tissue + '.csv','w')
            MoutFI.write("gene\tmean\tvar\n")            
        if mz == -1:    
            PoutFI = open('nci60prot_thresh/' + thresh_label + '/'
                          + tissue + '.csv','w')
            PoutFI.write("gene\tmean\tvar\n")
        MPoutFI = open('nci60prot_mRNA_thresh/' + thresh_label 
                       + '/' + tissue + '.csv','w')
        MPoutFI.write("gene\tmean\tvar\n")    
        for g in ModLPROTMRNA[clPROT]:
            mpEXP = np.array(ModLPROTMRNA[clPROT][g])
            mpEXP = mpEXP[np.logical_not(np.isnan(mpEXP))]
            if len(mpEXP) > 0:
                mpEXP = sum(mpEXP)
            else:
                mpEXP = NaN    
            mpEXP = str(mpEXP)
            outlist = [g, mpEXP, "1"]
            MPoutFI.write("\t".join(outlist)+"\n")
            if pz == -1:
                mEXP = NaN
                if ModLMICRD[clPROT].has_key(g):    
                    mEXP = np.array(ModLMICRD[clPROT][g])
                    mEXP = mEXP[np.logical_not(np.isnan(mEXP))]
                    if len(mEXP) > 0:
                        mEXP = sum(mEXP)
                    else:
                        mEXP = NaN
                mEXP = str(mEXP)
                outlist = [g, mEXP, "1"]
                MoutFI.write("\t".join(outlist)+"\n")
            if mz == -1:
                pEXP = NaN            
                if ModLPROTD[clPROT].has_key(g):
                    pEXP = np.array(ModLPROTD[clPROT][g])
                    pEXP = pEXP[np.logical_not(np.isnan(pEXP))]
                    if len(pEXP) > 0:
                        pEXP = sum(pEXP)
                    else:
                        pEXP = NaN    
                pEXP = str(pEXP)                
                outlist = [g, pEXP, "1"]
                PoutFI.write("\t".join(outlist)+"\n")
        if pz == -1:
            MoutFI.close()
        if mz == -1:
            PoutFI.close()
        MPoutFI.close()   

if True:        
    num_Pint = 50
    num_Mint = 100
    PZeroEnd = np.argmax(Rbott)
    # For now, it appears we have a reasonable range to search for zero
    # values in for protein intensities; this is not really the case
    # for microarray intensities. Once we have corresponding RNA-Seq data
    # for some cell lines, it will be interesting to know where these fall.
    print("Considering first " + str(PZeroEnd) + " protein pairs for zero"
          + " adjustment, up to max(R) = " + str(Rbott[PZeroEnd]) + ".")
    for i in range(-1, num_Mint):
        for j in range(-1, num_Pint):
            mz = -1
            pz = -1
            if i >= 0:    
                mzVAL = np.percentile(Mflat,int(100.0*i/num_Mint))
                mz = np.searchsorted(Mflat, mzVAL)
            if j >= 0:
                pzVAL = np.percentile(y_notDeep[0:PZeroEnd+1],int(100.0*j/num_Pint))
                pz = np.searchsorted(y_notDeep, pzVAL)                
            zeroAdjustExpression(OMRNAD, OPROTBoth, pz, mz)
