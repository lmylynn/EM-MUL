#! /usr/bin/env python
import sys, time, os, array, optparse, math
from collections import Counter
from decimal import *
import numpy as np
import pandas as pd
from pandas import Series,DataFrame
parser = optparse.OptionParser()

#score matrix
data1 = [[6,-18,-18,-18,-25],[-18,6,-18,3,-25],[-18,-18,6,-18,-25],[-18,-18,-18,3,-25],[-25,-25,-25,-25,25]]
data2 = [[3,-18,-18,-18,-25],[-18,6,-18,-18,-25],[3,-18,6,-18,-25],[-18,-18,-18,6,-25],[-25,-25,-25,-25,25]]
data3 = [[6,-18,-18,-18,-25],[-18,6,-18,3,-25],[-18,-18,6,-18,-25],[-18,3,-18,3,-25],[-25,-25,-25,-25,25]]
data4 = [[3,-18,3,-18,-25],[-18,6,-18,-18,-25],[3,-18,6,-18,-25],[-18,-18,-18,6,-25],[-25,-25,-25,-25,25]]
sc1={'SA':-48,'SC':-27,'SG':-48,'ST':-51}
sc2={'SA':-51,'SC':-48,'SG':-27,'ST':-48}
S=-349
S2=-328
gap_open = -100
gap_extension = -20
score_matrix1 =DataFrame(data1,columns=['A','C','G','T','N'],index=['A','C','G','T','N']) #index
score_matrix2 =DataFrame(data2,columns=['A','C','G','T','N'],index=['A','C','G','T','N']) #
score_matrix3 =DataFrame(data3,columns=['A','C','G','T','N'],index=['A','C','G','T','N']) #
score_matrix4 =DataFrame(data4,columns=['A','C','G','T','N'],index=['A','C','G','T','N']) #

parser.add_option("-r", "--ref", dest="reffile", metavar="FILE", default="")#ref
parser.add_option("-u", "--uni", dest="unifile", metavar="FILE", default="")#unique reads
parser.add_option("-m", "--multi", dest="multifile", metavar="FILE", default="")#multireads
parser.add_option("-o", "--over", dest="overfile", metavar="FILE", default="")#overfile

options, infiles = parser.parse_args()

#
ref1, ref2, cr, seq = {}, {}, '', ''
for line in open(options.reffile):
    if line[0] == '>': 
        if len(cr) > 0:
            ref1[cr] = seq.upper() #
            ref2[cr] = seq.upper()
        cr, seq = line[1:-1].split()[0], ''
    else: seq += line.strip()

ref1[cr] = seq.upper()
ref2[cr] = seq.upper()
del seq
print('reference file is already loaded: %s' % options.reffile)

cover = {}

#build arrays of ATCG，number，cr：
for cr in ref1:
    cover[cr] = array.array('H',[0])*len(ref1[cr])
print("building the array of coverage")


#get seq of unifile,ref，position
def get_info(line):
    if line[0] == '@': return []
    col = line.split('\t')
    flag = col[1]
    #ref_cr，flag, pos-1，cigar MID，seq
    cr, flag, pos, cigar, seq = col[2], col[1], int(col[3])-1, col[5], col[9]
    gap_pos, gap_size = 0, 0
    while 'I' in cigar or 'D' in cigar:
        for sep in 'MID':
            try: gap_size = int(cigar.split(sep, 1)[0]) #get number of M/I/D
            except ValueError: continue
            break
        if sep == 'M': gap_pos += gap_size
        elif sep == 'I': 
            seq = seq[:gap_pos] + seq[gap_pos+gap_size:]
        elif sep == 'D':
            seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
            gap_pos += gap_size
        cigar = cigar[cigar.index(sep)+1:]
    if pos + len(seq) >= len(ref1[cr]): return []
    return (seq, flag, cr, pos)

#get information of multireads
#read name，flag，ref_name，pos，cigar MID，n, seq
def get_multiinfo(line):
    if line[0] == '@': return []
    col = line.split('\t')
    read_name, flag, ref_name, pos, cigar, seq, qual = col[0], col[1], col[2], int(col[3])-1, col[5], col[9], col[10]
    gap_pos, gap_size = 0, 0
    while 'I' in cigar or 'D' in cigar:
        for sep in 'MID':
            try: gap_size = int(cigar.split(sep, 1)[0]) #
            except ValueError: continue
            break
        if sep == 'M': gap_pos += gap_size
        elif sep == 'I': 
            seq = seq[:gap_pos] + seq[gap_pos+gap_size:]
            qual = qual[:gap_pos] + qual[gap_pos+gap_size:]
        elif sep == 'D':
            seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
            qual = qual[:gap_pos] + '!' * gap_size + qual[gap_pos:]
            gap_pos += gap_size
        cigar = cigar[cigar.index(sep)+1:]
    if pos + len(seq) >= len(ref1[ref_name]): return []
    cigar = col[5]
    n = len(seq)
    return(read_name, flag, ref_name, pos, cigar, n, seq, qual)

def get_overinfo(cigar,seq,qual):
    gap_pos, gap_size = 0, 0
    while 'I' in cigar or 'D' in cigar:
        for sep in 'MID':
            try: gap_size = int(cigar.split(sep, 1)[0]) #
            except ValueError: continue
            break
        if sep == 'M': gap_pos += gap_size
        elif sep == 'I': 
            seq = seq[:gap_pos] + seq[gap_pos+gap_size:]
            qual = qual[:gap_pos] + qual[gap_pos+gap_size:]
        elif sep == 'D':
            seq = seq[:gap_pos] + '-' * gap_size + seq[gap_pos:]
            qual = qual[:gap_pos] + '!' * gap_size + qual[gap_pos:]
            gap_pos += gap_size
        cigar = cigar[cigar.index(sep)+1:]
    arr = phred_seq(qual)
    return (seq, arr)

#coverge array
before_cr = ''
for line in open(options.unifile):
    #print(line)
    getinfo = get_info(line)
    if len(getinfo) == 0:continue
    seq, flag, cr, pos = getinfo
    if before_cr != cr:
        cover_cr = cover[cr]
    pos2 = pos + len(seq)
    for i in range(pos, pos2):
        cover_cr[i] = cover_cr[i] + 1    
print("依次遍历unique reads，获得覆盖度数组")

tem_name = ''
a = []
b = []
list_p = {}
list_M_value = {}

def get_match_score1(s1,s2):                    #get score 
    score = score_matrix1[s1][s2]
    return score
def get_match_score2(s1,s2):                    #
    score = score_matrix2[s1][s2]
    return score
def get_match_score3(s1,s2):                    #
    score = score_matrix3[s1][s2]
    return score
def get_match_score4(s1,s2):                    #
    score = score_matrix4[s1][s2]
    return score    
pAC = 0.00025
pAT = 0.00025
pAG = 0.0005
pCA = 0.00025
pCT = 0.0005
pCG = 0.00025
pTC = 0.0005
pTA = 0.00025
pTG = 0.00025
pGA = 0.0005
pGT = 0.00025
pGC = 0.00025
pSNP = 0.001
CG = 0.85
CH = 0.02


def trASCII(a):
    Q = ord(a) - 33 
    p = math.pow(10, (-0.1*Q))
    p1 = '%.7f'%p
    return float(p1)

def phred_seq(qual):
    arr = [0.0 for i in range(len(qual))]
    for i in range(0,len(qual)):
        arr[i] = trASCII(qual[i])
    return arr
    
#0,256
def pro_ref_and_multi1(seq,s1,arr):
    pro_rm = [0.0 for i in range(len(seq))]
    gap_num = 0
    for i in range(0,len(seq)):
        if arr[i] == 1.0:
            gap_num = gap_num + 1
            if gap_num == 1:
                pro_rm[i] = 0.5*gap_open
            else:
                pro_rm[i] = 0.5*gap_extension
        elif arr[i] > 0:
            gap_num = 0 
            if s1[i] == 'A':
                if seq[i] == 'A':
                    pro_rm[i] = 0.5*((1-pSNP)*get_match_score1('A','A')+pSNP*(sc1['SA']-get_match_score1('A','A'))/3)
                elif seq[i] == 'C':
                    if ((i+1)<len(seq) and seq[i+1] == 'G'):
                        pro_rm[i] = 0.5*(pAC*CG*get_match_score1('C','A')+(1-pAC*CG)*(sc1['SA']-get_match_score1('C','A'))/3)
                    elif ((i+1)<len(seq) and s1[i+1] == 'G'):
                        pro_rm[i] = 0.5*(pAC*CG*get_match_score1('C','A')+(1-pAC*CG)*(sc1['SA']-get_match_score1('C','A'))/3)
                    else:
                        pro_rm[i] = 0.5*(pAC*CH*get_match_score1('C','A')+(1-pAC*CH)*(sc1['SA']-get_match_score1('C','A'))/3)
                elif seq[i] == 'G':
                    pro_rm[i] = 0.5*(pAG*get_match_score1('G','A')+(1-pAG)*(sc1['SA']-get_match_score1('G','A'))/3)
                elif seq[i] == 'T':
                    if ((i+1)<len(seq) and seq[i+1] == 'G'):
                        pro_rm[i] = 0.5*((pAT + pAC * (1-CG))*get_match_score1('T','A')+(1-(pAT + pAC * (1-CG)))*(sc1['SA']-get_match_score1('T','A'))/3)
                    elif ((i+1)<len(seq) and s1[i+1] == 'G'):
                        pro_rm[i] = 0.5*((pAT + pAC * (1-CG))*get_match_score1('T','A')+(1-(pAT + pAC * (1-CG)))*(sc1['SA']-get_match_score1('T','A'))/3)
                    else:
                        pro_rm[i] = 0.5*((pAT + pAC * (1-CH))*get_match_score1('T','A')+(1-(pAT + pAC * (1-CH)))*(sc1['SA']-get_match_score1('T','A'))/3)
                else:
                    pro_rm[i] = 0.5*get_match_score1('N','A')
            elif s1[i] == 'C':
                if seq[i] == 'A':
                    pro_rm[i] = 0.5*(pCA*get_match_score1('A','C')+(1-pCA)*(sc1['SC']-get_match_score1('A','C'))/3)
                elif seq[i] == 'C':
                    if ((i+1)<len(seq) and seq[i+1] == 'G'):
                        pro_rm[i] = 0.5*((1-pSNP)*CG*get_match_score1('C','C')+(1-(1-pSNP)*CG)*(sc1['SC']-get_match_score1('C','C'))/3)
                    elif ((i+1)<len(seq) and s1[i+1] == 'G'):
                        pro_rm[i] = 0.5*((1-pSNP)*CG*get_match_score1('C','C')+(1-(1-pSNP)*CG)*(sc1['SC']-get_match_score1('C','C'))/3)
                    else:
                        pro_rm[i] = 0.5*((1-pSNP)*CH*get_match_score1('C','C')+(1-(1-pSNP)*CH)*(sc1['SC']-get_match_score1('C','C'))/3)
                elif seq[i] == 'G':
                    pro_rm[i] = 0.5*(pCG*get_match_score1('G','C')+(1-pCG)*(sc1['SC']-get_match_score1('G','C'))/3)
                elif seq[i] == 'T':
                    if ((i+1)<len(seq) and seq[i+1] == 'G'):
                        pro_rm[i] = 0.5*((pCT+(1-pSNP)*(1-CG))*get_match_score1('T','C')+(1-pCT-(1-pSNP)*(1-CG))*(sc1['SC']-get_match_score1('T','C'))/3)
                    elif ((i+1)<len(seq) and s1[i+1] == 'G'):
                        pro_rm[i] = 0.5*((pCT+(1-pSNP)*(1-CG))*get_match_score1('T','C')+(1-pCT-(1-pSNP)*(1-CG))*(sc1['SC']-get_match_score1('T','C'))/3)
                    else:
                        pro_rm[i] = 0.5*((pCT+(1-pSNP)*(1-CH))*get_match_score1('T','C')+(1-pCT-(1-pSNP)*(1-CH))*(sc1['SC']-get_match_score1('T','C'))/3)
                else:
                    pro_rm[i] = 0.5*get_match_score1('N','C')
            elif s1[i] == 'G':
                if seq[i] == 'A':
                    pro_rm[i] = 0.5*(pGA*get_match_score1('A','G')+(1-pGA)*(sc1['SG']-get_match_score1('A','G'))/3)
                elif seq[i] == 'C':
                    if ((i+1)<len(seq) and seq[i+1] == 'G'):
                        pro_rm[i] = 0.5*(pGC*CG*get_match_score1('C','G')+(1-pGC*CG)*(sc1['SG']-get_match_score1('C','G'))/3)
                    elif ((i+1)<len(seq) and s1[i+1] == 'G'):
                        pro_rm[i] = 0.5*(pGC*CG*get_match_score1('C','G')+(1-pGC*CG)*(sc1['SG']-get_match_score1('C','G'))/3)
                    else:
                        pro_rm[i] = 0.5*(pGC*CH*get_match_score1('C','G')+(1-pGC*CH)*(sc1['SG']-get_match_score1('C','G'))/3)
                elif seq[i] == 'G':
                    pro_rm[i] = 0.5*((1-pSNP)*get_match_score1('G','G')+pSNP*(sc1['SG']-get_match_score1('G','G'))/3)
                elif seq[i] == 'T':
                    if ((i+1)<len(seq) and seq[i+1] == 'G'):
                        pro_rm[i] = 0.5*((pGT + pGC * (1-CG))*get_match_score1('T','G')+(1-(pGT + pGC * (1-CG)))*(sc1['SG']-get_match_score1('T','G'))/3)
                    elif ((i+1)<len(seq) and s1[i+1] == 'G'):
                        pro_rm[i] = 0.5*((pGT + pGC * (1-CG))*get_match_score1('T','G')+(1-(pGT + pGC * (1-CG)))*(sc1['SG']-get_match_score1('T','G'))/3)
                    else:
                        pro_rm[i] = 0.5*((pGT + pGC * (1-CH))*get_match_score1('T','G')+(1-(pGT + pGC * (1-CH)))*(sc1['SG']-get_match_score1('T','G'))/3)
                else:
                    pro_rm[i] = 0.5*get_match_score1('N','G')
            elif s1[i] == 'T':
                if seq[i] == 'A':
                    pro_rm[i] = 0.5*(pTA*get_match_score1('A','T')+(1-pTA)*(sc1['ST']-get_match_score1('A','T'))/3)
                elif seq[i] == 'C':
                    if ((i+1)<len(seq) and seq[i+1] == 'G'):
                        pro_rm[i] = 0.5*(pTC*CG*get_match_score1('C','T')+(1-pTC*CG)*(sc1['ST']-get_match_score1('C','T'))/3)
                    elif ((i+1)<len(seq) and s1[i+1] == 'G'):
                        pro_rm[i] = 0.5*(pTC*CG*get_match_score1('C','T')+(1-pTC*CG)*(sc1['ST']-get_match_score1('C','T'))/3)
                    else:
                        pro_rm[i] = 0.5*(pTC*CH*get_match_score1('C','T')+(1-pTC*CH)*(sc1['ST']-get_match_score1('C','T'))/3)
                elif seq[i] == 'G':
                    pro_rm[i] = 0.5*(pTG*get_match_score1('G','T')+(1-pTG)*(sc1['ST']-get_match_score1('G','T'))/3)
                elif seq[i] == 'T':
                    if ((i+1)<len(seq) and seq[i+1] == 'G'):
                        pro_rm[i] = 0.5*(((1-pSNP) + pTC * (1-CG))*get_match_score1('T','T')+(pSNP - pTC * (1-CG))*(sc1['ST']-get_match_score1('T','T'))/3)
                    elif ((i+1)<len(seq) and s1[i+1] == 'G'):
                        pro_rm[i] = 0.5*(((1-pSNP) + pTC * (1-CG))*get_match_score1('T','T')+(pSNP - pTC * (1-CG))*(sc1['ST']-get_match_score1('T','T'))/3)
                    else:
                        pro_rm[i] = 0.5*(((1-pSNP) + pTC * (1-CH))*get_match_score1('T','T')+(pSNP - pTC * (1-CH))*(sc1['ST']-get_match_score1('T','T'))/3)
                else:
                    pro_rm[i] = 0.5*get_match_score1('N','T')
            else:
                pro_rm[i] = 0.5*get_match_score1('N','N')
        else:
            gap_num = 0
            pro_rm[i] = 0.5*get_match_score1(seq[i], s1[i])
    #print(pro_rm)
    return (pro_rm,arr)
#16,272
def pro_ref_and_multi2(seq,s1,arr):
    pro_rm = [0.0 for i in range(len(seq))]
    gap_num = 0
    for i in range(0,len(seq)):
        if arr[i] == 1.0:
            gap_num = gap_num + 1
            if gap_num == 1:
                pro_rm[i] = 0.5*gap_open
            else:
                pro_rm[i] = 0.5*gap_extension
        elif arr[i] > 0:
            gap_num = 0 
            if s1[i] == 'T':
                if seq[i] == 'C':
                    pro_rm[i] = 0.5*(pTC*get_match_score2('C','T')+(1-pTC)*(sc2['ST']-get_match_score2('C','T'))/3)
                elif seq[i] == 'G':
                    if ((i-1)>=0 and seq[i-1] == 'C'):
                        pro_rm[i] = 0.5*(pTG*CG*get_match_score2('G','T')+(1-pTG*CG)*(sc2['ST']-get_match_score2('G','T'))/3)        
                    elif ((i-1)>=0 and s1[i-1] == 'C'):
                        pro_rm[i] = 0.5*(pTG*CG*get_match_score2('G','T')+(1-pTG*CG)*(sc2['ST']-get_match_score2('G','T'))/3)
                    else:
                        pro_rm[i] = 0.5*(pTG*CH*get_match_score2('G','T')+(1-pTG*CH)*(sc2['ST']-get_match_score2('G','T'))/3)
                elif seq[i] == 'T':
                    pro_rm[i] = 0.5*((1-pSNP)*get_match_score2('T','T')+ pSNP*(sc2['ST']-get_match_score2('T','T'))/3)
                elif seq[i] == 'A':
                    if ((i-1)>=0 and seq[i-1] == 'C'):
                        pro_rm[i] = 0.5*((pTA + pTG * (1-CG))*get_match_score2('A','T')+(1-(pTA + pTG * (1-CG)))*(sc2['ST']-get_match_score2('A','T'))/3)
                    elif ((i-1)>=0 and s1[i-1] == 'C'):
                        pro_rm[i] = 0.5*((pTA + pTG * (1-CG))*get_match_score2('A','T')+(1-(pTA + pTG * (1-CG)))*(sc2['ST']-get_match_score2('A','T'))/3)
                    else:
                        pro_rm[i] = 0.5*((pTA + pTG * (1-CH))*get_match_score2('A','T')+(1-(pTA + pTG * (1-CH)))*(sc2['ST']-get_match_score2('A','T'))/3)
                else:
                    pro_rm[i] = 0.5*get_match_score2('N','T')
            elif s1[i] == 'A':
                if seq[i] == 'C':
                    pro_rm[i] = 0.5*(pAC*get_match_score2('C','A')+(1-pAC)*(sc2['SA']-get_match_score2('C','A'))/3)
                elif seq[i] == 'G':
                    if ((i-1)>=0 and seq[i-1] == 'C'):
                        pro_rm[i] = 0.5*(AG*CG*get_match_score2('G','A')+(1-AG*CG)*(sc2['SA']-get_match_score2('G','A'))/3)
                    elif ((i-1)>=0 and s1[i-1] == 'C'):
                        pro_rm[i] = 0.5*(AG*CG*get_match_score2('G','A')+(1-AG*CG)*(sc2['SA']-get_match_score2('G','A'))/3)
                    else:
                        pro_rm[i] = 0.5*(AG*CH*get_match_score2('G','A')+(1-AG*CH)*(sc2['SA']-get_match_score2('G','A'))/3)
                elif seq[i] == 'T':
                    pro_rm[i] = 0.5*(pAT*get_match_score2('T','A')+(1-pAT)*(sc2['SA']-get_match_score2('T','A')))
                elif seq[i] == 'A':
                    if ((i-1)>=0 and seq[i-1] == 'C'):
                        pro_rm[i] = 0.5*((1-pSNP+pAG*(1-CG))*get_match_score2('A','A')+(pSNP-pAG*(1-CG))*(sc2['SA']-get_match_score2('A','A'))/3)
                    elif ((i-1)>=0 and s1[i-1] == 'C'):
                        pro_rm[i] = 0.5*((1-pSNP+pAG*(1-CG))*get_match_score2('A','A')+(pSNP-pAG*(1-CG))*(sc2['SA']-get_match_score2('A','A'))/3)
                    else:
                        pro_rm[i] = 0.5*((1-pSNP+pAG*(1-CH))*get_match_score2('A','A')+(pSNP-pAG*(1-CH))*(sc2['SA']-get_match_score2('A','A'))/3)
                else:
                    pro_rm[i] = 0.5*get_match_score2('N','A')
            elif s1[i] == 'C':
                if seq[i] == 'T':
                    pro_rm[i] = 0.5*(pCT*get_match_score2('T','C')+(1-pCT)*(sc2['SC']-get_match_score2('T','C')/3))
                elif seq[i] == 'G':
                    if ((i-1)>=0 and seq[i-1] == 'C'):
                        pro_rm[i] = 0.5*(pCG*CG*get_match_score2('G','C')+(1-pCG*CG)*(sc2['SC']-get_match_score2('G','C'))/3)
                    elif ((i-1)>=0 and s1[i-1] == 'C'):
                        pro_rm[i] = 0.5*(pCG*CG*get_match_score2('G','C')+(1-pCG*CG)*(sc2['SC']-get_match_score2('G','C'))/3)
                    else:
                        pro_rm[i] = 0.5*(pCG*CH*get_match_score2('G','C')+(1-pCG*CH)*(sc2['SC']-get_match_score2('G','C'))/3)
                elif seq[i] == 'C':
                    pro_rm[i] = 0.5*((1-pSNP)*get_match_score2('C','C')+pSNP*(sc2['SC']-get_match_score2('C','C'))/3)
                elif seq[i] == 'A':
                    if ((i-1)>=0 and seq[i-1] == 'C'):
                        pro_rm[i] = 0.5*((pCA + pCG * (1-CG))*get_match_score2('A','C')+(1-(pCA + pCG * (1-CG)))*(sc2['SC']-get_match_score2('A','C'))/3)
                    elif ((i-1)>=0 and s1[i-1] == 'C'):
                        pro_rm[i] = 0.5*((pCA + pCG * (1-CG))*get_match_score2('A','C')+(1-(pCA + pCG * (1-CG)))*(sc2['SC']-get_match_score2('A','C'))/3)
                    else:
                        pro_rm[i] = 0.5*((pCA + pCG * (1-CH))*get_match_score2('A','C')+(1-(pCA + pCG * (1-CH)))*(sc2['SC']-get_match_score2('A','C'))/3)
                else:
                    pro_rm[i] = 0.5*get_match_score2('N','C')
            elif s1[i] == 'G':
                if seq[i] == 'C':
                    pro_rm[i] = 0.5*(pGC*get_match_score2('C','G')+(1-pGC)*(sc2['SG']-get_match_score2('C','G'))/3)
                elif seq[i] == 'G':
                    if ((i-1)>=0 and seq[i-1] == 'C'):
                        pro_rm[i] = 0.5*((1-pSNP)*CG*get_match_score2('G','G')+(1-(1-pSNP)*CG)*(sc2['SG']-get_match_score2('G','G'))/3)
                    elif ((i-1)>=0 and s1[i-1] == 'C'):
                        pro_rm[i] = 0.5*((1-pSNP)*CG*get_match_score2('G','G')+(1-(1-pSNP)*CG)*(sc2['SG']-get_match_score2('G','G'))/3)
                    else:
                        pro_rm[i] = 0.5*((1-pSNP)*CH*get_match_score2('G','G')+(1-(1-pSNP)*CH)*(sc2['SG']-get_match_score2('G','G'))/3)
                elif seq[i] == 'T':
                    pro_rm[i] = 0.5*(pGT*get_match_score2('T','G')+(1-pGT)*(sc2['SG']-get_match_score2('T','G'))/3)
                elif seq[i] == 'A':
                    if ((i-1)>=0 and seq[i-1] == 'C'):
                        pro_rm[i] = 0.5*((pGA + (1-pSNP)* (1-CG))*get_match_score2('A','G')+(1-(pGA + (1-pSNP)* (1-CG)))*(sc2['SG']-get_match_score2('A','G'))/3)
                    elif ((i-1)>=0 and s1[i-1] == 'C'):
                        pro_rm[i] = 0.5*((pGA + (1-pSNP)* (1-CG))*get_match_score2('A','G')+(1-(pGA + (1-pSNP)* (1-CG)))*(sc2['SG']-get_match_score2('A','G'))/3)
                    else:
                        pro_rm[i] = 0.5*((pGA + (1-pSNP)* (1-CH))*get_match_score2('A','G')+(1-(pGA + (1-pSNP)* (1-CH)))*(sc2['SG']-get_match_score2('A','G'))/3)
                else:
                    pro_rm[i] = 0.5*get_match_score2('N','G')
            else:
                pro_rm[i] = 0.5*get_match_score2('N','N')
        else:
            gap_num = 0
            pro_rm[i] = 0.5*get_match_score2(seq[i],s1[i])
    #print(pro_rm)
    return (pro_rm,arr)

overinf= {}
for line in open(options.overfile):
    col = line.split('\t')
    if abs(int(col[4])-int(col[9]))==256 or abs(int(col[4])-int(col[9]))==0:
        name = col[3]+'\t'+col[0]+'\t'+str(int(col[1])-1)
        information = col[4]+'\t'+str(int(col[6])-1)+'\t'+col[11]+'\t'+col[15]+'\t'+col[16]
        if name in overinf.keys():
            overinf[name] = overinf[name]+'\t'+information
        else:
            overinf[name] = information
print("unique_overlap_read_file has been loaded.")

def rev_base(base):
    if base == 'A':
        return 'T'
    elif base == 'T':
        return 'A'
    elif base == 'C':
        return 'G'
    elif base == 'G':
        return 'C'
    elif base == 'N':
        return 'N'
    else:
        return '-'
    
def rev_seq(seq):
    seq = seq[::-1]
    rev = ''
    for i in range(0,len(seq)):
        rev=rev+rev_base(seq[i])
    return rev

def rev_qual(qual):
    rev = qual[::-1]
    return rev
def like(a):
    if a>0:
        return math.log(a,10)
    elif a<0:
        return -math.log(-a,10)
#multireads，seq=read(M and D)
for line in open(options.multifile):
    getm = get_multiinfo(line)
    if len(getm) == 0: continue
    read_name, flag, ref_name, pos, cigar, n, seq, qual = getm
    qual = qual.strip()
    likehood = 0.0
    maxbaseTotal = 0
    score_ref=0.0
    if flag == '0' or flag == '256':#forward
        s1 = ref1[ref_name][pos:(pos+n)]#ref
        (pro_rm,arr) = pro_ref_and_multi1(seq,s1,phred_seq(qual))
        for q in range(0,len(seq)):
            score_ref=score_ref+pro_rm[q]
        pro_S = [0.0 for i in range(len(seq))]
        pro_post = [0.0 for i in range(len(seq))]
        if (read_name+'\t'+ref_name+'\t'+str(pos)) in overinf.keys():
            col = overinf[read_name+'\t'+ref_name+'\t'+str(pos)].split('\t')
            number_over = int(len(col)/5)
            names=locals()
            for i in range(0,number_over):
                names['uni_seq'+str(i)], names['uni_err'+str(i)] = get_overinfo(col[2+i*5],col[3+i*5],col[4+i*5].strip())
                names['uni_pos'+str(i)] = int(col[1+i*5])
                names['gap_flag'+str(i)] = 0         
            for i in range(0,len(seq)):#position i of seq 
                baseTotal= 0#number of unique reads
                pro_observed = [0.0 for i in range(number_over)]
                for j in range(0,number_over):
                    uni_seq, uni_err, uni_pos = names['uni_seq'+str(j)], names['uni_err'+str(j)], names['uni_pos'+str(j)]
                    p = 0.0
                    if (i<len(qual) and uni_pos<=(pos+i) and (uni_pos+len(uni_seq)-1)>=(pos+i)):
                        if uni_seq[pos+i-uni_pos]==seq[i]:
                            if seq[i] in 'ACTGN' and uni_seq[pos+i-uni_pos] in 'ATCGN':
                                names['gap_flag'+str(j)] = 0
                                p = 1-uni_err[pos+i-uni_pos]-arr[i]+uni_err[pos+i-uni_pos]*arr[i]
                                pro_observed[j] = p*get_match_score3(seq[i],seq[i])+(1-p)*(S2-get_match_score3(seq[i],seq[i]))/24
                            else:
                                names['gap_flag'+str(j)] = names['gap_flag'+str(j)] + 1
                                if names['gap_flag'+str(j)] == 1:
                                    pro_observed[j] = gap_open
                                else:
                                    pro_observed[j] = gap_extension
                        else:
                            if seq[i] in 'ACTGN' and uni_seq[pos+i-uni_pos] in 'ATCGN':
                                names['gap_flag'+str(j)] = 0
                                p = uni_err[pos+i-uni_pos]+arr[i]-(uni_err[pos+i-uni_pos]*arr[i])
                                pro_observed[j] = p*get_match_score3(seq[i], uni_seq[pos+i-uni_pos])+(1-p)*(S2-get_match_score3(seq[i], uni_seq[pos+i-uni_pos]))/24
                            else:
                                names['gap_flag'+str(j)] = names['gap_flag'+str(j)] + 1
                                if names['gap_flag'+str(j)] == 1:
                                    pro_observed[j] = gap_open
                                else:
                                    pro_observed[j] = gap_extension
                        baseTotal = baseTotal + 1
                if baseTotal>maxbaseTotal:
                    maxbaseTotal = baseTotal
                if baseTotal!=0:
                    n = 0.0
                    for q in range(0,number_over):
                        n = n + pro_observed[q]
                    pro_post[i] = n/baseTotal
                    pro_S[i] = pro_rm[i]+0.5*pro_post[i]
                else:
                    pro_post[i] = 0
                    pro_S[i] = 2*pro_rm[i]   
            for f in range(0,len(seq)):#seq
                likehood =  likehood + like(pro_S[f])
        else:# no overllapped unique reads
            for i in range(0,len(seq)):
                pro_S[i] = 2*pro_rm[i]
            for f in range(0,len(seq)):#
                likehood =  likehood + like(pro_S[f])
                
    else:#flag=16/272
        s1 = ref2[ref_name][pos:(pos+n)]#ref
        (pro_rm,arr) = pro_ref_and_multi2(seq,s1,phred_seq(qual))     
        for q in range(0,len(seq)):
            score_ref=score_ref+pro_rm[q]        
        pro_S = [0.0 for i in range(len(seq))]
        pro_post = [0.0 for i in range(len(seq))]
        if (read_name+'\t'+ref_name+'\t'+str(pos)) in overinf.keys():
            col = overinf[read_name+'\t'+ref_name+'\t'+str(pos)].split('\t')
            number_over = int(len(col)/5)
            names=locals()
            for i in range(0,number_over):
                names['uni_seq'+str(i)], names['uni_err'+str(i)] = get_overinfo(col[2+i*5],col[3+i*5],col[4+i*5].strip())
                names['uni_pos'+str(i)] = int(col[1+i*5])
                names['gap_flag'+str(i)] = 0         
            for i in range(0,len(seq)):#
                baseTotal= 0#
                pro_observed = [0.0 for i in range(number_over)]
                for j in range(0,number_over):
                    uni_seq, uni_err, uni_pos = names['uni_seq'+str(j)], names['uni_err'+str(j)], names['uni_pos'+str(j)]
                    p = 0.0
                    if (i<len(qual) and uni_pos<=(pos+i) and (uni_pos+len(uni_seq)-1)>=(pos+i)):
                        if uni_seq[pos+i-uni_pos]==seq[i]:
                            if seq[i] in 'ACTGN' and uni_seq[pos+i-uni_pos] in 'ATCGN':
                                names['gap_flag'+str(j)] = 0
                                p = 1-uni_err[pos+i-uni_pos]-arr[i]+uni_err[pos+i-uni_pos]*arr[i]
                                pro_observed[j] = p*get_match_score4(seq[i],seq[i])+(1-p)*(S2-get_match_score4(seq[i],seq[i]))/24
                            else:
                                names['gap_flag'+str(j)] = names['gap_flag'+str(j)] + 1
                                if names['gap_flag'+str(j)] == 1:
                                    pro_observed[j] = gap_open
                                else:
                                    pro_observed[j] = gap_extension
                        else:
                            if seq[i] in 'ACTGN' and uni_seq[pos+i-uni_pos] in 'ATCGN':
                                names['gap_flag'+str(j)] = 0
                                p = uni_err[pos+i-uni_pos]+arr[i]-(uni_err[pos+i-uni_pos]*arr[i])
                                pro_observed[j] = p*get_match_score4(seq[i], uni_seq[pos+i-uni_pos])+(1-p)*(S2-get_match_score4(seq[i], uni_seq[pos+i-uni_pos]))/24
                            else:
                                names['gap_flag'+str(j)] = names['gap_flag'+str(j)] + 1
                                if names['gap_flag'+str(j)] == 1:
                                    pro_observed[j] = gap_open
                                else:
                                    pro_observed[j] = gap_extension
                        baseTotal = baseTotal + 1
                if baseTotal>maxbaseTotal:
                    maxbaseTotal = baseTotal                        
                if baseTotal!=0:
                    n = 0.0
                    for q in range(0,number_over):
                        n = n + pro_observed[q]
                    pro_post[i] = n/baseTotal
                    pro_S[i] = pro_rm[i]+0.5*pro_post[i]
                else:
                    pro_post[i] = 0
                    pro_S[i] = 2*pro_rm[i]   
            for f in range(0,len(seq)):#
                likehood =  likehood + like(pro_S[f])
        else:#
            for i in range(0,len(seq)):
                pro_S[i] = 2*pro_rm[i]
            for f in range(0,len(seq)):#
                likehood =  likehood + like(pro_S[f])
    
    data = '\t'.join([str(x) for x in [likehood, read_name, flag, ref_name, (pos+1), len(seq), maxbaseTotal,score_ref]])
    f = open('new_score_all.txt','a')
    f.write(data)
    f.write('\n')

read_name = ''
highscore = {}
secondscore = {}
highS = {}
secondS = {}
flagnum={}
a, b = [], []

def maxandsecond(b):
    b.sort()
    max=b[-1]
    second=b[0]
    flag=0
    for i in range(1,len(b)):
        if b[i] < max:
            second=b[i]
    if len(b)==1:
        flag=0
    elif b[-1]==b[-2]:#
        flag=1
    else:
        flag=0
    return (max,second,flag)
#highscore = {read name :high score}
f = open('new_score_all.txt')
lines = f.readlines()
last_line = lines[-1].strip()
for line in lines:
    line = line.strip()
    col = line.split('\t')
    if read_name.strip()=='':
        read_name = col[1]
        a.append(float(col[0]))
        b.append(float(col[7]))
    elif line == last_line:
        if col[1] == read_name:
            a.append(float(col[0]))
            b.append(float(col[7]))
            high,second,flag = maxandsecond(a)
            highscore[read_name] = high
            secondscore[read_name] = second  
            
            high,second,flag = maxandsecond(b)
            highS[read_name] = high
            secondS[read_name] = second
            flagnum[read_name] = flag
        else:
            high,second,flag = maxandsecond(a)
            highscore[read_name] = high
            secondscore[read_name] = second
            high,second,flag = maxandsecond(b)
            highS[read_name] = high
            secondS[read_name] = second
            flagnum[read_name] = flag
            
            highscore[col[1]] = col[0]
            secondscore[col[1]] = col[0]
            highS[col[1]] = col[7]
            secondS[col[1]] = col[7]
            flagnum[col[1]] = 0
    else:
        if col[1] == read_name:
            a.append(float(col[0]))
            b.append(float(col[7]))
        else:
            high,second,flag = maxandsecond(a)
            highscore[read_name] = high
            secondscore[read_name] = second
            high,second,flag = maxandsecond(b)
            highS[read_name] = high
            secondS[read_name] = second
            flagnum[read_name] = flag
            
            a,b=[],[]
            read_name = col[1]
            a.append(float(col[0]))
            b.append(float(col[7]))
f.close()


read_name = ''
a,b=[],[]
maxTotalnum = {}
Totalnum = {}
judge_if={}
def stroe_highS_and_highscore(line):
    col = line.split('\t')
    if float(col[0]) == highscore[col[1]] and float(col[7].strip()) == highS[col[1]]:
        judge_if[col[1]]='1'
#max number of unique reads
f = open('new_score_all.txt')
lines = f.readlines()
last_line = lines[-1]
for line in lines:
    col = line.split('\t')
    if read_name.strip()=='':
        read_name = col[1]
        if highS[read_name]==float(col[7].strip()):
            a.append(int(col[6]))
            b.append(int(col[6]))
            stroe_highS_and_highscore(line)
        elif secondS[read_name]==float(col[7].strip()):
            b.append(int(col[6]))
        else:
            continue       
    elif line == last_line:
        if col[1] == read_name:
            if highS[read_name]==float(col[7].strip()):
                a.append(int(col[6]))
                b.append(int(col[6]))
                stroe_highS_and_highscore(line)
            elif secondS[read_name]==float(col[7].strip()):
                b.append(int(col[6]))
            maxTotalnum[read_name] = a
            Totalnum[read_name] = b
        else:
            maxTotalnum[read_name] = a
            Totalnum[read_name] = b
            if highS[read_name]==float(col[7].strip()):
                maxTotalnum[col[1]] = [int(col[6])]
                Totalnum[col[1]] = [int(col[6])]
                stroe_highS_and_highscore(line)
            elif secondS[read_name]==float(col[7].strip()):
                Totalnum[col[1]] = [int(col[6])]        
    else:
        if col[1] == read_name:
            if highS[read_name]==float(col[7].strip()):
                a.append(int(col[6]))
                b.append(int(col[6]))
                stroe_highS_and_highscore(line)
            elif secondS[read_name]==float(col[7].strip()):
                b.append(int(col[6]))
            else:
                continue
        else:
            maxTotalnum[read_name] = a
            Totalnum[read_name] = b
            a,b = [],[]
            read_name = col[1]
            if highS[read_name]==float(col[7].strip()):
                a.append(int(col[6]))
                b.append(int(col[6]))
                stroe_highS_and_highscore(line)
            elif secondS[read_name]==float(col[7].strip()):
                b.append(int(col[6]))
            else:
                continue
f.close()

def sam_judge(li):
    if 0 in li:
        a = sum(li)
        if a == 0:
            return 1
        else:
            return 0
    else:
        return 1
def maxlist(a):
    a.sort()
    return a[-1]
def value(a):
    a.sort()
    if a[-1]-a[0]>=5:
        return 1
    
tem_dict = {}       
f = open('new_score_all.txt')
lines = f.readlines()
for line in lines:
    col = line.split('\t')
    inf = col[3]+'\t'+col[4]+'\t'+col[5] #read_name & pos & len(seq)
    if col[1] in highscore.keys():
        name=col[1]       
        if (float(highS[name])-float(secondS[name]))>3:##
            if flagnum[name]==0:#
                if (highscore[name]-secondscore[name])>pow(10,-5) and float(col[7].strip()) == highS[name] and float(col[0]) == highscore[name]:
                    f=open('highscore_after_rescole.txt','a')
                    f.write(line)
                elif float(col[7].strip()) == highS[name] and int(col[6])!=0 and int(col[6])==int(maxlist(Totalnum[name])):
                    f=open('highscore_after_rescole.txt','a')
                    f.write(line)
                else:
                    continue
            else:
                if float(col[7].strip()) == highS[name] and sam_judge(maxTotalnum[name])==1:
                    if float(col[0]) == highscore[name]:
                        f=open('highscore_after_rescole.txt','a')
                        f.write(line)
                elif float(col[7].strip()) == highS[name]:
                    if int(col[6])!=0:
                        f=open('highscore_after_rescole.txt','a')
                        f.write(line)                     
        else:
            num_list=Totalnum[name]
            if sam_judge(num_list)==1:#
                if (highscore[name]-secondscore[name])>0.2:
                    if value(num_list)==1:#
                        if int(col[6])==int(maxlist(num_list)):
                            f=open('highscore_after_rescole.txt','a')
                            f.write(line)
                    else:
                        if (highscore[name]-secondscore[name])>1:
                            if float(col[0]) == highscore[col[1]]:
                                f=open('highscore_after_rescole.txt','a')
                                f.write(line)
                        else:
                            if int(col[6])==int(maxlist(num_list)):
                                f=open('highscore_after_rescole.txt','a')
                                f.write(line)
            else:
                if int(col[6])!=0:
                    f=open('highscore_after_rescole.txt','a')
                    f.write(line)
f.close()

dict_m = {}
dict_uni = {}
name = ''
i = 0
f = open('highscore_after_rescole.txt')
lines = f.readlines()
last_line = lines[-1]
for line in lines:
    col = line.split('\t')
    if name.strip() == '':
        name = col[1]
        i = i + 1
    elif line == last_line:
        if name == col[1]:
            i = i + 1
            dict_m[name] = i
        else:
            dict_m[name] = i
            dict_m[col[1]] = 1
    else:
        if name == col[1]:
            i = i + 1
        else:
            dict_m[name] = i
            i = 0
            name = col[1]
            i = i + 1
f.close()

for key in dict_m:
    if dict_m[key] == 1:
        dict_uni[key] = 1
del dict_m

#tem_dict = {}
#f = open('highscore_after_rescole.txt')
#lines = f.readlines()
#for line in lines:
 #   col = line.split('\t')
  #  if col[1] in dict_uni.keys():
   #     f1 = open('uni_highscore_after_rescore.txt','a')
    #    f1.write(line)        
#f.close()   

#f = open('uni_highscore_after_rescore.txt')
#lines = f.readlines()
#for line in lines:
 #   col = line.split('\t')
  #  total_name[col[1]] = col[4]
total_name={}
f = open('highscore_after_rescole.txt')
lines = f.readlines()
for line in lines:
    col = line.split('\t')
    if col[1] in dict_uni.keys():
        total_name[col[1]] = col[4]

for line in open(options.multifile):
    col = line.split('\t')
    if col[0] in total_name.keys() and int(col[3])==int(total_name[col[0]]):
        f=open('Read_uni_highscore.sam','a')
        f.write(line)
    elif col[0] not in total_name.keys():
        f=open('coverage.txt','a')
        f.write(line)

tem_dict = {}
f=open('coverage.txt')
lines = f.readlines()
for line in lines:
    col = line.split('\t')
    inf = col[2]+'\t'+col[3]+'\t'+str(len(col[9])) #read_name & pos & len(seq)
    name = col[0]
    if name in tem_dict.keys():
        tem_dict[name] = tem_dict[name]+'\t'+inf
    else:
        tem_dict[name] = inf
    
def calcul(ref_name, pos, lens):#no ss1,have ss2
    cover_cr = cover[ref_name]
    sum = 0
    ss1, ss2 = 0.0, 0.0
    if (int(pos)-1-WinLa)>=0 and (int(pos)+int(lens)-1+WinLa)<len(cover_cr):
        for i in range(int(pos)-1-WinLa,int(pos)+int(lens)-1+WinLa):
            sum = sum + cover_cr[i]
        cov1 = sum/(4*WinLa)
        cov2 = (sum+int(lens))/(4*WinLa)
        for i in range(int(pos)-1-WinLa,int(pos)+int(lens)-1+WinLa):
            ss1 = ss1 + (cover_cr[i]-cov1)**2
            if i in range(int(pos)-1,int(pos)+int(lens)-1):
                ss2 = ss2 + (cover_cr[i]+1-cov2)**2
            else:
                ss2 = ss2 + (cover_cr[i]-cov2)**2
        ss1 = ss1/(4*WinLa)
        ss2 = ss2/(4*WinLa)
        return (ss1,ss2)
    else:
        if (int(pos)-1-WinLa)<0:
            leftpos, rightpos=0, 4*WinLa
        if (int(pos)+int(lens)-1+WinLa)>=len(cover_cr):
            rightpos, leftpos=len(cover_cr)-1, len(cover_cr)-1-4*WinLa
        for i in range(leftpos,rightpos):
            sum = sum + cover_cr[i]
        cov1 = sum/(4*WinLa)
        cov2 = (sum+int(lens))/(4*WinLa)
        for i in range(leftpos,rightpos):
            ss1 = ss1 + (cover_cr[i]-cov1)**2
            if i in range(int(pos)-1,int(pos)+int(lens)-1):
                ss2 = ss2 + (cover_cr[i]+1-cov2)**2
            else:
                ss2 = ss2 + (cover_cr[i]-cov2)**2
        ss1 = ss1/(4*WinLa)
        ss2 = ss2/(4*WinLa)
        return (ss1,ss2)  

uni_dict={}
WinLa = 50
for i in tem_dict.keys():
    col = tem_dict[i].split('\t')
    number_mul = int(len(col)/3)
    names=locals()
    for j in range(0,number_mul):
        names['cover-'+str(j)],names['cover+'+str(j)]=calcul(col[0+j*3], int(col[1+j*3]), int(col[2+j*3]))
    for j in range(0,1):
        if number_mul == 2:
            if (names['cover+'+str(j)]+names['cover-'+str(j+1)]) < (names['cover-'+str(j)]+names['cover+'+str(j+1)]):
                uni_dict[i]=col[1]
            elif (names['cover+'+str(j)]+names['cover-'+str(j+1)]) > (names['cover-'+str(j)]+names['cover+'+str(j+1)]):
                uni_dict[i]=col[4]
            else:
                continue
        elif number_mul == 3:
            if (names['cover+'+str(j)]+names['cover-'+str(j+1)]) < (names['cover-'+str(j)]+names['cover+'+str(j+1)]):
                if (names['cover+'+str(j)]+names['cover-'+str(j+2)]) < (names['cover-'+str(j)]+names['cover+'+str(j+2)]):
                    uni_dict[i]=col[1]
                elif (names['cover+'+str(j)]+names['cover-'+str(j+2)]) > (names['cover-'+str(j)]+names['cover+'+str(j+2)]):
                    uni_dict[i]=col[7]
                else:
                    continue
            elif (names['cover+'+str(j)]+names['cover-'+str(j+1)]) > (names['cover-'+str(j)]+names['cover+'+str(j+1)]):
                if (names['cover+'+str(j+1)]+names['cover-'+str(j+2)]) < (names['cover-'+str(j+1)]+names['cover+'+str(j+2)]):
                    uni_dict[i]=col[4]
                elif (names['cover+'+str(j+1)]+names['cover-'+str(j+2)]) > (names['cover-'+str(j+1)]+names['cover+'+str(j+2)]):
                    uni_dict[i]=col[7]
                else:
                    continue                
            else:
                continue
        else:
            continue
    
for line in open('coverage.txt'):
    col = line.split('\t')
    read_name, pos = col[0], int(col[3])
    if read_name in uni_dict.keys():
        if pos == int(uni_dict[read_name]):
            f=open('Read_uni_coverage.sam','a')
            f.write(line)

