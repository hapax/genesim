# genesim.py

# This program which simulates (X and autosomal) gene
# dropping through a pedigree, assuming no linkage.

# Map lengths (cM) and physical lengths (bp) of autosomes and X

physlengths = [247199719, 242751149, 199446827, 191263063, 180837866, 170896993, 158821424, 146274826, 140442298, 135374737, 134452384, 132289534, 114127980, 106360585, 100338915, 88822254, 78654742, 76117153, 63806651, 62435965, 46944323, 49528953, 154913754]
maplengths = [274.27,259.43,221.83,204.47,205.69,189.6,182.92,166.08,160.01,177.19,152.45,171.96,129.52,122.82,133.61,129.99,137.99,121.65,109.73,98.63,65.59,68.22,190.13]

import kinship as x
import numpy as n
import random as r

# Assuming there is no linkage, then crossover is
# spatial Poisson process. We simulate this by
# randomly sampling a Poisson rv whose parameter is the map
# length of the chromosome - this gives us the number of
# crossovers - then distribute these uniformly along interval.

def simcross(chrom_no):
    if chrom_no in range(1,24):
        k = n.random.poisson(maplengths[chrom_no - 1]/100)
        output = []
        for i in range(k):
            output.append(n.random.random())
        return sorted(output)
    else: print 'Bad input. Try again.'

# Get right labels for cross in meioses

def label(chrom, c):
    max_label = max(ints[0] for ints in chrom if ints[0] <= c) # Not sure this should be leq
    targ_int = [ints for ints in chrom if ints[0] == max_label]
    return targ_int[0][1]

def label_crs(chrom1, chrom2, cross):
    c1 = []
    c2 = []
    if cross:
        for k in range(len(cross)):
            if k % 2 == 0:
                c1.extend([[cross[k], label(chrom2, cross[k])]])
                c2.extend([[cross[k], label(chrom1, cross[k])]])
            if k % 2 == 1:
                c1.extend([[cross[k], label(chrom1, cross[k])]])
                c2.extend([[cross[k], label(chrom2, cross[k])]])
        return [sorted(c1), sorted(c2)]
    else: return cross

# Meiosis, chroms in format below

def meiosis(chrom1, chrom2, chrom_no, cross='r'):
    if cross == 'r': cross = simcross(chrom_no)
    elif cross:
        c1 = label_crs(chrom1,chrom2,cross)[0]
        c2 = label_crs(chrom1,chrom2,cross)[1]
        b11 = [ints for ints in chrom1 if ints[0] < cross[0]]
        b12 = []
        b21 = []
        b22 = [ints for ints in chrom2 if ints[0] < cross[0]]
        if len(cross) % 2 == 0:
            b11.extend([ints for ints in chrom1 if ints[0] > cross[-1]])
            b22.extend([ints for ints in chrom2 if ints[0] > cross[-1]])
            for k in range(len(cross)/2 - 2):
                b11.extend([ints for ints in chrom1 if cross[2*k+1] < ints[0] < cross[2*k+2]])
                b22.extend([ints for ints in chrom2 if cross[2*k+1] < ints[0] < cross[2*k+2]])
            for k in range(len(cross)/2 - 1):
                b12.extend([ints for ints in chrom1 if cross[2*k] < ints[0] < cross[2*k+1]])
                b21.extend([ints for ints in chrom2 if cross[2*k] < ints[0] < cross[2*k+1]])
            m1 = sorted(b11 + b21 + c1)
            m2 = sorted(b12 + b22 + c2)
            return [m1, m2]
        elif len(cross) % 2 == 1:
            b12.extend([ints for ints in chrom1 if ints[0] > cross[-1]])
            b21.extend([ints for ints in chrom2 if ints[0] > cross[-1]])
            for k in range(len(cross)/2 - 1):
                b11.extend([ints for ints in chrom1 if cross[2*k+1] < ints[0] < cross[2*k+2]])
                b12.extend([ints for ints in chrom1 if cross[2*k] < ints[0] < cross[2*k+1]])
                b21.extend([ints for ints in chrom2 if cross[2*k] < ints[0] < cross[2*k+1]])
                b22.extend([ints for ints in chrom2 if cross[2*k+1] < ints[0] < cross[2*k+2]])
            m1 = sorted(b11 + b21 + c1)
            m2 = sorted(b12 + b22 + c2)
        return [m1, m2]
    else: return [chrom1, chrom2]
                
# A chromosome is a list of genetic segments from founders.
# They are annoted by starting position (between 0 and 1)
# and the founder (and which chromosome from the founder).
# E.g., [[0, (1,1)], [0.5, (3, 2)]].

# Let's make a dictionary of chromosomes and update it
# as we go through pedigree.

def xdrop(pedigree):
    chroms = {}
    inds = x.inds(pedigree)
    founders = x.founders(pedigree)
    for f in founders:
        chroms[f] = {}
        if x.sex(pedigree, f) == 1:
            chroms[f][1] = 'y'
            chroms[f][2] = [[0,(f,2)]]
        else:
            chroms[f][1] = [[0,(f,1)]]
            chroms[f][2] = [[0,(f,2)]]
    while len(chroms) < len(inds):
        for i in inds:
            if i not in chroms.keys():
                chroms[i] = {}
                if x.mother(pedigree, i) in chroms.keys() and x.father(pedigree, i) in chroms.keys():
                    if x.sex(pedigree, i) == 1:
                        chroms[i][1] = 'y'
                        chroms[i][2] = r.choice(meiosis(chroms[x.mother(pedigree, i)][1], chroms[x.mother(pedigree, i)][2], 23, simcross(23)))
                    else:
                        chroms[i][1] = chroms[x.father(pedigree, i)][2]
                        chroms[i][2] = r.choice(meiosis(chroms[x.mother(pedigree, i)][1], chroms[x.mother(pedigree, i)][2], 23, simcross(23)))
    return chroms

def autodrop(pedigree, chrom_no):
    chroms = {}
    inds = x.inds(pedigree)
    founders = x.founders(pedigree)
    for f in founders:
        chroms[f] = {}
        chroms[f][1] = [[0,(f,1)]]
        chroms[f][2] = [[0,(f,2)]]
    while len(chroms) < len(inds):
        for i in inds:
            if i not in chroms.keys():
                chroms[i] = {}
                chroms[i][1] = r.choice(meiosis(chroms[x.father(pedigree, i)][1], chroms[x.father(pedigree, i)][2], chrom_no, simcross(chrom_no)))
                chroms[i][2] = r.choice(meiosis(chroms[x.mother(pedigree, i)][1], chroms[x.mother(pedigree, i)][2], chrom_no, simcross(chrom_no)))
    return chroms

def chrom_drop(ped, chrom_no):
    if chrom_no==23:
        return xdrop(ped)
    else: return autodrop(ped, chrom_no)

# Compare two individuals on x-chromosome - may become redundant

def emp_pihat(chroms, i, j, ped=0,chrom_no=23):
    if not chroms: chroms = chrom_drop(ped, chrom_no)
    c1 = chrom_to_list(chroms[i])
    c2 = chrom_to_list(chroms[j])
    return comp_chroms(c1,c2,2)+comp_chroms(c1,c2,1)/2

# Compare list of chroms

def ints(chrom):
    output = []
    for segs in chrom:
        output.append(segs[0])
    return output

def chrom_to_list(dic):
    lst = []
    if dic[1] == 'y':
        lst = [dic[2]]
    else: lst = [dic[1],dic[2]]
    return lst

def new_chrom(chrom_list1, chrom_list2):
    intvs = []
    for chrom in chrom_list1 + chrom_list2:
        intvs += ints(chrom)
    intvs = sorted(list(set(intvs)))
    new_chrom1 = [[i,set([])] for i in intvs]
    for i in range(len(intvs)):
        for chrom in chrom_list1+chrom_list2:
            new_chrom1[i][1].add(label(chrom, intvs[i]))
    return new_chrom1

def comp_chroms(chrom_list1, chrom_list2, ibd):
    new_chrom1 = new_chrom(chrom_list1, chrom_list2)
    if ibd == 0:
        nonshared = 0
        if len(chrom_list1) == len(chrom_list2) == 1:
            for k in range(len(new_chrom1)-1):
                if len(new_chrom1[k][1]) == 2:
                    nonshared += new_chrom1[k+1][0]-new_chrom1[k][0]
            if len(new_chrom1[-1][1]) == 2: nonshared += 1-new_chrom1[-1][0]
        elif len(chrom_list1) == 1 and len(chrom_list2)==2:
            for k in range(len(new_chrom1)-1):
                if len(new_chrom1[k][1]) == 2:
                    if label(chrom_list2[0],new_chrom1[k][0]) == label(chrom_list2[1],new_chrom1[k][0]):
                        nonshared += new_chrom1[k+1][0]-new_chrom1[k][0]
                if len(new_chrom1[k][1]) == 3: nonshared += new_chrom1[k+1][0]-new_chrom1[k][0]
            if len(new_chrom1[-1][1]) == 2:
                if label(chrom_list2[0],new_chrom1[-1][0]) == label(chrom_list2[1],new_chrom1[-1][0]):
                    nonshared += 1-new_chrom1[-1][0]
            if len(new_chrom1[-1][1]) == 3: nonshared += 1-new_chrom1[-1][0]
        elif len(chrom_list1) == 2 and len(chrom_list2)==1:
            for k in range(len(new_chrom1)-1):
                if len(new_chrom1[k][1]) == 2:
                    if label(chrom_list1[0],new_chrom1[k][0]) == label(chrom_list1[1],new_chro1m[k][0]):
                        nonshared += new_chrom1[k+1][0]-new_chrom1[k][0]
                if len(new_chrom1[k][1]) == 3: nonshared += new_chrom1[k+1][0]-new_chrom1[k][0]
            if len(new_chrom1[-1][1]) == 2:
                if label(chrom_list1[0],new_chrom1[-1][0]) == label(chrom_list1[1],new_chrom1[-1][0]):
                    nonshared += 1-new_chrom1[-1][0]
            if len(new_chrom1[-1][1]) == 3: nonshared += 1-new_chrom1[-1][0]
        elif len(chrom_list1) == len(chrom_list2) == 2:
            for k in range(len(new_chrom1)-1):
                if len(new_chrom1[k][1]) == 2:
                    if label(chrom_list1[0],new_chrom1[k][0]) == label(chrom_list1[1],new_chrom1[k][0]):
                        if label(chrom_list2[0],new_chrom1[k][0]) == label(chrom_list2[1],new_chrom1[k][0]):
                            nonshared += new_chrom1[k+1][0]-new_chrom1[k][0]
                elif len(new_chrom1[k][1]) == 3:
                    if (label(chrom_list1[0],new_chrom1[k][0]) == label(chrom_list1[1],new_chrom1[k][0])) != (label(chrom_list2[0],new_chrom1[k][0]) == label(chrom_list2[1],new_chrom1[k][0])):
                        nonshared += new_chrom1[k+1][0]-new_chrom[k][0]
                if len(new_chrom1[k][1]) == 4: nonshared += new_chrom1[k+1][0]-new_chrom1[k][0]
            if len(new_chrom1[-1][1]) == 2:
                if label(chrom_list1[0],new_chrom1[-1][0]) == label(chrom_list1[1],new_chrom1[-1][0]):
                    if label(chrom_list2[0],new_chrom1[-1][0]) == label(chrom_list2[1],new_chrom1[-1][0]):
                        nonshared += 1-new_chrom1[-1][0]
            if len(new_chrom1[-1][1]) == 3:
                if (label(chrom_list1[0],new_chrom1[-1][0]) == label(chrom_list1[1],new_chrom1[-1][0])) != (label(chrom_list2[0],new_chrom1[-1][0]) == label(chrom_list2[1],new_chrom1[-1][0])):
                    nonshared += 1-new_chrom1[-1][0]
            if len(new_chrom1[-1][1]) == 4: nonshared += 1-new_chrom1[-1][0]
        return nonshared
    elif ibd == 1:
        shared = 0
        if len(chrom_list1) == len(chrom_list2) == 1:
            for k in range(len(new_chrom1)-1):
                if len(new_chrom1[k][1]) == 1:
                    shared += new_chrom1[k+1][0]-new_chrom1[k][0]
            if len(new_chrom1[-1][1]) == 1: shared += 1-new_chrom1[-1][0]
        elif len(chrom_list1) == 1 and len(chrom_list2)==2:
            for k in range(len(new_chrom1)-1):
                if len(new_chrom1[k][1]) == 1: shared += new_chrom1[k+1][0]-new_chrom1[k][0]
                if len(new_chrom1[k][1]) == 2:
                    if label(chrom_list2[0],new_chrom1[k][0]) != label(chrom_list2[1],new_chrom1[k][0]):
                        shared += new_chrom1[k+1][0]-new_chrom1[k][0]
            if len(new_chrom1[-1][1]) == 1: shared += 1-new_chrom1[-1][0]
            if len(new_chrom1[-1][1]) == 2:
                if label(chrom_list2[0],new_chrom1[-1][0]) != label(chrom_list2[1],new_chrom1[-1][0]):
                    shared += 1-new_chrom1[-1][0]
        elif len(chrom_list1) == 2 and len(chrom_list2)==1:
            for k in range(len(new_chrom1)-1):
                if len(new_chrom1[k][1]) == 1: shared += new_chrom1[k+1][0]-new_chrom1[k][0]
                if len(new_chrom1[k][1]) == 2:
                    if label(chrom_list1[0],new_chrom1[k][0]) != label(chrom_list1[1],new_chrom1[k][0]):
                        shared += new_chrom1[k+1][0]-new_chrom1[k][0]
            if len(new_chrom1[-1][1]) == 1: shared += 1-new_chrom1[-1][0]
            if len(new_chrom1[-1][1]) == 2:
                if label(chrom_list1[0],new_chrom1[-1][0]) != label(chrom_list1[1],new_chrom1[-1][0]):
                    shared += 1-new_chrom1[-1][0]
        elif len(chrom_list1) == len(chrom_list2) == 2:
            for k in range(len(new_chrom1)-1):
                if len(new_chrom1[k][1]) == 2:
                    if label(chrom_list1[0],new_chrom1[k][0]) != label(chrom_list1[1],new_chrom1[k][0]):
                        if label(chrom_list2[0],new_chrom1[k][0]) != label(chrom_list2[1],new_chrom1[k][0]):
                            shared += new_chrom1[k+1][0]-new_chrom1[k][0]
                if len(new_chrom1[k][1]) == 3:
                    if label(chrom_list1[0],new_chrom1[k][0]) != label(chrom_list1[1],new_chrom1[k][0]) or label(chrom_list2[0],new_chrom1[k][0]) != label(chrom_list2[1],new_chrom1[k][0]):
                        shared += new_chrom1[k+1][0]-new_chrom1[k][0]
            if len(new_chrom1[-1][1]) == 2:
                if label(chrom_list1[0],new_chrom1[-1][0]) != label(chrom_list1[1],new_chrom1[-1][0]):
                    if label(chrom_list2[0],new_chrom1[-1][0]) != label(chrom_list2[1],new_chrom1[-1][0]):
                        shared += 1-new_chrom1[-1][0]
            if len(new_chrom1[-1][1]) == 3:
                if label(chrom_list1[0],new_chrom1[-1][0]) != label(chrom_list1[1],new_chrom1[-1][0]) or label(chrom_list2[0],new_chrom1[-1][0]) != label(chrom_list2[1],new_chrom1[-1][0]):
                    shared += 1-new_chrom1[-1][0]
        return shared
    elif ibd == 2:
        shared = 0
        if len(chrom_list1)*len(chrom_list2)<4:
            return shared
        for k in range(len(new_chrom1)-1):
            if len(new_chrom1[k][1]) == 1:
                shared += new_chrom1[k+1][0]-new_chrom1[k][0]
            if len(new_chrom1[k][1]) == 2:
                if label(chrom_list1[0],new_chrom1[k][0]) != label(chrom_list1[1],new_chrom1[k][0]):
                    if label(chrom_list2[0],new_chrom1[k][0]) != label(chrom_list2[1],new_chrom1[k][0]):
                        shared += new_chrom1[k+1][0]-new_chrom1[k][0]
        if len(new_chrom1[-1][1]) == 1: shared += 1-new_chrom1[-1][0]
        if len(new_chrom1[-1][1]) == 2:
                if label(chrom_list1[0],new_chrom1[-1][0]) != label(chrom_list1[1],new_chrom1[-1][0]):
                    if label(chrom_list2[0],new_chrom1[-1][0]) != label(chrom_list2[1],new_chrom1[-1][0]):
                        shared += 1-new_chrom1[-1][0]
        return shared

# Testing x ibd vs autosomal ibd content between two inds in a ped

def inner_prod(u, v):
    ip = 0
    for i in range(len(u)):
        ip += u[i]*v[i]
    return ip

def test(ped, ind_1, ind_2, no_reps):
    shared_content={}
    for i in range(1,24):
        shared_content[i] = []
    for i in range(no_reps):
        for j in range(1,24):
            shared_content[j].append(emp_pihat([], ind_1, ind_2, ped, j))
    for i in range(1,24):
        shared_content[i] = sum(shared_content[i])/len(shared_content[i])
    auto_content = inner_prod([v for (k,v) in shared_content.iteritems()][:22],maplengths[:22])
    return auto_content, shared_content[23]*maplengths[22]
