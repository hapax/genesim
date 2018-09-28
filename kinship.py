# kinship.py

# First, quick duplication remover:

def rem_dup(seq): 
    keys = {} 
    for e in seq: 
        keys[e] = 1 
    return keys.keys()

# We need a way of determining sex and parents of individuals.

def sex(pedigree, i):
    fs = open(pedigree, 'rb')
    text = fs.read().split('\n')
    for k in range(len(text)):
        if text[k]:
            lst = text[k].split()
            if lst[0] == str(i): return eval(lst[3])
    print "No individual with that id."

def father(pedigree, i):
    fs = open(pedigree, 'rb')
    text = fs.read().split('\n')
    for k in range(len(text)):
        if text[k]:
            lst = text[k].split()
            if lst[0] == str(i): return eval(lst[1])
    print "No individual with that id."

def mother(pedigree, i):
    fs = open(pedigree, 'rb')
    text = fs.read().split('\n')
    for k in range(len(text)):
        if text[k]:
            lst = text[k].split()
            if lst[0] == str(i): return eval(lst[2])
    print "No individual with that id."

def f_check(pedigree, i):
    return father(pedigree, i) == 0

def inds(pedigree):
    inds = []
    fs = open(pedigree, 'rb')
    text = fs.read().split('\n')
    for k in range(len(text)):
        if text[k]:
            lst = text[k].split()
            inds.append(eval(lst[0]))
    return inds

def founders(pedigree):
    return [x for x in inds(pedigree) if f_check(pedigree, x)]

# Easy stuff first. Relate X-linked IBD vector to Deltas, and Deltas to Psis.

def psi(pedigree, n, i, j):
    sex_i = sex(pedigree, i)
    sex_j = sex(pedigree, j)
    n = int(n)
    if [sex_i, sex_j] == [2,2]:
        if n == 1: return phi(pedigree, [[[i, 1], [i, 2], [j, 1], [j, 2]]])
        elif n == 2: return phi(pedigree, [[[i, 1], [i, 2]], [[j, 1], [j, 2]]])
        elif n == 3: return 2*phi(pedigree, [[[i, 1], [i, 2], [j, 1]], [[j, 2]]])
        elif n == 4: return phi(pedigree, [[[i, 1], [i, 2]], [[j, 1]], [[j, 2]]])
        elif n == 5: return 2*phi(pedigree, [[[i, 1], [j, 1], [j, 2]], [[i, 2]]])
        elif n == 6: return phi(pedigree, [[[j, 1], [j, 2]], [[i, 1]], [[i, 2]]])
        elif n == 7: return 2*phi(pedigree, [[[i, 1], [j, 1]], [[i, 2], [j, 2]]])
        elif n == 8: return 4*phi(pedigree, [[[i, 1], [j, 1]], [[i, 2]], [[j, 2]]])
        elif n == 9: return phi(pedigree, [[[i, 1]], [[j, 1]], [[i, 2]], [[j, 2]]])
        else: print "Please select valid n."
    elif [sex_i, sex_j] == [2,1]:
        if n == 1: return phi(pedigree, [[[i, 1], [i, 2], [j, 1]]])
        elif n == 2: return phi(pedigree, [[[i, 1], [i, 2]], [[j, 1]]])
        elif n == 3: return 2*phi(pedigree, [[[i, 1], [j, 1]], [[i, 2]]])
        elif n == 4: return phi(pedigree, [[[i, 1]], [[i, 2]], [[j, 1]]])
        else: print "Please select valid n."
    elif [sex_i, sex_j] == [1,2]:
        if n == 1: return phi(pedigree, [[[j, 1], [j, 2], [i, 1]]])
        elif n == 2: return phi(pedigree, [[[j, 1], [j, 2]], [[i, 1]]])
        elif n == 3: return 2*phi(pedigree, [[[j, 1], [i, 1]], [[j, 2]]])
        elif n == 4: return phi(pedigree, [[[j, 1]], [[j, 2]], [[i, 1]]])
        else: print "Please select valid n."
    elif [sex_i, sex_j] == [1,1]:
        if n == 1: return phi(pedigree, [[[i,1],[j,1]]])
        elif n == 2: return phi(pedigree, [[[i,1]],[[j,1]]])

def delta(pedigree, n, i, j):
    sex_i = sex(pedigree, i)
    sex_j = sex(pedigree, j)
    n = int(n)
    if [sex_i, sex_j] == [2,2]:
        if n == 1: return psi(pedigree, 1, i, j) - psi(pedigree, 3, i, j)/2 - psi(pedigree, 5, i, j)/2 + psi(pedigree, 7, i, j)/2 + psi(pedigree, 8, i, j)/4
        elif n == 2: return psi(pedigree, 2, i, j) - psi(pedigree, 3, i, j)/2 - psi(pedigree, 4, i, j) - psi(pedigree, 5, i, j)/2 - psi(pedigree, 6, i, j) + psi(pedigree, 7, i, j)/2 + 3*psi(pedigree, 8, i, j)/4 + psi(pedigree, 9, i, j)
        elif n == 3: return 2*psi(pedigree, 3, i, j) - 2*psi(pedigree, 7, i, j) - psi(pedigree, 8, i, j)
        elif n == 4: return 2*psi(pedigree, 4, i, j) - psi(pedigree, 8, i, j) - 2*psi(pedigree, 9, i, j)
        elif n == 5: return 2*psi(pedigree, 5, i, j) - 2*psi(pedigree, 7, i, j) - psi(pedigree, 8, i, j)
        elif n == 6: return 2*psi(pedigree, 6, i, j) - psi(pedigree, 8, i, j) - 2*psi(pedigree, 9, i, j)
        elif 0 < n < 10: return 4*psi(pedigree, n, i, j)
        else: print "Please pick a valid n."
    if [sex_i, sex_j] == [1,2] or [sex_i, sex_j] == [2, 1]:
        if n == 1: return psi(pedigree, 1, i, j) - psi(pedigree, 3, i, j)/2
        elif n == 2: return psi(pedigree, 2, i, j) - psi(pedigree, 3, i, j)/2 - psi(pedigree, 4, i, j)
        elif 0 < n < 5: return 2*psi(pedigree, n, i, j)
        else: print "Please pick a valid n."
    if [sex_i, sex_j] == [1,1]:
        if n == 1: return psi(pedigree, 1, i, j)
        elif n == 2: return psi(pedigree, 2, i, j)
        else: print "Please pick a valid n."

def inbreed(pedigree, i):
    if sex(pedigree, i) == 2: return delta(pedigree, 1, i, i)
    else: return 1

def xibdv(pedigree, i, j):
    sex_i = sex(pedigree, i)
    sex_j = sex(pedigree, j)
    if [sex_i, sex_j] == [2,2]:
        z0 = delta(pedigree, 2, i, j) + delta(pedigree, 4, i, j) + delta(pedigree, 6, i, j) + delta(pedigree, 9, i, j)
        z1 = delta(pedigree, 3, i, j) + delta(pedigree, 5, i, j) + delta(pedigree, 8, i, j)
        z2 = delta(pedigree, 1, i, j) + delta(pedigree, 7, i, j)
    elif [sex_i, sex_j] == [1,2] or [sex_i, sex_j] == [2, 1]:
        z0 = delta(pedigree, 2, i, j) + delta(pedigree, 4, i, j)
        z1 = delta(pedigree, 1, i, j) + delta(pedigree, 3, i, j)
        z2 = 0
    elif [sex_i, sex_j] == [1,1]:
        z0 = delta(pedigree, 2, i, j)
        z1 = delta(pedigree, 1, i, j)
        z2 = 0
    return [z0, z1, z2]

def xpihat(pedigree, i, j):
    if sex(pedigree, i)*sex(pedigree,j)!= 4:
        return xibdv(pedigree, i, j)[1]
    else: return xibdv(pedigree, i, j)[1]/2 + xibdv(pedigree, i, j)[2]

# OK, so here's the tricky part: calculating Phi.

# Let's first define a function which gets block inds with multiplicity

def block_inds(blocks):
    return [genes[0] for genes in blocks]

import copy

def phi(pedigree, partition):
    # The point here is to list the individuals involved, with multiplicity
    partition = copy.deepcopy(partition)
    sampled = []
    for blocks in partition:
        sampled.extend(block_inds(blocks))
    # B1
    for i in rem_dup(sampled):
        if sum(i in block_inds(blocks) for blocks in partition) > sex(pedigree, i): return 0
    # B2
    for blocks in partition:
        block_fs = [i for i in block_inds(blocks) if f_check(pedigree, i)]
        if len(rem_dup(block_fs)) > 1: return 0
    # B3
    if all(f_check(pedigree, i) for i in sampled):
        no_f_sampled = sum(sex(pedigree, i) == 2 for i in rem_dup(sampled))
        no_f_genes = sum(sex(pedigree, i) == 2 for i in sampled)
        return (1.0/2.0)**(no_f_genes - no_f_sampled)
    # R1
    for i in rem_dup(sampled):
        if not f_check(pedigree, i):
            if sum(i in block_inds(blocks) for blocks in partition) == 1:
                b_index = [i in block_inds(blocks) for blocks in partition].index(1)
                if sex(pedigree, i) == 1:
                    for genes in partition[b_index]:
                        if genes[0] == i:
                            partition[b_index].remove(genes)
                    for genes in partition[b_index]: # Try get rid of duplication
                        if genes[0] == i:            # at some point
                            partition[b_index].remove(genes)
                    partition1 = copy.deepcopy(partition)
                    partition1[b_index].append([mother(pedigree, i), 1])
                    return phi(pedigree, partition1)
                if sex(pedigree, i) == 2:
                    s = sum(genes[0] == i for genes in partition[b_index])
                    partition1 = copy.deepcopy(partition)
                    partition2 = copy.deepcopy(partition)
                    partition3 = copy.deepcopy(partition)
                    for genes in partition[b_index]:           #
                        if genes[0] == i:                      #
                            partition1[b_index].remove(genes)  #
                            partition2[b_index].remove(genes)  #
                            partition3[b_index].remove(genes)  #
                    for genes in partition1[b_index]:
                        if genes[0] == i:
                            partition1[b_index].remove(genes)
                    for genes in partition2[b_index]:
                        if genes[0] == i:
                            partition2[b_index].remove(genes)
                    for genes in partition3[b_index]:
                        if genes[0] == i:
                            partition3[b_index].remove(genes)
                    partition1[b_index].append([mother(pedigree, i), 1])
                    partition2[b_index].append([father(pedigree, i), 1])
                    partition3[b_index].extend([[mother(pedigree, i), 1], [father(pedigree, i), 1]])
                    return ((1.0/2.0)**s)*phi(pedigree, partition1) + ((1.0/2.0)**s)*phi(pedigree, partition2) + (1 - 2*((1.0/2.0)**s))*phi(pedigree, partition3)
    # R2
    for i in rem_dup(sampled):
        if not f_check(pedigree, i):
            if sum(i in block_inds(blocks) for blocks in partition) == 2:
                b_index1 = [i in block_inds(blocks) for blocks in partition].index(1)
                b_index2 = [i in block_inds(blocks) for blocks in partition[b_index1 + 1:]].index(1) + b_index1 + 1
                s = sum(genes[0] == i for genes in partition[b_index1])
                t = sum(genes[0] == i for genes in partition[b_index2])
                for genes in partition[b_index1]:
                    if genes[0] == i:
                        partition[b_index1].remove(genes)
                for genes in partition[b_index2]:
                    if genes[0] == i:
                        partition[b_index2].remove(genes)
                for genes in partition[b_index1]:           # WEIRD DUPLICATION
                    if genes[0] == i:                       # Try to get rid of this
                        partition[b_index1].remove(genes)   # at some point
                for genes in partition[b_index2]:           #
                    if genes[0] == i:                       #
                        partition[b_index2].remove(genes)   #
                partition1 = copy.deepcopy(partition)
                partition1[b_index1].append([mother(pedigree, i), 1])
                partition1[b_index2].append([father(pedigree, i), 1])
                partition2 = copy.deepcopy(partition)
                partition2[b_index1].append([father(pedigree, i), 1])
                partition2[b_index2].append([mother(pedigree, i), 1])
                return ((1.0/2.0)**(s + t))*phi(pedigree, partition1) + ((1.0/2.0)**(s + t))*phi(pedigree, partition2)
