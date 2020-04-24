import sys
import collections
import time


def menu(): # menu for user
    print("Welcome to the DNA Sequencing Data Analysis Tool!\n")
    print("Select an algorithm for DNA sequencing")
    print("1. Knuth-Morris-Pratt Algorithm")
    print("2. Boyer-Moore")
    print("3. EXIT")
    
    choice = input()
    if choice == '1':
        mykmp()
    elif choice == '2':
        mybm()
    elif choice == '3':
        sys.exit
    else:
        print("\nEnter choice 1 or 2 or 3")
        menu()    
 


def mybm():    
    print('\nBoyer-Moore Algorithm')
    print("Enter the Gene squence to search in the human genome(length > 1): ")
    gene=input()
    genome = readGenome('human.fasta')      #read genome file
    count = collections.Counter()
    for read in genome:
        count.update(read)
    print('\nStatistics of bases in the human genome are : ')
    print(count)  
    p_bm = BoyerMoore(gene, alphabet='ACGT')  #table creation
    start=time.time()
    occurrences=boyer_moore(gene, p_bm, genome)  #call BM algorithm
    end = time.time()
    if occurrences:
        print('\nGene sequence  entered by you is present at following offsets in the Human Genome:')
        print(occurrences)
        print('\nTotal number of occurrences of the gene are: %d' % len(occurrences))
        
    else:
        print('\nNo match for the entered gene found in the human genome')
        
    print('\nTime taken : %f\n\n'%(end -start))
    
    menu()
    
def mykmp():    
    print('\nKnuth-Morris-Pratt Algorithm')
    print("Enter the Gene squence to search in the human genome(length > 1): ")
    gene=input()
    genome = readGenome('human.fasta')     #read genome file
    count = collections.Counter()
    for read in genome:
        count.update(read)
    print('\nStatistics of bases in the human genome are : ')
    print(count)   
    
    start=time.time()
    occurrences=KMP(gene, genome)   #call KMP algorithm
    end = time.time()
    
    if occurrences:
        print('\nGene sequence  entered by you is present at following offsets in the Human Genome:')
        print(occurrences)
        print('\nTotal number of occurrences of the gene are: %d' % len(occurrences))
    else:
        print('\nNo match for the entered gene found in the human genome')
        
    print('\nTime taken : %f\n\n'%(end -start))
    menu()
    
    
   
    
def readGenome(filename):
    genome = ''
    try:     
       f = open(filename, 'r')
    except IOError:
        print("\nFile not found")
        
    for line in f:
        # ignore header line with genome information
        if not line[0] == '>':
            genome += line.rstrip()
    return genome




def z_array(string):
    # Z algorithm for preprocessing as given in Dan Gusfield's book
    assert len(string) > 1
    z = [len(string)] + [0] * (len(string)-1)
    for i in range(1, len(string)):
        if string[i] == string[i-1]:
            z[1] += 1
        else:
            break
    right, left = 0, 0
    if z[1] > 0:
        right, left = z[1], 1
    for k in range(2, len(string)):
        assert z[k] == 0
        if k > right:
            # Case 1
            for i in range(k, len(string)):
                if string[i] == string[i-k]:
                    z[k] += 1
                else:
                    break
            right, left = k + z[k] - 1, k
        else:
            # Case 2
            beta = right - k + 1
            zkp = z[k - left]
            if beta > zkp:
                # Case 2a
                z[k] = zkp
            else:
                # Case 2b
                match = 0
                for i in range(right+1, len(string)):
                    if string[i] == string[i - k]:
                        match += 1
                    else:
                        break
                left, right = k, right + match
                z[k] = right - k + 1
    return z


def n_array(string):
    return z_array(string[::-1])[::-1]


def big_l_prime_array(p, n):    #prime L array compilation
    lpa = [0] * len(p)
    for j in range(len(p)-1):
        i = len(p) - n[j]
        if i < len(p):
            lpa[i] = j + 1
    return lpa


def big_l_array(p, lpa): #big L array compilation
    la = [0] * len(p)
    la[1] = lpa[1]
    for i in range(2, len(p)):
        la[i] = max(la[i-1], lpa[i])
    return la


def small_l_prime_array(n): #LP array compilation
    small_lpa = [0] * len(n)
    for i in range(len(n)):
        if n[i] == i+1:  # prefix matching a suffix
            small_lpa[len(n)-i-1] = i+1
    for i in range(len(n)-2, -1, -1):  # shift to the left
        if small_lpa[i] == 0:
            small_lpa[i] = small_lpa[i+1]
    return small_lpa


def good_suffix_table(p): # table for good suffux rule
    n = n_array(p)
    lp = big_l_prime_array(p, n)
    return lp, big_l_array(p, lp), small_l_prime_array(n)


def good_suffix_mismatch(i, big_l_prime, sml_l_prime):  #calculate amount of shifts
    length = len(big_l_prime)
    assert i < length
    if i == length - 1:
        return 0
    i += 1  # i points to leftmost matching position of P
    if big_l_prime[i] > 0:
        return length - big_l_prime[i]
    return length - sml_l_prime[i]


def good_suffix_match(sml_l_prime):  #return the match
    return len(sml_l_prime) - sml_l_prime[1]


def dense_bad_char_tab(p, amap): #calculate the dense bad character
    tab = []
    nxtc = [0] * len(amap)
    for i in range(0, len(p)):
        c = p[i]
        assert c in amap
        tab.append(nxtc[:])
        nxtc[amap[c]] = i+1
    return tab


class BoyerMoore(object): # class for Boyer-Moore algo
 
    def __init__(self, p, alphabet='ACGT'):
        self.p = p
        self.alphabet = alphabet
        # Create map from alphabet characters to integers
        self.amap = {}
        for i in range(len(self.alphabet)):
            self.amap[self.alphabet[i]] = i
        # Make bad character rule table
        self.bad_char = dense_bad_char_tab(p, self.amap)
        # Create good suffix rule table
        _, self.big_l, self.small_l_prime = good_suffix_table(p)
    
    def bad_character_rule(self, i, char):  #bad char rule
        assert char in self.amap
        ci = self.amap[char]
        assert i > (self.bad_char[i][ci]-1)
        return i - (self.bad_char[i][ci]-1)
    
    def good_suffix_rule(self, j):      #good suffix rule
        length = len(self.big_l)
        assert j < length
        if j == length - 1:
            return 0
        j += 1  # point to leftmost position
        if self.big_l[j] > 0:
            return length - self.big_l[j]
        return length - self.small_l_prime[j]
    
    def match_skip(self):
        return len(self.small_l_prime) - self.small_l_prime[1]
    
  
    
    
def boyer_moore(p, p_bm, t):   #Boyer-Moore function
    i = 0
    occurrences = []
    while i < len(t) - len(p) + 1:
        shift = 1
        mismatched = False
        for j in range(len(p)-1, -1, -1):
            if p[j] != t[i+j]:
                skip_bc = p_bm.bad_character_rule(j, t[i+j])
                skip_gs = p_bm.good_suffix_rule(j)
                shift = max(shift, skip_bc, skip_gs)
                mismatched = True
                break
        if not mismatched:
            occurrences.append(i)
            skip_gs = p_bm.match_skip()
            shift = max(shift, skip_gs)
        i += shift
    return occurrences    
 

def KMP(gene, genome): #main KMP function
    n = len(genome)    #get the length of genome and gene
    m = len(gene)
    occurrences = []   #for storing offset of occurances of a gene
    prefix_as_suffix = [0]*m
    computePrefix(gene, m, prefix_as_suffix)   #find the prefix
    i=0
    j=0
    while i < n:
        if genome[i] == gene[j]: #compare gene and genome base by base
            i += 1
            j += 1
        else:               
            if j != 0:
                j = prefix_as_suffix[j-1]
            else:
                i += 1
        if j == m:
            occurrences.append((i-j))    #add the offset of match to occurrences
            j = prefix_as_suffix[j-1] 
    return occurrences               

def computePrefix(gene, m, prefix_as_suffix):   #find out the longest prefix of gene which is also suffix in genome 
    len = 0
    i = 1
    prefix_as_suffix[0] = 0
    while i < m:
        if gene[i] == gene[len]:
            prefix_as_suffix[i] = len + 1
            len += 1
            i += 1
        else:
            if len != 0:
                len = prefix_as_suffix[len-1]
            else:
                prefix_as_suffix[i] = 0
                i += 1
          

                       
    
    
if __name__=='__main__':
	menu()
