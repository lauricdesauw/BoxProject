## Kmer operators ## 

def reverse_complement(kmer)  :
    r = ""
    for i in range(len(kmer)) : 
        if kmer[i] == "A"  : 
            r = "T" + r
        elif kmer[i] =="T" :
            r = "A" +r
        elif kmer[i] == "C" :
            r = "G" + r
        elif kmer[i] == "G" :
            r = "C" + r
    return r

def get_canonical(kmer) : 
    if kmer < reverse_complement(kmer) : 
        return kmer 
    return reverse_complement(kmer)

def son(dbg,kmer): #Return a list of the following kmer that are in the dbg
    son = []
    part = kmer[1:]
    if(get_canonical(part + "A") in dbg) : 
        son.append(part + "A")

    if(get_canonical(part + "C") in dbg) : 
        son.append(part + "C")

    if( get_canonical(part + "T") in dbg) : 
        son.append(part + "T")

    if( get_canonical(part + "G") in dbg) : 
        son.append(part + "G")
    return son

## Dbg creation functions ## 

def add_kmer_dbg(dbg,kmer, freq) : #Add the canonical version of a k-mer to a dbg
    if (get_canonical(kmer) in freq.keys()) : 
        freq[get_canonical(kmer)] += 1
    else : 
        dbg.add( get_canonical(kmer))
        freq[get_canonical(kmer)] = 1

def add_kmer_from_a_read(dbg,read,k, freq) : #Add all the canonical representation of a read to a dbg
    for i in range(len(read) - k) : 
        add_kmer_dbg(dbg,read[i:i+k],freq) 

def create_dbg(read_file_name,k) :  #Create a dbg from a file full of read 
    read_file = open(read_file_name) 
    dbg = set()
    freq = {}
    while True : 
        line = read_file.readline()
        if not line : break 
        line = read_file.readline().strip()
        add_kmer_from_a_read(dbg,line,k,freq)
    print(len(dbg))
    return dbg,freq

def purge(dbg,freq) :  #We  keep only kmer that appear at least 4 times. Kmers that appears less than that may be errors of reading
    dbg_r = set()
    for kmer in dbg : 
        if freq[kmer] > 2 : 
            dbg_r.add(kmer)
    return dbg_r

## Utig functions ## 

def extend_right(dbg, kmer) :  #Extend a kmer on the right as long as the prolongement is unique
    curr = kmer
    r = ""
    while(True) : 
        k_son = son(dbg,curr)
        if (len(k_son) == 1) and (len(son(dbg,k_son[0])) == 1) : 
            r += k_son[0][-1]
            curr = k_son[0]
        else : 
            break;

    return r

def utig(dgb,kmer) : #Extend right and left
    return reverse_complement(extend_right(dgb,reverse_complement(kmer)) + kmer + extend_right(dbg,kmer))

def utig_son(dbg,ut,k) : #Same as the son function but with utig 
    kmer = ut[-k:-1]
    l = son(dbg,kmer) 
    for i in range(len(l)) : 
        l[i] = utig(l[i])
    return l

## Contig functions ##

def isbuble1(dbg,ut_list,k) : #Check if there is a buble of size one just after our utig
    visited  = set()
    common = True
    for ut in utig_son(dbg,ut_list[0],k) : 
        visited.add(ut)

    for ut in ut_list[1:] : 
        have_visited = True
        for utig in utig_son(dbg,ut,k) : 
            have_visited = have_visited and (utig in visited) 
        common = common and have_visited
    return common

def smooth_dbg(dbg) : #We are resolving some patern in the graph to make the deep parcour easier 
    smoothed_dbg = set()
    #TODO 
    return smoothed_dbg

def deep_parcour(dbg,source,target,k,visited) : #Find a way to the target from the source using a deep parcour of the graph
    sons_list = utig_son(dbg,source,k)
    print(sons_list)
    for next_ut in sons_list : 
        if next_ut == target : #We have reach the target
            return (True,source + target)
        elif next_ut in visited : #We are on an already visited node so he don't go to the target (because we are still running) 
            continue 
        else : 
            visited.add(next_ut)
            (a_way,path) = deep_parcour(dbg,next_ut,target,k,visited)
            if a_way : #There is a way to the target from that unitig
                return (true, source + path)
            continue

    return (False, "")

def contig(dbg,k,source,target) : #Generate contig from a dbg
    file = open("contig.txt", "w") 
    ut_source = utig(dbg,source)
    #dbg = smooth_dbg(dbg)
    visited = set()
    visited.add(ut_source)
    (found,ctg) = deep_parcour(dbg,ut_source,target,k,visited)
    file.write(ctg)
    file.write("\n")
    file.close()

## Program thread ## 
source = 'TTAAAATTTTATTGACTTAGGTCACTAAATA'
target = 'TCACCCCCGCACCATTACCCCCATCGCCCAG'

#Create the dbg from the read file
dbg,freq = create_dbg("Data/ecoli_2kb_perfect_reads.fasta",31)

#Purge the dbg to have only pertinent kmer
dbg = purge(dbg,freq)

print(len(dbg))
print(len(utig(dbg,source)))
print(utig_son(dbg,source,31))

#Generate the contig from the dbg (utig are hidden in the contig function)
contig(dbg,31, source,target)

