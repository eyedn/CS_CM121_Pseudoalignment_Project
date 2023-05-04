###############################################################################
#       CS CM121
#       Winter 2023
#       Project 2
#       Aydin Karayas
###############################################################################
import gzip
from sys import *

###############################################################################
#       DEFINE GLOBAL VARIABLE(S)
###############################################################################
chr_11_raw = open(argv[1], "r")
rna_seq_data = gzip.open(argv[2], "rt")
KMER_LENGTH = int(argv[3])

###############################################################################
#       FUNCTIONS
###############################################################################
# create a dictionary that pairs up complimentary bases
def create_base_pairs(): 
    base_pairs = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    return base_pairs

# create a 'transcipt':'seqeunce' dictionary given the raw transcriptome file
def create_transcriptome_dict(raw_transcriptome, bases_pairs):
    my_transcriptome = {}
    my_transcript = ""
    for my_line in raw_transcriptome:
        curr_line = my_line.strip()
        # ignore empty lines because they contain no useful information
        if len(curr_line) == 0:
            continue
        # lines starting with ">" identify transcipt names
        if curr_line[0] == ">":
            # every transcipt name is initialized as a key with
            my_transcript = curr_line[1:]
            my_transcriptome[my_transcript] = ["",""]
            continue
        # given a transcript: forward strand is [0], reverse compliment is [1]
        my_transcriptome[my_transcript][0] += curr_line.upper()
        my_transcriptome[my_transcript][1] += rev_com(curr_line.upper(), 
                                                      bases_pairs)
    return my_transcriptome

# create the reverse complitment of a given sequence
def rev_com(sequence, bases_pairs):
    compliment = ""
    for base in sequence:
        compliment += bases_pairs[base]
    # reverse the order of the complement to get the reverse complement
    reverse_compliment = compliment[::-1]
    return reverse_compliment

# create a 'kmer':'transcripts' dict given the transcriptome dict and kmer size
def create_kmer_dict(transcriptome, k):
    my_equiv_classes = {}
    for my_transcript in transcriptome:
        # generate kmers for the forward [0] and reverse compliment [1] strands
        fw_kmers = kmer_ize_sequence(transcriptome[my_transcript][0], k)
        rc_kmers = kmer_ize_sequence(transcriptome[my_transcript][1], k)
        # add kmers of this transcript to the 'kmer':'transcripts' dictionary
        my_equiv_classes = add_kmers(my_equiv_classes, fw_kmers, my_transcript)
        my_equiv_classes = add_kmers(my_equiv_classes, rc_kmers, my_transcript)
    return my_equiv_classes

# divide a given sequence into kmers of size k
def kmer_ize_sequence(sequence, k):
    kmers = []
    for i in range(len(sequence) - k + 1):
        kmers.append(sequence[i:i+k])
    return(kmers)

# add kmers from agiven transcript to a 'kmer':'transcripts' dictionary
def add_kmers(equiv_classes, kmers, transcript):
    for my_kmer in kmers:
        # new kmers are initalized as values; transcipts are stored in a set
        if my_kmer not in equiv_classes:
            equiv_classes[my_kmer] = set()
        # reason for using set: a transcipt appears a max of once per kmer key
        equiv_classes[my_kmer].add(transcript)
    return equiv_classes

# generate equivalence class counts for rna-seq reads
def generate_equiv_class_counts(rna_seq, transcriptome_kmers, k):
    my_output_table = []
    # my_label is the "read + true transcipt" given by the reads file
    my_label = ""
    for my_line in rna_seq:
        curr_line = my_line.strip()
        if len(curr_line) == 0:
            continue
        # lines that start with ">" identify read labels
        # the line also stores other information after a ";" character
        if curr_line[0] == ">":
            # we only want the read number and true transcript for this scipt
            my_label = curr_line[1:].split(";", 1)[0]
            continue
        # map a given read to the kmer dictionary; output is the equiv class
        my_equiv_class = map_read(transcriptome_kmers, curr_line.upper(), k)
        # format output csv: each row is a read from rna_seq with these columns
        #   string: label 
        #   int: number of transcipts in equiv class
        #   string: transcripts of equiv class
        #   bool: was the true transcipt found in equiv class (for report)
        my_ec_strings = f'"{",".join(sorted(my_equiv_class))}"'
        my_ec_len = str(len(my_equiv_class))
        correct = str(check_read_correctness(my_label, my_equiv_class))
        my_output_table.append([my_label, my_ec_len, my_ec_strings, correct])
    printable_results = [",".join(sublist) for sublist in my_output_table]
    return printable_results

# map read to kmer dict
def map_read(transcriptome_kmers, read, k):
    # first, generate the kmers for this read
    my_read_kmers = kmer_ize_sequence(read, k)
    # for every kmer found in the read, this will store its transcripts 
    my_kmer_tps = []
    for kmer in my_read_kmers:
        # the try block will run if the kmer is found in the transcriptome
        try:
            possible_transcripts = transcriptome_kmers[kmer]
            my_kmer_tps.append(possible_transcripts)
        # if the kmer is not found, it will be added as a null set
        except:
            null_set = set()
            my_kmer_tps.append(null_set)
    # to find the equivalence class, take the intersect of kmer transcripts
    # backtrack_errors does this and deals with sequence errors/null sets
    my_equiv_class, my_skips = backtrack_errors(my_kmer_tps, k)
    return my_equiv_class  

# takes the intersect of transcripts from kmers (deals with sequence errors)
def backtrack_errors(kmer_tps, k, skip = 0):
    # base case: if there are no kmer_tps the must read starts with a null set
    #            skip the first 31 kmers to final equiv class
    if len(kmer_tps) == 0:
        # NB: at this point skip == k + 1; see *** for why it was set to k + 1
        #     set skips left to the value of k to skip the first 31 kmers
        skips_left = skip - 1
        return set(), skip - 1
    # try to take the intersect of all kmer-specific transcipts 
    my_equiv_class = set.intersection(*kmer_tps)
    # a null set will recursively calls the function, dropping the last kmer
    if my_equiv_class == set():
        # *** NB: although we will be skipping k number of kmers at once a 
        #         valid equiv class is found, the skip counter will be initally
        #         set to k + 1. This is because the skip_errors function will
        #         begin decrementing the skip counter at the first instance of 
        #         a valid equiv class. We set skips to k + 1 ensuring all 
        #         possibly erroneous kmers are skipped.
        my_equiv_class, skips_left = backtrack_errors(kmer_tps[:-1], k, k + 1)
    else:
        skips_left = 0
    # once the shortest possible my_equiv_class is established, 
    # try to add more kmer-specific equiv classes by skipping errors
    # repeat this process until the end of the read
    my_equiv_class, skips_left = skip_errors(kmer_tps, my_equiv_class, 
                                             k, skips_left)
    return my_equiv_class, skips_left

# given a prelimiary equiv class, add more kmers by skipping erroneous kmers
def skip_errors(depth_specific_kmer_tps, curr_equiv_class, k, skips):
    # if skips == 0, we can try adding more kmers again
    if skips == 0:
        # current equiv class being null can arise if the first kmer was null
        # transcipts of the first kmer post-skips will be the new equiv class
        if curr_equiv_class == set():
            skip_intersect = set.intersection(depth_specific_kmer_tps[-1])
        else:
            # take intersect of current equiv class + current possibly valid kmer
            skip_intersect = set.intersection(curr_equiv_class, 
                                              depth_specific_kmer_tps[-1])
        # null sets indicate another error after skipping; restart skipping
        if skip_intersect == set():
            skips = k
        # the new kmer is added because its addition didn't cause a null set
        else:
            curr_equiv_class = skip_intersect
    # if told to skip, decrement skips and keep curr_equv_class as is
    else:
        skips -= skips
    return curr_equiv_class, skips

# when creating output data, also output if psuedoalignment was successful
def check_read_correctness(read_label, equiv_classes):
    # read_label looks like this: "read_number/true_transcipt"
    # extract the true transcipt and check if it's in the equiv class
    true_transcript = read_label.split('/', 1)[1]
    return true_transcript in equiv_classes

###############################################################################
#       MAIN CODE
###############################################################################
# define: A with T; G with C
DNA_base_pairs = create_base_pairs()

# reformating transcriptome data so all transcript sequences are 1 string
chr_11_transcriptome = create_transcriptome_dict(chr_11_raw, DNA_base_pairs)

# for each key (kmer of length k), we can find all transcripts it occurs in
chr_11_kmers = create_kmer_dict(chr_11_transcriptome, KMER_LENGTH)

# find all kmers exist each read and take the intersect of relevant transcripts
reads_against_chr11 = generate_equiv_class_counts(rna_seq_data, chr_11_kmers, KMER_LENGTH)

# print results
for line in reads_against_chr11:
    print(line)