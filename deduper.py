import argparse
# Cameron Kunstadt
# Bi624, University of Oregon
# Oct-19-2024


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file', help="Input must be a sam file sorted by chromosome")
parser.add_argument('-o', '--outfile')
parser.add_argument('-u', '--umi')
#parser.add_argument('-h', '--help')
args = parser.parse_args()

#TODO: Add help, test strand softclip adjustment, add RC UMI set

umi_set = {}
with open(args.umi) as umifile:
    umi_set = set(umifile.read().splitlines())

###########################################################

def get_line_info(line):
    '''Returns tuple of important info from provided SAM file line to use in duplicate removal. Position
       is softclip-corrected'''
    chrom: int =  get_chrom(line)
    strand: str = get_strand(line)
    position: int = get_corrected_position(line, strand)
    umi: str = get_umi(line)
    return (chrom, position, strand, umi)

def get_chrom(line):
    '''Returns chromosome from provided SAM file line'''
    groups = line.split('\t')
    chrom = groups[2]
    return int(chrom)

def get_strand(line):
    '''Returns strand from provided SAM file line'''
    groups = line.split('\t')
    flag = int(groups[1])
    if((flag & 16) != 16):
        strand = '+'
    else:
        strand = '-'
    return strand

def get_corrected_position(line, strand):
    '''Returns corrected position of SAM file line and strand'''
    groups = line.split('\t')
    position = int(groups[3])
    cigar =  groups[5]
    if strand == '+':
        return (position - plus_strand_softclip_adjustment(cigar))
    else:
        return (position + minus_strand_softclip_adjustment(cigar))

def plus_strand_softclip_adjustment(cigar):
    '''Calculates the needed softclip adjustment from cigar string, for + strand'''
    clip_num = ""
    for letter in cigar:
        if letter == 'S':
            return int(clip_num)
        elif letter == 'M':
            return 0
        else:
            clip_num += letter
    return int(clip_num)

def minus_strand_softclip_adjustment(cigar):
    '''Calculates the needed softclip adjustment from cigar string, for - strand'''
    clip_num = ""
    adjustment = 0
    past_first_clip = False

    for letter in cigar:
        if letter == 'S' and not past_first_clip:
            past_first_clip = True
        elif letter == 'S' and past_first_clip:
            adjustment += int(clip_num)
            clip_num = ""
        elif letter == 'M' or letter == 'D':
            past_first_clip = True
            adjustment += int(clip_num)
            clip_num = ""
        elif letter.isnumeric() == False:
            clip_num = ""
        else:
            clip_num += letter
    return adjustment

def get_umi(line):
    '''Returns umi from provided SAM file line'''
    groups = line.split('\t')
    qname = groups[0]
    qname_groups = qname.split(':')
    umi = qname_groups[7]
    return umi

def validate_umi(umi):
    '''Returns bool of if umi is in umi-set'''
    if umi in umi_set:
        return True
    else:
        return False
  
########################################################

with open(args.file, 'r') as infile, open(args.outfile, 'w') as outfile:
    dupset = set() # Set of unique read-info tuples to check if there are duplicates
    chr_num = 1
    while True:
        line = infile.readline()
        if line == "": # Break if EOF
            break
        elif line[0] == "@": # Write out all header lines regardless
            outfile.write(line)
        else:
            line_info = get_line_info(line) # line_info is a tuple used to compare read duplication
            chr = line_info[0]
            umi = line_info[3]

            if chr != chr_num:
                dupset = set() # Wipe the dupset every new chromosome to save on memory
                chr_num += 1
            
            if line_info not in dupset and validate_umi(umi): # if read is not duplicate, write out
                dupset.add(line_info)
                outfile.write(line)


