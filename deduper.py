import argparse
# Cameron Kunstadt
# Bi624, University of Oregon
# Oct-19-2024


parser = argparse.ArgumentParser()
parser.add_argument('-f', '--file')
parser.add_argument('-o', '--outfile')
parser.add_argument('-u', '--umi')
#parser.add_argument('-h', '--help')
args = parser.parse_args()

#TODO: Add help, go over rubric again

umi_set = {}
with open(args.umi) as umifile:
    umi_set = set(umifile.read().splitlines())

###########################################################
def get_line_info(line):
    '''Returns array of important info from provided SAM file line'''
    chrom = get_chrom(line)
    strand = get_strand(line)
    position = get_corrected_position(line, strand)
    umi = get_umi(line)
    return (chrom, position, strand, umi)

def get_chrom(line):
    '''Returns chromosome from provided SAM file line'''
    groups = line.split('\t')
    chrom = groups[2]
    return chrom

def get_position(line):
    '''Returns position from provided SAM file line'''
    groups = line.split('\t')
    position = groups[3]
    return position

def get_strand(line):
    '''Returns strand from provided SAM file line'''
    groups = line.split('\t')
    flag = int(groups[1])
    if((flag & 16) != 16):
        strand = '+'
    else:
        strand = '-'
    return strand

def get_umi(line):
    '''Returns umi from provided SAM file line'''
    groups = line.split('\t')
    qname = groups[0]
    qname_groups = qname.split(':')
    umi = qname_groups[7]
    return umi

def get_cigar(line):
    '''Returns CIGAR string from provided SAM file line'''
    groups = line.split('\t')
    cigar = groups[5]
    return cigar

def get_corrected_position(line, strand):
    '''Returns corrected position of SAM file line and strand'''
    position = int(get_position(line))
    cigar = get_cigar(line)
    if strand == '+':
        return (position - plus_softclip_adjustment(cigar))
    else:
        return (position - minus_softclip_adjustment(cigar))


def plus_softclip_adjustment(cigar):
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


def minus_softclip_adjustment(cigar):
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
            adjustment += int(clip_num)
            clip_num = ""
        elif letter.isnumeric() == False:
            clip_num = ""
        else:
            clip_num += letter
    return adjustment


def validate_umi(umi):
    '''Returns bool of if umi is in umi-set'''
    if umi in umi_set:
        return True
    else:
        return False
  
########################################################

with open(args.file, 'r') as infile, open(args.outfile, 'w') as outfile:
    dupset = set() # Set of read info to check if there are duplicates
    chr_num = 1
    while True:
        line = infile.readline()
        if line == "":
            break
        elif line[0] == "@": # Just write out the header regardless
            outfile.write(line)
        else:
            line_info = get_line_info(line) # Used to compare read duplication
            chr = int(line_info[0])
            umi = line_info[3]

            if chr != chr_num:
                dupset = set() # Wipe the dupset every new chromosome to save on memory
                chr_num += 1
            
            if line_info not in dupset and validate_umi(umi): # if read is not duplicate, write out
                dupset.add(line_info)
                outfile.write(line)


