__author__ = 'michaela'
import argparse
import csv
import os
from collections import Counter


parser = argparse.ArgumentParser(
    description='SNP type your genomes into introductions')
parser.add_argument('-f', dest='file',
                    help='text file with Gnumbers in one column. Can be the same input to the pipeline',
                    required=True)
parser.add_argument('-o', dest='output', help='output file requested')

parser.add_argument('-v', dest='version_pipeline', help='Whether the output of v1 or v2 should be used')
parser.add_argument('-m', dest='marker', help='which markers you want to use') # corresponds to Supplementary Table 2


args = parser.parse_args()

class AutoVivification(dict):
    """Implementation of autovivification in Python"""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

Shitikov_dic = AutoVivification()

with open(args.marker,"rU") as tsvfile:
    tsvreader = csv.DictReader(tsvfile, delimiter="\t")
    for row in tsvreader:
        Shitikov_dic[row['Import']][row['Position_ref']] = row['derived']

GNumber = []
outputfile = open(args.output,'w')
with open(args.file) as input_file:
    for line in input_file:
        line = line.strip()
        line = line.split()
        GNumber += [line[0]] # line[0] is column 1

for gnb in GNumber:
    #print gnb
    if args.version_pipeline == "v1":
        dir = '~/PipelineTB/v1/'
        if len(gnb) == 6:
            vcf = dir+gnb[0:3]+'/'+gnb[3:5]+'/'+gnb[5]+'/'+gnb+'.var.snp.vcf'
        else:
            vcf = dir+gnb[0:3]+'/'+gnb[3:5]+'/'+gnb[5:7]+'/'+gnb+'.var.snp.vcf'
    else:
        dir = '~/PipelineTB/v2/'

        if len(gnb) == 6:
            vcf = dir+gnb[0:3]+'/'+gnb[3:5]+'/'+gnb[5]+'/'+gnb+'.mutect2.filtered.allvariants.vcf'
        else:
            vcf = dir+gnb[0:3]+'/'+gnb[3:5]+'/'+gnb[5:7]+'/'+gnb+'.mutect2.filtered.allvariants.vcf'

    LIST_VCF = []
    if os.path.exists(vcf):

            with open(vcf, 'r') as VCF:  # open the VCF

                for row in VCF:  # Loop into the lines of the VCF
                    warning = False
                    if row[0] != '#':  # Ignore the header of the VCF

                        row = row.strip()
                        row = row.split()

                        POS = str(row[1])  # position is the second column of the VCF
                        ALT = str(row[4])  # ALT is the fifth column of the VCF
                        # sometimes there are two alt kept, since the het one was not removed (how annoying...)
                        # only add the position if the allele frequency is above 0.7
                        if args.version_pipeline == "v1":
                            if len(ALT.split(",")) > 1: # find the alternative with the higher allele frequency
                                FREQS = row[9].split(":")[6].split(",")
                                FREQS = [float(item.split("%")[0]) for item in FREQS]
                                highest = max((v,z) for z,v in enumerate(FREQS))[1]
                                ALT = ALT.split(",")[highest]
                                FREQ = FREQS[highest]
                            else:
                               FREQ = float(row[9].split(":")[6].split("%")[0])/100

                        elif args.version_pipeline == "v2":
                            if len(ALT.split(",")) > 1: # find the alternative with the higher allele frequency

                                FREQS = row[9].split(":")[2].split(",")
                                highest = max((v,z) for z,v in enumerate(FREQS))[1]
                                ALT = ALT.split(",")[highest]
                                FREQ = FREQS[highest]
                            else:
                                FREQ = row[9].split(":")[2]
                        try:
                            if float(FREQ) > 0.7:
                                LIST_VCF += [(POS,ALT)]  # Define a tuple variable which for each line of the VCF takes the value of (position,alternative base)
                        except ValueError:
                            print(gnb + POS)

            sublineage_shitikov = {}
            sublineage_shitikov_list = []
            sublineage_shitikov[gnb] = ''
            for sublineage in Shitikov_dic.items():
                lineage_name = sublineage[0]
                for i in sublineage[1].items():
                    if i in LIST_VCF:
                        sublineage_shitikov_list += [lineage_name]
            if len(set(sublineage_shitikov_list)) == 1:
                        if len(sublineage_shitikov_list) >= 7:
                            sublineage_shitikov[gnb] = list(set(sublineage_shitikov_list))[0]
                        else:
                            sublineage_shitikov[gnb] = list(set(sublineage_shitikov_list))[0] + " supported by that many SNPs:" + str(len(sublineage_shitikov_list))
            else: # to accept multiple at least 4 positions must be fulfilled
                wanted =[]
                counts = Counter(sublineage_shitikov_list)
                for y in range(len(set(sublineage_shitikov_list))):
                    if counts[list(set(sublineage_shitikov_list))[y]] >= 4:
                        wanted.append(list(set(sublineage_shitikov_list))[y])
                if len(wanted) == 1:
                    if counts[wanted[0]] >= 7:
                        sublineage_shitikov[gnb] = wanted[0]
                    else:
                        sublineage_shitikov[gnb] = wanted[0] + " supported by that many SNPs:" + str(counts[wanted[0]])
                else:   
                    sublineage_shitikov[gnb] = ';'.join(wanted)

                                


            outputfile.write(gnb+'\t'+sublineage_shitikov[gnb]+'\n')
    else:
        print('VCF not found for:', vcf)

outputfile.close()
