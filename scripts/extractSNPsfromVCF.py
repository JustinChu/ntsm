import argparse
from pyfaidx import Fasta
import os
import math
import sys

class VCFEntry:
    def __init__(self, chr, pos, wt, variant):
        self.chr = chr
        self.pos = int(pos)
        self.wt = wt
        self.variant = variant

class ExtractKmers:
    def __init__(self, vcf, fasta, k, prefix, ignore):
        self._vcf = vcf
        self._fasta = fasta
        self._k = k
        self._prefix = prefix
        self._ignore = ignore
        self._vcfEntries = {}
    
    #Returns True if bases are not purine to pyrimidine
    def _checkVariant(self, base1, base2):
        if(base1 == "A" or base1 == "T"):
            if(base2 == "A" or base2 == "T"):
                return True
            else:
                return False
        elif(base1 == "C" or base1 == "G"):
            if(base2 == "C" or base2 == "G"):
                return True
            else:
                return False
    
    #Return false if ordering needs to be reversed (output order not strand)
    def _orderVariant(self, base1, base2):
        if(base1 == "A" or base1 == "T"):
            return True
        else:
            return False
    
    def _parseVCF(self):
        vcfFH = open(self._vcf, 'r')
        #chr1    45508256    rs2275276    A    G    .    .    .  
        lines = vcfFH.readlines()
        # Strips the newline character
        idCounter = 0
        for line in lines:
            if line[0] != "#":
                tmpArr = line.rstrip().split("\t")
                snpID = tmpArr[2]
                if snpID == '.':
                    snpID = idCounter
                    idCounter += 1
                #throw an error if multiple alleles are listed in alternate
                if len(tmpArr[4]) > 1:
                    print("Error: Multiple alternate alleles found in VCF", file=sys.stderr)
                    exit(1)
                #store chr, pos, variant
                info = VCFEntry(tmpArr[0], tmpArr[1], tmpArr[3], tmpArr[4])
                self._vcfEntries[str(snpID)] = info
            
    def extract(self):
        self._parseVCF()
        #open fasta file
        fastaFile = Fasta(self._fasta)
        refFH = open(self._prefix + "_AT.fa", 'w')
        varFH = open(self._prefix + "_CG.fa", 'w')
        
        removeCount = 0
        #for each vcf entry extract wildtype and variant into a string
        for id in self._vcfEntries.keys():
            # print string to repective files
            chr = self._vcfEntries[id].chr
            offset = self._vcfEntries[id].pos - 1
            pos1 = math.ceil(offset - self._k / 2)
            pos2 = pos1 + self._k
            tmpStr = str(fastaFile[self._vcfEntries[id].chr][pos1:pos2])
            tmpStr = tmpStr.upper()
            if(self._vcfEntries[id].wt != tmpStr[int(self._k / 2)]):
                print("Wildtype allele does not match", file=sys.stderr)
                print("ref:" + self._vcfEntries[id].wt, file=sys.stderr)
                print("var:" + self._vcfEntries[id].variant, file=sys.stderr)
                print("fasta:" + str(fastaFile[self._vcfEntries[id].chr][offset]), file=sys.stderr)
                print("kmer:" + tmpStr, file=sys.stderr)
                removeCount += 1
                continue
            
            if(self._checkVariant(self._vcfEntries[id].wt, self._vcfEntries[id].variant)):
                print("Warning: " + id +" " +self._vcfEntries[id].variant + " " + self._vcfEntries[id].wt 
                      + " is not a purine <-> pyrimidine variant.", file=sys.stderr)
                removeCount += 1
                if(self._ignore):
                    continue
            modStr = tmpStr[0:int(self._k / 2)] + self._vcfEntries[id].variant + tmpStr[int(self._k / 2) + 1:]
            varFH.write(">" + id + "\n")
            refFH.write(">" + id + "\n")            
            if(self._orderVariant(self._vcfEntries[id].wt, self._vcfEntries[id].variant)):
                varFH.write(modStr + "\n")
                refFH.write(tmpStr + "\n")
            else:
                varFH.write(tmpStr + "\n")
                refFH.write(modStr + "\n")
        refFH.close()
        varFH.close()
        print("Removed " + str(removeCount) + " SNPs.", file=sys.stderr)
        
        
def main():
    parser = argparse.ArgumentParser(description='Extract k-mers listed from VCF file from reference. Produces 2 files (reference and variant)')
    parser.add_argument("-v", '--vcf', type=str, dest='vcf', help='vcf File')
    parser.add_argument("-f", '--fa', type=str, dest='fasta', help='fasta file with fai index')
    parser.add_argument("-k", '--kmer', type=int, dest='kmer', help='kmer size', default=25)
    parser.add_argument("-p", '--prefix', type=str, dest='prefix', help='output prefix', default = "")
    parser.add_argument("-i", '--ignoreReq', action='store_false', dest='ignore', help='ignore AT to CG conversion requirements', default = True)

    args = parser.parse_args()
    extractor = ExtractKmers(args.vcf, args.fasta, args.kmer, args.prefix, args.ignore)
    extractor.extract()
    

main()
    
