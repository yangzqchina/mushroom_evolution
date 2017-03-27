#!/usr/bin/env python3
#Short script to take output file from exonerate
#and parse out the gff output; convert it to GTF and concatenate into a single file.
#Julin Maloof
#April 19, 2013

import fileinput, sys, re

if len(sys.argv)  == 1 :
    print("This script takes a series of exonerate output files that include gff information,")
    print("extracts the gff information, converts it to gtf, and outputs it to standard out.")
    print("Use file names as arguments.")
    exit()

inGFF = False   
firstGene = True
transcripts = {}

#create regular expression that can find gene_id for each entry
regGene = re.compile(r"^.+\t.+\tgene(\t[0-9]+){3}\t[+-]\t.\tgene_id [0-9]+ ; sequence (?P<gene_id>\S+)")

#create regular expression to match first 8 fields of gff entry
regGff = re.compile(r"^(.+\t){8}")

for line in fileinput.input() :
    if line == "# --- START OF GFF DUMP ---\n" :
        inGFF = True
        gene_id = "" # reset gene name
        transcript_id = ""
        continue
    if line == "# --- END OF GFF DUMP ---\n" :
        firstGene = False
        inGFF = False
        continue
    if inGFF :    
        line = line.rstrip("\n")
        if line == "##gff-version 2" : 
            continue
        fields = line.split("\t")
        if (firstGene and line.startswith("#")) : #only print header info for first gene
            print(line)
        #search for line containing gene id information
        geneMatch = regGene.match(line)
        if geneMatch :
            gene_id = geneMatch.group("gene_id") 
            transcripts[gene_id] = transcripts.get(gene_id, 0) + 1
            transcript_id = gene_id + "." + str(transcripts[gene_id])
            fields[8] = fields[8] = 'gene_id "' + gene_id + '"'
            print("\t".join(fields))
            continue
        #reformat lines to proper gtf format, replacing 9th field with gene_id and transcript_if
        gffMatch = regGff.match(line)
        if gffMatch :
            fields[8] = 'gene_id "' + gene_id + '"; transcript_id "' + transcript_id + '"'
            print("\t".join(fields))
    

    