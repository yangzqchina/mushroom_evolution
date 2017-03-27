import sys
file=sys.argv[1]
def read_fasta(file):
    file2=open(file,'r')
    seq={}
    for line in file2:
        if line.startswith('>'):
            name=line.strip()
            seq[name]=''
        else:
            seq[name]+=line.upper().replace('\n','')
    file2.close()
    return(seq)
def write_fasta(file,seq):
    file2=open(file,'w')
    seq_heads=seq.keys()
    for head_i in seq_heads:
        file2.write(head_i+'\n'+seq[head_i]+'\n')
    file2.close()
outfile=sys.argv[2]
write_fasta(outfile,read_fasta(file))
