def count_fasta(fasta_file):
    '''count the number of sequences in a fasta file'''
    from Bio import SeqIO
    return sum(1 for record in SeqIO.parse(fasta_file, 'fasta'))

def check_seq_duplicates(fasta_file):
    '''check for duplicate sequences in a fasta file'''
    from Bio import SeqIO
    seqs = [str(record.seq) for record in SeqIO.parse(fasta_file, 'fasta')]
    return len(seqs) == len(set(seqs))

def fasta_to_dict(fasta_file):
    '''convert fasta file to dictionary'''
    from Bio import SeqIO
    return {record.description : record.seq for record in SeqIO.parse()}