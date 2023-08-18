from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import random, argparse, os
from utils import *

def arguments():
    parser = argparse.ArgumentParser(description="Tests the performance of Plasmid_reconstrutor using complete genomes")

    parser.add_argument('-i', '--indir', required=True, type=str, help="Input directory of complete genomes")
    parser.add_argument('-n', '--num_frag', required=True, type=int, help='Number of fragments to cut contigs into')
    parser.add_argument('-o', '--outdir', required=False, default='test_out',type=str, help='output directory')
    parser.add_argument('-t', '--threads', required=False, default=1,type=int, help='number of threads to use')

    args = parser.parse_args()
    return args


def fragment_contigs(seq_list, total_fragments):
    contigs = seq_list
    total_length = sum(len(record) for record in contigs)

    result = []
    current_fragments = 0

    for record in contigs:
        contig_fragments = int(len(record) / total_length * total_fragments)
        contig_fragments = min(contig_fragments, total_fragments - current_fragments)
        if contig_fragments == 0:
            result.append(record)
            continue

        cut_positions = sorted(random.sample(range(1, len(record)), contig_fragments - 1))
        cut_positions = [0] + cut_positions + [len(record)]
        for i in range(len(cut_positions) - 1):
            start = cut_positions[i]
            end = cut_positions[i + 1]
            fragment_seq = record.seq[start:end]
            fragment = SeqRecord(fragment_seq, id=record.id, description="")
            result.append(fragment)
        current_fragments += contig_fragments
    return result

def main():
    args=arguments()
    indir=args.indir
    frags=args.num_frag
    outdir=args.outdir

    makedir(outdir)

    for file in os.listdir(indir):

        file_path=f"{indir}/{file}"
        seq_list=read_fasta(file_path)

        plasmid_num = 1
        for record in seq_list:
            length=len(record)
            if length<2000000:
                record.id=f"plasmid_{plasmid_num}"
                plasmid_num += 1
            else:
                record.id="chromosome"



        fragments = fragment_contigs(seq_list, frags)
        output_file=f"{outdir}/{file}"
        write_fasta(output_file, fragments)

main()