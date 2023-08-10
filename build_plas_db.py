import argparse, shutil
import subprocess as sp
from codecs import decode
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
import numpy as np
from utils import *


def arguments():
    parser = argparse.ArgumentParser(description="Creates a plasmid database for homoplas using skani and makeblastdb")
    parser.add_argument('-i', '--input', type=str, help="input fasta", required=True)
    parser.add_argument('-o', '--outdir', type=str, help="output directory", required=True)
    parser.add_argument('-t', '--threads', type=int, default=1, help='Number of threads to use')
    args = parser.parse_args()
    return args

def makeblastdb(input, outdir):
    out_path=f"{outdir}/clust_map"
    cmd=f"makeblastdb -in {input} -dbtype nucl -out {out_path}"
    sp.run(cmd, shell=True)
    return out_path

def skani_triangle(input, output, threads):
    cmd=f'skani triangle -i {input} -t {threads} --distance --slow -m 300 -o {output}'
    sp.run(cmd, shell=True)

def rename_seq(db, clus_arr, output):
    fasta_seq = read_fasta(db)
    for i in range(len(clus_arr[0])):
        id=''
        for cluster in clus_arr:
            id+=f"{str(cluster[i])}_"
            fasta_seq[i].id = id
    write_fasta(output, fasta_seq)

def cluster_id_fasta(distance_file, input, output):
    with open(distance_file, 'r') as file:
        lines = file.readlines()
    lines=lines[1:]
    dimension_sq=len(lines[1:])
    arr = np.zeros((dimension_sq, dimension_sq), dtype=float)

    line_num=0

    for line in lines:
        line_arr=line.split("\t")[1:]
        entry_num=0

        for entry in line_arr:
            try:
                arr[line_num, entry_num]=entry
                try:
                    arr[entry_num, line_num]=entry
                except:
                    pass
                entry_num+=1
            except:
                pass
        line_num+=1

    condensed_distance = squareform(arr)

    # Calculate the linkage matrix using single linkage
    linkage_matrix = linkage(condensed_distance, method='single')


    start_val=1
    cluster_arr=[]
    for i in range(99):
        num=(i+1)*start_val
        cluster=fcluster(linkage_matrix, num, criterion='distance').tolist()
        cluster_arr.append(cluster)


    rename_seq(input, cluster_arr, output)

def main():
    args = arguments()
    input=args.input
    outdir=args.outdir
    threads=args.threads

    makedir(outdir)

    if __name__ == '__main__':
        print('\n\n\nRunning build_plas_db.py')

    #make a blast database

    #copy database fasta into blast database directory
    db_fasta_path=f"{outdir}/original_seq.fasta"
    shutil.copy(input, db_fasta_path)

    #make a distance matrix with skani
    print("\n\nBuilding distance array with skani\n\n")
    output_dist=f'{outdir}/dist.tsv'
    skani_triangle(db_fasta_path, output_dist, threads)

    output_fasta=f"{outdir}/clust_map.fasta"

    print("\n\nWriting Cluster IDs")
    cluster_id_fasta(output_dist, db_fasta_path, output_fasta)

    makeblastdb(output_fasta, outdir)

    print("\n\nDone running build_plas_db.py")
main()