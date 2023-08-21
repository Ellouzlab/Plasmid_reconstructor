import argparse, shutil
from utils import *


def arguments():
    parser = argparse.ArgumentParser(description="Finds Pantoea plasmids based on homology")

    parser.add_argument('-i', '--infile', required=True, type=str, help="Fasta with contigs to be classified fasta")
    parser.add_argument('-p', '--plasdb', default="plasmid_db/clust_map", type=str, help="DB containing plasmid fasta")
    parser.add_argument('-c', '--chromdb', default="chromdb/chromdb", type=str, help="DB containing chromosome fasta")
    parser.add_argument('-r', '--repdb', default='rep_db/rep_dna', type=str, help="DB containing ORI fasta")
    parser.add_argument('--tempdir', default="tmp", type=str, help="Directory to contain temporary files.")
    parser.add_argument('-t', '--threads', default=1, type=int, help="Number of threads for blast searches.")
    parser.add_argument('-o', '--outdir', default="output", type=str, help="output directory")

    args = parser.parse_args()
    return args

def isplasmid(plas_df, chrom_df):
    max_plas = plas_df["query_coverage"].max()
    if len(plas_df) == 0:
        max_plas = 0

    max_chrom = chrom_df["query_coverage"].max()
    if len(chrom_df) == 0:
        max_chrom = 0

    if max_chrom > max_plas:
        return False
    elif max_plas > max_chrom:
        return True
    else:
        return False

def contig_class(contigs_loc, blast_plas_loc, blast_chrom_loc, blast_rep_loc, plasdb, chromdb, repdb ,threads):
    contig_data_arr=[]
    for contig in os.listdir(contigs_loc):
        contig_path=f"{contigs_loc}/{contig}"

        blast_plas_path=f"{blast_plas_loc}/{contig}.csv"
        blast_chrom_path=f"{blast_chrom_loc}/{contig}.csv"
        blast_rep_path = f"{blast_rep_loc}/{contig}.csv"

        plas_df=blast_run(contig_path, plasdb, threads=threads)
        plas_df.to_csv(blast_plas_path)

        chrom_df=blast_run(contig_path, chromdb, threads=threads)
        chrom_df.to_csv(blast_chrom_path)

        rep_df=blast_run(contig_path, repdb, threads=threads)
        rep_df.to_csv(blast_rep_path)

        plas_result=isplasmid(plas_df, chrom_df)

        #print(contig_path)
        try:
            rep=True
            rep_id=rep_df.iloc[0]['qtitle'].split('|')[-1]
        except:
            rep=False
            rep_id=None

        try:
            idx_max_value = plas_df['Match_base_pid'].idxmax()

            # Retrieve the corresponding title
            cluster = plas_df.loc[idx_max_value, 'qtitle']
        except:
            cluster=None
        #print(cluster)

        contig_data_arr.append((contig_path, chrom_df, plas_df, plas_result, rep, cluster, rep_id))
    return contig_data_arr

def write_chrom_fasta(contig_arr, location):
    chromosome_records=[]
    for contig in contig_arr:
        if not contig[3]:
            record=read_fasta_single(contig[0])
            chromosome_records.append(record)
    if not len(chromosome_records)==0:
        SeqIO.write(chromosome_records, location, "fasta")


def write_plas_fasta(contig_info, outdir):
    plasmid_arr=[]

    for contig in contig_info:
        if contig[3]:
            plasmid_arr.append(contig)
    #print(contig[0])

    prim_clust_dict={}
    unclustered=[]
    for contig in plasmid_arr:
        clust=contig[5]

        if contig[4]:
            if clust not in prim_clust_dict:
                prim_clust_dict[clust]=[contig]
            else:
                prim_clust_dict[clust].append(contig)
        else:
            unclustered.append(contig)


    final_unclustered=[]
    for contig in unclustered:
        i=0
        clusters_unc=contig[5].split('_')
        match = False
        while i<99:

            for key in prim_clust_dict:
                clusters_cl=key.split('_')
                if clusters_cl[i]==clusters_unc[i]:
                    match=True
                    prim_clust_dict[key].append(contig)
                    break
            if match:
                break
            i+=1
        if not match:
            final_unclustered.append(contig)


    for key in prim_clust_dict:
        plasclust2fasta(prim_clust_dict[key], outdir)

    if len(final_unclustered)==0:
        pass
    else:
        unclust2fasta(final_unclustered, outdir)

def main():
    args=arguments()
    chromdb=args.chromdb
    plasdb=args.plasdb
    repdb=args.repdb
    tempdir = args.tempdir
    threads=args.threads
    outdir=args.outdir

    makedir(tempdir)
    makedir(outdir)

    contig_location = makedir(f"{tempdir}/contigs")
    blast_plas_location = makedir(f"{tempdir}/blast_plas_res")
    blast_chrom_location = makedir(f"{tempdir}/blast_chrom_res")
    blast_rep_location = makedir(f"{tempdir}/blast_rep_res")

    fasta_breaker(args.infile, contig_location)

    contig_info = contig_class(
        contig_location, blast_plas_location, blast_chrom_location, blast_rep_location, plasdb, chromdb, repdb, threads
    )
    chrom_path = f"{outdir}/chromosome.fasta"
    write_chrom_fasta(contig_info, chrom_path)
    write_plas_fasta(contig_info, outdir)
    shutil.rmtree(tempdir)


main()