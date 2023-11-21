#!/usr/bin/python3
import argparse, shutil, traceback
from utils import *


def arguments():
    parser = argparse.ArgumentParser(description="Finds Pantoea plasmids based on homology")

    parser.add_argument('-i', '--infile', required=True, type=str, help="Fasta with contigs to be classified fasta")
    parser.add_argument('-p', '--plasdb', default="plasmid_db/clust_map", type=str, help="DB containing plasmid fasta")
    parser.add_argument('-c', '--chromdb', default="chromdb/chromdb", type=str, help="DB containing chromosome fasta")
    parser.add_argument('-r', '--repdb', default='rep_db/rep_dna', type=str, help="DB containing ORI fasta")
    parser.add_argument('--tempdir', type=str, help="Directory to contain temporary files.")
    parser.add_argument('-t', '--threads', default=1, type=int, help="Number of threads for blast searches.")
    parser.add_argument('-o', '--outdir', default="output", type=str, help="output directory")
    parser.add_argument('--keep', action='store_true',help="If used, temporary files not deleted")
    args = parser.parse_args()
    if args.tempdir is None:
        args.tempdir = f"{args.outdir}/tmp"
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

def remove_overlap(df):
    df = df.sort_values(by=['qstart', 'bitscore'], ascending=[True, False])
    to_remove = set()
    last_qend = -1
    for i, row in df.iterrows():
        if row['qstart'] <= last_qend:
            to_remove.add(i)
        else:
            last_qend = row['qend']
    df = df.loc[~df.index.isin(to_remove)]

    return df
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

        plas_result=isplasmid(plas_df, chrom_df)

        #print(contig_path)
        try:
            if not len(rep_df)==0:
                rep = True
                filtered_rep_df= rep_df[rep_df['sseqcov']>0.2]
                filtered_rep_df.loc[filtered_rep_df["sseqcov"] < 0.8, 'sseqid'] = filtered_rep_df['sseqid'] + '_like'
                filtered_rep_df=remove_overlap(filtered_rep_df)
                filtered_rep_df.to_csv(blast_rep_path)
                unique_values = filtered_rep_df['sseqid'].str.split('|').apply(lambda x: x[-1]).unique()
                rep_id='&'.join(unique_values)
            else:
                rep=False
                rep_id=False
        except:
            rep=False
            rep_id=None

        try:
            idx_max_value = plas_df['bitscore'].idxmax()
            cluster = plas_df.loc[idx_max_value, 'sseqid']
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

def similar_bitscore(score1, score2, threshold):
    return abs(score1 - score2) / score1 < threshold

def get_similar_hits(df, threshold=0.01):

    highest_score = df["Match_base_pid"].max()
    similar_seq_ids = df[df['Match_base_pid'].apply(lambda x: similar_bitscore(x, highest_score, threshold))]['sseqid'].tolist()
    return [similar_seq_ids]


def write_plas_fasta(contig_info, outdir):
    plasmid_arr=[]

    for contig in contig_info:
        if contig[3]:
            plasmid_arr.append(contig)
    #print(contig[0])

    prim_clust_dict={}
    unclustered=[]
    for contig in plasmid_arr:
        plas_ids=get_similar_hits(contig[2])[0]
        if len(plas_ids)==1:
            clust=plas_ids[0]
            if contig[4]:
                if clust not in prim_clust_dict:
                    prim_clust_dict[clust]=[contig]
                else:
                    prim_clust_dict[clust].append(contig)
            else:
                unclustered.append(contig)

        else:
            match=False
            for clust in plas_ids:
                if contig[4]:
                    if clust in prim_clust_dict:
                        prim_clust_dict[clust].append(contig)
                        match=True
                        break

            if match == False:
                if contig[4]:
                    prim_clust_dict[contig[5]] = [contig]
                else:
                    unclustered.append(contig)


    final_unclustered=[]
    for contig in unclustered:
        i=0
        clusters_unc=contig[5].split('_')
        match = False
        while i<9:

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
def except_action(tempdir, outdir, keep_tmp):
    if not keep_tmp:
        shutil.rmtree(tempdir)
    shutil.rmtree(outdir)

def main():
    args=arguments()
    chromdb=args.chromdb
    plasdb=args.plasdb
    repdb=args.repdb
    tempdir = args.tempdir
    threads=args.threads
    outdir=args.outdir
    keep_tmp=args.keep

    makedir(outdir)
    tempdir = makedir(tempdir)
    try:
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
        if not keep_tmp:
            shutil.rmtree(tempdir)
    except KeyboardInterrupt:
        print("\nprogram was interrupted")
        except_action(tempdir, outdir, keep_tmp)
    except Exception:
        print(f"\nAn error occured:\n")
        traceback.print_exc()
        except_action(tempdir, outdir, keep_tmp)

main()

