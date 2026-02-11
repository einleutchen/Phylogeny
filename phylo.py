import os
import sys
import subprocess
import pandas as pd
import glob
import numpy as np
from Bio import SeqIO, Phylo
import matplotlib.pyplot as plt
from collections import defaultdict

# specify threads
threads = 4 

# define input and output paths
genomes = "finalproject/genomes"
protein = "finalproject/prodigal/proteins"
genes = "finalproject/prodigal/genes"
stats = "finalproject/prodigal/stats"
database_prot = "finalproject/database/proteins"
database_nucl = "finalproject/database/genes"
blastn = "finalproject/blast/blastn"
blastp = "finalproject/blast/blastp"

dir1 = [genomes, protein, genes, stats, database_nucl, database_prot, blastn, blastp]

# make input and output directories if they don't exist
for items in dir1:
    os.makedirs(items, exist_ok=True)

# create a run_command function
def run_command(command, filename=None, logfile=False):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if logfile:
        stdout, stderr = process.communicate()
        with open("logfile", "a") as log_file: # the file is created if it does not exist, if it does, append to logfile
            log_file.write(f"Output for {filename}:\n\n") # file name for logging
            log_file.write(f"Stdout:\n{stdout.decode('utf-8')}\n\n") # log: capture standard output
            log_file.write(f"Log:\n{stderr.decode('utf-8')}\n\n") # capture any errors or warnings, decode bytes b'' to string
            log_file.write(f"Return code: {process.returncode}\n\n") # 0 means success
    return process.returncode

# automatically download genomes from NCBI
ncbi_cmd = ["./datasets", "download", "genome", "accession",
            "--inputfile", "./accessions.txt",
            "--include", "genome",
            "--filename", "genomes_only.zip"]
# run the ncbi command
if not os.path.exists("genomes_only.zip"):
    print("Downloading genomes from NCBI...")
    run_command(ncbi_cmd, filename="NCBI_Download", logfile=True)
    # unzip the downloaded genomes
    print("Unzipping downloaded genomes...")
    os.system("unzip genomes_only.zip")
    # move the fasta files to the genomes directory
    print("Moving genome files to the genomes directory...")
    os.system("mv ncbi_dataset/data/GCA_*/* finalproject/genomes")

# remove trailing newlines from fasta files
for fasta_file in os.listdir(genomes):
    if fasta_file.endswith(".fna"):
        print(f"Remove trailing newlines from file: {fasta_file}")
        with open(os.path.join(genomes, fasta_file), "r") as f:
            for line in f:
                if not line.startswith(">"):
                    line.rstrip()
    else:
        print(f"Skipping non-fasta file: {fasta_file}")

# predict ORFs using prodigal
# run prodigal for each fasta file in input directory
for fasta_file in os.listdir(genomes):
    print(f"Running prodigal on file: {fasta_file}") 
    print(os.path.join(genomes, fasta_file))      
    prodigal_command = [
    "./prodigal-gv",
    "-i", os.path.join(genomes, fasta_file),
    "-a", os.path.join(protein, fasta_file.replace(".fna", ".faa")),
    "-d", os.path.join(genes, fasta_file)
    ]
    # execute command only if output file does not exist
    if not os.path.exists(os.path.join(protein, fasta_file.replace(".fna", ".faa"))):
        run_command(prodigal_command, filename="prodigal-protein", logfile=True)

# define a function to create blast databases
def create_blast_db(input_dir, db_type, ext, output_path, extensions:list):
    for protein_file in os.listdir(input_dir):
        print(f"Creating BLAST database for file: {protein_file}")
        db_name = protein_file.replace(ext, "_db")
        db_files = [os.path.join(output_path, db_name + ext) for ext in extensions]
        makeblastdb_command = [
            "makeblastdb",
            "-in", os.path.join(input_dir, protein_file),
            "-dbtype", db_type,
            "-out", os.path.join(output_path, db_name),
            "-parse_seqids"
        ]
        # execute command only if output file does not exist
        if not all(os.path.exists(db_file) for db_file in db_files):
            run_command(makeblastdb_command, filename="makeblastdb-protein", logfile=True)

# create blast database from predicted proteins
create_blast_db(protein, "prot", ".faa", database_prot, [".phr", ".pin", ".psq"])
# create blast database from predicted genes
create_blast_db(genes, "nucl", ".fna", database_nucl, [".nhr", ".nin", ".nsq"])

# generate a summary tsv file for all genomes
os.system("touch finalproject/genomes_summary.tsv")

# select a reference genome based on the completeness and contamination results
# Pass the whole command as one long string
summary_command= (
    "./datasets summary genome accession --inputfile ./accessions.txt "
    "--assembly-source genbank --assembly-level complete --as-json-lines | "
    "./dataformat tsv genome "
    "--fields organism-name,accession,checkm-completeness,checkm-contamination,assmstats-number-of-contigs,assmstats-atgc-count,ani-best-ani-match-ani"
    "> finalproject/genomes_summary.tsv"
)
# You must set shell=True for the pipe '|' to work
subprocess.run(summary_command, shell=True, check=True)

# read and print the summary tsv file
df = pd.read_csv("finalproject/genomes_summary.tsv", sep="\t")

# hiererarchically filter genomes with single contig > length > completeness > contamination
single_contig_genomes = df[df['Assembly Stats Number of Contigs'] == 1].sort_values(by='CheckM completeness', ascending=False)
min_completeness = single_contig_genomes['CheckM completeness'].min()
max_contamination = single_contig_genomes['CheckM contamination'].max()
filtered_genome = single_contig_genomes[(single_contig_genomes['Assembly Stats ATGC Count'] == single_contig_genomes['Assembly Stats ATGC Count'].max()) &
                                        (single_contig_genomes['CheckM completeness'] > min_completeness) &
                                        (single_contig_genomes['CheckM contamination'] < max_contamination)]

pattern = filtered_genome.iloc[0]['Assembly Accession']
search_file = os.path.join(genomes, f"{pattern}*.fna")
found_files = glob.glob(search_file)
reference_genome = os.path.basename(found_files[0]) if found_files else None

print("High quality genomes after filtering:", reference_genome)

# create a function to perform blastp and blastn
def perform_blast(which_blast, query, db, extensions, output_dir, evalue, outfmt, threads=4): 
    list_db = [item.replace(extensions, '') for item in os.listdir(db) if item.endswith(extensions)]
    for db_file in list_db:
        if db_file.replace("_db", '') != reference_genome.replace(".fna", ''):
            output_filename = db_file + ".blastout"
            output_path = os.path.join(output_dir, output_filename)
            print(f"Running {which_blast} of {query} against database {db_file}")
            blast_command = [
                which_blast,
                "-query", query,
                "-db", os.path.join(db, db_file),
                "-out", output_path,
                "-evalue", evalue,
                "-outfmt", outfmt,
                "-num_threads", str(threads)
            ]
        # execute command only if output file does not exist
            if not os.path.exists(output_path):
                run_command(blast_command, filename="blast", logfile=True)

# define blast output format
outformat = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"

#  perform blastp 
perform_blast("blastp", os.path.join(protein, reference_genome.replace(".fna", ".faa")), database_prot, ".pdb", blastp, "1e-5", outformat)

#  perform blastn
perform_blast("blastn", os.path.join(genes, reference_genome), database_nucl, ".ndb", blastn, "1e-6", outformat)

# find the gene sets overlap between blastp and blastn results
def overlap_sets(input_dir, start=1, end=1000, pident=70, querycov=70, matched_key=None, match_pairs=None):
    if matched_key is None:
        matched_key = []
    if match_pairs is None:
        match_pairs = []
    list_keys = ["CP043998.1_" + str(i) for i in range(start, end+1)]
    #print("List keys: ", list_keys[0:10])  # print first 10 keys for verification
    dict_ids = {key: 0 for key in list_keys}
    all_pairs = []
    for file in os.listdir(input_dir):
        print(f"Processing {input_dir} output file: {file}, coverage: {querycov}%, identity: {pident}%, range: {start}-{end}")
        df = pd.read_csv(os.path.join(input_dir, file), sep="\t", header=None)
        #print(df.head())
        aligncov = (df.iloc[:, 7] - df.iloc[:, 6]) / df.iloc[:, 12] * 100 # calculate alignment coverage
        pairs = df[(df.iloc[:, 2] > pident) & (aligncov > querycov)].iloc[:,0:2].values.tolist()
        single_pairs = list(set(tuple(pair) for pair in pairs))  # convert to set of tuples to remove duplicates
        all_pairs.extend(single_pairs)
        #print("Single pairs:", single_pairs[0:10])  # print first 10 pairs for verification
        #print("Last 10 values in pairs", type(pairs), pairs[-10:])
        for key in dict_ids:
            if key in [pair[0] for pair in single_pairs]:
                dict_ids[key] += 1
            else:
                continue
        #print("Dictionary: ", dict_ids)
    if any(value == 10 for value in dict_ids.values()):
        new_match = [key for key in dict_ids if dict_ids[key] == 10]
        new_pairs = [pair for pair in all_pairs if pair[0] in new_match]
        match_pairs.extend(new_pairs)
        matched_key.extend(new_match)
        #print("Matched Key:", matched_key, len(matched_key), len(dict_ids))
    if len(matched_key) >= 5:
        #print(f"Found matches: {matched_key}")
        return match_pairs
    else:
        if end >= 5000:
            #print(f"Found matches: {matched_key}")
            return match_pairs
        else:
            print(f"rerunning with extended range...")
            return overlap_sets(input_dir, start + 1000, end + 1000, pident, querycov, matched_key, match_pairs)

overlap_proteins = overlap_sets(blastp, querycov=90)
#print("Overlapping proteins in all BLAST results:", overlap_proteins)
overlap_genes = overlap_sets(blastn, querycov=50)
#print("Overlapping genes in all BLAST results:", overlap_genes)

# define a function to extract the sets
def selected_groups(array:list, colnames=["Query", "Subject"], num=5):
    selected_proteins = pd.DataFrame(array, columns=colnames)
    grouped_df = selected_proteins.groupby("Query")["Subject"].apply(list)
    set_ids = []
    group_10 = grouped_df[grouped_df.apply(len) == 10]
    if num > len(grouped_df):
        num = len(grouped_df)
        for index in range(0, num):
            set_ids.extend([grouped_df.iloc[index] + [grouped_df.index[index]]])
    else:
        num = num
        for index in range(0, num):
            set_ids.extend([group_10.iloc[index] + [group_10.index[index]]])
    return set_ids

# collects all the ids to lists
prot_ids = selected_groups(overlap_proteins)
nucl_ids = selected_groups(overlap_genes)

# clean the lists to length of 11
cleaned_prot = []
cleaned_nucl = []

for sublist in prot_ids:
    unique_items = []
    seen_prefixes = set()
    
    for item in sublist:
        if len(unique_items) == 11:
            break
        prefix = item.rsplit('_', 1)[0]
        if prefix not in seen_prefixes:
            unique_items.append(item)
            seen_prefixes.add(prefix)
            
    cleaned_prot.append(unique_items)

for sublist in nucl_ids:
    unique_items = []
    seen_prefixes = set()
    
    for item in sublist:
        if len(unique_items) == 11:
            break
        prefix = item.rsplit('_', 1)[0]
        if prefix not in seen_prefixes:
            unique_items.append(item)
            seen_prefixes.add(prefix)
            
    cleaned_nucl.append(unique_items)

# name the sets
prot1, prot2, prot3, prot4, prot5 = cleaned_prot[0], cleaned_prot[1], cleaned_prot[2], cleaned_prot[3], cleaned_prot[4]
nucl1, nucl2, nucl3 = cleaned_nucl[0], cleaned_nucl[1], cleaned_nucl[2]

# # list of all sets
protein_list = [prot1, prot2, prot3, prot4, prot5]
gene_list = [nucl1, nucl2, nucl3]

os.system(f"rm -r finalproject/MSA/ finalproject/iqtree finalproject/tree_plots")

# directories for MSA and iqtree
MSA_protein = "finalproject/MSA/proteins"
MSA_gene = "finalproject/MSA/genes"
MSA_alignprot = "finalproject/MSA/alignment/proteins"
MSA_aligngene = "finalproject/MSA/alignment/genes"
iqtree_protein = "finalproject/iqtree/protein"
iqtree_gene = "finalproject/iqtree/genes"
plots_p = "finalproject/tree_plots/protein"
plots_n = "finalproject/tree_plots/gene"

# create directories and subdirectories
dir2 = [MSA_aligngene, MSA_alignprot, MSA_gene, MSA_protein, iqtree_gene, iqtree_protein, plots_n, plots_p]

for dir in dir2:
    os.makedirs(dir, exist_ok=True)

# define a function to extract fasta sequences from id sets
def extract_seq(ids_list:list, filename=None, extension=None, input_dir=None, output_dir=None):
    output_file = os.path.join(output_dir, f"{filename}.fasta")
    if not os.path.exists(output_file):
        for id in ids_list:
            seqkit_command = f"seqkit grep -p {id} {input_dir}/*.{extension} >> {output_file}" 
            subprocess.run(seqkit_command, shell=True, check=True)
idx1 = 0
idx2 = 0
# extract protein sequences for concatenation
for item in protein_list:
    idx1 += 1
    print(f"Extracting protein sequences: {item}")
    outname = f"prot_{str(idx1)}"
    extract_seq(item, filename=outname, extension="faa", input_dir=protein, output_dir=MSA_protein)

# extract nucleotide sequences for concatenation
for item in gene_list:
    idx2 += 1
    print(f"Extracting nucleotide sequences: {item}")
    outname = f"prot_{idx2}"
    extract_seq(item, filename=outname, extension="fna", input_dir=genes, output_dir=MSA_gene)

# modify header and remove * at the end of each sequences in FASTA files
def oneline(input_dir):
    for file in os.listdir(input_dir):
        print("Modifying ", file)
        file_path = os.path.join(input_dir, file)
        if file.endswith(".fasta"):
            basename = os.path.splitext(file)[0]
            output_file = os.path.join(input_dir, f"{basename}_cleaned.fasta")
            if os.path.exists(output_file):
                os.system(f"rm {output_file}")
            with open(output_file, "w") as out_handle:
                for record in SeqIO.parse(file_path, "fasta"):
                    new_id = record.id.split("_")[0]
                    seq = str(record.seq).replace("\n", "").rstrip("*").strip()
                    out_handle.write(f">{new_id}\n{seq}\n")

oneline(MSA_protein)
oneline(MSA_gene)

# define a function for MUSCLE process
def msa(input_dir, output_dir):
    for file in os.listdir(input_dir):
        if file.endswith("_cleaned.fasta"):
            print(f"Running MUSCLE on {file}")
            infile = os.path.join(input_dir, file)
            basename = os.path.splitext(file)[0]
            output_path = os.path.join(output_dir, f"{basename}_MSA.fasta" )
            if not os.path.exists(output_path):
                muscle_command = f"muscle -in {infile} -out {output_path}"
                subprocess.run(muscle_command, shell=True, check=True)

msa(MSA_protein, MSA_alignprot)
msa(MSA_gene, MSA_aligngene)

# define a function to call seqkit concat
def concatenate(input_dir):
    output_file = os.path.join(input_dir, "all_cleaned_MSA.fasta")
    if not os.path.exists(output_file):
        seqkit_concat = f"seqkit concat {input_dir}/*cleaned_MSA.fasta > {output_file}" 
        subprocess.run(seqkit_concat, shell=True, check=True)

concatenate(MSA_alignprot)
concatenate(MSA_aligngene)

# define a function to rename duplicated headers in MSA files
def rename_dup(input_dir):
    for file in os.listdir(input_dir):
        print("Renaming duplicated headers ", file)
        file_path = os.path.join(input_dir, file)
        if file.endswith("_MSA.fasta"):
            basename = os.path.splitext(file)[0]
            output_file = os.path.join(input_dir, f"{basename}_renamed.fasta")
            if os.path.exists(output_file):
                print(f"Skipping {file}, cleaned version already exists.")
                continue
            record_count = defaultdict(int)
            with open(output_file, "w") as out_handle:
                for record in SeqIO.parse(file_path, "fasta"):
                    initial_id = record.id
                    record_count[initial_id] += 1
                    if record_count[initial_id] > 1: 
                        new_id = f"{initial_id}_{record_count[initial_id]}" 
                    else:
                        new_id = initial_id                  
                    seq = str(record.seq).replace("\n", "").rstrip("*").strip()
                    out_handle.write(f">{new_id}\n{seq}\n")

rename_dup(MSA_alignprot)
rename_dup(MSA_aligngene)


# define a function for IQTREE2 process
def tree(input_dir, output_dir):
    for file in os.listdir(input_dir):
        if file.endswith("_renamed.fasta"):
            print(f"Running IQTREE on {file}")
            base_name = os.path.splitext(file)[0]
            file_prefix = os.path.join(output_dir, base_name)
            check = file_prefix + ".treefile"
            if not os.path.exists(check):
                infile = os.path.join(input_dir, file)
                iqtree_command = f"iqtree2 -s {infile} --prefix {file_prefix}"
                subprocess.run(iqtree_command, shell=True, check=True)

tree(MSA_alignprot, iqtree_protein)
tree(MSA_aligngene, iqtree_gene)

# visualize contree tree

mapping = {
    "CP009933.1": "Clostridium scatologenes",
    "CP009225.1": "Clostridium sporogenes",
    "CP009687.1": "Clostridium aceticum",
    "CP009688.1": "Clostridium aceticum",
    "CP020559.1": "Clostridium formicaceticum",
    "CP020953.1": "Clostridium drakei",
    "CP023671.1": "Clostridium septicum",
    "CP023672.1": "Clostridium septicum",
    "CP043998.1": "Clostridium diolis",
    "CP040626.1": "Clostridium butyricum",
    "CP040627.1": "Clostridium butyricum",
    "CP040628.1": "Clostridium butyricum",
    "CP040629.1": "Clostridium butyricum",
    "CP065681.1": "Clostridium perfringens",
    "CP065680.1": "Clostridium perfringens",
    "CP073653.1": "Clostridium beijerinckii",
    "CP073654.1": "Clostridium beijerinckii",
    "LT906477.1": "Clostridium cochlearium"
}

def plot_tree(input_dir, output_dir, id_map):
    os.makedirs(output_dir, exist_ok=True)
    for file in os.listdir(input_dir):
        if file.endswith(".treefile"):
            input_file = os.path.join(input_dir, file)
            output_path = os.path.join(output_dir, f"{file}.png")
            if not os.path.exists(output_path):
                tree = Phylo.read(input_file, "newick")
                for clade in tree.get_terminals():
                    if clade.name and clade.name in id_map:
                        clade.name = id_map[clade.name]
                fig = plt.figure(figsize=(10,8))
                Phylo.draw(tree, do_show=False)
                # add axes representing branch length
                ax = plt.gca()
                print(f"Saved tree plot: {output_path}")
                ax.spines['top'].set_visible(False)
                ax.spines['right'].set_visible(False)
                ax.spines['left'].set_visible(False)
                ax.get_yaxis().set_visible(False)
                ax.spines['bottom'].set_visible(True) 
                ax.set_xlabel("Branch Length (Substitutions per site)") 
                # save the pltots
                plt.savefig(output_path)
                plt.close()

plot_tree(iqtree_protein, plots_p, mapping)
plot_tree(iqtree_gene, plots_n, mapping)

print(f"-------------------------------THE END------------------------------------")