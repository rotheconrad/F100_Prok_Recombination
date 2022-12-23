# F100 as a signal for prokaryote recombination

Contains the code and workflow for the F100 recombination project.

This workflow will assess the frequency of genes with 100% sequence similarity (F100) between a pair genomes compared to the total reciprocal best matches (RBMs) of genes shared between them for all genome pairs within a species. A collection of genomes in fasta format belonging to the same species is all that is required as input.

Assuming that a high F100 score for a given average nucleotide identity (ANI) indicates recent homologous recombination between a pair of genomes, this workflow will identify genomes from a set of species genomes with the greatest amount of recent recombination (i.e. which genome pairs share the most identical genes).

The users genomes may be compared to a model of 330+ species genomes from NCBI that have at least 10 Complete genomes per species, a model of simulated random neutral evolution genomes along a 95%-100% ANI gradient with zero recombination, or a model built from the users own set of genomes.

The workflow will yeild many intermediate files and several publication ready figures.

1. All vs All genome pair fastANI results
1. Shared fraction vs ANI scatter plot PDF
1. CDC gene prediction fasta files from Prodigal
1. RBM sequence similarities for each genome pair
1. F100 score for each genome pair in .tsv file
1. (OPTIONAL) Histogram of RBM sequence similarity for each genome pair (PDFs)
1. F100 vs ANI with GAM model PDF
1. 

See 01a_Building_Complete_Genomes_Model.txt for detailed notes/methods used to build the model with NCBI Complete genomes. Code referenced in these notes can be found in 00b_Complete_Genomes_Model directory.

See 01b_Building_Simulated_Genomes_Model.txt for detailed notes/methods used to build the model with simulated genomes. Code referenced in these notes can be found in 00c_Simulated_Genomes_Model directory.

# Require dependencies

- [Prodigal](https://github.com/hyattpd/Prodigal)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [enveomics](http://enve-omics.ce.gatech.edu/enveomics/docs) (for aai.rb)
- [ruby](https://www.ruby-lang.org/en/) (for enveomics)
- [python](https://www.python.org/)

### References

1. Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics. 2010 Dec;11(1):1-1.
1. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC bioinformatics. 2009 Dec;10(1):1-9.
1. Rodriguez-R LM, Konstantinidis KT. The enveomics collection: a toolbox for specialized analyses of microbial genomes and metagenomes. PeerJ Preprints; 2016 Mar 27.
1. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology. 2017 Nov;35(11):1026-8.
1. Flanagan D, Matsumoto Y. The Ruby Programming Language: Everything You Need to Know. " O'Reilly Media, Inc."; 2008 Jan 25.
1. Sanner MF. Python: a programming language for software integration and development. J Mol Graph Model. 1999 Feb 1;17(1):57-61.

# Required packages for Python

- [pandas](https://pandas.pydata.org/) 
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- [datashader](https://datashader.org/)
- [pygam](https://pygam.readthedocs.io/)

*Python and all packages can be easily installed with conda or pip. Prodigal, BLAST+ and MMseqs2 can also be installed with conda.*

### References

1. Van Rossum G, Drake FL. Python 3 Reference Manual. Scotts Valley, CA: CreateSpace; 2009.
1. McKinney W, others. Data structures for statistical computing in python. In: Proceedings of the 9th Python in Science Conference. 2010. p. 51–6.
1. Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, et al. Array programming with NumPy. Nature. 2020;585:357–62.
1. Virtanen P, Gommers R, Oliphant TE, Haberland M, Reddy T, Cournapeau D, et al. SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods. 2020;17:261–72.
1. Hunter JD. Matplotlib: A 2D graphics environment. Computing in science & engineering. 2007;9(3):90–5.
1. Waskom ML. Seaborn: statistical data visualization. Journal of Open Source Software. 2021 Apr 6;6(60):3021.
1. James A. Bednar, Jim Crist, Joseph Cottam, and Peter Wang (2016). "Datashader: Revealing the Structure of Genuinely Big Data", 15th Python in Science Conference (SciPy 2016).
1. Servén D., Brummitt C. (2018). pyGAM: Generalized Additive Models in Python. Zenodo. DOI: 10.5281/zenodo.1208723

# STEP 01: Prepare your genomes

This workflow is intended for a collection of genomes belonging to the same species. Start with your genome files in fasta format in their own directory. We will refer to this directory as ${genomes_dir}.

1. Rename the fasta deflines of your genome files. This is necessary to ensure all contigs (or chromosomes/plasmids) in your genome files follow the same format for downstream processing. 

Because genomes downloaded from NCBI follow a typical naming convention of, "GCF_000007105.1_ASM710v1_genomic.fna," the default behavior of this script is cut the third underscore position ("\_") and use it as a prefix for renaming the fasta deflines in numeric consecutive order.

So for a fasta file with three contigs or chromosomes the script cuts
"ASM710v1" from filename and renames fasta deflines as:

>ASM710v1_1
>ASM710v1_2
>ASM710v1_n

```bash
for f in ${genomes_dir}/*; do python 00d/Workflow_Scripts/01a_rename_fasta.py -i $f; done
```

Alternatively, the user can input their own desired prefix.

```bash
for f in ${genomes_dir}/*; do n=`echo basename $f | cut -d_ -f3`; python 00d/Workflow_Scripts/01a_rename_fasta.py -i $f -p $n; done
```

2. (OPTIONAL) Check shared fraction vs ANI of your genomes (fastANI)

It is a good idea to know how similar your genomes are to eachother. Sometimes you may split your genomes into two or more groups, or remove a few outlier genomes. This can be done manually as needed by removing the desired fasta files from ${genomes_dir} after reviewing the fastANI plot(s).

The end result we want is a single file with all vs. all genome pairs. There are many ways to achieve and depending if you're running fastANI locally or on a cluster there are multiple parallization options. Here is a simple example using the many to many option (consult the fastANI documentation and/or your clusters best practices for other options).

```bash
# make a list of your fasta files with the path
for f in ${genomes_dir}/*; do echo $f; done > genome_file_list.txt

# run fastANI in many to many mode
./fastANI --ql genome_file_list.txt --rl genome_file_list.txt -o fastANI_allV.ani
```

Plot results

```bash
# look at script details and options
python 00d_Workflow_Scripts/01b_fastANI_scatter_pyGAM.py -h
# basic plot
python 00d_Workflow_Scripts/01b_fastANI_scatter_pyGAM.py -i fastANI_allV.ani -s test_species -o fastANI_allV_sharedfrac_ANI.pdf -m True
```

# Step 02: F100 frequency model

Brief description.

1. Predict genes using Prodigal

```bash
# make new directory for gene CDS fastas we'll refer to this as ${genes_dir}
mkdir ${genes_dir}
# loop through genome fasta files and predict genes for each with prodigal
for f in ${genomes_dir}/*; do n=`basename $f | cut -d. -f1`; prodigal -q -d ${genes_dir}/${n}.fnn -i $f; echo $f; done
# in our experience Prodigal tends to predict a lot of short genes. While some of these may be real, we think a lot of them are likely noise.
# Filter by gene length to remove tiny genes
for f in ${gene_dir}/*; do python 00d_Workflow_Scripts/02a_len_filter_genes.py -i $f; done
```

2. All vs all aai.rb in nucleotide mode

```bash
# make new directory for aai.rb reciprocal best match results we'll refer to this as ${rbm_dir}
mkdir ${rbm_dir}

# generate all vs all gene file name combinations list
f=(${genes_dir}/*)
for ((i = 0; i < ${#f[@]}; i++)); do for ((j = i + 1; j < ${#f[@]}; j++)); do echo ${f[i]} ${f[j]}; done; done > genes_filenames_allV.txt

# run aai.rb in nucleotide mode
# loop over genes_filenames_allV.txt file and run each genome pair combination
# x and y are the full file paths to genome 1 and genome 2 from genes_filenames_allV.txt. m and n are used to create the new filenames and should be adjusted to fit your particular naming scheme if needed.
while read p;
  do
    	x=`echo $p | cut -d' ' -f1`;
        m=`basename $x | cut -d. -f1`;
        y=`echo $p | cut -d' ' -f2`;
        n=`basename $y | cut -d. -f1`;
        aai.rb -1 ${x} -2 ${y} -N -R ${rbm_dir}/${m}-${n}.rbm;
  done < genes_filenames_allV.txt

# concatenate rbm results to a single file.
cat ${rbm_dir}/* > RBMs_allV.rbm

# clean up individual rbm files
# CAUTION - double the concatenated file is good before removing data
rm -r ${rbm_dir}/
```

*The all vs all genome pair computation with aai.rb can take a while. There are many ways to parallalize this step if you are working on a cluster and we recommend following the best practices for your particular cluster. One quick method, you can use the bash split command to break up the genes_filenames_allV.txt file into many smaller files, say like 50 lines per file and then run several instances of the while read loop as separate jobs.*

```bash
split -l 50 -a -d genes_filenames_allV.txt genes_filenames_allV_
for f in genes_filenames_allV_*; do (run pbs or sbatch script with the while read loop); done
```

3. Compute F100s

*replace ${my_species_name} with your species name (use underscore instead of space) and replace ${my_species} with the output prefix of your choice. This will create an important file ${my_species}_F100.tsv needed for the next steps.*

```bash
# Calculate F100 scores and plot RBM sequence similarity histograms
# to view script info/params
python 00d_Workflow_Scripts/02b_AAI_RBM_F100.py -h
# run with default
python 00d_Workflow_Scripts/02b_AAI_RBM_F100.py -i RBMs_allV.rbm -s ${my_species_name} -o ${my_species}
```

(OPTIONAL) set -p True to Build histograms of RBM sequence similarities for each genome pair

```bash
# this will create a pdf for each genome pair
python 00d_Workflow_Scripts/02b_AAI_RBM_F100.py -i RBMs_allV.rbm -s ${my_species_name} -o ${my_species} -p True

# create a directory to store the PDFs
mkdir 04_rbm_pdf
mv ${my_species}*.pdf 04_rbm_pdf
```

4. Compare your genomes to a model

To view model info/options:

```bash
python 00d_Workflow_Scripts/02c_f100_scatter_pyGAM.py -h
```

This step will create a PDF of your genomes on top of the models and write out a sig-pairs.tsv file that can be used in Excel, R or etc. to sort by the xx column and select interesting genome pairs. The sig-pairs file labels any genome pairs with F100 and ANI values outside the 95% confidence interval of the model as "Recombining," and anything within the confidence interval as Non-recombining. But please be aware this is just a simplified terminology. This method is ranking genome pairs with higher to lower F100 scores based on the ANI between the genome pair. Genes with 100% sequence similarity are a proxy for recent homologous recombination or strict evolutionary conservation. Any genome pair, even the points at the bottom of the confidence interval may still have some 100% similar genes and as such could possibly have undergone a small amount of recent homologous recombination. 

### Complete genomes model

*This model comes from 330 species with at least 10 complete genomes in the NCBI RefSeq database as of April 20, 2022. Replace ${my_species_complete_model} with a  file name prefix of your choosing.*

```bash
# Uncompress the Complete_Genome_Model_Data.tsv.zip file
unzip Complete_Genome_Model_Data.tsv.zip
# this input file is part of the git repository
python 00d_Workflow_Scripts/02c_f100_scatter_pyGAM.py -i Complete_Genome_Model_Data.tsv -i2 ${my_species}_F100.tsv -o ${my_species_complete_model}
```

### Simulated Neutral model

*This model comes from [SpecSim](link-to-spec-sim-github) introducing random single point mutations across genes to fit a gamma distribution of RBMs along an ANI gradient from 95%-100% ANI with a step size of 0.01 ANI and 10 genomes per step. Replace ${my_species_simulated_model} with a  file name prefix of your choosing.*

```bash
# Uncompress the Simulated_Neutral_Model_Data.tsv.zip file
unzip Simulated_Neutral_Model_Data.tsv.zip
# this input file is part of the git repository
python 00d_Workflow_Scripts/02c_f100_scatter_pyGAM.py -i Simulated_Neutral_Model_Data.tsv -i2 ${my_species}_F100.tsv -o ${my_species_simulated_model}
```

### Create a model from your own genomes

*You can also just use your genomes to build the model and then feed that output back into the script with your genomes on top. Replace ${my_species_custom_model} with a  file name prefix of your choosing.*

```bash
# For this step we use the input file generate as part of this STEP 02
python 00d_Workflow_Scripts/02c_f100_scatter_pyGAM.py -i ${my_species}_F100.tsv -i2 ${my_species}_F100.tsv -o ${my_species_custom_model}
```

# STEP 03: Investigate recombinant positions in specific genome pairs

STEP 02 generates output containing F100 data for all genome pairs. In this step we select specific genome pairs of interest from STEP 02 and take a closer look.

### Setup

1. Concatenate all gene CDS

```bash
cat ${genes_dir}/* > all_genes_CDS.fnn
```

2. create some directories that we'll need

*replace ${mmseqs_dir} with your own directory name*

```bash
mkdir ${mmseqs_dir}
```
### mmseqs

1. Create an mmseqs database using all_genes_CDS.fnn

*replace ${my_db} with your own database name - just pick a name, whatever you want to call it*

```bash
mmseqs createdb all_genes_CDS.fnn ${mmseqs_dir}/${my_db}
```

2. Cluster at 90% nucleotide ID

```bash
mmseqs cluster ${mmseqs_dir}/${my_db} ${mmseqs_dir}/DBclustered tempfiles --min-seq-id 0.90 --cov-mode 1 -c 0.5 --cluster-mode 2 --cluster-reassign
```

3. Add sequence identity and alignment information to clustering result

*Super high -e prevents sequences getting dropped from clusters. If you notice you are missing sequencing downstream try increasing -e even higher. See [MMseqs2 GitHub Issue](https://github.com/soedinglab/MMseqs2/issues/598)*

```bash
mmseqs align ${mmseqs_dir}/${my_db} ${mmseqs_dir}/${my_db} ${mmseqs_dir}/DBclustered ${mmseqs_dir}/DBaligned -e 1.0E30
```

3. Write mmseqs database to TSV format

```bash
mmseqs createtsv ${mmseqs_dir}/${my_db} ${mmseqs_dir}/${my_db} {mmseqs_dir}/DBaligned all_genes_CDS_aligned.tsv

# You can remove the tempfiles directory at this point
rm -r tempfiles
```

4. Create a binary matrix of genomes and gene clusters

Each genome is a column. Each gene cluster is row. 1 if the genome has a gene in the cluster else 0. This is used to identify core genes or accessory genes. Conserved genes are the core gene clusters with the least average sequence difference. You can proceed to STEP 05 after you have the binary matrix file. 

```bash
python 00d_Workflow_Scripts/03a_MMSeqsTSV-to-BinaryMatrix.py -i all_genes_CDS_aligned.tsv -o pangenome_matrix.tsv
```

*side quest: create two optional plots just for fun because you can once you have the binary matrix file.*

5. create pangenome model

```bash
# for script info/options
python 00d_Workflow_Scripts/03b_Pangenome_Calculate_Model_Plot.py -h

# with default settings
python 00d_Workflow_Scripts/03b_Pangenome_Calculate_Model_Plot.py -b pangenome_matrix.tsv -o pangenome_model -n my_species_name
```

6. create clustermap

```bash
# for script info/options
python 00d_Workflow_Scripts/03c_Clustermap_fromBinary.py -h

# with default settings
python 00d_Workflow_Scripts/03c_Clustermap_fromBinary.py -b pangenome_matrix.tsv -o pangenome_clustermap.pdf
```

# STEP 05: Investiage specific genome pairs recombinant positions

1. Create list of genes with pangenome category (conserved, core, accessory)

Writes a tsv file needed for next step plust a Histogram PDF of average gene distance within core gene clsuters. The bottom 5% are designated as conserved genes.

```bash
python 00d_Workflow_Scripts/04a_Get_Genes_Clusters_PanCat.py -b pangenome_matrix.tsv -c all_genes_CDS_aligned.tsv -o pancat_file.tsv
```

2. Analyze interesting genome pairs

Repeat this step as many times as you have genome pairs you're interested in.

```bash
# for script info/option
python 00d_Workflow_Scripts/04b_F100_distance_analysis.py -h

# with default settings
python 00d_Workflow_Scripts/04b_F100_distance_analysis.py -rbm RBMs_allV.rbm -PC pancat_file.tsv -cA ${genes_dir}/genomeA.fnn -cB ${genes_dir}/genomeB.fnn -gA ${genomes_dir}/genomeA.fna -gB ${genomes_dir}/genomeB.fna -o genomeA-genomeB -draft True
```
