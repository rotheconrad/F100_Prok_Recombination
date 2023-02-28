# F100 as a signal for prokaryote recombination

![F100 Recent Recombination Model with Complete level Genomes from NCBI RefSeq for 330 Species](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/01_Complete_Genomes_Model.png)

Contains the code and workflow for the F100 recombination project.

This workflow uses 100% sequence similarity between reciprocal best matched genes (RBMs) from a pair of genomes as a proxy for recent homologous recombination events. It assesses the frequency of 100% similar RBMs to total RBMs (F100) for each genome pair input by the user and compares that to the expected value based on aggregate-genome average nucleotide identity (ANI) of the genome pair compared to 1) a model of 330+ species genomes from NCBI that have at least 10 Complete genomes per species, 2) a model of simulated random neutral evolution genomes along a 95%-100% ANI gradient with zero recombination, or 3) an easy to build custom model from the users own colleciton of genomes. This workflow also evaluates the genomic positions of recent recombinant events for genome pairs and for one genome against many genomes. It can also create clusters of genomes with high recombination activity and it can perform hypthosesis testing of broad level functional gene annotations of recombinant vs. non-recombinant genes.

A collection of genomes in fasta format is all that is required as input to begin. This workflow was designed focused on genome collections from the same species (≥95% ANI) but it will work at broader or finer genome similarity groupings as long as some 100% RBMs exist between the genomes.

This workflow will yeild many intermediate files and several publication ready figures.

*Figures are publication quality PDF files, DataTables are tab separated value (tsv) files*

1. DataTable: All vs All genome pair fastANI results (PART 01, Step 02)
1. Figure: Shared fraction vs ANI scatter plot (PART 01, Step 02)
1. Figure: ANI distance hierarchical clustered heatmap (PART 01, Step 02)
1. Fasta: CDC gene predictions from Prodigal (PART 02, Step 01)
1. Figures: Histogram of gene length distributions (PART 02, Step 01)
1. DataTable: RBM sequence similarities for each genome pair (PART 02, Step 02)
1. DataTable: F100 score for each genome pair (PART 02, Step 03)
1. Figure: Histogram of RBM sequence similarity for each genome pair (PART 02, Step 03)
1. Figure: F100 vs ANI with various GAM models (PART 02, Step 04)
1. DataTable: F100 vs ANI data, confidence interval and p value for each genome pair (PART 02, Step 04)
1. Figure: F100 distance hierarchical clustered heatmap (PART 02, Step 05)
1. DataTable: Gene clusters from MMSeqs2 clustered and aligned (PART 03, Step 06)
1. DataTable: Presence/Absence binary matrix of genomes and gene cluster. (PART 03, Step 07)
1. Figure: Pangenome model figures (PART 03, Step 08-09)
1. DataTable: Genes assigned to pangenome classes: Conserved, Core, Accessory (PART 04, Step 01)
1. Figure: Recombinant gene class and location between genome pairs (PART 04, Step 02)
1. Figure: Stats test recombinant gene spacing vs. Poisson distribution (PART 04, Step 02)

See 01a_Building_Complete_Genomes_Model.txt for detailed notes/methods used to build the model with NCBI Complete genomes. Code referenced in these notes can be found in 00b_Complete_Genomes_Model directory.

See 01b_Building_Simulated_Genomes_Model.txt for detailed notes/methods used to build the model with simulated genomes. Code referenced in these notes can be found in 00c_Simulated_Genomes_Model directory.

## Required dependencies

- [FastANI](https://github.com/ParBLiSS/FastANI)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [EggNog Mapper](https://github.com/eggnogdb/eggnog-mapper) (optional for functional annotation hypothesis testing)
- [COGclassifier](https://github.com/moshi4/COGclassifier/) (optional alternative to EggNog Mapper)
- [Enveomics Collection](http://enve-omics.ce.gatech.edu/enveomics/docs) (for aai.rb)
- [Ruby](https://www.ruby-lang.org/en/) (for enveomics)
- [Python](https://www.python.org/) (for all custom code in this workflow)

#### References

1. Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nature communications. 2018 Nov 30;9(1):1-8.
1. Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics. 2010 Dec;11(1):1-1.
1. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC bioinformatics. 2009 Dec;10(1):1-9.
1. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology. 2017 Nov;35(11):1026-8.
1. Rodriguez-R LM, Konstantinidis KT. The enveomics collection: a toolbox for specialized analyses of microbial genomes and metagenomes. PeerJ Preprints; 2016 Mar 27.
1. Flanagan D, Matsumoto Y. The Ruby Programming Language: Everything You Need to Know. " O'Reilly Media, Inc."; 2008 Jan 25.
1. Sanner MF. Python: a programming language for software integration and development. J Mol Graph Model. 1999 Feb 1;17(1):57-61.

## Required packages for Python

- [pandas](https://pandas.pydata.org/) 
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- [datashader](https://datashader.org/)
- [pygam](https://pygam.readthedocs.io/)

*Python and all packages can be easily installed with conda or pip. Prodigal, BLAST+ and MMseqs2 can also be installed easily with [Conda](https://docs.conda.io/en/latest/miniconda.html). Just search "conda install name"*

#### References

1. Van Rossum G, Drake FL. Python 3 Reference Manual. Scotts Valley, CA: CreateSpace; 2009.
1. McKinney W, others. Data structures for statistical computing in python. In: Proceedings of the 9th Python in Science Conference. 2010. p. 51–6.
1. Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, et al. Array programming with NumPy. Nature. 2020;585:357–62.
1. Virtanen P, Gommers R, Oliphant TE, Haberland M, Reddy T, Cournapeau D, et al. SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods. 2020;17:261–72.
1. Hunter JD. Matplotlib: A 2D graphics environment. Computing in science & engineering. 2007;9(3):90–5.
1. Waskom ML. Seaborn: statistical data visualization. Journal of Open Source Software. 2021 Apr 6;6(60):3021.
1. James A. Bednar, Jim Crist, Joseph Cottam, and Peter Wang (2016). "Datashader: Revealing the Structure of Genuinely Big Data", 15th Python in Science Conference (SciPy 2016).
1. Servén D., Brummitt C. (2018). pyGAM: Generalized Additive Models in Python. Zenodo. DOI: 10.5281/zenodo.1208723

# PART 01: Prepare your genomes

This workflow is intended for a collection of genomes belonging to the same species. Start with your genome files in fasta format in their own directory. We will refer to this directory as ${genomes_dir}.

### Step 01: Rename fasta deflines

Rename the fasta deflines of your genome files. This is necessary to ensure all contigs (or chromosomes/plasmids) in your genome files follow the same naming format for downstream processing. 

Because genomes downloaded from NCBI follow a typical naming convention of, "GCF_000007105.1_ASM710v1_genomic.fna," the default behavior of this script is to cut the third underscore position ("\_") and use it as a prefix for renaming the fasta deflines in numeric consecutive order.

So with default settings the script will cut "ASM710v1" from filename "GCF_000007105.1_ASM710v1_genomic.fna" and rename the fasta deflines (Contigs/Scaffolds/Chromosomes) as:

>\>ASM710v1_1  
>AATGGATCAGTCCGCCGACCGCGCCTGGAACGAATGTCTCGACATCATCCGGGACAATGT...  
>\>ASM710v1_2  
>GAGCCGCCAGAGCTTCACGACCTGGTTTGAGCCGCTGGAGGCCCACTCCTTGGAGGACGA...  
>\>ASM710v1_n  
>GGACGACCTGCGCAAGCTGACGATCCAACTTCCGAGCCGGTTTTACTACGAGTGGATTGA...  

This step requires Python.

Input: genome fasta files in ${my_genomes} directory

Output: overwrites genome fasta files with new names 

To use the renaming script on all files in a directory with default setting:

```bash
for f in ${genomes_dir}/*; do python 00d/Workflow_Scripts/01a_rename_fasta.py -i $f; done
```

Alternatively, the user can input their own desired prefix using the "-p" flag in which case the input filename is ignored. Replace "${name}" with anything you want:

```bash
for f in ${genomes_dir}/*; do name=`echo basename $f | cut -d_ -f3`; python 00d/Workflow_Scripts/01a_rename_fasta.py -i $f -p ${name}; done
```

So with -p my_genome the script will output:

>\>my_genome_1  
>AATGGATCAGTCCGCCGACCGCGCCTGGAACGAATGTCTCGACATCATCCGGGACAATGT...  
>\>my_genome_2  
>GAGCCGCCAGAGCTTCACGACCTGGTTTGAGCCGCTGGAGGCCCACTCCTTGGAGGACGA...  
>\>my_genome_n  
>GGACGACCTGCGCAAGCTGACGATCCAACTTCCGAGCCGGTTTTACTACGAGTGGATTGA...  

### Step 02: Inspect genome similarity

(OPTIONAL) Check shared genome fraction vs ANI of your genomes using fastANI.

It is a good idea to know how similar your genomes are to each other. Sometimes you may split your genomes into two or more groups, or remove a few outlier genomes based on their ANI distribution. This can be done manually as needed by removing the desired fasta files from ${genomes_dir} or partition files into separate directories after reviewing the fastANI plot(s).

To generate the fastANI plots, we want a single file with all vs. all genome pair results from fastANI. There are many ways to achieve this depending on how many genomes you have and if you're running fastANI locally or on a cluster. Here is a simple example using the many to many option (consult the fastANI documentation and/or your clusters best practices for other options - I frequently use the one to many mode and distribute the work across multiple instances and node. At the end I simply concatenate the outputs for the same results as the many to many mode).

This step requires fastANI.

Input: genome fasta files in ${my_genomes} directory

Output: all vs. all fastANI tsv outfile file saved as fastANI_allV.ani

```bash
# make a list of your fasta files with the path
for f in ${genomes_dir}/*; do echo $f; done > genome_file_list.txt

# run fastANI in many to many mode
fastANI --ql genome_file_list.txt --rl genome_file_list.txt -o fastANI_allV.ani
```

Once you have the all vs. all fastANI output, here is a script to easily create a scatter plot of the results. The shared genome fraction is the number of shared fragments between two genomes divided by the total number of fragments from the larger genome as determined by fastANI. ANI is the genome-aggregate average nucleotide identity as determined by fastANI. Species-like genome clusters are generally ≥95% ANI between them but it may be interesting to investigate the differences between groupings of genomes you find at any ANI level.

```bash
# look at script details and options
python 00d_Workflow_Scripts/01b_fastANI_scatter_pyGAM.py -h
# run with default settings
python 00d_Workflow_Scripts/01b_fastANI_scatter_pyGAM.py -i fastANI_allV.ani -s test_species -o fastANI_allV_sharedfrac_ANI.pdf -m True
```

![Shared fraction vs. ANI](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/fastANI_allV_sharedfrac_ANI.png)

And here is another script to easily create a hierarchical clustered heatmap of the results. A square ANI distance (100 - ANI) matrix is generated from the all vs. all fastANI output. The result shows which genomes are the most similar to each other. (future addition: clustering algorithm to automatically partition the distance matrix into clusters)

```bash
# look at script details and options
python 00d_Workflow_Scripts/01c_fastANI_clustermap.py -h
# run with default settings
python 00d_Workflow_Scripts/01c_fastANI_clustermap.py -i fastANI_allV.ani -o fastANI_allV.pdf
```

![ANI Distance Clustered Heatmap](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/fastANI_allV_heatmap.png)

# PART 02: Recombinant genomes analysis

In this section we identify which genome pairs from your species have the most or least amount of recent horizontal gene transfer. We use sequence similarity from reciprocal best blast matches (RBMs) of the genes between two genomes to calculate the frequency of 100% identical genes in the genome (F100). F100 is the number of RBMs with 100% sequence similarity divided by the total RBMs between two genomes. Using the F100 as a signal for recent recombination events, we fit a generalized additive model (GAM) to a set of data (1. Complete genomes from NCBI, 2. Simulated neutral evoltuion genomes without recombination, or 3. a custom genome set provided by the user) to show the expected F100 per ANI of a genome pair. We also identify clusters of frequently recombining genome pairs by creating a hierarchical clustered heatmap from F100 scores as a distance metric matrix and use the HDBSCAN algorithm to partition this into clusters of highly recombining genomes.

### Step 01: Predict genes using Prodigal

We use the predicted CDS to calculate reciprocal best matches (RBMs) from which we derive the frequency of genes with 100% RBM sequence similarity (F100) and we also use the genes for clustering in the pangenome analysis and to get functional annotation information. We need the nucleotide sequence in fasta format to calculate RBMs and F100 and also for gene clustering. We need the amino acid sequence in fasta format for functional annotation. Prodigal also writes a gbk or gff file as its main output. We won't use the gbk/gff output in this workflow, but it is sometimes useful to have on hand so we'll place them in their own directory.

Prodigal has a habit of sometimes predicting very short genes and/or also very long genes which are suspicious and likely have high error rate potential. Since we are more interested in the average gene behavior for this analysis, which is about 1000bp in length, we have a filter script to remove predicted genes that are very short or long. You can skip this step if you don't think it is necessary for your genomes. The script will write a histogram of the gene length distribution with quantile markers, and output how many genes it removed, which genes, and their lengths. The default to remove genes ≤ 100 base pairs and ≥ 8000 base pairs. Use the -aa True flag for amino acid sequence files. The default min and max lengths can be changed with the -min and max flags.

This step requires Prodigal and Python with Numpy and Matplotlib packages.

Input: Genomes fasta files

Output: CDS fasta files of nucleotide and amino acid sequence. Histogram gene length distribution PDFs.

```bash
# make new directory for gene CDS fastas we'll refer to this as ${genes_dir}
mkdir ${genes_dir_fnn} ${genes_dir_faa} ${genes_dir_gff}
# loop through genome fasta files and predict genes for each with prodigal
for f in ${genomes_dir}/*; do name=`basename $f | cut -d. -f1`; prodigal -q -f gff -o ${genes_dir_gff}/${name}.gff -d ${genes_dir_fnn}/${name}.fnn -a ${genes_dir_faa}/${name}.faa -i $f; echo $f; done
# in our experience Prodigal tends to predict a lot of short genes. While some of these may be real, we think a lot of them are likely noise.
# filter sequences. use -h for options
python 00d_Workflow_Scripts/02a_len_filter_genes.py -h
# Filter nucleotide sequences
for f in ${genes_dir_fnn}/*.fnn; do python 00d_Workflow_Scripts/02a_len_filter_genes.py -i $f; done
# Filter amino acid sequences
for f in ${genes_dir_faa}/*.faa; do python 00d_Workflow_Scripts/02a_len_filter_genes.py -i $f -aa True; done
```

![Gene Length Distribution](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/geneLenDist.png)

### Step 02: All vs all aai.rb in nucleotide mode

This step calculates the RBMs for each genome pair that we need to compute the F100 scores. We us the enveomics collections' aai.rb script in nucleotide mode for this step. aai.rb uses blast-plus to compute the one-way and two-way gene alignments and the best reciprocal best match hits for gene pairs (RBMs). We filter the RBM alignments to remove spurious short high identity sequence alignments (alignments alignment length / gene length) ≥ 0.50 (e.g. an RBM gene pair need to have at least a 50% alignment across the length of the shorter gene sequence).

This step requires the aai.rb script from the enveomics collection, Ruby, and Blast-plus.

Input: CDS fasta files of nucleotide and amino acid sequence

Output: Tabular blast table with filtered RBMs (RBMs_alV.rbm)

```bash
# make new directory for aai.rb reciprocal best match results we'll refer to this as ${rbm_dir}
mkdir ${rbm_dir}

# generate all vs all gene file name combinations list
f=(${genes_dir_fnn}/*)
for ((i = 0; i < ${#f[@]}; i++)); do for ((j = i + 1; j < ${#f[@]}; j++)); do echo ${f[i]} ${f[j]}; done; done > genes_filenames_allV.txt

# run aai.rb in nucleotide mode
# loop over genes_filenames_allV.txt file and run each genome pair combination
# x and y are the full file paths to genome 1 and genome 2 from genes_filenames_allV.txt. m and n are used to create the new filenames and should be adjusted to fit your particular naming scheme if needed.
while read p;
  do
    	x=`echo $p | cut -d' ' -f1`;
        m=`basename $x | cut -d. -f1`;
        y=`echo $p | cut -d' ' -f2`;
        name=`basename $y | cut -d. -f1`;
        aai.rb -1 ${x} -2 ${y} -N -R ${rbm_dir}/${m}-${name}.rbm -L 0.5;
  done < genes_filenames_allV.txt

# concatenate rbm results to a single file.
cat ${rbm_dir}/* > RBMs_allV.rbm

# clean up individual rbm files
# CAUTION - double check the concatenated file is good before removing data
rm -r ${rbm_dir}/
```

*The all vs all genome pair computation with aai.rb can take a while. There are many ways to parallalize this step if you are working on a cluster and we recommend following the best practices for your particular cluster. One quick method, you can use the bash split command to break up the genes_filenames_allV.txt file into many smaller files, say like 50 lines per file and then run several instances of the while read loop as separate jobs.*

```bash
split -l 50 -a -d genes_filenames_allV.txt genes_filenames_allV_
for f in genes_filenames_allV_*; do (run pbs or sbatch script with the while read loop); done
```

### Step 03: Compute F100s

This step calculates the F100 for each genome pair using the outputs from aai.rb in the previous step. It can write out histograms of the RBM sequence identity distribution.

*replace ${my_species} with the output prefix of your choice. This will create an important file ${my_species}_F100.tsv needed for the next steps.*

This step requires Python with Matplotlib, and Seaborn packages.

Input: RBMs_alV.rbm

Output: ${my_species}\_F100.tsv

```bash
# Calculate F100 scores and plot RBM sequence similarity histograms
# to view script info/params
python 00d_Workflow_Scripts/02b_AAI_RBM_F100.py -h
# run with default
python 00d_Workflow_Scripts/02b_AAI_RBM_F100.py -i RBMs_allV.rbm -o ${my_species}
```

(OPTIONAL) set -p True to Build histograms of RBM sequence similarities for each genome pair. The x-axis is automatically calculated to fit the data so if the values extend low down into the 70s etc it means at least one RBM has a low sequence similarity but these typically aren't visible in the plot due to the high counts above 95.

```bash
# this will create a pdf for each genome pair
python 00d_Workflow_Scripts/02b_AAI_RBM_F100.py -i RBMs_allV.rbm -s ${my_species_name} -o ${my_species} -p True

# create a directory to store the PDFs
mkdir 04_rbm_pdf
mv ${my_species}*.pdf 04_rbm_pdf
```

![Histogram of gene RBMs](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/test_genome01-genome02.png)

### Step 04: Compare your genomes to a model

To view model info/options:

```bash
python 00d_Workflow_Scripts/02c_f100_scatter_pyGAM.py -h
```

This step shows which genome pairs have a higher frequency of recently recombining genes (F100) expected for their ANI range compared to models built from other collections of other genomes.

This step will create a PDF of your genomes on top of the models and write out a sig-pairs.tsv file that can be used in Excel, R or etc. to sort by the xx column and select interesting genome pairs. The sig-pairs file labels any genome pairs with F100 and ANI values outside the 95% confidence interval of the model as "Recombining," and anything within the confidence interval as Non-recombining. But please be aware this is just a simplified terminology. This method is ranking genome pairs with higher to lower F100 scores based on the ANI between the genome pair. Genes with 100% sequence similarity are a proxy for recent homologous recombination or strict evolutionary conservation. Any genome pair, even the points at the bottom of the confidence interval may still have some 100% similar genes and as such could possibly have undergone a small amount of recent homologous recombination. 

This step requires Python with Numpy, Pandas, Matplotlib, Seaborn, Datashader, and PyGAM packages.

Input: ${my_species}\_F100.tsv from Step 02 and one of the three model tsv files from the options below.

Output: 2 files: ${my_species}\_complete_model_sig-pairs.tsv and ${my_species}\_complete_model_GAMplot.pdf

#### Option 01: Complete genomes model

*This model comes from 330 species with at least 10 complete genomes in the NCBI RefSeq database as of April 20, 2022. Replace ${my_species_complete_model} with a  file name prefix of your choosing.*

```bash
# Uncompress the Complete_Genome_Model_Data.tsv.zip file
unzip Complete_Genome_Model_Data.tsv.zip
# this input file is part of the git repository
python 00d_Workflow_Scripts/02c_f100_scatter_pyGAM.py -i Complete_Genome_Model_Data.tsv -i2 ${my_species}_F100.tsv -o ${my_species}_complete_model
```

![Your data on top of the complete genomes model](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_complete_model_GAMplot.png)

#### Option 02: Simulated Neutral model

*This model comes from [SpecSim](link-to-spec-sim-github) introducing random single point mutations across genes to fit a gamma distribution of RBMs along an ANI gradient from 95%-100% ANI with a step size of 0.01 ANI and 10 genomes per step. Replace ${my_species_simulated_model} with a  file name prefix of your choosing. I plan to make a github for this code and users can generate their own sets of simulated genomes by tweaking various parameters. Then follow Option 03 for building a custom model.*

```bash
# Uncompress the Simulated_Neutral_Model_Data.tsv.zip file
unzip Simulated_Neutral_Model_Data.tsv.zip
# this input file is part of the git repository
python 00d_Workflow_Scripts/02c_f100_scatter_pyGAM.py -i Simulated_Neutral_Model_Data.tsv -i2 ${my_species}_F100.tsv -o ${my_species_simulated_model}
```

![Your data on top of the simulated neutral genomes model](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_simulated_model_GAMplot.png)

#### Option 03: Custom model

*You can build the model using your species genomes and compare them to themselves, or build any other custom model from any set of genomes by repeating this pipeline with that set of genomes. Replace ${my_species_custom_model} with a  file name prefix of your choosing.*

```bash
# For this step we use the input file generate in PART 02, Step 03 twice, or you can generate two diffent F100.tsv files for different genome sets to create your own custom genome model.
python 00d_Workflow_Scripts/02c_f100_scatter_pyGAM.py -i ${my_species}_F100.tsv -i2 ${my_species}_F100.tsv -o ${my_species_custom_model}
```

![Your data on top of a model build from your own data](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_custom_model_GAMplot.png)

### Step 05: Clusters of the most frequently recombining genomes

This step reads in the ${my_species}\_F100.tsv file from Step 02 and generates a square F100 distance (1 - F100) matrix then creates a hierarchical clustered heatmap figure saved as a PDF. (future addition: clustering algorithm to automatically partition the distance matrix into clusters)

This step requires Python with Pandas, Matplotlib, and Seaborn packages.

Input: ${my_species}\_F100.tsv from Step 02.

Output: Clustered heatmap PDF 

```bash
python  00d_Workflow_Scripts/02d_F100_clustermap.py -i ${my_species}_F100.tsv -o ${my_species}_F100.pdf
```

![F100 Distance Heatmap](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_F100.png)

# PART 03: Recombinant Genes Analysis

STEP 02 generates output containing F100 data for all genome pairs. In this step we investigate recombinant gene positions and annotions from specific genome pairs of interest, and we investigate recombinant gene positions and annotations for 1 genome to many genomes.

### Step 01: Concatenate all gene CDS

In this step we prepare to cluster and annotate our genes with MMSeqs2 and EggNog mapper (or COGclassifier) by concatenate predicted CDS from all genomes into a single file for nucleotide and a single file for amino acid sequence. 

```bash
cat ${fnn_genes_dir}/* > all_genes_CDS.fnn
cat ${faa_genes_dir}/* > all_genes_CDS.faa
```

### Step 02: create some directories that we'll need

*replace ${mmseqs_dir} with your own directory name*

This step requires

Input:

Output:

```bash
mkdir ${mmseqs_dir}
```

### Step 03: Create an mmseqs database using all_genes_CDS.fnn

*replace ${my_db} with your own database name - just pick a name, whatever you want to call it*

This step requires

Input:

Output:

```bash
mmseqs createdb all_genes_CDS.fnn ${mmseqs_dir}/${my_db}
```

### Step 04: Cluster at 90% nucleotide ID

This step requires

Input:

Output:

```bash
mmseqs cluster ${mmseqs_dir}/${my_db} ${mmseqs_dir}/DBclustered tempfiles --min-seq-id 0.90 --cov-mode 1 -c 0.5 --cluster-mode 2 --cluster-reassign
```

### Step 05: Add sequence identity and alignment information to clustering result

*Super high -e prevents sequences getting dropped from clusters. If you notice you are missing sequencing downstream try increasing -e even higher. See [MMseqs2 GitHub Issue](https://github.com/soedinglab/MMseqs2/issues/598)*

This step requires

Input:

Output:

```bash
mmseqs align ${mmseqs_dir}/${my_db} ${mmseqs_dir}/${my_db} ${mmseqs_dir}/DBclustered ${mmseqs_dir}/DBaligned -e 1.0E60
```

### Step 06: Write mmseqs database to TSV format

This step requires

Input:

Output:

```bash
mmseqs createtsv ${mmseqs_dir}/${my_db} ${mmseqs_dir}/${my_db} ${mmseqs_dir}/DBaligned all_genes_CDS_aligned.tsv

# You can remove the tempfiles directory at this point
rm -r tempfiles
```

### Step 07: Write out cluster representative fasta file

This step requires

Input:

Output:

```bash
mmseqs createsubdb ${mmseqs_dir}/DBaligned ${mmseqs_dir}/${my_db} ${mmseqs_dir}/my_rep_seqs
mmseqs convert2fasta ${mmseqs_dir}/my_rep_seqs my_rep_seqs.fna
```

### Step 08: Annotate representative sequences

This step requires

Input:

Output:

```bash
# First we have to grab the amino acid sequence for the representative sequences
```

### Step 09: Create a binary matrix of genomes and gene clusters

Each genome is a column. Each gene cluster is row. 1 if the genome has a gene in the cluster else 0. This is used to identify core genes or accessory genes. Conserved genes are the core gene clusters with the least average sequence difference. You can proceed to STEP 05 after you have the binary matrix file. 

This step requires

Input:

Output:

```bash
python 00d_Workflow_Scripts/03a_MMSeqsTSV-to-BinaryMatrix.py -i all_genes_CDS_aligned.tsv -o pangenome_matrix.tsv
```

*side quest: create two optional plots just for fun because you can once you have the binary matrix file.*

### Step 10 (OPTIONAL): create pangenome model

This step requires

Input:

Output:

```bash
# for script info/options
python 00d_Workflow_Scripts/03b_Pangenome_Calculate_Model_Plot.py -h

# with default settings
python 00d_Workflow_Scripts/03b_Pangenome_Calculate_Model_Plot.py -b pangenome_matrix.tsv -o pangenome_model -n my_species_name
```

![Pangenome curve model of your genomes](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/pangenome_model_pangenome_curves.png)

### Step 11 (OPTIONAL): create clustermap

This step requires

Input:

Output:

```bash
# for script info/options
python 00d_Workflow_Scripts/03c_Clustermap_fromBinary.py -h

# with default settings
python 00d_Workflow_Scripts/03c_Clustermap_fromBinary.py -b pangenome_matrix.tsv -o pangenome_clustermap.pdf
```

![Gene clusters vs genomes from your data](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/pangenome_clustermap.png)

# PART 04: Investigate specific genome pairs recombinant positions

In this section, we tie all previous steps together and look at the types and distribution of recent plausible recombination events between specific genome pairs.

### Step 01: Assign pangenome class to genes

Create list of genes with pangenome category (conserved, core, accessory). This step writes a tsv file needed for next step plust a Histogram PDF of average gene distance within core gene clsuters. The bottom 5% are designated as conserved genes.

This step requires

Input:

Output:

```bash
python 00d_Workflow_Scripts/04a_Get_Genes_Clusters_PanCat.py -b pangenome_matrix.tsv -c all_genes_CDS_aligned.tsv -o pancat_file.tsv
```

![Average sequence distance within gene clusters](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/pancat_file.png)


### OPTIONAL STEP: Reorder contigs of draft genomes and MAGs

Contigs in draft genome and MAG assemblies are not typically aligned to any particular order. It can be helpful for the figure created in Step 02 to align the genome pair to each other or to a common reference genome. One way to do that is with [Mauve](https://darlinglab.org/mauve/mauve.html). See "Reordering contigs" section of the User Guide here: [https://darlinglab.org/mauve/user-guide/reordering.html](https://darlinglab.org/mauve/user-guide/reordering.html).

### Step 02: Run analysis for each genome pair of interest

Repeat this step as many times as you have genome pairs you're interested in.

cA and cB flags denote the predicted CDS in nucleotides fasta file from prodigal (using .fnn here) and gA and gB flags are for the genome fasta files (using .fna here).

This step requires

Input:

Output:

```bash
# for script info/option
python 00d_Workflow_Scripts/04b_F100_distance_analysis.py -h

# with default settings
python 00d_Workflow_Scripts/04b_F100_distance_analysis.py -rbm RBMs_allV.rbm -PC pancat_file.tsv -cA ${genes_dir}/genomeA.fnn -cB ${genes_dir}/genomeB.fnn -gA ${genomes_dir}/genomeA.fna -gB ${genomes_dir}/genomeB.fna -o genomeA-genomeB
```

This script writes a tsv file for each genome with columns: Genome, Gene, F100, PanCat, Start, Stop, Strand, Width.
The gene name should contain the genome identifier_contigNumberFromAssembly_geneNumberOnContig. The F100 column is assigned a 0 or 1. 1 indicates the gene sequence has 100% identity with its corresponding RBM in the other genome and thus a candidate for recent homoloug recombination. A 0 indicates the gene does not have 100% sequence identity with its RBM. The PanCat column indicates the pangenome class assigned to the gene. Start and stop positions are relative to the genome with all contigs concatenated in fasta file order. Strand indicates the strand sense (1) or antisense (-1) the gene is predicted on. Width indicates the gene length or distance between the start and stop positions.

The first figure labeled as \_genomes.pdf shows the location of recombinant genes on the two genomes labeled by pangenome class (Conserved, Core, Accessory, or non-recombinant). In this instance, non-recombinant indicates less than 100% sequence similarity between two genes and thus a recent recombination event involving the gene pair in question is unlikely.

![Recombinant gene positions in genome by pangenome class](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/genome01-genome02_genomes.png)

The distribution of recombinant gene location is also assessed and compared to a Null Model using the Poisson distribution as a proxy for evenly distributed recombination events across the genome. If the p-value is low and the and the k value is close to 1, the spacing of genes in that category does not fit a Poisson model for evenly spaced events. If the p-value is high, but the k value is still close 1 (and not to 0) this indicates a majority of the data falls inside the Poisson distribution but the overal shape of the distribution is still not a great fit. See the Q-Q plot to visualy why. The mean of the data is shown as a dashed line, a poisson model based on this mean is shown as a red curve, and the emperical data (the number of genes between events) is plotted as a histogram.

![Distance between events vs. Poisson](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/genome01-genome02_B_distance.png)

The Q-Q plot shows how the quantiles from your emperical data align with quantiles from the Poisson distribution fit to your data.

![Q-Q plot of your data to a Poisson fit](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/genome01-genome02_B_distance-qq.png)

Example of what a Q-Q plot of a good fit looks like:
![Good Q-Q plot example](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/Good_QQ_example.png)

### Step 03: Hypothesis test gene annotation bins

After the mmseq clustering step we end up with a fasta file containing representative genes for each gene cluster. Here we will annotate this file using either EggNog Mapper or COGclassifier, summarize the gene annotations by gene pair and by recombinant pangenome category, and then perform a parametric and non-parametric hypothesis test for recombinant genes against non-recombinant genes.

This step requires

Input:

Output:

(make sure we output the rep gene fasta from mmseq step)

- add EggNog Mapper to depency list (optional)
- add COGclassifier to depency list (optional)

- breif commands for option 1 or option 2

- write code to parse inputs

