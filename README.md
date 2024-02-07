# Investigating the role of homologous recombination in driving sequence-discrete units at the species and intra-species levels.

![F100 Recent Recombination Model with Complete level Genomes from NCBI RefSeq for 330 Species](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/01_Complete_Genomes_Model.png)

Contains the code and workflow for the F100 recombination project.

This workflow uses 100% sequence similarity between reciprocal best matched genes (RBMs) from a pair of genomes as a proxy for recent homologous recombination events. It assesses the frequency of 100% similar RBMs to total RBMs (F100) for each genome pair input by the user and compares that to the expected value based on aggregate-genome average nucleotide identity (ANI) of the genome pair compared to 1) a model of 330+ species genomes from NCBI that have at least 10 Complete genomes per species, 2) a model of simulated random neutral evolution genomes along a 95%-100% ANI gradient with zero recombination, or 3) an easy to build custom model from the users own colleciton of genomes. This workflow also evaluates the genomic positions of recent recombinant events for genome pairs and for one genome against many genomes. It can also create clusters of genomes with high recombination activity and it can perform hypthosesis testing of broad level functional gene annotations of recombinant vs. non-recombinant genes.

A collection of genomes in fasta format is all that is required as input to begin. This workflow was designed focused on genome collections from the same species (≥95% ANI) but it will work at broader or finer genome similarity groupings as long as some 100% RBMs exist between the genomes.

The steps are left separately so the user can more easily follow the workflow, and so individual steps can be more efficiently parallelized depending on the users system.

# Table of Contents

1. [Part 01: Genome Preparation](#part-01-genome-preparation)
	1. [Step 01: Rename fasta deflines](#step-01-rename-fasta-deflines)
	1. [Step 02: Inspect genome similarity](#step-02-inspect-genome-similarity)
	1. [Step 03: Assign clades, phylogroups, and genomovars](#step-03-assign-clades-phylogroups-and-genomovars)

1. [Part 02: Genome Analysis](#part-02-genome-analysis)
	1. [Step 01: Predict genes with Prodigal](#step-01-predict-genes-with-Prodigal)
	1. [Step 02: Compute Reciprocal Best Match Genes](#step-02-compute-reciprocal-best-match-genes)
	1. [Step 03: Compute F<sub>100</sub> scores](#step-03-compute-f100-scores)
	1. [Step 04: Compare User Genomes to GAM Models](#step-04-compare-user-genomes-to-gam-models)
	1. [Step 05: Identify Significant Outliers](#step-05-identify-significant-outliers)
	1. [Step 06: F<sub>100</sub> score Clustered Heatmap](#step-06-f100-score-clustered-heatmap)
	1. [Step 07: Identical gene fractions by groupings](#step-07-identical-gene-fractions-by-groupings)

1. [Part 03: Gene Analysis](#part-03-gene-analysis)
	1. [Step 01: Generate gene clusters with MMSeqs2](#step-01-generate-gene-clusters-with-mmseqs2)
	1. [Step 02: Annotate representative genes with EggNog Mapper or COGclassifier](#Step-02-annotate-representative-genes-with-eggNog-mapper-or-cogclassifier)
	1. [Step 03: Assign pangenome class to genes](#step-03-assign-pangenome-class-to-genes)
	1. [Step 04: Reorder-align contigs for MAGs, SAGs, and draft genomes](#step-04-reorder-align-contigs-for-MAGs-SAGs-and-draft-genomes)
	1. [Step 05: Explore genome pairs of interest](#step-05-explore-genome-pairs-of-interest)
	1. [Step 06: Explore one vs. many genome groups of interest](#step-06-explore-one-vs-many-genome-groups-of-interest)

1. [Software Dependencies](#software-dependencies)
1. [How to Cite](#how-to-cite)
1. [Future Improvements](#future-improvements)

## Data table and Figure Outputs

This workflow will yeild many intermediate files and several publication ready figures.

*Figures are publication quality PDF files, DataTables are tab separated value (tsv) files*

#### Part 01:
1. DataTable: All vs All genome pair fastANI results
1. Figure: Shared fraction vs ANI scatter plot
1. Figure: ANI distance hierarchical clustered heatmap

#### Part 02:
1. Fasta: Predicted CDS gene sequences from Prodigal
1. Figures: Histograms of gene length distributions
1. DataTable: RBM sequence similarities for each genome pair
1. DataTable: F100 score for each genome pair
1. Figure: Histogram of RBM sequence similarity for each genome pair
1. Figure: F100 vs ANI with various GAM models
1. DataTable: F100 vs ANI data, confidence interval and p value for each genome pair
1. Figure: F100 distance hierarchical clustered heatmap

#### Part 03:
1. DataTable: Gene clusters from MMSeqs2
1. Fasta: Representative sequence fasta for each gene cluster
1. DataTable: Presence/Absence binary matrix of genomes and gene cluster
1. DataTable: tsv gene list for Coinfinder input
1. Figure: Pangenome curve model
1. Figure: Pangenome clustered heatmap
1. DataTable: Gene annotations
1. DataTable: Genes assigned to pangenome classes: Conserved, Core, Accessory, Specific
1. Figure: Histogram of average within gene cluster sequence distance for Core genes
- Genome pairs
1. Figure: genome pairs: Recombinant gene position by pangenome class
1. Figure: genome pairs: Distance between recombination events distribution test
1. Figure: genome pairs: Recombinant vs. Non-recombinant gene annotation test
1. Figure: genome pairs: Sequence identity of RBMs vs. genome position
1. DataTable: Gene RBM info, position info, annoation info
- One genome to many genomes
1. Figure: genome group: Recombinant gene position by pangenome class
1. Figure: genome group: Distance between recombination events distribution test
1. Figure: genome group: Recombinant vs. Non-recombinant gene annotation test
1. Figure: genome group: Sequence identity of RBMs vs. genome position
1. Figure: genome group: Core vs total recombinant positions rarefaction curve
1. DataTable: Gene RBM info, position info, annoation info
1. DataTable: RBM Matrix of gene/genome recombinant sites
1. Figure: Rarefaction curve of recombinant sites per genome



# PART 01: Genome Preparation

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

([Return to Table of Contents](#table-of-contents))

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


As you can see in the figures above, there is a small cluster of genomes between 3-5% different than the other genomes. In other words we have two fairly distinct clusters. This indicates we may want to split them into two separate analyses. For now, we will leave them together.

([Return to Table of Contents](#table-of-contents))

### Step 03: Assign clades, phylogroups, and genomovars

A Note on preparing clade, phylogroup, genomovar assignments.

The shared fraction vs. ANI plot can be used to look for values between genome clusters such as around 99.2-99.8% ANI for the genomovar gap or below for phylogroup or clade values. Tools such as MiGA Clades or dREP can be used to sort genomes into clades and genomovars, or the ANI Distance Clustered Heatmap Genomovar version can be used to make manual assignments.

Additional, many well studied species already phylogroup or sequencing typing tools such clermontyping etc. that can be used to assign genomes to groups.

Prepare a metadata file with the genome name in column 1 and additional metadata values in the remaining columns. Any category may be used such as location, sample, niche, sequence type, phylogroup, clade, and genomovar.

The Heatmap figure with metadata is useful to determine genomes and genome groups of interest to investigate in more detail in step 3.

INCLUDE METADATA FILE EXAMPLE TABLE HERE

([Return to Table of Contents](#table-of-contents))

# PART 02: Genome Analysis

In this section we identify which genome pairs from your species have the most or least amount of recent horizontal gene transfer. We use sequence similarity from reciprocal best blast matches (RBMs) of the genes between two genomes to calculate the frequency of 100% identical genes in the genome (F100). F100 is the number of RBMs with 100% sequence similarity divided by the total RBMs between two genomes. Using the F100 as a signal for recent recombination events, we fit a generalized additive model (GAM) to a set of data (1. Complete genomes from NCBI, 2. Simulated neutral evoltuion genomes without recombination, or 3. a custom genome set provided by the user) to show the expected F100 per ANI of a genome pair. We also identify clusters of frequently recombining genome pairs by creating a hierarchical clustered heatmap from F100 scores as a distance metric matrix and use the HDBSCAN algorithm to partition this into clusters of highly recombining genomes.

### Step 01: Predict genes with Prodigal

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

([Return to Table of Contents](#table-of-contents))

### Step 02: Compute Reciprocal Best Match Genes

In this step we calculate all vs. all genomes reciprocal best match (RBM) genes for each genome pair. RBM genes are used to compute the F<sub>100</sub> score. We use a basic two way BLAST search to identify matching best hits between the genes of two genomes. The script filters the RBM alignments to remove spurious short high identity sequence alignments (alignments alignment length / gene length) ≥ 0.50 (e.g. an RBM gene pair need to have at least a 50% alignment across the length of the shorter gene sequence). This script retains ties for RBMs.

This step requires Python and BLAST+.

Input: CDS fasta files of nucleotide sequences

Output: Tabular blast table with filtered RBMs (RBMs_alV.rbm)

```bash
# make new directory for reciprocal best match results we'll refer to this as ${rbm_dir}
mkdir ${rbm_dir}

# generate all vs all gene file name combinations list
f=(${genes_dir_fnn}/*)
for ((i = 0; i < ${#f[@]}; i++)); do for ((j = i + 1; j < ${#f[@]}; j++)); do echo ${f[i]} ${f[j]}; done; done > genes_filenames_allV.txt

# loop over genes_filenames_allV.txt file and run each genome pair combination
# x and y are the full file paths to genome 1 and genome 2 from genes_filenames_allV.txt. m and n are used to create the new filenames and should be adjusted to fit your particular naming scheme if needed.
while read p;
  do
    	x=`echo $p | cut -d' ' -f1`;
        m=`basename $x | cut -d. -f1`;
        y=`echo $p | cut -d' ' -f2`;
        name=`basename $y | cut -d. -f1`;
        python 00d_Workflow_Scripts/02b_get_RBMs.py -g1 ${x} -g2 ${y} -o ${rbm_dir}/${m}-${name}.rbm;
  done < genes_filenames_allV.txt

# concatenate rbm results to a single file.
cat ${rbm_dir}/* > RBMs_allV.rbm

# clean up individual rbm files
# CAUTION - double check the concatenated file is good before removing data
rm -r ${rbm_dir}/
```

*The all vs all genome pair RBM computations can take a while. There are many ways to parallalize this step if you are working on a cluster and we recommend following the best practices for your particular cluster. One quick method, you can use the bash split command to break up the genes_filenames_allV.txt file into many smaller files, say like 50 lines per file and then run several instances of the while read loop as separate jobs.*

```bash
split -l 50 -a -d genes_filenames_allV.txt genes_filenames_allV_
for f in genes_filenames_allV_*; do (run pbs or sbatch script with the while read loop); done
```
([Return to Table of Contents](#table-of-contents))

### Step 03: Compute F<sub>100</sub> scores

This step calculates the F100 for each genome pair using the outputs from 02b_get_RBMs.py in the previous step. It can write out histograms of the RBM sequence identity distribution.

*replace ${my_species} with the output prefix of your choice. This will create an important file ${my_species}_F100.tsv needed for the next steps.*

This step requires Python with Matplotlib and Seaborn packages.

Input: RBMs_alV.rbm

Output: ${my_species}\_F100.tsv

```bash
# Calculate F100 scores and plot RBM sequence similarity histograms
# to view script info/params
python 00d_Workflow_Scripts/02c_get_F100.py -h
# run with default
python 00d_Workflow_Scripts/02c_get_F100.py -i RBMs_allV.rbm -o ${my_species}
```

(OPTIONAL) set -p True to Build histograms of RBM sequence similarities for each genome pair. The x-axis is automatically calculated to fit the data so if the values extend low down into the 70s etc it means at least one RBM has a low sequence similarity but these typically aren't visible in the plot due to the high counts above 95.

```bash
# this will create a pdf for each genome pair
python 00d_Workflow_Scripts/02c_get_F100.py -i RBMs_allV.rbm -s ${my_species_name} -o ${my_species} -p True

# create a directory to store the PDFs
mkdir 04_rbm_pdf
mv ${my_species}*.pdf 04_rbm_pdf
```

![Histogram of gene RBMs](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/test_genome01-genome02.png)

([Return to Table of Contents](#table-of-contents))

### Step 04: Compare User Genomes to GAM Models

To view model info/options:

```bash
python 00d_Workflow_Scripts/02d_f100_scatter_pyGAM.py -h
```

This step shows which genome pairs have a higher frequency of recently recombining genes (F100) expected for their ANI range compared to models built from other collections of other genomes.

This step will create a PDF of your genomes on top of the models and write out a sig-pairs.tsv file that can be used in Excel, R or etc. to sort by the xx column and select interesting genome pairs. The sig-pairs file labels any genome pairs with F100 and ANI values outside the 95% confidence interval of the model as "Recombining," and anything within the confidence interval as Non-recombining. But please be aware this is just a simplified terminology. This method is ranking genome pairs with higher to lower F100 scores based on the ANI between the genome pair. Genes with 100% sequence similarity are a proxy for recent homologous recombination or strict evolutionary conservation. Any genome pair, even the points at the bottom of the confidence interval may still have some 100% similar genes and as such could possibly have undergone a small amount of recent homologous recombination. 

This step requires Python with Numpy, Pandas, Matplotlib, Seaborn, Datashader, and PyGAM packages.

Input: ${my_species}\_F100.tsv from Step 02 and one of the three model tsv files from the options below.

Output:
      1) ${my_species}\_complete_model_sig-pairs.tsv
      2) ${my_species}\_complete_model_GAMplot.pdf

#### Option 01: Complete genomes model (Subsampled model also available)

*This model comes from 330 species with at least 10 complete genomes in the NCBI RefSeq database as of April 20, 2022. Replace ${my_species_complete_model} with a  file name prefix of your choosing.*

See 01a_Building_Complete_Genomes_Model.txt for detailed notes/methods used to build the model with NCBI Complete genomes. Code referenced in these notes can be found in 00b_Complete_Genomes_Model directory.

The subsampled model filters the complete genome set for species with n genome pairs above 95% ANI with n ≥ 1000 and it randomly downsamples species with ≥ 5000 genome pairs above 95% ANI to n = 5000. The effect is a wider 95% confidence.

```bash
# Uncompress the Complete_Genome_Model_Data.tsv.zip file
unzip Complete_Genome_Model_Data.tsv.zip
# this input file is part of the git repository
python 00d_Workflow_Scripts/02d_f100_scatter_pyGAM.py -i Complete_Genome_Model_Data.tsv -i2 ${my_species}_F100.tsv -o ${my_species}_complete_model
```

![Your data on top of the complete genomes model](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_complete_model_GAMplot.png)

#### Option 02: Simulated Neutral model

*This model comes from [Population-Genome-Simulator](https://github.com/rotheconrad/Population-Genome-Simulator) which simulates a population of genomes by introducing random single point mutations across genes to fit a gamma distribution of RBMs along an ANI gradient from 95%-100% ANI with a step size of 0.01 ANI and 10 genomes per step. Replace ${my_species_simulated_model} with a  file name prefix of your choosing. Users can generate their own sets of simulated genomes by tweaking various parameters. Then follow Option 03 for building a custom model.*

```bash
# Uncompress the Simulated_Neutral_Model_Data.tsv.zip file
unzip Simulated_Neutral_Model_Data.tsv.zip
# this input file is part of the git repository
python 00d_Workflow_Scripts/02d_f100_scatter_pyGAM.py -i Simulated_Neutral_Model_Data.tsv -i2 ${my_species}_F100.tsv -o ${my_species_simulated_model}
```

![Your data on top of the simulated neutral genomes model](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_simulated_model_GAMplot.png)

#### Option 03: Custom model

*You can build the model using your species genomes and compare them to themselves, or build any other custom model from any set of genomes by repeating this pipeline with that set of genomes. Replace ${my_species_custom_model} with a  file name prefix of your choosing.*

```bash
# For this step we use the input file generate in PART 02, Step 03 twice, or you can generate two diffent F100.tsv files for different genome sets to create your own custom genome model.
python 00d_Workflow_Scripts/02d_f100_scatter_pyGAM.py -i ${my_species}_F100.tsv -i2 ${my_species}_F100.tsv -o ${my_species_custom_model}
```

![Your data on top of a model build from your own data](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_custom_model_GAMplot.png)

([Return to Table of Contents](#table-of-contents))

### Step 05: Identify Significant Outliers

In this step we'll use the sigpairs.tsv output files from any of the models in step 04 to build some bar plots looking at which groups from the meta data file are most represented above or below the GAM modeled F<sub>100</sub> score confidence intervals.

NEED TO ADD THE CORRECT SCRIPT AND COMMANDS AND A NEW FIGURE

```bash
python  00d_Workflow_Scripts/02e_F100_clustermap.py -i ${my_species}_F100.tsv -o ${my_species}_F100.pdf
```

![F100 Distance Heatmap](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_F100.png)

([Return to Table of Contents](#table-of-contents))

### Step 06: F<sub>100</sub> score Clustered Heatmap

This step reads in the ${my_species}\_F100.tsv file from Step 02 and generates a square F100 distance (1 - F100) matrix then creates a hierarchical clustered heatmap figure saved as a PDF. (future addition: clustering algorithm to automatically partition the distance matrix into clusters)

This step requires Python with Pandas, Matplotlib, and Seaborn packages.

Input: ${my_species}\_F100.tsv from Step 02.

Output: Clustered heatmap PDF 

```bash
python  00d_Workflow_Scripts/02e_F100_clustermap.py -i ${my_species}_F100.tsv -o ${my_species}_F100.pdf
```

![F100 Distance Heatmap](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_F100.png)

([Return to Table of Contents](#table-of-contents))

### Step 07: Identical gene fractions by groupings

This step is for the violin plots showing the identical gene fractions on the y-axis for various different one vs. many groupings on the x-axis. based on the metadata file.

NEED TO ADD THE CORRECT SCRIPT AND COMMANDS AND A NEW FIGURE

```bash
python  00d_Workflow_Scripts/02e_F100_clustermap.py -i ${my_species}_F100.tsv -o ${my_species}_F100.pdf
```

![F100 Distance Heatmap](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/my_species_F100.png)

([Return to Table of Contents](#table-of-contents))

# PART 03: Gene Analysis

STEP 02 generates output containing F100 data for all genome pairs. In this step we investigate recombinant gene positions and annotions from specific genome pairs of interest, and we investigate recombinant gene positions and annotations for 1 genome to many genomes.

### Step 01: Generate gene clusters with MMSeqs2

For this part of the analysis, we want to identify highly conserved genes, core genes, accessory genes, and genome specific genes to make our recent recombination position analysis more informative. To do this, we will use MMSeqs2 to cluster the nucleotide sequences for the predicted CDS into gene clusters based on 90% sequence similarity and 50% sequence alignment overlap. As a bonus, we can plot the permutational pangenome curve, and a clustermap of shared genomes across the pangenome.

This step requires MMSeqs2 and Python with Numpy, Pandas, LMFIT, , Matplotlib, and Seaborn packages.

Input: all_genes_CDS.fnn

Output:
      1) MMSeqs2 cluster file
      2) Cluster representative sequence fasta file
      3) Binary matrix for gene to genome presence/absence
      4) Pangenome model figure
      5) Pangenome clustermap

#### Concatenate all gene CDS to single file

In this step we prepare to cluster our representative genes with MMSeqs2 by concatenating the nucleotide sequences of the predicted CDS from all genomes into a single fasta file.

```bash
cat ${genes_dir_fnn}/*.fnn > all_genes_CDS.fnn
```

#### Create a directory that we'll use for MMSeq2 intermediates

*replace ${mmseqs_dir} with your own directory name*

```bash
mkdir ${mmseqs_dir}
```

#### Create an mmseqs database using all_genes_CDS.fnn

*replace ${my_db} with your own database name - just pick a name, whatever you want to call it*

```bash
mmseqs createdb all_genes_CDS.fnn ${mmseqs_dir}/${my_db}
```

#### Cluster at 90% nucleotide identity

```bash
mmseqs cluster ${mmseqs_dir}/${my_db} ${mmseqs_dir}/DBclustered tempfiles --min-seq-id 0.90 --cov-mode 1 -c 0.5 --cluster-mode 2 --cluster-reassign
```

#### Write mmseqs database to TSV format

```bash
mmseqs createtsv ${mmseqs_dir}/${my_db} ${mmseqs_dir}/${my_db} ${mmseqs_dir}/DBclustered all_genes_CDS_mmseqs_clusters.tsv
```

#### Write out cluster representative fasta file

```bash
mmseqs createsubdb ${mmseqs_dir}/DBclustered ${mmseqs_dir}/${my_db} ${mmseqs_dir}/my_rep_seqs
mmseqs convert2fasta ${mmseqs_dir}/my_rep_seqs my_rep_seqs.fnn
```

#### Cleanup tempory files

```bash
# You can remove the tempfiles directory at this point
rm -r tempfiles
```

#### Create a binary matrix of genomes and gene clusters

Each genome is a column. Each gene cluster is row. 1 if the genome has a gene in the cluster else 0. This is used to identify core genes or accessory genes. Conserved genes are the core gene clusters with the least average sequence difference. You can proceed to STEP 05 after you have the binary matrix file. 

```bash
python 00d_Workflow_Scripts/03a_MMSeqsTSV-to-BinaryMatrix.py -i all_genes_CDS_mmseqs_clusters.tsv -o pangenome_matrix.tsv
```

*side quest: create two plots just for fun because we can once we have the binary matrix file. The code for these side quest figures was developed for a [previous publication](https://doi.org/10.1038/s41396-021-01149-9).*

#### (OPTIONAL): create Coinfinder input file

```bash
python 00d_Workflow_Scripts/03a_coinfinder_format.py -i pangenome_matrix.tsv -o my_coinfinder_input.tsv
```

#### (OPTIONAL): create pangenome model

```bash
# for script info/options
python 00d_Workflow_Scripts/03b_Pangenome_Calculate_Model_Plot.py -h

# with default settings
python 00d_Workflow_Scripts/03b_Pangenome_Calculate_Model_Plot.py -b pangenome_matrix.tsv -o pangenome_model -n my_species_name
```

![Pangenome curve model of your genomes](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/pangenome_model_pangenome_curves.png)

#### (OPTIONAL): create clustermap

Genomes are the  columns and genes are the rows. A 1 (dark gray) indicates the gene cluster is present and a 0 (light gray) indicates a gene is absent in a genome.

```bash
# for script info/options
python 00d_Workflow_Scripts/03c_Clustermap_fromBinary.py -h

# with default settings
python 00d_Workflow_Scripts/03c_Clustermap_fromBinary.py -b pangenome_matrix.tsv -o pangenome_clustermap.pdf
```

![Gene clusters vs genomes from your data](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/pangenome_clustermap.png)

([Return to Table of Contents](#table-of-contents))

### Step 02: Annotate representative genes with EggNog Mapper or COGclassifier

In this step we prepare to annotate our representative genes with EggNog mapper (or COGclassifier) by concatenating the amino acid sequences of the predicted CDS from all genomes into a single fasta file. 

This step requires EggNog mapper or COGclassifier

Input: Representative gene fasta from mmseqs (my_rep_seqs.fnn)

Output: tsv file of gene annotations from EggNog or COGclassifier

#### Concatenate all amino acid sequence predicted CDS

```bash
cat ${genes_dir_faa}/*.faa > all_genes_CDS.faa
```

#### Retrieve amino acid sequence for representative genes

```bash
python 00d_Workflow_Scripts/03d_get_AA_reps_fasta.py -r my_rep_seqs.fnn -a all_genes_CDS.faa -o my_rep_seqs.faa
```

#### Annotate genes with EggNog Mapper

For EggNog we followed the installation instructions on their [GitHub](https://github.com/eggnogdb/eggnog-mapper) page to create a bacteria database for diamond. However, since the EggNog mapper output format is the same regardless of the database you use, you can create and use whatever database you'd like. Annotation can take several hours.

```bash
# create directory for EggNog Mapper
mkdir EggNog
# specify paths to EggNog databases on your system
db1="/path/to/eggnog-mapper/data"
db2="${db1}/bacteria.dmnd"
# run emapper
emapper.py --data_dir $db1 --dmnd_db $db2 -i my_rep_seqs.faa -o $outpre --output_dir EggNog/${my_annotations}
```

#### Annotate genes with COGclassifier

For COGclassifier we follow the installation instructions on their [GitHub](https://github.com/moshi4/COGclassifier/) page. Annotation can take several hours.

```bash
# create output directory
mkdir COGclass
# run COGclassifier
COGclassifier -i my_rep_seqs.faa -o COGclass
```

([Return to Table of Contents](#table-of-contents))

### Step 03: Assign pangenome class to genes

In this step we assign pangenome gene categories to our gene clusters. This step generates a .tsv file with columns: Gene_Name, Cluster_Name, Pangenome_Category, n/N, Distance, and a PDF hitsogram of average within cluster sequence distance of core gene clusters.

 - Gene_Name is the name of all individual genes in the pangenome
 - Cluster_Name is the representative gene of the assigned cluster
 - Pangenome_Category is one of: (Conserved, Core, Accessory, Specific)
 - n/N is number of genomes with gene category over total genomes in analysis
 - Distance is the average within cluster sequence distance from the representative gene.
 - Conserved genes are a subset (10%) of core genes with the least within cluster sequence distance.
 - Core genes are found in ≥ 90% of genomes in the analysis.
 - Specific genes are found in only 1 genome in the analysis.
 - Accessory genes are all other genes.

The histogram shows the distribution of average within cluster sequence distance. Sequence distance is calculated for each gene as 1 - sequence identity [0:1] compared to the cluster representative sequence. The vertical dashed red line shows the 0.10 quantile.

Input: pangenome_matrix.tsv, all_genes_CDS_mmseqs_clusters.tsv

Output: 1 tsv file, 1 pdf file

```bash
python 00d_Workflow_Scripts/03e_Get_Genes_Clusters_PanCat.py -b pangenome_matrix.tsv -m all_genes_CDS_mmseqs_clusters.tsv -r RBMs_allV.rbm -o pancat_file.tsv
```

![Average sequence distance within gene clusters](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/pancat_file.png)

([Return to Table of Contents](#table-of-contents))

### Step 04: Reorder-align contigs for MAGs, SAGs, and draft genomes

This step is optional.

Contigs in draft genome and MAG assemblies are not typically aligned to any particular order. It can be helpful for the figure created in Step 02 to align the genome pair to each other or to a common reference genome. One way to do that is with [Mauve](https://darlinglab.org/mauve/mauve.html). See "Reordering contigs" section of the User Guide here: [https://darlinglab.org/mauve/user-guide/reordering.html](https://darlinglab.org/mauve/user-guide/reordering.html).

([Return to Table of Contents](#table-of-contents))

### Step 05: Explore genome pairs of interest

In this section, we tie all previous steps together and look at the types and distribution of recent plausible recombination events between specific genome pairs.

This step has several outputs. It plots the positions of recombinant and non-recombinant genes by their genome coordinates separated by their pangenome class. It does some statitistical testing on the distance between recombinant events, and it plots recombinant vs. non-recombinant gene annotations by COG categories and performs Chi-square hypothesis testing.

This step requires Python with Numpy, Pandas, Scipy, StatsModels, MatPlotLib, and Seaborn packages.

#### Run analysis for each genome pair of interest

Repeat this step as many times as you have genome pairs you're interested in.

Input: 1) RBMs_allV.rbm 2) pancat_file.tsv 4) annotation file 3) genome pair fastas cA, cB, gA, gB

cA and cB flags denote the predicted CDS in nucleotides fasta file from prodigal (using .fnn here) and gA and gB flags are for the genome fasta files (using .fna here). For input 4) annotation file, use -ano EggNog/my_annotations.emapper.annotations for EggNog Mapper results or -ano COGclassifier/classifier_result.tsv for COGclassifier results. You only need one or the other annotation file, not both.

Output:

	1. tsv file
	1. recombinant gene position plot
	1. recombinant distribution plot
	1. recombinant annotations plot
	1. RBM sequence identity vs. genome position plot

The tsv file contians the file columns and data:

	- Genome: A or B as input by the user
	- Gene: Sequence identifier contains the genome_identifier_contigNumberFromAssembly_geneNumberOnContig.
	- PanCat:  column indicates the pangenome class assigned to the gene.
	- pID: Nucleotide sequence identity of pairwise gene alignment of the RBM.
	- REC: Assigned a 0 or 1. 1 indicates the gene sequence has ≥ REC% identity with its corresponding RBM in the other genome and thus a candidate for recent homologous recombination. A 0 indicates the gene does not. The REC threshold is controlled by the -rec parameter (default 99.8) and sets the threshold for an RBM to be considered recombinant or not. It affects the results of the gene annotations plot and chi-square hypothesis tests and it affects the recombinant positions plot and Poisson and Geometric distribution tests.
	- Recombinant: Recombinant or Non-Recombinant classification based on "REC"
	- Start: Start of gene in genome coordinates
	- Stop: End of gene in genome coordinates
	- Strand: 1 for positive or sense and -1 for negative or anti-sense strand of predicted CDS.
	- COG Category: Higher level COG category assignement based on COG.
	- COG: Single letter COG assignment.
	- Gene Annotation: Short gene name of the assigned gene annotation
	- Annotation Description: Long form gene name or description of the assigned gene annotation
	- Mismatch: Number of mismatches in the blast alignment for RBMs.
	- Width: Gene length. Distance between stop and start.

```bash
# for script info/option
python 00d_Workflow_Scripts/03f_Recombinant_pair_analysis.py -h

# with default settings
python 00d_Workflow_Scripts/03f_Recombinant_pair_analysis.py -rbm RBMs_allV.rbm -pc pancat_file.tsv -ano annotation_file -cA ${genes_dir}/genomeA.fnn -cB ${genes_dir}/genomeB.fnn -gA ${genomes_dir}/genomeA.fna -gB ${genomes_dir}/genomeB.fna -o genomeA-genomeB
```

#### Recombinant gene position by pangenome class

The first figure labeled as "\_genomes.pdf" shows the location of recombinant genes on the two genomes labeled by pangenome class (Conserved, Core, Accessory, or non-recombinant). In this instance, non-recombinant indicates less than 100% sequence similarity between two genes and thus a recent recombination event involving the gene pair in question is unlikely.

The percentages on the right for each gene class are 1) the percent of genes in that class out of the total genes (top number). They should add to 100%. And 2) number of mismatched base pairs / total base pairs * 100 where mismatched base pairs are calculated from the RBM gene alignments and total base pairs is the length of genes.

![Recombinant gene positions in genome by pangenome class](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/g01-g02_genomes.png)

#### Distance between recombination events distribution test

The distribution of recombinant gene location is also assessed and compared to a Null Model using the Poisson distribution (top left panel) as a proxy for evenly distributed recombination events across the genome or a Geometric distribution (top right panel). If the p-value is low and the k value is close to 1, the spacing of genes in that category does not fit the named distribution well. If the p-value is high, but the k value is still close 1 (and not to 0) this indicates a majority of the data falls inside the named distribution but the overal shape of the distribution is still not a great fit. See also the Q-Q plots (two lower panels) for a different perspective on distribution fit. The mean of the data is shown as a dashed line, the named model (Poisson left; Geometric right) fit to the data is shown as a red curve, and the emperical data (the number of genes between events) is plotted as a histogram in grey.

![Distance between events vs. Poisson or Geometric distribution](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/g01-g02_A_distance.png)

Here is an example of what a Q-Q plot of a good fit looks like:

![Good Q-Q plot example](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/Good_QQ_example.png)

#### Recombinant vs. Non-recombinant gene annotation test

This step also bins annotations by COG category and plots recombinant vs non-recombinant gene annotations. Performs Chi-square test of independence with post hoc chi-square tests if initial significance is found. Asterisks in the plot denote significant difference between a category after Benjamini/Hochberg multiple test correction. Detailed Chi-square results are printed to the screen.

Annotation categories are based on COG categories with specific labeling as follows:

	- X: Mobile
	- C: Metabolism 1
	- G: Metabolism 1
	- E: Metabolism 1
	- F: Metabolism 1
	- H: Metabolism 1
	- I: Metabolism 1
	- P: Metabolism 2
	- Q: Metabolism 2
	- J: Ribosomal (Translation ribosomal structure and biogenesis)
	- A: Information (Information Storage and Processing) 
	- K: Information (Information Storage and Processing) 
	- L: Information (Information Storage and Processing) 
	- B: Information (Information Storage and Processing) 
	- D: Cellular (Cellular Processing and Signaling)
	- Y: Cellular (Cellular Processing and Signaling)
	- V: Cellular (Cellular Processing and Signaling)
	- T: Cellular (Cellular Processing and Signaling)
	- M: Cellular (Cellular Processing and Signaling)
	- N: Cellular (Cellular Processing and Signaling)
	- Z: Cellular (Cellular Processing and Signaling)
	- W: Cellular (Cellular Processing and Signaling)
	- U: Cellular (Cellular Processing and Signaling)
	- O: Cellular (Cellular Processing and Signaling)
	- S: Conserved Hypothetical
	- R: Conserved Hypothetical

At the time of writing, the current EggNog 5.5 database does not include an X category and so Mobile genes are assigned based on a keyword search for the following: transposase, phage, integrase, viral, plasmid, integron, or transposon. Here is a [list of COG categories](http://www.sbg.bio.ic.ac.uk/~phunkee/html/old/COG_classes.html) and [here is another](https://ecoliwiki.org/colipedia/index.php/Clusters_of_Orthologous_Groups_(COGs)). Also see the [NCBI COG database](https://www.ncbi.nlm.nih.gov/research/cog/).

Any gene without an assigned anotation is labeled as Hypothetical.

![Gene Annotation COG Categories](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/g01-g02_annotations_bar.png)

Detail script output to screen for Chi-square tests:

```bash
Running Script...


	Reading genome fasta files ...

	Reading CDS fasta files ...

	Reading RBM file ...

	Reading PanCat file ...

	Reading EggNog annotation file ...

	Building recombinant position plots for genome A ...

	Building recombinant position plots for genome B ...

	Building annotation plot and tests ...

Initial Chi2 test contigency table:

Recombinant             Recombinant  Non-recombinant  Total
Annotation                                                 
Hypothetical                    216              996   1212
Conserved Hypothetical          168              852   1020
Ribosomal                        32              272    304
Information                      76              484    560
Cellular                        170              975   1145
Metabolism 1                    274             1227   1501
Metabolism 2                     96              281    377
Mobile                            6              129    135
Total                          1038             5216   6254

Chi2 expected frequency table:

Recombinant             Recombinant  Non-recombinant   Total
Annotation                                                  
Hypothetical             201.160217      1010.839783  1212.0
Conserved Hypothetical   169.293252       850.706748  1020.0
Ribosomal                 50.456028       253.543972   304.0
Information               92.945315       467.054685   560.0
Cellular                 190.039974       954.960026  1145.0
Metabolism 1             249.126639      1251.873361  1501.0
Metabolism 2              62.572114       314.427886   377.0
Mobile                    22.406460       112.593540   135.0
Total                   1038.000000      5216.000000  6254.0

chi2 statistic: 54.4502, dof: 7, chi2 pvalue: 0.000000

Post hoc Chi2 test contingency table for Hypothetical:

Recombinant   Recombinant  Non-recombinant  Total
Annotation                                       
OTHERs                822             4220   5042
Hypothetical          216              996   1212
Total                1038             5216   6254

Chi2 expected frequency table:

Recombinant   Recombinant  Non-recombinant   Total
Annotation                                        
OTHERs         836.839783      4205.160217  5042.0
Hypothetical   201.160217      1010.839783  1212.0
Total         1038.000000      5216.000000  6254.0

chi2 statistic: 1.5203, dof: 1, chi2 pvalue: 0.217580

Post hoc Chi2 test contingency table for Conserved Hypothetical:

Recombinant             Recombinant  Non-recombinant  Total
Annotation                                                 
OTHERs                          870             4364   5234
Conserved Hypothetical          168              852   1020
Total                          1038             5216   6254

Chi2 expected frequency table:

Recombinant             Recombinant  Non-recombinant   Total
Annotation                                                  
OTHERs                   868.706748      4365.293252  5234.0
Conserved Hypothetical   169.293252       850.706748  1020.0
Total                   1038.000000      5216.000000  6254.0

chi2 statistic: 0.0053, dof: 1, chi2 pvalue: 0.941827

Post hoc Chi2 test contingency table for Ribosomal:

Recombinant  Recombinant  Non-recombinant  Total
Annotation                                      
OTHERs              1006             4944   5950
Ribosomal             32              272    304
Total               1038             5216   6254

Chi2 expected frequency table:

Recombinant  Recombinant  Non-recombinant   Total
Annotation                                       
OTHERs        987.543972      4962.456028  5950.0
Ribosomal      50.456028       253.543972   304.0
Total        1038.000000      5216.000000  6254.0

chi2 statistic: 8.0532, dof: 1, chi2 pvalue: 0.004542

Post hoc Chi2 test contingency table for Information:

Recombinant  Recombinant  Non-recombinant  Total
Annotation                                      
OTHERs               962             4732   5694
Information           76              484    560
Total               1038             5216   6254

Chi2 expected frequency table:

Recombinant  Recombinant  Non-recombinant   Total
Annotation                                       
OTHERs        945.054685      4748.945315  5694.0
Information    92.945315       467.054685   560.0
Total        1038.000000      5216.000000  6254.0

chi2 statistic: 3.8319, dof: 1, chi2 pvalue: 0.050285

Post hoc Chi2 test contingency table for Cellular:

Recombinant  Recombinant  Non-recombinant  Total
Annotation                                      
OTHERs               868             4241   5109
Cellular             170              975   1145
Total               1038             5216   6254

Chi2 expected frequency table:

Recombinant  Recombinant  Non-recombinant   Total
Annotation                                       
OTHERs        847.960026      4261.039974  5109.0
Cellular      190.039974       954.960026  1145.0
Total        1038.000000      5216.000000  6254.0

chi2 statistic: 2.9488, dof: 1, chi2 pvalue: 0.085941

Post hoc Chi2 test contingency table for Metabolism 1:

Recombinant   Recombinant  Non-recombinant  Total
Annotation                                       
OTHERs                764             3989   4753
Metabolism 1          274             1227   1501
Total                1038             5216   6254

Chi2 expected frequency table:

Recombinant   Recombinant  Non-recombinant   Total
Annotation                                        
OTHERs         788.873361      3964.126639  4753.0
Metabolism 1   249.126639      1251.873361  1501.0
Total         1038.000000      5216.000000  6254.0

chi2 statistic: 3.7620, dof: 1, chi2 pvalue: 0.052429

Post hoc Chi2 test contingency table for Metabolism 2:

Recombinant   Recombinant  Non-recombinant  Total
Annotation                                       
OTHERs                942             4935   5877
Metabolism 2           96              281    377
Total                1038             5216   6254

Chi2 expected frequency table:

Recombinant   Recombinant  Non-recombinant   Total
Annotation                                        
OTHERs         975.427886      4901.572114  5877.0
Metabolism 2    62.572114       314.427886   377.0
Total         1038.000000      5216.000000  6254.0

chi2 statistic: 22.1090, dof: 1, chi2 pvalue: 0.000003

Post hoc Chi2 test contingency table for Mobile:

Recombinant  Recombinant  Non-recombinant  Total
Annotation                                      
OTHERs              1032             5087   6119
Mobile                 6              129    135
Total               1038             5216   6254

Chi2 expected frequency table:

Recombinant  Recombinant  Non-recombinant   Total
Annotation                                       
OTHERs        1015.59354       5103.40646  6119.0
Mobile          22.40646        112.59354   135.0
Total         1038.00000       5216.00000  6254.0

chi2 statistic: 13.8379, dof: 1, chi2 pvalue: 0.000199

Post hoc p values:
 [0.21757988944901951, 0.9418273956553266, 0.00454232505440978, 0.0502849911270642, 0.0859413619504977, 0.05242922236818852, 2.5759744269878553e-06, 0.00019927122428812193]

Benjamini/Hochberg corrected p values:
 [2.48662731e-01 9.41827396e-01 1.21128668e-02 8.38867558e-02
 1.14588483e-01 8.38867558e-02 2.06077954e-05 7.97084897e-04]


Complete success space cadet!! Finished without errors.
```

#### Sequence identity of RBMs vs. genome position

And finally, we have a figure with the sequence identity of RBMs vs. genome position. Each tik mark is a gene. The color represents the pangenome color Conserved (yellow), Core (pink), Accessory (green), or Specific (red) only. The colors in this figure do not indicate the recombinant classification. The y-axis sequence identities of the RBM alignments and the REC threshold (default 99.8) indicate this. The genome is divided into panels to fit the range of the genome sequence in base pairs (or nucleotides) along the x-axis. Each panel shows the genome range indicated in the lower rate. The -subs parameter allows to zoom in or out on the genome from the default of 10 panels. Only the first two panels are shown below.

![RBM sequence identity vs. genome position](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/g01-g02_A_posline_out.png)

([Return to Table of Contents](#table-of-contents))

### Step 06: Explore one vs. many genome groups of interest

This script is similar to the 03f script for genome pairs except it looks at one genome compared to many genomes.

The input is a two column tsv file with paths to the genome fastas in the first column and paths to the gene fastas in the second column. See the example input file group_1_list.tsv. The first genome in the list will be compared to the remaining genomes in the list.

And then it takes the same rbm, annotation, and pancate files as above.

This script returns data, graphics and statistics concerning the genomic positions between genes categorized by highly conserved core genes, recombining core genes, recombing accessory genes, and non-recombining genes. A gene is said to be recombining if the RBM value equals 100 and if it is not highly conserved.

This script performs some statistical tests on the distance between recombinant genes, and on the distributions of gene annotations.

The -rec parameter (default 99.8) sets the threshold for an RBM to be considered recombinant or not. It affects the results of the gene annotations plot and chi-square hypothesis tests and it affects the recombinant positions plot and Poisson and Geometric distribution tests.

Input:

	Input is a two column tab separated file with genome fasta paths in first column and gene fasta paths in the second column. The first genome/gene file in the list becomes the main genome the others are compared to.

	genomes_dir/x_x_genome01_x.fna	genes_dir_fnn/x_x_genome01_x.fnn
	genomes_dir/x_x_genome02_x.fna	genes_dir_fnn/x_x_genome02_x.fnn
	genomes_dir/x_x_genome03_x.fna	genes_dir_fnn/x_x_genome03_x.fnn
	genomes_dir/x_x_genome04_x.fna	genes_dir_fnn/x_x_genome04_x.fnn
	genomes_dir/x_x_genome05_x.fna	genes_dir_fnn/x_x_genome05_x.fnn
	genomes_dir/x_x_genome06_x.fna	genes_dir_fnn/x_x_genome06_x.fnn
	genomes_dir/x_x_genome07_x.fna	genes_dir_fnn/x_x_genome07_x.fnn
	genomes_dir/x_x_genome08_x.fna	genes_dir_fnn/x_x_genome08_x.fnn
	genomes_dir/x_x_genome09_x.fna	genes_dir_fnn/x_x_genome09_x.fnn
	genomes_dir/x_x_genome10_x.fna	genes_dir_fnn/x_x_genome10_x.fnn
	genomes_dir/x_x_genome11_x.fna	genes_dir_fnn/x_x_genome11_x.fnn
	genomes_dir/x_x_genome12_x.fna	genes_dir_fnn/x_x_genome12_x.fnn
	genomes_dir/x_x_genome13_x.fna	genes_dir_fnn/x_x_genome13_x.fnn
	genomes_dir/x_x_genome14_x.fna	genes_dir_fnn/x_x_genome14_x.fnn
	genomes_dir/x_x_genome15_x.fna	genes_dir_fnn/x_x_genome15_x.fnn
	genomes_dir/x_x_genome16_x.fna	genes_dir_fnn/x_x_genome16_x.fnn
	genomes_dir/x_x_genome17_x.fna	genes_dir_fnn/x_x_genome17_x.fnn
	genomes_dir/x_x_genome18_x.fna	genes_dir_fnn/x_x_genome18_x.fnn
	genomes_dir/x_x_genome19_x.fna	genes_dir_fnn/x_x_genome19_x.fnn
	genomes_dir/x_x_genome20_x.fna	genes_dir_fnn/x_x_genome20_x.fnn
	genomes_dir/x_x_genome21_x.fna	genes_dir_fnn/x_x_genome21_x.fnn
	genomes_dir/x_x_genome22_x.fna	genes_dir_fnn/x_x_genome22_x.fnn
	genomes_dir/x_x_genome23_x.fna	genes_dir_fnn/x_x_genome23_x.fnn
	genomes_dir/x_x_genome24_x.fna	genes_dir_fnn/x_x_genome24_x.fnn
	genomes_dir/x_x_genome25_x.fna	genes_dir_fnn/x_x_genome25_x.fnn
	genomes_dir/x_x_genome26_x.fna	genes_dir_fnn/x_x_genome26_x.fnn
	genomes_dir/x_x_genome27_x.fna	genes_dir_fnn/x_x_genome27_x.fnn
	genomes_dir/x_x_genome28_x.fna	genes_dir_fnn/x_x_genome28_x.fnn
	genomes_dir/x_x_genome29_x.fna	genes_dir_fnn/x_x_genome29_x.fnn
	genomes_dir/x_x_genome30_x.fna	genes_dir_fnn/x_x_genome30_x.fnn

Output:

	1. tsv file: Gene RBM table
	1. recombinant gene position plot
	1. recombinant distribution plot
	1. recombinant annotations plot
	1. RBM sequence identity vs. genome position plot
	1. tsv file: Binary matrix for recombinant RBM positions
	1. recombinant rRBM curve plot
	1. tsv file: rRBM data
	1. recombinant RBM clustermap

The Gene table tsv file contains the following columns and data:

	- Genome: The one genome that faced many.
	- Contig: Each contig from the one genome that faced many.
	- Gene: Sequence identifier contains the genome_identifier_contigNumberFromAssembly_geneNumberOnContig.
	- Match Genome: Genome of the RBM
	- Match Conitg: Contig for the RBM
	- Match Gene: Gene of the RBM
	- PanCat:  column indicates the pangenome class assigned to the gene.
	- pID: Nucleotide sequence identity of pairwise gene alignment of the RBM.
	- REC: Assigned a 0 or 1. 1 indicates the gene sequence has ≥ REC% identity with its corresponding RBM in the other genome and thus a candidate for recent homologous recombination. A 0 indicates the gene does not. The REC threshold is controlled by the -rec parameter (default 99.8) and sets the threshold for an RBM to be considered recombinant or not. It affects the results of the gene annotations plot and chi-square hypothesis tests and it affects the recombinant positions plot and Poisson and Geometric distribution tests.
	- Recombinant: Recombinant or Non-Recombinant classification based on "REC"
	- Start: Start of gene in genome coordinates
	- Stop: End of gene in genome coordinates
	- Strand: 1 for positive or sense and -1 for negative or anti-sense strand of predicted CDS.
	- COG Category: Higher level COG category assignement based on COG.
	- COG: Single letter COG assignment.
	- Gene Annotation: Short gene name of the assigned gene annotation
	- Annotation Description: Long form gene name or description of the assigned gene annotation
	- Mismatch: Number of mismatches in the blast alignment for RBMs.
	- Width: Gene length. Distance between stop and start.

```bash
# for script info/option
python 00d_Workflow_Scripts/03g_Recombinant_group_analysis.py -h

# with default settings
python 00d_Workflow_Scripts/03g_Recombinant_group_analysis.py -i group_1_list.tsv -o group_g1 -rbm RBMs_allV.rbm -pc pancat_file.tsv -ano annotation_file
```

#### Recombinant gene position by pangenome class

![Recombinant gene position by pangenome class](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/group_g1_posbar.png)

#### Distance between recombination events distribution test

![Distance between recombination events distribution test](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/group_g1_gene_distribution.png)

#### Recombinant vs. Non-recombinant gene annotation test

![Recombinant vs. Non-recombinant gene annotation test](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/group_g1_annotations_bar.png)

#### Sequence identity of RBMs vs. genome position

![Sequence identity of RBMs vs. genome position](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/group_g1_posline.png)

#### Recombinant RBM curve plot

##### this is an old plot adapted from my pangenome work. Replacing it with the rarefaction plot below.

```bash
# Recombinant RBM (rRBMs) curve
python 03h_RBM-Curve_Calculate_Model_Plot.py -b group_g1_rbm_matrix.tsv -o group_g1_rRBMcurve.pdf
```

![Recombinant RBM curve plot](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/group_g1_rRBMcurve.png)

#### Recombinant RBM gene clustermap

```bash
# gene RBM clustermap
python 03i_RBM_Clustermap.py -b group_g1_rbm_matrix.tsv -o group_g1_rbmclustermap.pdf
```

![Recombinant RBM clustermap](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/group_g1_rbmclustermap.png)

#### Recombinant RBM gene rarefaction plot

This script takes 3 input files.

The first file is the tab separated RBM matrix generated by the
03g_Recombinant_group_analysis.py script. Each column represents a
genome and each row represents an RBM gene. Genome IDs should be
on the first line (row). Each row after that is a RBM gene with a value
for each genome of 0 (absent), 1 (present), or 2 (conserved).

The second file is generated by the user and is a two column, comma-
separated list of genomeID,userAssignedClade where genomeID matches the
genome IDs from line 1 of the rmb matrix file and the userAssignedClade
can be anything the user assigns. Each genome in the rbm matrix should
have a clade assignment. Example:

	genome1,clade1
	genome2,clade1
	genome3,clade2
	genome4,clade3
	genome5,clade1
	genome6,clade3
	genome7,clade3
	genome8,clade2

The third file is a two column, comma-separated list of clade,color.

	clade1,#e41a1c
	clade2,#377eb8
	clade3,#4daf4a

Color guide: https://colorbrewer2.org/#type=qualitative&scheme=Set1&n=9

*conserved sites are masked (not included)*

*If you don't have multiple clades, label all genomes the same clade and color.*

```bash
python 03j_RBM-Clade_Rarefaction.py -r group_g1_rbm_matrix.tsv -l clade_list.txt -c color_list.txt -o group_g1_rbm_rarefaction
```

![Recombinant RBM clustermap](https://github.com/rotheconrad/F100_Prok_Recombination/blob/main/00a_example_figures/rbm_rarefaction_test_plot.png)

([Return to Table of Contents](#table-of-contents))

# Software Dependencies

## External dependencies

- [FastANI](https://github.com/ParBLiSS/FastANI)
- [Prodigal](https://github.com/hyattpd/Prodigal)
- [BLAST+](https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- [Diamond]()
- [MMseqs2](https://github.com/soedinglab/MMseqs2)
- [EggNog-mapper](https://github.com/eggnogdb/eggnog-mapper)
- [COGclassifier](https://github.com/moshi4/COGclassifier/) (optional alternative to EggNog Mapper)
- [Python](https://www.python.org/) version 3.6+ (for all custom code in this workflow)

#### References

1. Jain C, Rodriguez-R LM, Phillippy AM, Konstantinidis KT, Aluru S. High throughput ANI analysis of 90K prokaryotic genomes reveals clear species boundaries. Nature communications. 2018 Nov 30;9(1):1-8.
1. Hyatt D, Chen GL, LoCascio PF, Land ML, Larimer FW, Hauser LJ. Prodigal: prokaryotic gene recognition and translation initiation site identification. BMC bioinformatics. 2010 Dec;11(1):1-1.
1. Camacho C, Coulouris G, Avagyan V, Ma N, Papadopoulos J, Bealer K, Madden TL. BLAST+: architecture and applications. BMC bioinformatics. 2009 Dec;10(1):1-9.
1. Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. Nature methods. 2021 Apr;18(4):366-8.
1. Steinegger M, Söding J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nature biotechnology. 2017 Nov;35(11):1026-8.
1. Cantalapiedra CP, Hernández-Plaza A, Letunic I, Bork P, Huerta-Cepas J. eggNOG-mapper v2: functional annotation, orthology assignments, and domain prediction at the metagenomic scale. Molecular biology and evolution. 2021 Dec 1;38(12):5825-9.
1. Huerta-Cepas J, Szklarczyk D, Heller D, Hernández-Plaza A, Forslund SK, Cook H, Mende DR, Letunic I, Rattei T, Jensen LJ, von Mering C. eggNOG 5.0: a hierarchical, functionally and phylogenetically annotated orthology resource based on 5090 organisms and 2502 viruses. Nucleic acids research. 2019 Jan 8;47(D1):D309-14.
1. Sanner MF. Python: a programming language for software integration and development. J Mol Graph Model. 1999 Feb 1;17(1):57-61.

## Required packages for Python

- [pandas](https://pandas.pydata.org/) 
- [numpy](https://numpy.org/)
- [scipy](https://scipy.org/)
- [statsmodels](https://www.statsmodels.org)
- [matplotlib](https://matplotlib.org/)
- [seaborn](https://seaborn.pydata.org/)
- [lmfit](https://lmfit.github.io/lmfit-py/)
- [datashader](https://datashader.org/)
- [pygam](https://pygam.readthedocs.io/)

*Python and all packages can be easily installed with conda or pip. Prodigal, BLAST+ and MMseqs2 can also be installed easily with [Conda](https://docs.conda.io/en/latest/miniconda.html). Just search "conda install name"*

#### References

1. Van Rossum G, Drake FL. Python 3 Reference Manual. Scotts Valley, CA: CreateSpace; 2009.
1. McKinney W, others. Data structures for statistical computing in python. In: Proceedings of the 9th Python in Science Conference. 2010. p. 51–6.
1. Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D, et al. Array programming with NumPy. Nature. 2020;585:357–62.
1. Virtanen P, Gommers R, Oliphant TE, Haberland M, Reddy T, Cournapeau D, et al. SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods. 2020;17:261–72.
1. Seabold S, Perktold J. Statsmodels: Econometric and statistical modeling with python. InProceedings of the 9th Python in Science Conference 2010 Jun 28 (Vol. 57, No. 61, pp. 10-25080).
1. Hunter JD. Matplotlib: A 2D graphics environment. Computing in science & engineering. 2007;9(3):90–5.
1. Waskom ML. Seaborn: statistical data visualization. Journal of Open Source Software. 2021 Apr 6;6(60):3021.
1. Newville M, Stensitzki T, Allen DB, Rawlik M, Ingargiola A, Nelson A. LMFIT: Non-linear least-square minimization and curve-fitting for Python. Astrophysics Source Code Library. 2016 Jun:ascl-1606.
1. James A. Bednar, Jim Crist, Joseph Cottam, and Peter Wang (2016). "Datashader: Revealing the Structure of Genuinely Big Data", 15th Python in Science Conference (SciPy 2016).
1. Servén D., Brummitt C. (2018). pyGAM: Generalized Additive Models in Python. Zenodo. DOI: 10.5281/zenodo.1208723

([Return to Table of Contents](#table-of-contents))

# How to Cite

If you use any part of this workflow in your research please cite the following manuscript:

PLACEHOLDER FOR MANUSCRIPT CITATION AND LINK

# Future Improvements

1. More efficient All vs. All RBM with Diamond or other.
1. Filter MMSeqs2 Gene clusters for improved true ortholog assignments.
1. Compute number of unique alleles for each gene cluster.
1. Compute pN/pS for each gene cluster.
