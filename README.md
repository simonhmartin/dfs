# The *D* Frequency Spectrum

* [Introduction](#Introduction)
* [Code in this repository](#Code-in-this-repository)
* [Explore *D*<sub>FS</sub> using a local R SHINY app](#Explore-DFS-using-a-local-R-SHINY-app)
* [Compute *D*<sub>FS</sub> from empirical data](#Compute-DFS-from-empirical-data)

## Introduction

*D*<sub>FS</sub> is a simple extension to the classic D statistic (AKA the ABBA BABA test). *D*<sub>FS</sub> partitions *D* according to the frequencies of derived alleles. Thus, it reveals how the overall *D* value is broken down across different allele frequency classes (or bins). *D*<sub>FS</sub> therefore carries information about the timing and direction of introgression.

Check out this [SHINY app](https://shmartin.shinyapps.io/shiny_plot_dfs_moments/) that allows you to explore how *D*<sub>FS</sub> behavies over a range of model parameters.

## Code in this repository

`DFS.R` provides functions to compute *D*<sub>FS</sub>, as well as the standard *D*, and some related statistics.

`plot_DFS_from_SFS.R` gives code to compute and plot *D*<sub>FS</sub> given an input frequency spectrum (see below for details on the input). Example data for six taxa are provided in the `empirical_data` directory.

`sim_DFS_moments.R` provides code to run simulations locally (either batch or singlular). The simulations are run by sending commands to the python script `sim_SFS_moments.py`, which requires the [`moments`](https://bitbucket.org/simongravel/moments/src/master/) package and `numpy` to be installed on your system.

`SHINY_sim_DFS_moments/` is a directory containing an R SHINY app, that you can run locally to simulate and visualise *D*<sub>FS</sub> for any model you choose.

## Explore *D*<sub>FS</sub> using a local R SHINY app

The easiest way to simulate *D*<sub>FS</sub> for any model you like is to run the provided R SHINY app locally. You will need `Rstudio` with the `shiny` package installed, as well as the [`moments`](https://bitbucket.org/simongravel/moments/src/master/) and `numpy` python packages installed.

To run the app, lauch `Rstudio` and set the working directory to `SHINY_sim_DFS_moments/`.

Then, in the console run:

```R
library(shiny)
runApp()
```


## Compute *D*<sub>FS</sub> from empirical data

The script `plot_DFS_from_SFS.R` gives code to compute and plot *D*<sub>FS</sub> given an input frequency spectrum.

#### Tabular frquency spectrum

These scripts make use of 3D or 4D frequency spectra. However, because these can have very many entries (sometimes more than the number of SNPs in the genome if sample sizes are large), and many of these will be zero, the scripts do not make use of standard SFS formats. Instead the they use of a sparse tabular format, which only records non-zero entries of the SFS. This is simply a table in which the first three columns (or 4 in the case of a 4D SFS) give the allele count in each population (equivalent to the indices of the multidimensional SFS). The final column gives the corresponding number of sites.

Within the script, there are lines that specify which columns correspond to each population and which column corresponds to the number of sites.

#### Compute SFS from a VCF file

Given a VCF file with three focal populations and an outgroup, you can compute a 3D SFS using a few steps. **Note that *D*<sub>FS</sub> is only useful if your sample size for populations P1 and P2 is at least 10 haploid genomes (5 diploids)**.

Before doing so, it is necessary to download the `genomics_general` reporitory.

```bash
#download package
wget https://github.com/simonhmartin/genomics_general/archive/v0.4.tar.gz
#extract files from zipped archive
tar -xzf v0.4.tar.gz
#delete zipped file
rm v0.4.tar.gz
```
To ensure that the libraries are recognisable by python, add the `genomics_general' directory to the Python path

```bash
export PYTHONPATH=$PYTHONPATH:genomics_general-0.4
```

Then, convert the VCF to `.geno` format for downstream analysis. Here we add a filter to retain only genotypes suppoted by a read depth of 8.

```bash
python genomics_general-0.4/VCF_processing/parseVCFs.py --threads 6 \
-i mydata.vcf.gz --skipIndels --gtf flag=DP min=8 -o mydata.geno.gz
```

Now we can compute the SFS. This command has two parts. We first run `freq.py` to compute allele frequencies at each site in each population. This requires that you also provide a populations file, that gives all samples names in the first column and the population they belong to in the second column (you can use any names for your populations). While it is possibel to compute the derived allele frequency, or the minor allele frequency, here we will not request either of these, so instead it will just return the counts for all four bases at each site in each population.


``` bash
python genomics_general-0.4/freq.py --threads 6 -g mydata.geno.gz \
-p P1 -p P2 -p P3 -p OG --popsFile populations.txt | gzip > mydata.basecounts.tsv.gz
```

Next we run `sfs.py`, which produces one or more SFS files. By including the `--outgroup OG` flag, we are telling it to use the outgroup population ("OG") to polarise allele frequencies for the ingroups. We then tell it to output just one 3D SFS, representing populations P1, P2 and P3. We also specify that it should subsample the data down to 10 haplotypes for P1 and P2 and 2 for P3. This is useful if you have larger sample sizes but lots of missing data. At each site, the script will try to use the available data to make up the number of genomes required. **Note that your SFS must have the same number of samples for P1 and P2**.

```
python genomics_general-0.4/sfs.py -i mydata.basecounts.tsv.gz --inputType baseCounts \
--outgroup OG  --FSpops P1 P2 P3 --subsample 10 10 2 --pref mydata. --suff .subsample10.sfs
```

For ease of use, it is also possible to pipe the output of `freq.py` directly to `sfs.py` (in which case you should not include the `-i` option).

The resulting SFS can be loaded into the script `plot_DFS_from_SFS.R`. Be sure to set the number of (haploid) samples correctly.

