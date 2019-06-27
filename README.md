# Genome-wide association (GWA) tutorial quantitative analyss

## fork to rafalcode
* some scripts tried to install packages: that is not the job of scripts.
* packages should all be pre-installed.
* takes about 3 hours to run, even when num.thread is set to 4 (although Lpruning() only ever seems to run 1 thread)
* allinone.R takes the three scripts and puts them into one.
* there is a prepath variable saying where the data is, you may need to change this.
* output0.txt can be useful, it gives the typical console output for a successful run.

## Additional files

For this tutorial you will additionally need the files

- 117Malay_282lipids.txt
- 120Indian_282lipids.txt
- 122Chinese_282lipids.txt
- 105Indian_2527458snps.bed, .bim, .fam
- 108Malay_2527458snps.bed, .bim, .fam
- 110Chinese_2527458snps.bed, .bim, .fam

stored in the folders 'Lipidomic' and 'Genomics' contained in the following .gz:
http://phg.nus.edu.sg/StatGen/public_html/Iomics/downloads/iOmics_data.tar.gz

## Instructions

1. Combine the folders 'Lipidomic' and 'Genomics' and all files from this repo in your working directory.
2. Install all packages listed on top of the scripts. `snpStats` and `SNPRelate` are deposited in BioConductor, all other packages in CRAN.
3. Run the scripts in their exact numbered order.

## Acknowledgements

This work was largely based on the following publications:

- *Establishing multiple omics baselines for three Southeast Asian populations in the Singapore Integrative Omics Study*, Saw et al. (2017), Nat. Comm. (data source)
- *A guide to genome-wide association analysis and post-analytic interrogation*, Reed et al. (2015), Stats. in Med. (method source)

Enjoy, any feedback is welcome!

Francisco
