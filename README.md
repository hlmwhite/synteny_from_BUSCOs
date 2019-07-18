# synteny_from_BUSCOs
Generate draft genome synteny between related organisms using BUSCO (.tsv) output.

Currently requires:
- BUSCO (https://busco.ezlab.org/)
- exonerate utility 'fastalength' (https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate)
- R (or Rstudio for plot viewing, rather than through html) 
- R package BioCircos (https://cran.r-project.org/web/packages/BioCircos/vignettes/BioCircos.html)

steps to run:

0. optional step - rename fasta headers (e.g. >chr1, >chr2, etc). Stops R getting stuck on funny characters

1. run BUSCO over each genome of interest. Make sure to use the same database between genomes!!!

> ./run_BUSCO.py -i SEQUENCE_FILE -o OUTPUT_NAME -l LINEAGE -m geno 

2. provide full.tsv files and fasta files of both genomes to study. Currently requires running with palette.txt and transparent_palette.txt in working directory.

> ./generate_Rscript.sh full_busco_genome1.tsv full_busco_genome2.tsv genome1.fasta genome2.fasta

3. run the Rscript "busco_dot.R" in R/Rstudio.

> source("busco_dot.R", echo = TRUE)


The output should resemble something like this...

![Image of canu_scaf](https://github.com/hlmwhite/synteny_from_BUSCOs/blob/master/synteny_from_busco.out.png)

### why do this?

The tool should hopefully provide a quick and simple method to observe relationships and conserved synteny between complete/near complete genomes of related organisms. 

### need to do...

1. maybe change fasta inputs to chrom.sizes files
2. how to denote genome 1 and 2 in the plot
3. output chromosome/scaffold names and sizes to R terminal (any other important info? e.g  num BUSCOs)
4. rename chromosome/scaffold names in fasta and in .tsv files? easier for R!
5. direct output to .png/.jpg
6. probably need to clean up the bash script 
7. output other useful info - e.g. % of genome containing complete buscos, and therefore % missing from plot, etc
8. add an option for colourblind colour palette


