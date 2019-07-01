#!/bin/bash

# BUSCO dot plots

# how to run: generate_R_script.sh genome1.tsv genome2.tsv genome1.fa genome2.fa

if [ "$#" != 4 ]; then
echo ''
echo 'Usage = ./generate_R_script.sh genome1.tsv genome2.tsv genome1.fa genome2.fa'
echo ''
        exit 0
fi


if [ -d "work_dir_synteny" ]; then
echo ''
echo 'remove previous work_dir_synteny first!!'
echo ''
        exit 0
fi


mkdir work_dir_synteny
cp $1 work_dir_synteny/full_1.tsv
cp $2 work_dir_synteny/full_2.tsv
cp $3 work_dir_synteny/1.fa
cp $4 work_dir_synteny/2.fa

cd work_dir_synteny

cp ../shell_R .
cp ../pallete.txt .
cp ../transparent_pallete.txt .

# generate chorm sizes 
cat 1.fa | fastalength | sort -nr > 1.sizes

cat 2.fa | fastalength | sort -nr > 2.sizes


# complete buscos for genome 2

awk '$2 == "Complete"' full_2.tsv | cut -f 1,3,4 | sort -k1,1 > 2_comp.tsv

cat 2_comp.tsv | cut -f 2 | sort | uniq > 2_busco.chroms

# complete buscos for genome 1

awk '$2 == "Complete"' full_1.tsv | cut -f 1,3,4 | sort -k1,1 > 1_comp.tsv

cat 1_comp.tsv | cut -f 2 | sort | uniq > 1_busco.chroms




### generate links.table


cat <(grep -f 1_busco.chroms 1.sizes | sort -nr ) <(grep -f 2_busco.chroms 2.sizes | sort -n ) \
        | awk '{print $2" = "$1}' | sed -e 's/^/"/' -e 's/ = /" = /' | sed 's/$/,/' | sed '$s/,$/\)/' > genome.list


grep -f <(grep -f <(cat 1_comp.tsv | cut -f 1 ) 2_comp.tsv | cut -f 1 ) 1_comp.tsv > shared_1.comp

grep -f <(cat shared_1.comp | cut -f 1 ) 2_comp.tsv > shared_2.comp

paste shared_1.comp shared_2.comp > links.table

cat links.table | cut -f 2 | sort | uniq > gen1_buscos_only

# number of scaffold colours to use

num_colours=$(wc -l gen1_buscos_only | awk '{print $1}' )

head -n $num_colours transparent_palette.txt > colours_to_use.txt

#paste gen1_buscos_only colours_to_use.txt > busco_chrom_colours.tbl
paste <(grep -f gen1_buscos_only 1.sizes | cut -f 2 ) colours_to_use.txt > busco_chrom_colours.tbl

# generate links.table for each scaffold with busco on it

while read scaff pattern ; do
        
        grep $scaff links.table > "$scaff".table

echo $scaff'links_chromosomes_1 = c(' > links_chromosomes_1
echo $scaff'links_chromosomes_2 = c(' > links_chromosomes_2

echo $scaff'links_pos_1 = c(' > links_pos_1
echo $scaff'links_pos_2 = c(' > links_pos_2

echo $scaff'links_labels = c(' > links_scaff_table

echo ')' > bracket.txt
        
cat "$scaff".table | cut -f 2 | sed -e 's/^/\'\''/' -e 's/$/\'\'',/' | sed '$s/,$//' | cat links_chromosomes_1 - bracket.txt > links_"$scaff".1

cat "$scaff".table | cut -f 5 | sed -e 's/^/\'\''/' -e 's/$/\'\'',/' | sed '$s/,$//' | cat links_chromosomes_2 - bracket.txt > links_"$scaff".2


cat "$scaff".table | cut -f 3 | sed -e 's/$/,/' | sed '$s/,$//' | cat links_pos_1 - bracket.txt > links_pos_"$scaff".1

cat "$scaff".table | cut -f 6 | sed -e 's/$/,/' | sed '$s/,$//' | cat links_pos_2 - bracket.txt > links_pos_"$scaff".2


cat "$scaff".table | cut -f 1 | sed -e 's/^/"/' -e 's/$/",/' | sed '$s/,$//' | cat links_scaff_table - bracket.txt > links_"$scaff".labels

echo 'tracklist = tracklist + BioCircosLinkTrack("myLinkTrack'$scaff'", '$scaff'links_chromosomes_1, '$scaff'links_pos_1, '$scaff'links_pos_1 , '$scaff'links_chromosomes_2, '$scaff'links_pos_2, '$scaff'links_pos_2 
  
cat links_"$scaff".1 links_"$scaff".2 links_pos_"$scaff".1 links_pos_"$scaff".2 links_"$scaff".labels >> body1

cat "$scaff"_tracklist >> body3

done < busco_chrom_colours.tbl

# generate colours for scaffolds/chromosomes

scaff_colours=$(cat 1_busco.chroms 2_busco.chroms | wc -l )

scaff_col_pal=$(head -n $scaff_colours pallete.txt | sed -e 's/^/"/' -e 's/$/",/' | sed '$s/,$//' ) # > scaff_col_pal.txt

head -n $scaff_colours pallete.txt | sed -e 's/^/"/' -e 's/$/",/' | sed '$s/,$//'  > scaff_col_pal.txt

echo 'tracklist = BioCircosBackgroundTrack("myBackgroundTrack", minRadius = 0, maxRadius = 0.9, borderSize = 0, fillColors = "#EEFFEE")' > body2 

echo 'BioCircos(tracklist, genome = myGenome, genomeFillColor = c('$scaff_col_pal'), chrPad = 0.02, displayGenomeBorder = FALSE, yChr =  FALSE, genomeTicksDisplay = FALSE,  genomeLabelTextSize = 0, genomeLabelDy = 

cat <(head -n 5 shell_R.txt ) genome.list body1 body2 body3 body4 > ../busco_dot.R


