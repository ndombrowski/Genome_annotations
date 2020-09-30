Annotation pipeline
================
Nina Dombrowski
2020-09-30

This example is run on a set of \~20 Crenarchaeal Genomes and the
location dirs are kept for convenience. To run on this example please
copy the genomes to your own dir and *do not run this in the original
folder*

# 
# General
# 


## Version programs

# 
1.  prokka 1.14-dev
2.  Python 2.7.5
3.  perl v5.16.3
4.  HMMER 3.1b2
5.  blastp (Protein-Protein BLAST 2.7.1+)
6.  diamond 0.9.22
7.  interproscan: InterProScan-5.29-68.0

## Version databases
# 

1.  KAAS: Automatic Annotation Server Ver. 2.1
2.  TIGRFAMs\_15.0\_HMM.LIB: downloaded Sept 2018
3.  PFams: RELEASE 31.0, downloaded Sept 2018
4.  CAZymes: dbCAN-HMMdb-V7, dbCAN v2 on 21 Sept 2018
5.  Merops: downloaded from Merops server in Nov18
6.  HydDB: downloaded form HydDB webserver in Nov18
7.  COGs: downloaded from NCBI Nov18
8.  TransporterDB: downloaded from TCDB Nov 18
9.  ncbi\_nr (maintained by Alejandro): files from Nov18
10. ArCOGs. [link to
    ftp](ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/zip.aliar14.tgz)
    (update date: 2018)
11. KOhmms: downloaded 90419 from [link
    here](https://www.genome.jp/tools/kofamkoala/)

# 
# Setup of genomes
# 


## Set working directory
# 

``` bash
cd /export/lv1/user/spang_team/Projects/BlackSea18/Database_Annotations/v2
```

## Location of relevant data
# 

Notice: These genomes can be used for testing the pipeline

``` bash
#1. BS18 v2 genomes
cd /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/*/*faa
```

We have 5001 genomes to analyse

## Make a file to loop through new bins

# 

``` bash
mkdir FileLists

cat /export/lv1/user/spang_team/Projects/BlackSea18/Contaminant_cleaning/*/FileLists/Bins_to_keep > FileLists/FileLists/Files.txt
```

# 
# Do the Annotations
# 

## Get list to link old with new contig name, test name: NIOZ136\_mx\_b338\_temp

# 

``` bash
mkdir contig_maping

#get contig length for controlling things
for sample in `cat FileLists/Files.txt`; do perl ~/../spang_team/Scripts/Others/length+GC.pl ../../Contaminant_cleaning/*/fna/v2_cleaned/renamed/${sample}_2.fna > contig_maping/${sample}_2_temp1; done

#add number of contigs as separate column
for sample in `cat FileLists/Files.txt`; do awk '$1=(FNR FS $1 FILENAME)' contig_maping/${sample}_2_temp1 > contig_maping/${sample}_temp2; done

#add in binIDs
for sample in `cat FileLists/Files.txt`; do awk 'BEGIN{OFS="\t"}{split($2,a,"-"); print a[2],$1,$3,$4}' contig_maping/${sample}_temp2  > contig_maping/${sample}_temp3; done

for sample in `cat FileLists/Files.txt`; do awk 'BEGIN{OFS="\t"}{split($1,a,"/"); print a[1],a[2], $2,$3,$4}' contig_maping/${sample}_temp3 | sed 's/contig_maping//g' | sed 's/_temp1//g'  > contig_maping/${sample}_temp4; done

#combine and cleanup
cat contig_maping/*temp4 > contig_maping/Contig_Old_mapping.txt

rm contig_maping/*temp*

#create file for merging with prokka IDs
awk 'BEGIN{OFS="\t"}{print $2"_contig_"$3,$0}' contig_maping/Contig_Old_mapping.txt > contig_maping/Contig_Old_mapping_for_merging.txt
```

## Prepare list that links proteins with binIDs

# 

``` bash
#get list of protein accession nrs
grep "^>" /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa  > temp
cut -f1 -d " " temp > temp2
sed 's/>//g' temp2 > AllProteins_list.txt

rm temp*

#Modify protein list to add in a column with binID
awk -F'\t' -v OFS='\t' '{split($1,a,"-"); print $1, a[1]}' AllProteins_list.txt | LC_ALL=C sort > 1_Bins_to_protein_list.txt

wc -l 1_Bins_to_protein_list.txt
```

We have 9,785,373 proteins

## Prepare a file with prokka annotations $ Lengths $ other genome info

# 

``` bash
#match contig names with protein IDs
gzip /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/*/*gbk

#make list for looping
for f in /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/*/*gbk.gz; do echo $f >> list_of_gzips; done

#control that the nr of files we have is correct
wc -l list_of_gzips

'''
5001
'''

#convert the gbk file to get the nr of proteins per contig
python ~/../spang_team/Scripts/Others/Parse_prokka_for_MAGs_from_gbk-file.py -i list_of_gzips -t Contig_Protein_list.txt

#make prokka to binId links
awk 'BEGIN{OFS="\t"}{split($1,a,"-"); print a[2],$2}' 1_Bins_to_protein_list.txt | awk 'BEGIN{OFS="\t"}{split($1,a,"_"); print a[1],$2}' | sort | uniq > Prokka_to_BinID.txt

#link old to new binIDs
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' Prokka_to_BinID.txt Contig_Protein_list.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $1,$7,$2,$2,$3,$4,$5}'  > temp_Bin_Contig_Protein_list.txt

#add in extra column for contig nr
awk 'BEGIN{FS="\t";OFS="\t"}{split($4,a,"_"); print $2"_contig_"a[2],$1,$2,$3,a[2],$5,$6,$7}' temp_Bin_Contig_Protein_list.txt > Bin_Contig_Protein_list.txt

#merge with old contig IDs
#headers: accession, BinID, newContigID, oldContigID, mergeContigID,ContigLengthNew, LengthContigOld, GC,ProteinID, prokka
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' contig_maping/Contig_Old_mapping_for_merging.txt Bin_Contig_Protein_list.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $3"-"$7,$3,$4,$10,$1,$6,$14,$13,$7,$8}' >  Bin_Contig_Protein_list_merged.txt

#prepare a file with protein length and gc to add protein length into our file
perl ~/../spang_team/Scripts/Others/length+GC.pl /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa > temp

#merge with contig/protein info
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' temp Bin_Contig_Protein_list_merged.txt | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$3,$4,$5,$6,$7,$8,$9,$12,$13,$10}' > temp4

#add in taxon string
#Comment: Bin_to_tax.txt info was copied from bin_stats
#1st column = binID, 2nd column = tax string from gtdb (can be exchange by other tax info of course)
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' /export/lv1/user/spang_team/Projects/BlackSea18/Database_Annotations/v1/Bin_to_tax.txt temp4 | awk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$14,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12}' > temp5

#add a header
echo -e "accession\tBinID\tTaxString\tNewContigID\tOldContigId\tContigIdMerge\tContigNewLength\tContigOldLength\tGC\tProteinID\tProteinGC\tProteinLength\tProkka" | cat - temp5 > B_GenomeInfo.txt

#control that all is ok
wc -l B_GenomeInfo.txt
#9,785,374

#remove temp files
rm temp*
```

## ArCOGs search –\> done

# 

[download files from here](ftp://ftp.ncbi.nih.gov/pub/wolf/COGs/arCOG/)

``` bash
#modify file for merging
#sed 's/ /_/g' ar14.arCOGdef.tab > ar14.arCOGdef.nina.tab
#newer version:  sed 's/ /_/g' ar14.arCOGdef18.tab > ar14.arCOGdef18.nina.tab

#modifications to the arcog alignment
#cd /export/lv1/user/spang_team/Databases/arCOGs2014/ar14.ali

#1A. rename and add in filename/arcog (old version 2014)
#for i in *sr; do
#awk '/^/{sub("^","&"FILENAME"_");sub(/\.sr/,x)}1' $i > renamed/$i
#done

#1B. rename and add in filename/arcog (new version 2018, done 06092019)
#cd /export/lv1/user/spang_team/Databases/arCOGs2019/ali.ar14_v2018
#mkdir renamed
#for i in *sr; do
#awk '/^/{sub("^","&"FILENAME"_");sub(/\.sr/,x)}1' $i > renamed/$i
#done

#make list of files to use
#ls *sr > arcog_list

mkdir arcogs

#make db out of concat genomes
makeblastdb -in /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa  -out /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean   -dbtype prot -parse_seqids

#run arcog search
~/../spang_team/Scripts/Hmmer/hmmsearchTable /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa  ~/../spang_team/Databases/arCOGs2019/All_Arcogs_2018.hmm 40 -E 1e-5 > arcogs/All_arcogs_hmm.txt

#separate header with arcog and ID
awk -F'\t' -v OFS='\t' '{split($3,a,"."); print $1, a[1], $5,$6}' arcogs/All_arcogs_hmm.txt | LC_ALL=C sort > arcogs/temp3

#merge with contig list
LC_ALL=C join -a1  -j1 -e'-' -t $'\t' -o 0,2.2,2.3 <(LC_ALL=C sort 1_Bins_to_protein_list.txt) <(LC_ALL=C sort arcogs/temp3) | LC_ALL=C sort  > arcogs/temp4

#merge with arcog names
LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\t' -o1.1,0,2.3,2.4,2.2,1.3 <(LC_ALL=C sort -k2 arcogs/temp4) <(LC_ALL=C sort -k1 ~/../spang_team/Databases/arCOGs2019/ar14.arCOGdef18.nina.tab) | LC_ALL=C  sort > arcogs/temp5

#add in headers
echo -e "accession\tarcogs\tarcogs_geneID\tarcogs_Description\tPathway\tarcogs_evalue" | cat - arcogs/temp5 > arcogs/C_arcogs.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' arcogs/C_arcogs.tsv B_GenomeInfo.txt > merge/temp_BC.tsv

#clean-up
rm arcogs/temp*
```

## KO search using Hmm’s –\> done

# 

``` bash
#1. Setup directories
mkdir KO_hmm

#2. run hmmsearch against all pfams
~/../spang_team/Scripts/Hmmer/hmmsearchTable /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa  ~/../spang_team/Databases/KO_terms/All_KOs.hmm 40 -E 1e-5 > KO_hmm/All_KO_hmm.txt

#merge with protein list
LC_ALL=C join -a1 -t $'\t' -j1 -o 0,2.3,2.5,2.6 <(LC_ALL=C sort 1_Bins_to_protein_list.txt) <(LC_ALL=C sort KO_hmm/All_KO_hmm.txt) | sort -u -k1 > KO_hmm/temp

#get rid of empty space
awk 'BEGIN {FS = OFS = "\t"} {for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "-" }; 1' KO_hmm/temp > KO_hmm/temp2

#merge with KO_hmm names
LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\t' -o1.1,1.2,1.3,1.4,2.2,2.12  <(LC_ALL=C sort -k2 KO_hmm/temp2) <(LC_ALL=C sort -k1 /export/lv1/user/spang_team/Databases/KO_terms/ko_list_for_mapping) | LC_ALL=C  sort > KO_hmm/temp3

#add in an extra column that lists whether hits have a high confidence score
awk  -v OFS='\t' '{ if ($4 > $5){ $7="high_score" }else{ $7="-" } print } ' KO_hmm/temp3 > KO_hmm/temp4

#add header
echo -e "accession\tKO_hmm\te_value\tbit_score\tbit_score_cutoff\tDefinition\tconfidence" | cat - KO_hmm/temp4 > KO_hmm/M_KO_hmm.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' KO_hmm/M_KO_hmm.tsv merge/temp_BC.tsv > merge/temp_BCM.tsv

#control lines
wc -l merge/*

#clean up
rm KO_hmm/temp*
```

## Pfam search –\> done

# 

``` bash
#1. Setup directories
mkdir pfam

#2. run hmmsearch against all pfams
~/../spang_team/Scripts/Hmmer/hmmsearchTable /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa  ~/../spang_team/Databases/Pfam/Pfam-A.hmm 40 -E 1e-5 > pfam/All_pfam.txt

#remove header to avoid issues with join
sed '1d' pfam/All_pfam.txt > pfam/temp

#cp 4th column into first and remove everything after the dot.
awk -F "\t" -v OFS="\t" '{gsub(/\./, "\t", $4); print}'  pfam/temp > pfam/temp2

#merge with protein list
LC_ALL=C join -a1 -t $'\t' -j1 -o 0,1.2,2.6,2.4 <(LC_ALL=C sort 1_Bins_to_protein_list.txt) <(LC_ALL=C sort pfam/temp2) | sort -u -k1 > pfam/temp3

#get rid of empty space
awk 'BEGIN {FS = OFS = "\t"} {for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "-" }; 1' pfam/temp3  > pfam/temp4

#merge with pfam names
LC_ALL=C join -a1 -1 4 -2 1 -t $'\t' -e'-' -o1.1,0,2.5,1.3 <(LC_ALL=C sort -k4 pfam/temp4) <(LC_ALL=C sort -k1 /export/lv1/user/spang_team/Databases/Pfam/Pfam-A.clans.cleaned.tsv) | LC_ALL=C sort > pfam/temp5

#add in header
echo -e "accession\tPFAM_hmm\tPFAM_description\tPfam_Evalue" | cat - pfam/temp5 > pfam/E_Pfam.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' pfam/E_Pfam.tsv merge/temp_BCM.tsv > merge/temp_BCME.tsv

#control lines
wc -l merge/*

#clean up
rm pfam/temp*
```

## TIRG search –\> done

# 

``` bash
#1. Setup directories
mkdir TIGRs

#2. run hmmsearch
~/../spang_team/Scripts/Hmmer/hmmsearchTable /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa  ~/../spang_team/Databases/TIGRPFAM/TIGRFAMs_15.0_HMM.LIB 40 -E 1e-5 > TIGRs/All_TIGR.txt

#merge with protein list
LC_ALL=C join -a1 -t $'\t' -j1 -o 0,1.2,2.3,2.5 <(LC_ALL=C sort 1_Bins_to_protein_list.txt) <(LC_ALL=C sort TIGRs/All_TIGR.txt) | sort -u -k1 > TIGRs/temp

#get rid of empty space
awk 'BEGIN {FS = OFS = "\t"} {for(i=1; i<=NF; i++) if($i ~ /^ *$/) $i = "-" }; 1' TIGRs/temp > TIGRs/temp2

#merge with pfam names
LC_ALL=C join -a1 -1 3 -2 1 -e'-' -t $'\t' -o1.1,0,2.3,2.4,1.4  <(LC_ALL=C sort -k3 TIGRs/temp2) <(LC_ALL=C sort -k1 /export/lv1/user/spang_team/Databases/TIGRPFAM/TIGR_Info_clean2.txt)| LC_ALL=C  sort > TIGRs/temp3

#add header
echo -e "accession\tTIRGR\tTIGR_description\tEC\tTIGR_Evalue" | cat - TIGRs/temp3 > TIGRs/F_TIGR.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' TIGRs/F_TIGR.tsv merge/temp_BCME.tsv > merge/temp_BCMEF.tsv

#control lines
wc -l merge/*

#clean up
rm TIGRs/temp*
```

## CazyDB search –\> done

# 

``` bash
#modify mapping file
#sed 's/ /_/g' CAZY_mapping.txt > CAZY_mapping_2.txt

#make folder
mkdir CAZYmes

#2. run hmmsearch against all pfams
~/../spang_team/Scripts/Hmmer/hmmsearchTable /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa /export/lv1/user/spang_team/Databases/CAZymes/Database/dbCAN-HMMdb-V7.txt 20 -E 1e-20 > CAZYmes/All_vs_CAZYmes.txt

#remove .hmm
sed s'/\.hmm//g'  CAZYmes/All_vs_CAZYmes.txt > CAZYmes/temp

#remove header to avoid issues with join
sed '1d' CAZYmes/temp > CAZYmes/temp2

#merge with contig list
LC_ALL=C join -a1  -j1 -e'-' -t $'\t' -o 0,2.3,2.5 <(LC_ALL=C sort 1_Bins_to_protein_list.txt) <(LC_ALL=C sort CAZYmes/temp2)  | LC_ALL=C sort  > CAZYmes/temp3

#merge with mapping file
LC_ALL=C join -a1  -1 2 -2 1 -e'-' -t $'\t' -o 1.1,1.2,1.3,2.2 <(LC_ALL=C sort -k2 CAZYmes/temp3) <(LC_ALL=C sort /export/lv1/user/spang_team/Databases/CAZymes/Database/CAZY_mapping_2.txt) | LC_ALL=C sort  > CAZYmes/temp4

#add in headers
echo -e "accession\tCAZy\tCAZy_evalue\tDescription" | cat - CAZYmes/temp4 > CAZYmes/G_CAZy.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' CAZYmes/G_CAZy.tsv merge/temp_BCMEF.tsv > merge/temp_BCMEFG.tsv

#control lines
wc -l merge/*

#clean up
rm CAZYmes/temp*
```

## Merops search

# 

``` bash
#prep db
#cd /export/lv1/user/spang_team/Databases/Merops

#grep ">" Merops_DB.faa
#split -l 600000  Names.txt segment
#scp to deskop and modify in exel (one column Merops ID, one colum description)
#back to ada

#cat segment1_.txt segment2_.txt > Merops_Description.txt
#cut -f1 -d " " Merops_DB.faa > Merops_DB_short.faa
#sed 's/"//g' Merops_Description.txt > Merops_Description2.txt

#make blastdb
#makeblastdb -in Merops_DB_short.faa -out Merops_DB_short   -dbtype prot -parse_seqids

#make folder
mkdir Merops

#2. run blast against all MeropsDB
blastp -num_threads 40 -outfmt 6 -query /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa -db /export/lv1/user/spang_team/Databases/Merops/Merops_DB_short -out Merops/All_vs_Merops.tsv -evalue 1e-20

#find best hit
perl ~/../spang_team/Scripts/Others/best_blast.pl Merops/All_vs_Merops.tsv Merops/temp

'''
Total # records = 766728
Best only # records = 7844
'''

#merge with contig list
LC_ALL=C join -a1  -j1 -e'-' -t $'\t' -o 0,2.2,2.11 <(LC_ALL=C sort 1_Bins_to_protein_list.txt) <(LC_ALL=C sort Merops/temp) -t $'\t' | LC_ALL=C sort  > Merops/temp2

#merge with Merops names
LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\t' -o1.1,0,2.2,1.3  <(LC_ALL=C sort -k2  Merops/temp2) <(LC_ALL=C sort -k1 /export/lv1/user/spang_team/Databases/Merops/Merops_Description2.txt) | LC_ALL=C  sort > Merops/temp3

#add in headers
echo -e "accession\tHMerops\tMeropsDescr\tMerops_evalue" | cat - Merops/temp3 > Merops/H_Merops.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' Merops/H_Merops.tsv merge/temp_BCMEFG.tsv > merge/temp_BCMEFGH.tsv

#control stuff
wc -l merge/*tsv

#clean up
rm Merops/temp*
```

## Transporter DB search –\> done

# 

[wget files from here](http://www.tcdb.org)

``` bash
#cd /export/lv1/user/spang_team/Databases/TransporterDB

#sed 's/gnl|TC-DB|//g' tcdb.faa > tcdb_renamed.faa
#cut -f1 -d " " tcdb_renamed.faa > tcdb_renamed_short.faa

#not used
#sed 's/>[0-9]*|/>/' tcdb_renamed.faa > tcdb_renamed2.faa

#not used
#join -a1  -1 1 -2 1 -e'-' -o 1.2,0,1.3,2.2 <(LC_ALL=C sort TCDB_Desc_1.txt) <(LC_ALL=C sort TCDB_Desc_2.txt) -t $'\t' | LC_ALL=C sort  > TCDB_Desc_final.txt

#make blastdb
#makeblastdb -in tcdb_renamed_short.faa -out tcdb_renamed_short -dbtype prot -parse_seqids

#make folder
mkdir TransporterDB

#2. run blast against all MeropsDB
blastp -num_threads 10 -outfmt 6 -query /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa -db /export/lv1/user/spang_team/Databases/TransporterDB/tcdb_renamed_short -out TransporterDB/All_vs_TPDB.tsv -evalue 1e-20

#find best hit
perl ~/../spang_team/Scripts/Others/best_blast.pl TransporterDB/All_vs_TPDB.tsv TransporterDB/temp

'''
Total # records = 26770500
Best only # records = 979506
'''

#split ids
awk -F'\t' -v OFS='\t' '{split($2,a,"|"); print $1, a[1], a[2], $11}' TransporterDB/temp| LC_ALL=C sort > TransporterDB/temp2

#merge with contig list
LC_ALL=C join -a1  -j1 -e'-' -t $'\t' -o 0,2.3,2.4 <(LC_ALL=C sort 1_Bins_to_protein_list.txt) <(LC_ALL=C sort TransporterDB/temp2)  | LC_ALL=C sort  > TransporterDB/temp3

#merge with TPDB names
#join -a1 -1 3 -2 2 -e'-' -o1.1,0,2.3,1.4  <(LC_ALL=C sort -k1 TransporterDB/temp3) <(LC_ALL=C sort -k2 /export/lv1/user/spang_team/Databases/TransporterDB/TCDB_Desc_final.txt) -t $'\t' | column -t | LC_ALL=C  sort > TransporterDB/temp4

#add in headers
echo -e "accession\tTDBD_ID\tTPDB_evalue" | cat - TransporterDB/temp3 > TransporterDB/I_TPDB.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' TransporterDB/I_TPDB.tsv merge/temp_BCMEFG.tsv > merge/temp_BCMEFGHI.tsv

#control lines
wc -l merge/*

#clean up
rm TransporterDB/temp*
```

## HydDB search –\> done

# 

``` bash

#preparation
#done in: /export/lv1/user/spang_team/Databases/HydDB

#remove duplicate headers
#awk 'BEGIN{RS=">"}NR>1{sub("\n","\t"); gsub("\n",""); print RS$0}' HydDB.faa | awk '!seen[$1]++' | awk -v OFS="\n" '{print $1,$2}' > HydDB_uniq.faa

#make blastdb
#makeblastdb -in HydDB_uniq.faa -out HydDB_uniq   -dbtype prot -parse_seqids

mkdir HydDB

blastp -num_threads 20 -outfmt 6 -query /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa -db /export/lv1/user/spang_team/Databases/HydDB/HydDB_uniq -out HydDB/All_vs_HydDB.tsv -evalue 1e-20

#find best hit
perl ~/../spang_team/Scripts/Others/best_blast.pl HydDB/All_vs_HydDB.tsv HydDB/temp

'''
Total # records = 6751109
Best only # records = 109635
'''

#merge with contig list
LC_ALL=C join -a1  -j1 -e'-'  -t $'\t' -o 0,2.2,2.3,2.11,2.12 <(LC_ALL=C sort 1_Bins_to_protein_list.txt) <(LC_ALL=C sort HydDB/temp) | LC_ALL=C sort  > HydDB/temp2

#merge with HydDB names
LC_ALL=C join -a1 -1 2 -2 1 -e'-' -t $'\t' -o1.1,0,2.2,1.4  <(LC_ALL=C sort -k2  HydDB/temp2) <(LC_ALL=C sort -k1 /export/lv1/user/spang_team/Databases/HydDB/HydDB_mapping)  | LC_ALL=C  sort > HydDB/temp3

#add in headers
echo -e "accession\tHydDB\tDescription\tHydDB_evalue" | cat - HydDB/temp3 > HydDB/J_HydDB.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' HydDB/J_HydDB.tsv merge/temp_BCMEFGHI.tsv > merge/temp_BCMEFGHIJ.tsv

#control lines
wc -l merge/*

#clean up
rm HydDB/temp*
```

## IPRscan –\> done

# 

``` bash
mkdir IPRscan

#needed modifications for prokka files
cd ../../Prokka/V2
mkdir split_faa
cd split_faa

python ~/../spang_team/Scripts/Others/Split_Multifasta.py -m ../All_Genomes_clean.faa -n 5000
#9785373 split into 1958 files

mkdir split_faa
mv split_faa/ ../../Database_Annotations/v2
cd ../../Database_Annotations/v2

#setup parallel once before running IPRscan

#iprscan: (interproscan-5.31-70.0)
#"cite:Tange (2011): GNU Parallel - The Command-Line Power Tool, ;login: The USENIX Magazine, February 2011:42-47."

#test if all files are  grapped
parallel -j10 'i={}; echo $i' ::: split_faa/File*.faa

#run search on laplace
nano interproscan.sh

'''
#!/bin/sh
#SBATCH --partition=LONG1
#SBATCH --nodelist=no71
#SBATCH --nodes=1   # require 1 nodes
#SBATCH --ntasks-per-node=80  # (by default, "ntasks"="cpus")
#SBATCH --error=job.%J.err
#SBATCH --output=job.%J.out
# The above two lines reflect file locations for standard error and output
# Executable commands :

parallel -j20 'i={}; /export/lv1/user/spang_team/Scripts/Interproscan/interproscan-5.31-70.0/interproscan.sh -i $i -d IPRscan/ -T IPRscan/temp --iprlookup --goterms' ::: /export/lv1/user/spang_team/Projects/BlackSea18/Database_Annotations/v2/split_faa/File*.faa

'''

#start job
sbatch interproscan.sh

#Concat result files
cd IPRscan
mkdir single_files
mv File* single_files
mkdir Concat_results/

cat single_files/File*.faa.xml > Concat_results/All_bins_iprscan-results.xml
cat single_files/File*.faa.gff3 > Concat_results/All_bins_bins_iprscan-results.gff3
cat single_files/File*.faa.tsv > Concat_results/All_bins_bins_iprscan-results.tsv

#cleanup
rm single_files/File*

cd ..

#clean up fasta header so that it is exactly the same as the accession ID in the interproscan results
python ~/../spang_team/Scripts/Others/parse_IPRdomains_vs2_GO_2.py -s /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa -i IPRscan/Concat_results/All_bins_bins_iprscan-results.tsv -o IPRscan/All_bins_bins_iprscan-results_parsed.tsv

#remove issue with spaces
sed 's/ /_/g' IPRscan/All_bins_bins_iprscan-results_parsed.tsv | LC_ALL=C  sort > IPRscan/temp.txt

#print only columns of interest
awk -F'\t' -v OFS="\t"  '{print $1, $2,$3,$4, $5}' IPRscan/temp.txt > IPRscan/temp2

#add header
echo -e "accession\tPFAM\tPFAMdescription\tIPR\tIPRdescription" | cat - IPRscan/temp2 > IPRscan/K_IPR.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' IPRscan/K_IPR.tsv merge/temp_BCMEFGHIJ.tsv > merge/temp_BCMEFGHIJK.tsv

#control lines
wc -l merge/*

#clean up
rm IPRscan/temp*
```

## Diamond against NCBI NR

# 

``` bash
mkdir NCBI_NR

#run diamond against NR database
diamond blastp -q /export/lv1/user/spang_team/Projects/BlackSea18/Prokka/V2/All_Genomes_clean.faa --more-sensitive --evalue 1e-5 --threads 20 --seq 50 --no-self-hits --db /export/data01/databases/ncbi_nr/diamond/nr.dmnd --taxonmap /export/data01/databases/taxmapfiles/ncbi_nr/prot.accession2taxid.gz --outfmt 6 qseqid qtitle qlen sseqid salltitles slen qstart qend sstart send evalue bitscore length pident staxids -o NCBI_NR/All_NCBInr.tsv

#Select columns of interest in diamond output file
awk -F'\t' -v OFS="\t" '{ print $1, $5, $6, $11, $12, $14, $15 }'  NCBI_NR/All_NCBInr.tsv > NCBI_NR/temp

#Parse Diamnond Results
python ~/../spang_team/Scripts/Others/parse_diamond_blast_results_id_taxid.py -i AllProteins_list.txt -d NCBI_NR/temp -o NCBI_NR/temp2

#rm header
sed 1d NCBI_NR/temp2 > NCBI_NR/temp3

#add an '-' into empty columns
awk -F"\t" '{for(i=2;i<=NF;i+=2)gsub(/[[:blank:]]/,"_",$i)}1' OFS="\t" NCBI_NR/temp3 > NCBI_NR/temp4

#the python script above sometimes leaves an empty 7th column, this gets rid of that issue
awk -F'\t' -v OFS="\t"  '{if (!$7) {print $1,$2, $4 , $6, "-"} else {print $1, $2, $4, $6, $7}}' NCBI_NR/temp4 | LC_ALL=C sort > NCBI_NR/temp5

#split columns with two tax ids
awk -F'\t' -v OFS='\t' '{split($5,a,";"); print $1, $2, $3, $4, a[1]}' NCBI_NR/temp5 > NCBI_NR/temp6

#in column 2 remove everything after < (otherwise the name can get too long)
awk -F'\t' -v OFS='\t' '{split($2,a,"<"); print $1, a[1], $3, $4, $5}' NCBI_NR/temp6 > NCBI_NR/temp6a

#merge with tax names
LC_ALL=C join -a1 -1 5 -2 1 -e'-' -t $'\t'  -o1.1,1.2,1.3,1.4,1.5,2.2  <(LC_ALL=C sort -k5  NCBI_NR/temp6a) <(LC_ALL=C sort -k1 ~/../spang_team/Databases/NCBI_taxonomy/Jul2019/taxonomy5.txt ) | LC_ALL=C  sort > NCBI_NR/temp7

#add in header
echo -e "accession\tTopHit\tE_value\tPecID\tTaxID\tTaxString" | cat - NCBI_NR/temp7 > NCBI_NR/L_Diamond.tsv

#combine with previous dataframe
awk 'BEGIN{FS="\t";OFS="\t"}FNR==NR{a[$1]=$0;next}{print $0,a[$1]}' NCBI_NR/L_Diamond.tsv merge/temp_BCMEFGHIJK.tsv > merge/temp_BCMEFGHIJKL.tsv

#control lines
wc -l merge/*

#clean up
rm NCBI_NR/temp*
```

# 

# 

# Modify final dataframe

# 

# 

``` bash
#rm redundant headers and reconstruct accession ID
awk 'BEGIN{FS="\t";OFS="\t"}NR==1{for(x=2;x<=NF;x++)if($x!="accession")l[x]++;}{for(i=1;i<=NF;i++)if(i in l)printf (i==NF)?$i"":$i"\t";printf "\n"}' merge/temp_BCMEFGHIJKL.tsv > merge/temp_1.tsv

#reconstruct one accession ID
awk 'BEGIN{FS="\t";OFS="\t"}{print $1"-"$9,$0}' merge/temp_1.tsv > merge/temp_2.tsv

#change column name
awk 'BEGIN{FS="\t";FS="\t"; OFS="\t"}{if(NR==1) $1="accession"} {print $0 }' merge/temp_2.tsv > merge/Annotations_BS18_v2.txt

#control for a bin
#grep "NIOZ128_mb_b154_1" merge/Annotations_BS18_IlluminaBins_v1.txt > merge/NIOZ128_mb_b154_1.txt

#control first 50 columns
head -50 merge/Annotations_BS18_v2.txt > merge/test.txt

#clean up
rm merge/temp_[12].tsv
gzip merge/*tsv
```

# 

# 

# Subset final dataframe

# 

# 

## Extract hits for archaea only

``` bash
mkdir FileList

nano FileLists/DPANN_list
nano FileLists/Archaea_list

#subset annotation list
mkdir merge/subset

grep -f FileLists/Archaea_list merge/Annotations_BS18_v2.txt > merge/subset/Annotations_archaea.txt

```
