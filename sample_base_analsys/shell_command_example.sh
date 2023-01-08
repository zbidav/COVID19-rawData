#/bin/sh
# Here I will provide the code that we used for analsyis of COVID-19 samples. I will give an example for one sample. The samp steps were repeted for all the samples. 

# download the samples to the server
sra_list="SRR_Acc_List.txt"
out_dir="raw/"
ngc_file="prj_30451.ngc" #obtained from dbgap
cat $list_SRR | while read -r SRR_RUN; do 
	echo $run
	prefetch.2.10.8 -C yes --ngc $NGC_FILE -p $SRR_RUN -O $out_dir > /dev/null 2>&1
	fasterq-dump.2.10.8 --split-files -p -e 20 --ngc $NGC_FILE -O $out_dir ${SRR_RUN}.sra > /dev/null 2>&1
	rm ${SRR_RUN}.sra
done

#FASTQC 
cat SRR13250339_1.fastq > SRR13250339_1_temp.fastq
fastqc -q --extract -o SRR13250339 SRR13250339_1_temp.fastq
cat SRR13250339_2.fastq > SRR13250339_2_temp.fastq
fastqc -q --extract -o SRR13250339 SRR13250339_2_temp.fastq
rm SRR13250339/SRR13250339_*_temp.fastq

#trimmomatic
java -jar trimmomatic-0.33.jar PE -phred33 raw/SRR13250339_1.fastq raw/SRR13250339_2.fastq SRR13250339_1.fastq SRR13250339_1.unpaired.fastq SRR13250339_2.fastq SRR13250339_2.unpaired.fastq ILLUMINACLIP:TruSeq-All.fa:2:30:10:2 LEADING:3 TRAILING:3 MINLEN:40
lzma -7 -T 30 SRR13250339*.fastq

#primseq
lzcat SRR13250339_1.fastq.lzma > SRR13250339_1.tmp.fastq
lzcat SRR13250339_2.fastq.lzma > SRR13250339_2.tmp.fastq 
perl516 prinseq-lite.pl -fastq SRR13250339_1.tmp.fastq -fastq2 SRR13250339_2.tmp.fastq -out_bad null -derep 14 -trim_to_len 100000 -trim_left 0 -trim_right 0 -min_qual_score 0 -min_qual_mean 0  -out_good SRR13250339
rm /SRR13250339/SRR13250376_1.tmp.fastq
rm /SRR13250339/SRR13250376_2.tmp.fastq
find SRR13250339* -type f|xargs cat
find SRR13250339 -empty -type d -delete

#STAR
STAR-2.7.3a --genomeDir hg38.STAR.7.ReadsLn150.gencode28/  --genomeLoad LoadAndExit
STAR-2.7.3a --alignSJoverhangMin 8 --alignIntronMax 1000000 --alignMatesGapMax 600000 --outFilterMismatchNoverLmax 0.3 --outFilterMismatchNoverReadLmax 1 --outFilterMatchNminOverLread  0.66 --outFilterMultimapNmax 1 --outReadsUnmapped Fastx --outSAMattributes All --outSAMtype BAM Unsorted --quantMode GeneCounts  --genomeLoad LoadAndKeep --limitBAMsortRAM 37580963840 --outBAMsortingThreadN 10 --runThreadN 40 --genomeDir /hg38.STAR.7.ReadsLn150.gencode28/ --readFilesCommand cat --readFilesIn SRR13250339_1.fastq SRR13250339_2.fastq --outFileNamePrefix SRR13250339/SRR13250339
samtools-1.9 sort -l 9 -m 1G -@ 30 -o SRR13250339/SRR13250339Aligned.sortedByCoord.out.bam SRR13250339/SRR13250339Aligned.out.bam
rm SRR13250339/SRR13250339Aligned.out.bam
samtools-1.9 index SRR13250339/SRR13250339Aligned.sortedByCoord.out.bam
find SRR13250339 -name '*Unmapped.out.mate?' -or -name '*(STAR_genes_counts_suffix)s' | xargs -I magicstring bash -c 'chmod 775 magicstring'
lzma -7 -T 30 /SRR13250339/SRR13250339ReadsPerGene.out.tab
find SRR13250339 -name '*Unmapped.out.mate?'|xargs -I magicstring bash -c 'lzma -7 -T 30 magicstring'
rm -r SRR13250339/SRR13250339_STARtmp

```
python "STARstatsV2.py" -r $output_directory -o $summeryFolder --plot -po $summeryFolder --force >> $cmdPath
```

#SALMON
mkdir -p SALMON/ResultFiles
bash /home/alu/fulther/PipelineConfigs/BashSupport/runSalmon_1.4.0.sh -1 SRR13250339_1.fastq -2 SRR13250339_2.fastq -o SALMON/SRR13250339 -uc cat -l A -i /SALMON/index/hg38 -t 8 -g /SALMON/index/hg38/gencode_v32.transcriptToGeneID.tab --extra ""
mv SALMON/SRR13250339/quant.sf SALMON/SRR13250339/../ResultFiles/SRR13250339.quant.sf
mv SALMON/SRR13250339/quant.quant.genes.sf SALMON/SRR13250339/../ResultFiles/SRR13250339.quant.genes.sf
rm -r /private3/COVID_dbgap_Nasal/raw/Group1/SALMON_1.4.0/SRR13250339
 
```
IN_DIR="SALMON/ResultFiles/";
OUT_DIR="SALMON/";
mkdir -p $OUT_DIR;
Rscript "uniteSalmonFIles.R" -i $IN_DIR -o $OUT_DIR --wide_genes --wide_samples 

Rscript ZMAD_ScoreCalc_V2.R -i "SALMON/SalmonTPM.wideSamples.csv" -o "SALMON/" -p ArticleCalcLog2 --log2 
```
  
#RNAEditingIndex
## we applied the tool using basal human paramaeters on all Alu elements
IN_DIR="STAR/"
OUT_DIR="RNAEditingIndex/Alu/"
mkdir -p $OUT_DIR
nohup RNAEditingIndex1.1 -d $IN_DIR -o $OUT_DIR -l $OUT_DIR -os $OUT_DIR --genome hg38 --paired_end --keep_cmpileup > ${OUT_DIR}/runAluIndex_Log.txt &

## we applied the tool using basal human paramaeters on 3'UTR Alu elements
IN_DIR="STAR/"
OUT_DIR="RNAEditingIndex/Alu3pUTR/"
REGIOUNS="hg38.Alu3pUTR_minLen200_17022021.RepeatsInOppositeOrientationAtRegions.sorted.merged.bed";
mkdir -p $OUT_DIR
RNAEditingIndex1.1 -d $IN_DIR -o $OUT_DIR -l $OUT_DIR -os $OUT_DIR -rb $REGIOUNS --genome hg38 --paired_end --per_region_output --per_sample_output > ${OUT_DIR}/runAluIndex_Log.txt

#REDI-Known
cp Orshay_AG_list.bed known/SRR13250339/SRR13250339/SRR13250339known_sites_bed.txt

python REDITOOLS/REDItools-1.0.4/reditools/REDItoolKnownHillel.py -i STAR/SRR13250339/SRR13250339Aligned.sortedByCoord.out.bam -f hg38.fa -l known/SRR13250339/SRR13250339/SRR13250339_sites_bed.txt -v 1 -n 0.001 -r 1 -T5-5 -c 20 -Q 33 -o known/SRR13250339/SRR13250339/SRR13250339_Folder

mv known/SRR13250339/SRR13250339/SRR13250339_Folder/*/outTable_* known/SRR13250339/SRR13250339/SRR13250339.REDItoolKnown.out.tab
rm /known/SRR13250339/SRR13250339/SRR13250339known_sites_bed.txt*
rm -r known/SRR13250339/SRR13250339/SRR13250339_Folder
