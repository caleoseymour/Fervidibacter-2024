## Cale Seymour
## University of Nevada, Las Vegas
## July 2023
## Fervidibacteria sacchari PD1 transcriptomics + carbon metabolism pipeline.

## Output directory.
out_dir=FSac-sealed-genome_transcriptome

## Input files.
genome_dir=FSac-sealed-genome
contig_file=FSac-sealed-genome/FSac-sealed-genome.fasta
feature_file=FSac-sealed-genome/GCF_030520105.1.gff3

## Grab read mappings from LIBRARIES file.
grep -A5000 -m1 -e ' SAMPLE SUMMARY' Fersaccriptomics/Transcriptome\ Analysis/LIBRARIES.txt | tail -n +3 > read_mappings.txt.xls

## Build STAR ref genome.
include/STAR --runMode genomeGenerate --runThreadN 1 --genomeDir $genome_dir --genomeSAsparseD 2 --genomeSAindexNbases 9 --genomeFastaFiles $contig_file


mkdir $out_dir
## Iterate over read mappings files and run STAR.
fields=(`head read_mappings.txt.xls -n 1 | sed -e 's/[0-9]*=//g'`)
while IFS=$'\t' read "${fields[@]}"; do
    
    echo "mapping library '$libraryName' using readfile '$fileUsed'..."
    
    ## Unzip raw fastq to temp file.
    gunzip -ckq Fersaccriptomics/Filtered\ Raw\ Data/$fileUsed > .temp.fastq
    
    ## Run STAR to align.
    include/STAR \
        --genomeDir $genome_dir \
        --alignIntronMax 1 \
        --runThreadN 1 \
        --readFilesIn .temp.fastq \
        --outFileNamePrefix $out_dir/$libraryName/$libraryName.
    
    samtools view -S -b $out_dir/$libraryName/$libraryName.Aligned.out.sam > $out_dir/$libraryName/$libraryName.Aligned.out.bam
    
    ## Cleanup. Don't need to keep the old sam-file and temp file around.
    rm $out_dir/$libraryName/$libraryName.Aligned.out.sam
    rm .temp.fastq
    
done <<< $(tail read_mappings.txt.xls -n +2)

## Run featureCounts.
featureCounts -F GFF -t CDS -g ID -a $feature_file -o $out_dir/feature-counts.tsv.xls $(find ./$out_dir/ -type f -name '*.Aligned.out.bam')

