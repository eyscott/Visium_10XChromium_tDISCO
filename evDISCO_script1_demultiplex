#####run fastq header parser###

module load nixpkgs/16.09
module load python/3.7.0

source ./python_env/bin/activate
for f in $(ls *R2_001.fastq | awk -F'[_]' '{print $1 "_" $2 "_" $3 "_" $4}');
do python ./sequence_fixer.py --fq1 ${f}_R1_001.fastq --fq2 ${f}_R2_001.fastq --barcodes barcodes_v3.yaml; done 

##trimming DISCO reads with fastp
module load StdEnv/2020
module load nixpkgs/16.09
module load fastp/0.20.1

for f in $(ls *R2_001_demultiplexed.fastq | awk -F'[-_.]' '{print $1 "_" $2 "_" $3 "_" $4}');
do fastp -i ${f}_R2_001_demultiplexed.fastq -x -p  \
--stdout ${f} -o ${f}_R2_demultiplexed_trimmed.fastq;
done

##Align with Star###
for f in $(ls *R2_demultiplexed_trimmed.fastq | awk -F'[-_.]' '{print $1 "_" $2 "_" $3 "_" $4}');
do STAR --runThreadN 8 --genomeDir ../mouse_genome \
--readFilesIn ${f}_R2_demultiplexed_trimmed.fastq \
--outFileNamePrefix ${f}_starA \
--outSAMtype BAM SortedByCoordinate --outFilterScoreMinOverLread 0 --outFilterMatchNminOverLread 0 \
--outFilterMatchNmin 0 --outFilterMismatchNmax 2; done

########convert these .bam files into .sam files########
for file in ./*.bam
do
    echo $file 
    samtools view -h $file > ${file/.bam/.sam}
done

####run UMI and Sam parser ####
module load nixpkgs/16.09
module load python/3.7.0

source ./python_env/bin/activate
for f in $(ls *.sam | awk -F'[-_.]' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5 "." $6 "." $7 }');
do echo python umi_parser.py --sam_file ${f}.sam ;done

for f in $(ls *_unique_only.sam | awk -F'[-_.]' '{print $1 "_" $2 "_" $3 "_" $4 "_" $5 "." $6 "." $7 }');
python sam_barcode_demultiplexer.py --sam_file ${f}_unique_only.sam;done

##Make gene tables with FeatureCounts

module load StdEnv/2020
module load r/4.1.2    

R CMD BATCH featureCounts.R 

