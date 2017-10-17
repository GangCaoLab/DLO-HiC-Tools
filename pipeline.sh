#!/bin/bash


set -e

help=$(cat << EOF
Pipeline for analyze DLO HiC data.

usage: pipeline -i INPUT [-o OUTPUT] [-h]

optional arguments:
~~~~~~~~~~~~~~~~~~~
-i      input 
-o      output directory, store the result and intermediate results
-h      show help

Requirements:
~~~~~~~~~~~~~
UNIX like OS
GNU parallel
Fastx-toolkit (http://hannonlab.cshl.edu/fastx_toolkit/download.html)
bedtools 
samtools 
BWA

Python packages:
    numpy, matplotlib
    regex
    HTSeq (http://www-huber.embl.de/HTSeq/doc/overview.html)
    hiclib/iced (https://github.com/hiclib/iced)

EOF
)

bash check_requirement.sh # check requirements

## configs:
# ~~~~~~~~
# data directories. 
# NOTE: all data file must in fastq format, if it's compressed please extract it firstly.
# NOTE: please edit "./linkers.json", it provide the linker info to the "./extract_hookers.py"
# FASTQ ASCII quality standard
QUA=33
# P7 adapter sequence
adapter="GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGCCGTCTTCTGCTTG"
# information about preprocessing
prep=1 # if need to preprocessing 1, else 0
# information about trim barcode
trim_barcode=0      # if need to trim barcode 1, else 0
len_barcode_left=4
len_barcode_right=3
# restriction sites bed file
rest_sites_file="./data/ref_genome/rest_site.position.bed"
# BWA index file
BWA_index="./data/BWA_index/genome.fa"
# linkers config file
linkers="./data/linkers.json"
# chromosome infomation for creat hic result matrix
export chromosome_length="./data/ref_genome/chromosome_length.txt"
# thread num
thread_num=2
# sub part size
sub_part_size=10000
##

## parse arguments
while getopts "h?i:o:" opt; do
    case "$opt" in
        h|\?)
            echo "$help"
            exit 0
            ;;
        i)
            input=$OPTARG
            ;;
        o)
            tmp_dir=$OPTARG
            ;;
    esac
done

if [ -z $tmp_dir ];then
    echo "-o (output directory) required"
    exit 1
fi
if [ -z $input ];then
    echo "-i (input fastq file) required"
    exit 1
fi

## make intermediate directory
if [ -d $tmp_dir ];then
    rm -r $tmp_dir
fi
mkdir $tmp_dir


## STEP1 ##
## preproccessing: remove adapter

# get identifer of the target
if [[ "$input" =~ .*fq$ ]]; then # get basename of file
    id=$(basename "$input" .fq)
else
    id=$(basename "$input" .fastq)
fi

# clipper adapter
if [ $prep -gt 0 ];then
    echo "[info][$id] preproccessing: remove adapter"
    cat $input \
        | parallel --no-notice -k --jobs "$thread_num" -l "$sub_part_size" \
            --recstart @ --pipe \
            fastx_clipper -Q"$QUA" -a "$adapter" -i - \
    > $tmp_dir/clean_"$id".fq
    input=$tmp_dir/clean_"$id".fq
fi


## STEP2 ##
## extract hookers(sub part of reads, for do alignment.)
echo "[info][$id] extract hookers"
# extract hookers
mkdir $tmp_dir/sub_parts_tmp
split -l "$sub_part_size" "$input" "$tmp_dir"/sub_parts_tmp/sub_
parallel --no-notice --jobs "$thread_num"\
    "extract_hookers.py -i {} --linker $linkers\
        --left  $tmp_dir/sub_parts_tmp/left{/.}\
        --right $tmp_dir/sub_parts_tmp/right{/.}\
        --check --outtype fastq  --len 50 > /dev/null"\
    ::: "$tmp_dir"/sub_parts_tmp/sub_*
cat "$tmp_dir"/sub_parts_tmp/left*  > "$tmp_dir"/hookers_"$id"_left.fq
cat "$tmp_dir"/sub_parts_tmp/right* > "$tmp_dir"/hookers_"$id"_right.fq
rm -r "$tmp_dir"/sub_parts_tmp

# trim barcode
if [ $trim_barcode -gt 0 ]; then

    cat $tmp_dir/hookers_"$basename"_left.fq  \
        | parallel --no-notice --jobs "$thread_num" --pipe -l "$sub_part_size" -k  \
            fastx_trimmer -Q"$QUA" -f "$len_barcode_left" \
    > $tmp_dir/tmpfile_l

    cat $tmp_dir/hookers_"$basename"_left.fq  \
        | parallel --no-notice --jobs "$thread_num" --pipe -l "$sub_part_size" -k  \
            fastx_trimmer -Q"$QUA" -f "$len_barcode_right" \
    > $tmp_dir/tmpfile_r

    mv $tmp_dir/tmpfile_l $tmp_dir/hookers_"$id"_left.fq
    mv $tmp_dir/tmpfile_r $tmp_dir/hookers_"$id"_right.fq
fi


## STEP3 ##
## mapping hookers use BWA
echo "[info][$id] mapping hookers use BWA"
bwa aln -t "$thread_num" "$BWA_index" $tmp_dir/hookers_"$id"_left.fq > $tmp_dir/"$id"_left.bwa
bwa aln -t "$thread_num" "$BWA_index" $tmp_dir/hookers_"$id"_right.fq > $tmp_dir/"$id"_right.bwa
#parallel --no-notice "bwa samse ${BWA_index} $tmp_dir/${id}_{}.bwa \
    #$tmp_dir/hookers_${id}_{}.fq > $tmp_dir/mapping_${id}_{}.sam" ::: left right
bwa samse ${BWA_index} $tmp_dir/${id}_left.bwa \
    $tmp_dir/hookers_${id}_left.fq > $tmp_dir/mapping_${id}_left.sam
bwa samse ${BWA_index} $tmp_dir/${id}_right.bwa \
    $tmp_dir/hookers_${id}_right.fq > $tmp_dir/mapping_${id}_right.sam
rm $tmp_dir/*.bwa


## STEP4 ##
## filter mapping results
echo "[info][$id] filter mapping results"
# filter out low quality mapping
parallel --no-notice "samtools view -bSq 20\
    ${tmp_dir}/mapping_${id}_{}.sam > ${tmp_dir}/filtered_${id}_{}.bam" ::: left right
# from bam to bed file and sort it
for j in {"left","right"};
do
    bedtools bamtobed -i "$tmp_dir"/filtered_"$id"_"$j".bam \
        | parallel --no-notice --jobs "$thread_num" --pipe -l "$sub_part_size" -k \
            'cut -f 1-4,6' \
        | sort -k 4,4 \
        > "$tmp_dir"/"$id"_"$j".bed
done

# join left and right mapping results to a pair file in bedpe format
join -j 4 "$tmp_dir"/"$id"_left.bed "$tmp_dir"/"$id"_right.bed \
    | parallel --no-notice --jobs "$thread_num" --pipe -l "$sub_part_size" -k \
        "awk 'OFS=\"\t\"{print \$2,\$3,\$4,\$6,\$7,\$8,\$5,\$9,\$1}'"\
    > "$tmp_dir"/pair_"$id".bedpe

# Removing redundancy
cat "$tmp_dir"/pair_"$id".bedpe \
    | parallel --no-notice --jobs "$thread_num" --pipe -l "$sub_part_size" -k \
        'remove_redundancy.py '\
> "$tmp_dir"/nr_pair_"$id".bedpe

# Noise reduction
cat "$tmp_dir"/nr_pair_"$id".bedpe \
    | parallel --no-notice --jobs "$thread_num" --pipe -l "$sub_part_size" -k \
        noise_reduce.py --restriction "$rest_sites_file" \
> "$tmp_dir"/clean_pair_"$id".bedpe


## STEP5 ##
## creat matrix
echo "[info][$id] creat matrix"
creat_matrix.py --chr_len "$chromosome_length" -i "$tmp_dir"/clean_pair_"$id".bedpe\
    -o "$tmp_dir"/origin_"$id" --figure "$tmp_dir"/origin_"$id".png


## STEP6 ##
## ICE adjustment
echo "[info][$id] ICE adjustment"
ice_adjust.py -i "$tmp_dir"/origin_"$id" -o "$tmp_dir"/adj_"$id"\
    --figure "$tmp_dir"/adj_"$id".png

