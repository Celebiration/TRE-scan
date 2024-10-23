#!/usr/bin/bash
#source activate bam-readcount
set -e
Usage()
{
    echo "bam-readcountm.sh -i {input.bam} -q {min_mapping_quality} -b {min_base_quality} -f {reference.fa} -l {site_file} -m {threads} -o {output}"
}
while getopts ':i:q:b:f:l:m:o:' OPT; do
    case $OPT in
        i) input="$OPTARG";;
        q) min_mapping_quality="$OPTARG";;
        b) min_base_quality="$OPTARG";;
        f) reference="$OPTARG";;
        l) site_file="$OPTARG";;
        m) thread="$OPTARG";;
        o) output="$OPTARG";;
        *) Usage; exit 1;;
    esac
done
if [ -z $input ];then Usage; exit 1; fi
if [ -z $min_mapping_quality ];then Usage; exit 1; fi
if [ -z $min_base_quality ];then Usage; exit 1; fi
if [ -z $reference ];then Usage; exit 1; fi
if [ -z $site_file ];then Usage; exit 1; fi
if [ -z $thread ];then Usage; exit 1; fi
if [ -z $output ];then Usage; exit 1; fi

mkdir -p .tmp_${input}.${site_file}
num=`wc -l $site_file|awk '{print $1}'`
ll=`expr $num / $thread`
if [ `expr $ll % $thread` -gt 0 ];then ((ll++));fi
split --numeric-suffixes -l $ll $site_file .tmp_${input}.${site_file}/${site_file}
file_list=`ls .tmp_${input}.${site_file}/${site_file}*`
for ii in $file_list
do
    {
        bam-readcount -q $min_mapping_quality -b $min_base_quality -f $reference $input -l $ii -w 0 > ${ii}_out
        rm $ii
    }&
done
wait
status=0
for ii in $file_list
do
    if [ ! -e ${ii}_out ];then echo "Error: output ${ii}_out doesn't exist.";status=1;fi
done
if [ $status == 0 ];then echo "merging files...";else exit 1;fi

ls .tmp_${input}.${site_file}/${site_file}*_out|sort|xargs cat > $output
rm -rf .tmp_${input}.${site_file}