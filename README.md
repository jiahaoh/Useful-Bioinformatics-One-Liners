# Useful Bioinformatics One Liner Collection

Here are useful bash one liners collected from various sources shown down below, I also added some of my own tricks. Please start a pull request if you would like to add yours on the list. Thank you! - Jiahao

[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/jiahaoh/Useful-Bioinformatics-One-Liners/graphs/commit-activity)
[![HitCount](http://hits.dwyl.io/jiahaoh/Useful-Bioinformatics-One-Liners.svg)](http://hits.dwyl.io/jiahaoh/Useful-Bioinformatics-One-Liners)

## Source

* [The Williams lab @ Cornell University](http://williamslab.bscb.cornell.edu/?page_id=235)
* [Stephen Turner's Github](https://github.com/stephenturner/oneliners)
* [Ming Tang's Github](https://github.com/crazyhottommy/bioinformatics-one-liners)
* [Zhigang Lu's Blog](https://zhiganglu.com/post/one-liners-collection/)
* [竹子-博客](http://www.cnblogs.com/peida/tag/%E6%AF%8F%E6%97%A5%E4%B8%80linux%E5%91%BD%E4%BB%A4/)
* [Seqtk](https://github.com/lh3/seqtk)
* [Bioawk](https://github.com/lh3/bioawk)
* [Seqkit](https://bioinf.shenwei.me/seqkit/)
* [生信技能树 BiliBili](https://www.bilibili.com/video/av28813815/?p=1)
* [UCSC Data File Formats](<https://genome.ucsc.edu/FAQ/FAQformat.html#format1>)

## Content

|   **Sortware**    |      **Format**     |
| :---------------: | :-----------------: |
|    [awk](#awk)    | [FASTA/Q](#fasta-q) |
|    [sed](#sed)    | [SAM/BAM](#sam-bam) |
|   [grep](#grep)   |     [VCF](#vcf)     |
|    [tar](#tar)    | [GFF/GTF](#gff-gtf) |
|   [perl](#perl)   |     [BED](#bed)     |
| [Bioawk](#bioawk) |     [PSL](#psl)     |
|  [Seqtk](#seqtk)  |     [WIG](#wig)     |
| [SeqKit](#seqkit) |    :construction:   |

|     **Other**     |
| :---------------: |
|  [alias](#alias)  |
| [tricks](#tricks) |


## awk

Print line xx content
```bash
awk '(NR==xx){print $0}' input_file
```

Calculate **sum** of column 1
```bash
awk '{sum+=$1} END {print sum}' input_file
```

Calculate **average** of column 1
```bash
awk '{sum+=$1}END{print sum/NR}' input_file
```

Combine multiple rows based on the same column value
```bash
awk '$1!=p{if(p)print s; p=$1; s=$0; next}{sub(p,x); s=s $0} END{print s}' input_file
```

Split multiple columns into two columns with common first column
```bash
awk -v OFS='\t' '{for (i=2;i<=NF;i++) print $1,$i}' input_file
```

Get a range of text
```bash
awk '/START-WORD/, /END-WORD/' input_file > output_file
```
<br>

## sed

Trim leading and tailing white space
```bash
sed 's/^[ \t]*//;s/[ \t]*$//' input_file
```

Delete blank line
```bash
sed '/^$/d' input_file
```

Delete everything after target
```bash
sed -n '/target/,$!p' input_file
```

Delete everything before target
```bash
sed -n '/target/,$p' input_file
```

Get a range of text
```bash
sed -n "/START-WORD-HERE/,/END-WORD-HERE/p" input_file > output_file
sed "/START-WORD-HERE/,/END-WORD-HERE/!d" input_file > output_file
```

Print line xx content
```bash
sed -n xxp input_file
```
<br>

## grep

Search with a file containing queries
```bash
grep -f query_file input_file
```

Search for any string in all txt files
```bash
grep -r 'STRING1\|STRING2\|STRING3' input_file
```
<br>

## tar

Create tar, gzip tar, bz2 Tar
```bash
tar -cvf tar_name.tar input_file
tar -zcvf tar_name.tar.gz input_file
tar -jcvf tar_name.tar.bz2 input_file
```

Extract
```bash
tar -xvf tar_name.tar file_name
tar -zxvf tar_name.tar.gz file_name
tar -jxvf tar_name.tar.bz2 file_name
```

List
```bash
tar -tvf tar_name.tar
tar -ztvf tar_name.tar.gz
tar -jtvf tar_name.tar.bz2
```
<br>

## perl

Reverse complement of seq
```bash
echo <SEQUENCE> | perl -nle 'print map{$_ =~ tr/ACGT/TGCA/; $_} reverse split("",$_)'
```
<br>

## Bioawk

:construction: ...

<br>

## Seqtk

Reverse complement of fasta/q
```bash
seqtk seq -r input.fq > output.fq
```

Extract sequences with names in file name.lst, one sequence name per line
```bash
seqtk subseq input.fq name.lst > output.fq
```

Extract sequences in regions contained in file reg.bed
```bash
seqtk subseq input.fa reg.bed > output.fa
```

Trim low-quality bases from both ends using the Phred algorithm
```bash
seqtk trimfq input.fq > output.fq
```

Trim 5bp from the left end of each read and 10bp from the right end
```bash
seqtk trimfq -b 5 -e 10 input.fa > output.fa
```

Convert fastq to fasta
```bash
seqtk seq -a input.fq.gz > output.fa
```
<br>

## Seqkit
Simple statistics of fasta/q
```bash
seqkit stats *.f{a,q}.gz
```

Calculate GC content
```bash
seqkit fx2tab -ignlH *.f{a,q}.gz
```

<br>

## FASTA/Q <a name="fasta-q"/>

Convert fastq to fasta
```bash
sed -n '1~4s/^@/>/p;2~4p' input.fq > output.fa

seqtk seq -a input.fq.gz > output.fa
```

Grep fastq reads containing a pattern but maintain the fastq format
```bash
grep -A 2 -B 1 'PATTERN' input.fq | sed '/^--$/d' > output.fq
```

Output sequence name and length
```bash
cat input.fa | awk '$0 ~ ">" {print c; c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }'
```

Get the sequence length distribution from fastq
```bash
zcat input.fastq.gz | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'
```

Split a multi-fasta file into individual fasta files
```bash
awk '/^>/{s=++d".fa"} {print > s}' input.fa
```

Calculate GC content for each seq
```bash
awk ' \
BEGIN { \
    FS=""; \
    cg=0; \
    t=0; \
    print "Header""\t""GC#""\t""Total#""\t""GC%"; \
} \
{ \
    if ($1 != ">") { \
        for (i = 1; i <= NF; i++) { \
            if ($i ~ /[ACTGactg]/) { \
                t++;
            } \
            if ($i ~ /[CGcg]/) { \
                cg++;
            } \
        } \
    } \
    else { \
        if (t > 0) { \
            print h"\t"cg"\t"t"\t"(cg/t); \
            cg = 0; \
            t = 0; \
        } \
        h = substr($0,2); \
    } \
} \
END { \
    print h"\t"cg"\t"t"\t"(cg/t); \
}' input.fa
```
<br>

## SAM/BAM <a name="sam-bam"/>

Convert bam to fastq
```bash
samtools view input.bam | awk 'BEGIN {FS="\t"} {print "@" $1 "\n" $10 "\n+\n" $11}' > output.fq
```

Convert bam to bed
```bash
samtools view input.bam | perl -F'\t' \
-ane '$strand=($F[1]&16)?"-":"+";$length=1;$tmp=$F[5];$tmp =~ s/(\d+)[MD]/$length+=$1/eg;print "$F[2]\t$F[3]\t".($F[3]+$length)."\t$F[0]\t0\t$strand\n";' > output.bed
```

Convert bam to wig
```bash
samtools mpileup -BQ0 input.sorted.bam | \
perl -pe '($c, $start, undef, $depth) = split;if ($c ne $lastC || $start != $lastStart+1) {print "fixedStep chrom=$c start=$start step=1 span=1\n";}$_ = $depth."\n";($lastC, $lastStart) = ($c, $start);' | \
gzip -c > output.wig.gz
```

Index your bam files in parallel (GNU parallel required)
```bash
find *.bam | parallel 'samtools index {}'
```
<br>

## VCF

Extract PASS calls from vcf file
```bash
cat input.vcf | awk -F '\t' '{if($0 ~ /\#/) print; else if($7 == "PASS") print}' > output_PASS.vcf
```

Sort vcf file with header
```bash
cat input.vcf | awk '$0~"^#" { print $0; next } { print $0 | "sort -k1,1V -k2,2n" }'
```

Convert vcf to bed
```bash
sed -e 's/chr//' input.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' > output.bed
```
<br>

## GFF/GTF <a name="gff-gtf"/>

Extract all gene IDs from a GFF3 file
```bash
cat input.gff3 | grep $'\tgene\t' | perl -ne '/ID=([^;]+)/ and printf("%s\n", $1)'
```

Print all sequences annotated in a GFF3 file
```bash
cut -s -f 1,9 input.gff3 | grep $'\t' | cut -f 1 | sort | uniq
```

Determine all feature types annotated in a GFF3 file
```bash
grep -v '^#' input.gff3 | cut -s -f 3 | sort | uniq
```

Determine the number of genes annotated in a GFF3 file
```bash
grep -c $'\tgene\t' input.gff3
```

Print length of each gene in a GFF3 file
```bash
grep $'\tgene\t' input.gff3 | cut -s -f 4,5 | perl -ne '@v = split(/\t/); printf("%d\n", $v[1] - $v[0] + 1)'
```
<br>

## BED

Split a bed file by chromosome
```bash
cat input.bed | \
sort -k1,1 -k2,2n | \
sed 's/^chr//' | \
awk '{close(f);f=$1}{print > f".bed"}'
```
<br>

## PSL

:construction: ...

<br>

## WIG

:construction: ...

<br>

## Alias

Show the PATH
```bash
alias path='echo -e ${PATH//:/\\n}'
```

CDs
```bash
alias ..='cd ..'
alias ...='cd ../../'
alias ....='cd ../../../'
alias .....='cd ../../../../'
alias ......='cd ../../../../../'
alias u="cd ..;ls"
```

Ask before execute
```bash
alias mv="mv -i"
alias cp="cp -i"
alias rm="rm -i"
```

Refresh/Edit .bashrc/.bash_profile
```bash
alias refresh="source ~/.bashrc"
alias eb="vi ~/.bashrc"


alias refresh="source ~/.bash_profile"
alias eb="vi ~/.bash_profile"
```

Clear
```bash
alias c="clear"
```

Count seq number for FASTA
```bash
alias countfa="grep -c '^>'"
```

Count seq number for FASTQ
```bash
alias countfq="bioawk -cfastx 'END{print NR}'"
```
<br>

## Tricks

Extract function as suggested by Mendel Cooper in "Advanced Bash Scripting Guide"
```bash
extract () {
   if [ -f $1 ] ; then
       case $1 in
        *.tar.bz2)      tar xvjf $1 ;;
        *.tar.gz)       tar xvzf $1 ;;
        *.tar.xz)       tar Jxvf $1 ;;
        *.bz2)          bunzip2 $1 ;;
        *.rar)          unrar x $1 ;;
        *.gz)           gunzip $1 ;;
        *.tar)          tar xvf $1 ;;
        *.tbz2)         tar xvjf $1 ;;
        *.tgz)          tar xvzf $1 ;;
        *.zip)          unzip $1 ;;
        *.Z)            uncompress $1 ;;
        *.7z)           7z x $1 ;;
        *)              echo "don't know how to extract '$1'..." ;;
       esac
   else
       echo "'$1' is not a valid file!"
   fi
}
```

mkdir then cd
```bash
function mcd { mkdir -p "$1" && cd "$1";}
```

Reverse complement
```bash
echo <SEQUENCE> | rev | tr 'ACTG' 'TGAC'
```

Split a bed file by chromosome
```bash
cat nexterarapidcapture_exome_targetedregions_v1.2.bed | \
sort -k1,1 -k2,2n | \
sed 's/^chr//' | \
awk '{close(f);f=$1}{print > f".bed"}'
```

Rename files
```bash
for file in *gz
do zcat $file > ${file/bed.gz/bed} done
```

Make dir with current date
```bash
mkdir $(date +%F)
```

Loop through all chromosomes
```bash
for i in {1..22} X Y
do
  echo $i
done
```

For loop usage
```bash
for chr in {1..22}; do ./command data_chr$chr.txt; done
for ext in bed bim fam; do ./command data.$ext; done
for file in *.txt; do ./command $file; done
```

Rename all .txt files to .bak
```bash
find . -name "*.txt" | sed "s/\.txt$//" | xargs -i echo mv {}.txt {}.bak
```

Run FASTQC in parallel 12 jobs at a time
```bash
find *.fq | parallel -j 12 "fastqc {} --outdir ."
```
