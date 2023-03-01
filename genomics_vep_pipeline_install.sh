#!/bin/bash

# This file should be sourced
# Please run the below command in the shell
# printf "y\ny\n466\ny\n123\n" | . maf2maf.sh

REFDIR=$HOME/genomics_project

echo "creating working directory at $REFDIR"

mkdir -p $REFDIR

echo "change directory to $REFDIR"

cd $REFDIR

echo "downloading and installing samtools-1.15 into $REFDIR"

wget https://github.com/samtools/samtools/releases/download/1.15/samtools-1.15.tar.bz2
tar -xf samtools-1.15.tar.bz2
cd samtools-1.15/
./configure --prefix=$REFDIR
make
make install
cd ..
rm samtools-1.15.tar.bz2

echo "installation done: samtools-1.15"


echo "downloading and installing htslib-1.15 into $REFDIR"

wget https://github.com/samtools/htslib/releases/download/1.15/htslib-1.15.tar.bz2
tar -xf htslib-1.15.tar.bz2
cd htslib-1.15/
./configure --prefix=$REFDIR
make
make install
cd ..
rm htslib-1.15.tar.bz2

echo "installation done: htslib-1.15"


echo "downloading and installing bcftools-1.15 into $REFDIR"


wget https://github.com/samtools/bcftools/releases/download/1.15/bcftools-1.15.tar.bz2
tar -xf bcftools-1.15.tar.bz2
cd bcftools-1.15/
./configure --prefix=$REFDIR
make
make install
cd ..
rm bcftools-1.15.tar.bz2

echo "installation done: bcftools-1.15"

echo "installing dependencies (Archive::Zip/DBI/Module::Build/Try::Tiny) for Variant Effects Predictor (VEP)"


echo "installing cpanm package which will further help in installing dependencies for VEP"
sudo apt-get -y update
sudo apt-get install -y cpanminus
sudo apt-get install -y tabix

echo "installation done: cpanminus"

echo "creating directory $REFDIR/cpanm"
mkdir -p $REFDIR/cpanm

echo "exporting path to $PERL5LIB"

export PERL5LIB=$PERL5LIB:$REFDIR/cpanm/lib/perl5

echo "installing dependencies (Archive::Zip/DBI/Module::Build/Try::Tiny/Bio::DB::HTS) for Variant Effects Predictor (VEP)"
cpanm -l $REFDIR/cpanm Archive::Zip
cpanm -l $REFDIR/cpanm DBI
cpanm -l $REFDIR/cpanm Module::Build
cpanm -l $REFDIR/cpanm Try::Tiny
cpanm -l $REFDIR/cpanm Bio::DB::HTS


echo "installation of dependencies for VEP done"


echo "installing bioperl"
sudo apt-get -y update
sudo apt-get -y install bioperl


echo "installation done: bioperl"


echo "git cloning and installing VEP"

git clone https://github.com/Ensembl/ensembl-vep.git


cd ensembl-vep/

perl INSTALL.pl

cd ..

echo "installation done: VEP"


echo "downloading reference genome (GRCh38) from gencode database"

mkdir $REFDIR/fasta
cd fasta/
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_39/GRCh38.p13.genome.fa.gz
gunzip GRCh38.p13.genome.fa.gz; bgzip -c GRCh38.p13.genome.fa > GRCh38.p13.genome.bgzip.gz
cd ..
echo "done: GRCh38"


echo "exporting the $REFDIR/bin to the path"

export PATH=$REFDIR/bin:$PATH

echo "Done."



#echo "downloading vcf2maf from GitHub"


#wget https://github.com/mskcc/vcf2maf/archive/refs/tags/v1.6.21.tar.gz
#tar -xvzf v1.6.21.tar.gz
#rm v1.6.21.tar.gz

#echo "changing directory to vcf2maf"
#cd vcf2maf-1.6.21/

#commands to run annotator using vcf2maf package

#echo "running maf2vcf.pl to convert maf file to vcf file"
#perl maf2vcf.pl --input-maf $REFDIR/TCGA-CHOL-mutect2.maf --output-dir $REFDIR/ --ref-fasta $REFDIR/GRCh38.p13.genome.bgzip.gz

#echo "done. please find the output vcf file in $REFDIR"


#echo "running maf2maf.pl to convert maf file to maf via vcf file annotated with VEP package"
#perl maf2maf.pl --input-maf $REFDIR/TCGA-CHOL-mutect2.maf --output-maf $REFDIR/TCGA-CHOL-mutect2.vep.maf --ref-fasta $REFDIR/GRCh38.p13.genome.fa --vep-path $REFDIR/ensembl-vep --ncbi-build GRCh38

#echo "done. please find the output vcf file (*.vep.maf) in $REFDIR"


# command to run vep tool directly on vcf files

#./vep --species homo_sapiens --assembly GRCh38 --offline --no_progress --no_stats --buffer_size 5000 --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --canonical --protein --biotype --uniprot --tsl --pubmed --variant_class --shift_hgvs 1 --check_existing --total_length --allele_number --no_escape --xref_refseq --failed 1 --vcf --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length --dir ../../.vep --fasta ../fasta/GRCh38.p13.genome.fa --format vcf --input_file ../data/0059990c-d7df-47b4-a0e1-844746e97c1e.wxs.aliquot_ensemble_masked.vcf --output_file ../data/0059990c-d7df-47b4-a0e1-844746e97c1e.wxs.aliquot_ensemble_masked.vep.vcf