# P.patens_geno_design
Python and R scripts for genome design as a supplement to GenoDesigner (https://github.com/SynMoss/GenoDesigner) to mitigate GenoDesigner’s potential drawbacks like PCRmark design. These scripts can design PCRmarks for most genes and generate drafts of synthetic sequences of chromosome arms in a command-line manner. 
The example directory contains example data for GenoDesigner and scripts provided here.
# Installation
---
1. Clone the repo
```
git clone https://github.com/WenfeiY/P.patens_geno_design.git
cd P.patens_geno_design
```
2. Create Anaconda environment
```
conda env create -f environment.yml
```
3. Install OligoArrayAux and primer3
```
wget http://www.unafold.org/download/oligoarrayaux-3.8.1.tar.gz
tar -zxvf oligoarrayaux-3.8.1.tar.gz
sudo apt-get install -y build-essential g++ cmake git-all
git clone https://github.com/primer3-org/primer3.git primer3
cd primer3/src
make
make test
```
4. Add PATH to .bashrc
```
echo 'export PATH="<path_to_software>/oligoarrayaux-3.8.1/src:$PATH"' >> ~/.bashrc
echo 'export PATH="<path_to_software>/primer3-2.6.1/src:$PATH"' >> ~/.bashrc
source ~/.bashrc
```
5. Activate environment

```
conda activate geno
```
# Usage & Example
---
PCRmark to distinguish wild-type and synthetic sequences can be generated and verified by PCRmark_design.r:
$ Rscript PCRmarker_design.R -h
usage: PCRmarker_design.R -g <Genome_file> -f <gff_file> -c <Cen_file> -n <Chr_number(INT)> -a <L/R> -t <thread> -m <tm_location> -o <output_path>
Options:
-g GENOME_FILE, --genome_file=GENOME_FILE   genome file in fasta format
-f GFF_FILE, --gff_file=GFF_FILE    genome annotation file in gff format
-c CEN_FILE, --cen_file=CEN_FILE    centromere file in csv format
-n CHR, --chr=CHR   chromosome number (1, 2...)
-a ARM, --arm=ARM   specify chromosome arm to design (L or R)
-t THREAD, --thread=THREAD  thread count
-m TM_PATH, --tm_path=TM_PATH   the path of oligotm
-o OUTPUT_PATH, --output_path=OUTPUT_PATH   the output path for the generated PCRmark information
-h, --help  Show this help message and exit
```
The genomic information of the Physcomitrium patens can be obtained from our recent work (https://doi.org/10.6084/m9.figshare.22975925.v1).
```
Example:
```
Rscript PCRmarker_design.R -g Physcomitrium_patens_V4_genome.fasta -f Physcomitrium_patens_V4_rename.gff3 -c centromere.csv -n 16 -a L -t 20 -m ~/primer3/src/oligotm -o ~/synMoss_design/
```
The output file Pp.Chr<chromosome number><L/R>_PCRmark.txt contains information on designed wtPCRmarks and relevant synPCRmarks for each gene, see Pp.Chr16L_PCRmark.txt in the example directory.
```
Python script Pp.syn.generate.py can be used to generate designed chromosome sequences. You can get usage by:
$ python Pp.syn.generate.py -u
Usage: python Pp.syn.generate.py -g <Genbank_file_path> -n <Chr_number(INT)> -a <L/R> -p <PCRmark_file_path>
Options and arguments:
-g, --gb_file		chromosome genbank file
-n, --chr		chromosome number (1, 2...)
-a, --arm		specify chromosome arm to design (L or R)
-p, --primer_file	PCRmark information for the chromosome arm
--clean			not generate semi-designed sequence files
-u, --usage		print usage
```
Example:
```
python Pp.syn.generate.py -g Pp.Chr16L.gb -n 16 -a L --clean
```
The output file SynMoss.Chr16L.gb in the syn directory contains sequence and annotations of designed Chr16L. You can add or delete elements and split the synthetic sequence into small fragments for further synthesis and construction.
Note: SynMoss.Chr16L.gb in the example directory is the example output file of GenoDesigner, which is a little bit different from the output file. The difference results from the "Replace Stop Codons" step. The Python script here swaps all the stop codons to TAA and might cause amino acid changes in other CDSs. The GenoDesigner will not swap stop codons that overlap with other coding regions and provide a report (See Help-GenoDesigner at https://github.com/SynMoss/GenoDesigner).
