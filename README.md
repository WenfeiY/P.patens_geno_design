# P.patens_geno_design
------
Python and R script for Physcomitrium patens genome design
# Installation
---
1. Clone the repo
```
git clone https://github.com/WenfeiY/P.patens_geno_design.git
cd P.patens_geno_design
```
2. Create Anaconda enviroment
```
conda env create -f environment.yml
```
3. Install OligoArrayAux
```
wget http://www.unafold.org/download/oligoarrayaux-3.8.1.tar.gz
tar -zxvf oligoarrayaux-3.8.1.tar.gz
```
4. Add PATH to .bashrc
```
echo 'export PATH="<path_to_software>/oligoarrayaux-3.8.1/src:$PATH"' >> ~/.bashrc
source ~/.bashrc
```
5. Activate environment

```
conda activate geno
```
# Usage & Example
---
PCRmark to distinguish wild-type and synthetic sequences can be generated and verified by PCRmark_design.r:

The output file <Chromosome arm>_PCRmark_info.bed contains information of designed wtPCRmarks and relevent synPCRmarks for each gene, see Chr16L_PCRmark_info.bed in example_output directory.
```
Python script Pp.syn.generate.py can be used to generate designed chromosome sequences. You can get usage by:
$ python Pp.syn.generate.py -u
Usage: python Pp.syn.generate.py -g <Genbank_file_path> -n <Chr_number(INT)> -a <L/R> -p <PCRmark_file_path>
Options and arguments:
-g, --gb_file		chromosome genbank file
-n, --chr		chromosome number (1, 2...)
-a, --arm		specify chromosome arm to design (L or R)
-p, --primer_file	PCRmark information for the chromosome arm
-c, --cen		centromere file in bed format, default as Centromere.bed
--clean			not generate semi-designed sequence files
-u, --usage		print usage
```
Example:
```
python Pp.syn.generate.py -g wt/Chr16.gb -n 16 -a L --clean
```
The output file synMoss.Chr16L.gb in syn directory (same as example_output/synMoss.Chr16L.gb) contain sequence and annotations of designed Chr16L. Based on this, you can add or delete elements and split into small fragments for further synthesis and construction.



