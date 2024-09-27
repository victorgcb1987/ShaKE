
## Requirements

**smudgeplot v0.2.5** : We recommend to use conda (`conda install smudgeplot==0.2.5`) to install this module.  

kmc https://github.com/refresh-bio/KMC. You will need to add kmc binaries to PATH variable for example using `export PATH=$PATH:kmc/install/bin/`

clone KATULU respository: `git clone https://github.com/victorgcb1987/KATULU.git`

## How to use 

### Omics diversity calculator
Omics diversity calculator is a pipeline that will calculate Shannon diversity index for files of your choice. These file can be fasta/fastq files or gene expression files like the ones provided by stringtie.

In order to runt this program you will need to prepare a table similar to this one:

| onegroup | R | transcriptome | 0 | 9999999999 | /data/users/carpinterov/doctorado/test/test1_r1.fastq.gz,/data/users/carpinterov/doctorado/test/test1_r2.fastq.gz |
|---|---|---|---|---|---|
| onegroup | R | transcriptome | 0 | 9999999999 | /data/users/carpinterov/doctorado/test/test2_r1.fastq.gz,/data/users/carpinterov/doctorado/test/test2_r2.fastq.gz |
| onegroup | G | genome | 0 | 9999999999 | /data/users/carpinterov/doctorado/test/test5.fasta.gz |
| othergroup G transcriptome | 0 | 999999999 | /data/users/carpinterov/doctorado/test/test4_r1.fastq.gz,/data/users/carpinterov/doctorado/test/test4_r2.fastq.gz |
| othergroup | E | expression | 0 | 000000000 | /data/users/carpinterov/doctorado/test/Petal.guided.abund.tsv |
| othergroup |E | expression | 0 | 000000000 | /data/users/carpinterov/doctorado/test/adult_vascular_leaf.guided.abund.tsv |
| anothergroup | G  |genome | 0 | 9999999999 | /data/users/carpinterov/doctorado/test/test5.fasta.gz |
| anothergroup | G |genome | 0 | 9999999999 | /data/users/carpinterov/doctorado/test/test5.fasta.gz |
