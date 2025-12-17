# 00 concatenate the fastqs per read to make it easier the processing
zcat tinygex_S1_L001_R1_001.fastq.gz tinygex_S1_L002_R1_001.fastq.gz | gzip > tinygex_S2_L001_R1_001.fastq.gz
zcat tinygex_S1_L001_R2_001.fastq.gz tinygex_S1_L002_R2_001.fastq.gz | gzip > tinygex_S2_L001_R2_001.fastq.gz

# 01 subset the reads 
seqtk sample -s42 tinygex_S2_L001_R1_001.fastq.gz 0.5 | gzip > tinygex_S3_L001_R1_001.fastq.gz
seqtk sample -s42 tinygex_S2_L001_R2_001.fastq.gz 0.5 | gzip > tinygex_S3_L001_R2_001.fastq.gz

# 02 identify the reads in the subset S3
seqtk seq -A tinygex_S3_L001_R1_001.fastq.gz | grep '^>' | sed 's/^>//' > read_names_subset_R1.txt
seqtk seq -A tinygex_S3_L001_R2_001.fastq.gz | grep '^>' | sed 's/^>//' > read_names_subset_R2.txt

# 03 get the full reads
seqtk seq -A tinygex_S2_L001_R1_001.fastq.gz | grep '^>' | sed 's/^>//' > read_names_full_R1.txt
seqtk seq -A tinygex_S2_L001_R2_001.fastq.gz | grep '^>' | sed 's/^>//' > read_names_full_R2.txt

# 04 get the complementary reads
cat read_names_full_R1.txt | grep -v -F -f read_names_subset_R1.txt > read_names_complementary_R1.txt
cat read_names_full_R2.txt | grep -v -F -f read_names_subset_R2.txt > read_names_complementary_R2.txt

# 05 extract the complementary reads
zcat tinygex_S2_L001_R1_001.fastq | grep -A 3 -F -f read_names_complementary_R1.txt --no-group-separator | gzip > tinygex_S4_L001_R1_001.fastq.gz
zcat tinygex_S2_L001_R2_001.fastq | grep -A 3 -F -f read_names_complementary_R2.txt --no-group-separator | gzip > tinygex_S4_L001_R2_001.fastq.gz

