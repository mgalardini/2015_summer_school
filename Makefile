# Input files
# The genome file should be put in the same directory as the makefile
GENOME = reference.gbk
RFASTA = reference.fasta
# Start of the genome fragment
START = 48972142
# Scripts directory
SRCDIR = $(CURDIR)/src
# Gene to be mutated
GENE = ENSG00000125084
# Targeted mutations
MUTATIONS = C218G R235W
# Random mutations
RMUT = 10
# Sample name
SAMPLE = patient1

# Software
ARTDIR = $(SOFTDIR)art_bin_ChocolateCherryCake

MUTFASTA = mutated.fasta
$(MUTFASTA): $(GENOME)
	$(SRCDIR)/mutate_gbk $(GENOME) $(GENE) $(MUTFASTA) $(MUTATIONS)

READ1 = $(SAMPLE)1.fq
READ2 = $(SAMPLE)2.fq
$(READ1): $(MUTFASTA)
	$(ARTDIR)/art_illumina -p -i $(MUTFASTA) -l 100 -f 100 -m 200 -s 10 -o $(SAMPLE) -1 $(ARTDIR)/Illumina_profiles/Emp100R1.txt -2 $(ARTDIR)/Illumina_profiles/Emp100R2.txt

TREAD1 = $(SAMPLE)1.trim.fq
TREAD2 = $(SAMPLE)2.trim.fq
$(TREAD1): $(READ1)
	interleave_pairs $(READ1) $(READ2) | trim_quality -q 20 -w 5 --paired_reads | deinterleave_pairs -o $(TREAD1) $(TREAD2)

VCF = $(SAMPLE).vcf
$(VCF): $(TREAD1) $(RFASTA)
	bwa index $(RFASTA)
	bwa mem $(RFASTA) $(TREAD1) $(TREAD2) > aln.sam && \
		samtools view -bS aln.sam -q 20 -f 2 -F 256 -o aln.bam && \
		samtools sort aln.bam aln.sorted && \
		samtools index aln.sorted.bam && \
		freebayes -f $(RFASTA) --ploidy 1 --genotype-qualities --standard-filters aln.sorted.bam > $(VCF).tmp && \
		$(SRCDIR)/correct_vcf $(VCF).tmp $(START) $(VCF) && \
		rm $(VCF).tmp		

all: $(VCF)

.PHONY: all
