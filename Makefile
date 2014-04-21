#all:	bam_rs.100
#all:	KO_H3K4me2_vs_input KO_H3K9_vs_input KO_H3K27me3_vs_input WT_H3K4me2_vs_input WT_H3K9_vs_input WT_H3K27me3_vs_input KO_vs_WT_H3K4me2 KO_vs_WT_H3K9 KO_vs_WT_H3K27me3
all:	MACS2.RS0.summary/macs2.rs0_broad_peaks_mtx.txt
#all:	MEME.RS0/macs2.rs0_merged

bam_rs.100:
	GDIR="/Volumes/hts_core/Shared/repeat_soaker/rmsk.mm9.bed" ; \
	if [ ! -d $@ ]; then \
		mkdir $@ ; \
	fi ; \
	for file in `find pepe_chip/bam -type f -name '*.bam'`; do \
		fname=`basename $$file` ; \
		fname=$${fname%????} ; \
		qsub -b y -l h_vmem=8G -V -j y -m e -cwd -q all.q -N job-$$fname-rs "repeat-soaker -r $$GDIR -p 1 -o $@/$$fname.bam $$file" ; \
	done

KO_H3K4me2_vs_input:
	DIR="bam_rs.0/" ; \
	T=$$DIR"LB_ngs-7a4a651415_AGTCAA.bam "$$DIR"LB_ngs-3c9897e4tk_CCGTCC.bam"; \
	C=$$DIR"LB_ngs-cd4xdr6dt1_GTGGCC.bam "$$DIR"LB_ngs-a4t53e0yy5_GTTTCG.bam" ; \
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$@ "macs2 callpeak -t $$T -c $$C -g mm -B --broad --name $@" 

KO_H3K9_vs_input:
	DIR="bam_rs.0/" ; \
	T=$$DIR"LB_ngs-27d40x839x_AGTTCC.bam"; \
	C=$$DIR"LB_ngs-cd4xdr6dt1_GTGGCC.bam "$$DIR"LB_ngs-a4t53e0yy5_GTTTCG.bam" ; \
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$@ "macs2 callpeak -t $$T -c $$C -g mm -B --broad --name $@"

KO_H3K27me3_vs_input:
	DIR="bam_rs.0/" ; \
	T=$$DIR"LB_ngs-d8xx4ca9c9_ACTTGA.bam "$$DIR"LB_ngs-61652449kc_ATGTCA.bam"; \
	C=$$DIR"LB_ngs-cd4xdr6dt1_GTGGCC.bam "$$DIR"LB_ngs-a4t53e0yy5_GTTTCG.bam" ; \
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$@ "macs2 callpeak -t $$T -c $$C -g mm -B --broad --name $@"

WT_H3K4me2_vs_input:
	DIR="bam_rs.0/" ; \
	T=$$DIR"LB_ngs-e08axk9ky4_TTAGGC.bam "$$DIR"LB_ngs-x498eyt45k_GCCAAT.bam"; \
	C=$$DIR"LB_ngs-erdak51t7y_GTCCGC.bam "$$DIR"LB_ngs-97tt5r5282_GTGAAA.bam" ; \
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$@ "macs2 callpeak -t $$T -c $$C -g mm -B --broad --name $@"

WT_H3K9_vs_input:
	DIR="bam_rs.0/" ; \
	T=$$DIR"LB_ngs-179022204r_ATCACG.bam "$$DIR"LB_ngs-1da684dc9r_TGACCA.bam"; \
	C=$$DIR"LB_ngs-erdak51t7y_GTCCGC.bam "$$DIR"LB_ngs-97tt5r5282_GTGAAA.bam" ; \
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$@ "macs2 callpeak -t $$T -c $$C -g mm -B --broad --name $@"

WT_H3K27me3_vs_input:
	DIR="bam_rs.0/" ; \
	T=$$DIR"LB_ngs-kddt6k043t_CGATGT.bam "$$DIR"LB_ngs-48r54dkakr_ACAGTG.bam"; \
	C=$$DIR"LB_ngs-erdak51t7y_GTCCGC.bam "$$DIR"LB_ngs-97tt5r5282_GTGAAA.bam" ; \
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$@ "macs2 callpeak -t $$T -c $$C -g mm -B --broad --name $@"

KO_vs_WT_H3K4me2:
	DIR="bam_rs.0/" ; \
	T=$$DIR"LB_ngs-7a4a651415_AGTCAA.bam "$$DIR"LB_ngs-3c9897e4tk_CCGTCC.bam"; \
	C=$$DIR"LB_ngs-e08axk9ky4_TTAGGC.bam "$$DIR"LB_ngs-x498eyt45k_GCCAAT.bam" ; \
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$@ "macs2 callpeak -t $$T -c $$C -g mm -B --broad --name $@"

KO_vs_WT_H3K9:
	DIR="bam_rs.0/" ; \
	T=$$DIR"LB_ngs-27d40x839x_AGTTCC.bam"; \
	C=$$DIR"LB_ngs-179022204r_ATCACG.bam "$$DIR"LB_ngs-1da684dc9r_TGACCA.bam" ; \
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$@ "macs2 callpeak -t $$T -c $$C -g mm -B --broad --name $@"

KO_vs_WT_H3K27me3:
	DIR="bam_rs.0/" ; \
	T=$$DIR"LB_ngs-d8xx4ca9c9_ACTTGA.bam "$$DIR"LB_ngs-61652449kc_ATGTCA.bam"; \
	C=$$DIR"LB_ngs-kddt6k043t_CGATGT.bam "$$DIR"LB_ngs-48r54dkakr_ACAGTG.bam" ; \
	qsub -b y -l h_vmem=10G -V -j y -m e -M dozmorovm@omrf.org -cwd -q all.q -N macs2_$@ "macs2 callpeak -t $$T -c $$C -g mm -B --broad --name $@"

MACS2.RS0.summary/macs2.rs0_merged.bed:
	for file in `find MACS2.RS0 -type f -name '*_broad_peaks.bed'`; do \
		cat $$file >> macs2_all.bed ; \
	done && \
	mergeBed -i macs2_all.bed > $@ && \
	rm macs2_all.bed

MACS2.RS0.summary/macs2.rs0_broad_peaks_mtx.txt:	MACS2.RS0.summary/macs2.rs0_merged.bed
	k="" ; \
	for file in `find bam_rs.0 -type f -name '*.bam' | sort`; do \
		k=$$k" "$$file ; \
		fname=`basename $$file` ; \
		fname=$${fname%????} ; \
		echo $$fname >> MACS2.RS0.summary/macs2.rs0_broad_peaks_head.txt ; \
	done && \
	bedtools multicov -bams $$k -bed $<  >> $@

MEME.RS0/macs2.rs0_merged:	MACS2.RS0.summary/macs2.rs0_merged.bed
	if [ ! -d $@ ]; then \
		mkdir $@ ; \
	fi ; \
	bedtools getfasta -fi mm9.fasta -bed $< -fo $@.fa && \
	qsub -V -cwd -b y -j y -l h_vmem=30G -m e -pe openmpi 8 meme-chip -noecho -oc $@ -time 300 -db /Volumes/hts_core/Shared/meme_db/JASPAR_CORE_2014_vertebrates.meme -db /Volumes/hts_core/Shared/meme_db/uniprobe_mouse.meme -meme-mod zoops -meme-minw 6 -meme-maxw 30 -meme-nmotifs 3 -dreme-e 0.05 -centrimo-score 5 -centrimo-ethresh 10 $@.fa 