import os

### funziona con python2.7. Richiede bowtie2
### se le read forward e reverese sono indicate con _1 e _2, usa questi invece di _1 e _2


file_name = 'trimmed_a_1_'#### nome del file fastq senza estensione. Per esempio, se il tuo file fastq si chiama WT_A_1.fastq.gz, il file name è WT_A_
contaminants = "/mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/bowtie2_analysis/bowtie_index/mouse_rrnas" ### path all'index di bowtie2 (da creare a partire dal file dei contaminanti)
paired_1 = '/mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/cutadapt_out/trimmed_a_1_1.fq.gz' ### path alle read 1 after trimming. Contiene il file_name
paired_2 = '/mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/cutadapt_out/trimmed_a_1_2.fq.gz' ### path alle read 2 after trimming. Contiene il file_name
contaminants_folder = '/mnt/large/jfiorentino/Cerase_Data/bulkRNA_seq/raw_data/a_1/cutadapt_out/bowtie2/'### folder di output

#if os.path.isdir(contaminants_folder)==False:#
#	os.mkdir(contaminants_folder)

### setta il parametro -p di bowtie sul numero di core che ti servono
### bowtie allinea le read su contaminant, scrive le read che allineano su contaminants_folder+file_name+"_1 e quelle che non allineano su contaminants_folder+file_name+"_1_unmapped"
#os.system("/home/jonny/miniconda3/bin/bowtie2 -p 20 -x "+contaminants+" -U "+paired_1+" -S "+contaminants_folder+file_name+"_1 --un "+contaminants_folder+file_name+"_1_unmapped")
#os.system("/home/jonny/miniconda3/bin/bowtie2 -p 20 -x "+contaminants+" -U "+paired_2+" -S "+contaminants_folder+file_name+"_2 --un "+contaminants_folder+file_name+"_2_unmapped")
#print "sorting files..."

#ordino le read non mappate
#os.system( "cat "+contaminants_folder+file_name+"_1_unmapped | paste - - - - | sort -k1,1 -t ' ' | tr '\t' '\n' > "+contaminants_folder+file_name+"_1_unmapped_sorted")
#os.system( "cat "+contaminants_folder+file_name+"_2_unmapped | paste - - - - | sort -k1,1 -t ' ' | tr '\t' '\n' > "+contaminants_folder+file_name+"_2_unmapped_sorted")

#rimuovo le read non mappate non ordinate
#os.system("rm "+contaminants_folder+file_name+"_1_unmapped")
#os.system("rm "+contaminants_folder+file_name+"_2_unmapped")

#rimuovo le read mappate
#os.system("rm "+contaminants_folder+file_name+"_1")
#os.system("rm "+contaminants_folder+file_name+"_2")

### fastqcombinepairedend.py prende in input i file fastq di read non mappate sortate e rimette in fase le read.
### @SEQ sono le prime 4 lettere di ogni riga del fastq che comincia con @

os.system("python fastqcombinepairedend.py '@A00' ' ' "+contaminants_folder+file_name+"_1_unmapped_sorted "+contaminants_folder+file_name+"_2_unmapped_sorted")

#le read decontaminate conterranno stillpaired nel nome. Rimuovi _unmapped_sorted dopo aver controllato che tutti sia OK

