'''
#Download SRR fastqs into data/fastq/
bash scripts/01_download_fastq_ncbi.sh --srr_list config/srr_ids.txt --out_dir data/fastq --rename_casava

#Build manifest from that FASTQ directory
Rscript scripts/02_make_manifest.R --fastq_dir data/fastq --read_type paired --out manifests/manifest_paired.txt

#Install QIIME2 2026.1 (one-shot)
bash scripts/03_install_qiime2_2026.1.sh

#Run QIIME2 taxa barplot pipeline + export level-5/6/7.csv
bash scripts/04_qiime2_taxa_barplot_pipeline.sh \
  --manifest manifests/manifest_paired.txt \
  --metadata metadata.txt \
  --classifier classifier.qza \
  --out_dir qiime2_out

#Run mbX on level-7.csv
Rscript scripts/05_run_mbx.R --metadata metadata.txt --group_var BMIClass

#Diversity (YOU set sampling depth)
bash scripts/06_diversity.sh \
  --table qiime2_out/feature_table.qza \
  --rep_seqs qiime2_out/representative_sequences.qza \
  --metadata metadata.txt \
  --sampling_depth 10000 \
  --out_dir diversity_out

#ANCOMBC + MaAsLin2 (separate dirs)
Rscript scripts/07_ancombc.R --microbiome level-7.csv --metadata metadata.txt --group_var BMIClass --out_dir ancombc_out
Rscript scripts/08_maaslin2.R --microbiome level-7.csv --metadata metadata.txt --group_var BMIClass --out_dir maaslin2_out
'''
