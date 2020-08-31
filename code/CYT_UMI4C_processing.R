library(umi4cPackage)
p4cLoadConfFiles(conf_dir="/home/labs/lplab/mramos/Projects/CYT_hg19/data/UMI4C/merged/conf/")

p4cCreate4CseqTrack(sample_ids=1:36, overwrite.if.exist=TRUE)

