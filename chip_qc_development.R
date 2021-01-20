#!/usr/bin/env Rscript

##Based on Joe Boyd's FRIP QC script
##with modifications by Mike Mariani 
##for general chip QC. This is a work 
##in progress. 

##Frietze Lab UVM 2021

library(seqsetvis)
library(data.table)
library(GenomicRanges)
library(ggplot2)

##config_dt = data.table(
##    mark_keys = c("ATAD", 
##				  "[0-9]ac_", 
##				  "ATAD2_", 
##				  "ATAD2B_", 
##				  "BRD4_", 
##				  "H3K4me3", 
##				  "H4K12ac", 
##				  "H4K16ac", 
##				  "H4K20ac", 
##				  "H4K5ac", 
##				  "H4K8ac"),
##    file_tags = c("all_ATAD", 
##				  "all_acetyl", 
##				  "ATAD2", 
##				  "ATAD2B", 
##				  "BRD4", 
##				  "H3K4me3", 
##				  "H4K12ac", 
##				  "H4K16ac", 
##				  "H4K20ac", 
##				  "H4K5ac", 
##				  "H4K8ac"),
##    min_fractions = c(.5),
##    min_number = 2,
##    view_sizes = c(2e3),
##    modes = c("seacr", 
##			  "broad", 
##			  rep("seacr", 3), 
##			  rep("broad", 6))
##)

##mode_vars = list(
##    broad = list(
##        peak_dir_FUN = function(search_dirs){
##            file.path(search_dirs, "macs2.broad.all.frag.aug18")
##        },
##        peak_file_key = ".+peaks.broadPeak$",
##        peak_load_FUN = easyLoad_broadPeak
##    ),
##    narrow = list(
##        peak_dir_FUN = function(search_dirs){
##            file.path(search_dirs, "macs2.narrow.all.frag.aug18")
##        },
##        peak_file_key = ".+peaks.narrowPeak$",
##        peak_load_FUN = easyLoad_narrowPeak
##    ), 
##    seacr = list(
##        peak_dir_FUN = function(search_dirs){
##            file.path(search_dirs, "seacr.aug12.all.frag")
##        },
##        peak_file_key = ".stringent.bed$",
##        peak_load_FUN = easyLoad_seacr
##    )
##)

make_dt = function(files, group_depth = ifelse(grepl("bam$", files[1]), 3, 2)){
  p_dt = data.table(file = files)
  p_dt[, name := sub("\\..+", "", basename(file))]
  p_dt[, name := sub("(?<=rep[0-9])_.+", "", name, perl = TRUE)]
  if(group_depth == 1){
    p_dt[, group := basename(dirname(file))]   
  }else if(group_depth == 2){
    p_dt[, group := basename(dirname(dirname(file)))]    
  }else if(group_depth == 3){
    p_dt[, group := basename(dirname(dirname(dirname(file))))]
  }
  p_dt[, batch := LETTERS[as.numeric(factor(group, levels = group_lev))]]
  p_dt[, name := paste0(name, ".", batch)]
  dupe_num = 1
  if(any(duplicated(p_dt$name))){
    p_dt[, dupe_num := seq(.N)-1, .(name)]   
    p_dt[dupe_num > 0, name := paste0(name, ".", dupe_num) ]
    p_dt$dupe_num = NULL
  }
  # p_dt[, condition := ifelse(grepl("DF", name), "DF", "WT")]
  # p_dt[, xlink := ifelse(grepl("[T4]x_", name), "xlinked", "not")]
  p_dt[]
}

peaks.files <- list.files(path="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_macs2_output/hg38/narrow",
                          patter=".narrowPeak",
                          full.names = TRUE)
bws.files   <- list.files(path='/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws',
                          pattern=".bw$",
                          full.names = TRUE)
bams.files  <- list.files(path="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38",
                          pattern=".bam$",
                          full.names = TRUE)
fastq.files <- list.files(path="/slipstream/home/mmariani/projects/vzv_interactions/vzv_ctcf_cohrs_chip_data/trimmomatic_output",
                          pattern=".fastq.gz",
                          full.names=TRUE)

length(peaks.files)
length(bws.files)
length(bams.files)
length(fastq.files)

base_dir = "/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_macs2_output/hg38/narrow"
group_lev = factor(basename(base_dir))

##qc_dir = "~/R/qc_cutnrun/qc_output_KQ.high_consensus"
##search_dirs = c("~/cutruntools/workdir_102920_KQ_MCF7.hg38")

## ".neg" indicates input file(s)

##bam_dirs    <- file.path(search_dirs, "aligned.aug10/sorted")
##bw_dirs     <- file.path(search_dirs, "macs2.narrow.all.frag.aug18")
##bw_dirs.neg <- file.path(search_dirs, "macs2.narrow.all.frag.aug18")

##group_lev = factor(basename(search_dirs))
##add fastq reads
##add total peaks

##stopifnot(config_dt$modes %in% names(mode_vars))
##negative_key = "IgG_"
##negative_tag = "IgG"

##i = 1
##for(i in seq(7, nrow(config_dt))){
##    mark_key     <- config_dt$mark_keys[i]
##    file_tag     <- config_dt$file_tags[i]
##    min_fraction <- config_dt$min_fractions[i]
##    min_number   <- config_dt$min_number[i]
##    view_size    <- config_dt$view_sizes[i]
##    mode_var     <- mode_vars[[config_dt$modes[i]]]
##    peak_dirs    <- mode_var$peak_dir_FUN(search_dirs)
##    peak_key     <- mode_var$peak_file_key
##    message(file_tag)
##    o_dir = file.path(qc_dir, paste0("qc_", file_tag))
##    dir.create(o_dir, showWarnings = FALSE, recursive = TRUE)
##    setwd(o_dir)
##    ## peak_dirs = file.path(search_dirs, "macs2.narrow.all.frag.aug18")
    
message("find files...")
    
    ## peak_files = dir(peak_dirs, pattern = "[kK]27[aA][cC].+peaks.narrowPeak", full.names = TRUE)
    ##peak_files  = dir(peak_dirs, pattern = paste0(mark_key, ".+", peak_key), full.names = TRUE)
    ##bam_files   = dir(bam_dirs, pattern = paste0(mark_key, ".+bam$"), full.names = TRUE)
    ##bw_files    = dir(bw_dirs, pattern = paste0(mark_key, ".+bw$"), full.names = TRUE)
    ##fastq_files = c(
    ##    dir(search_dirs, paste0(mark_key, '.+R1_001.fastq.gz$'), full.names = TRUE),
    ##    dir(file.path(search_dirs, 'fastq_done'), paste0(mark_key, '.+R1_001.fastq.gz$'), full.names = TRUE)
    ##)

peak_files      <- peaks.files[c(2,5,6)]
bw_files        <- bws.files[c(1,3,5,6,9,10)]  
bam_files       <- bams.files[c(1,3,5,6,9,10)]
fastq_files     <- fastq.files[c(1,3,5,6,9,10)]

peak_files.neg  <- peaks.files[c(1,3,4)]
bam_files.neg   <- bams.files[c(2,4,7,8,11,12)]  
bw_files.neg    <- bws.files[c(2,4,7,8,11,12)]
fastq_files.neg <- fastq.files[c(2,4,7,8,11,12)]

##peak_files  <- bam_dirs   
##bam_files   <- bw_dirs 
##bw_files    <- bw_dirs
##fastq_files <- 
      
    ##peak_files.neg = dir(peak_dirs, pattern = paste0(negative_key,  ".+", peak_key), full.names = TRUE)
    ##bam_files.neg = dir(bam_dirs, pattern = paste0(negative_key, ".+bam$"), full.names = TRUE)
    ##bw_files.neg = dir(bw_dirs, pattern = paste0(negative_key, ".+bw$"), full.names = TRUE)
    ##fastq_files.neg = c(
    ##    dir(search_dirs, paste0(negative_key, '.+R1_001.fastq.gz$'), full.names = TRUE),
    ##    dir(file.path(search_dirs, 'fastq_done'), paste0(negative_key, '.+R1_001.fastq.gz$'), full.names = TRUE)
    ##)
    
    ##list the fastq files:
    system("ls -d /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/* | grep '.bam$'")
    bam.files.check <- system("ls -d /slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/* | grep '.bam$'", intern=TRUE)
    
    message("count fastq reads...")
    fq_dt = rbindlist(
                    pbmcapply::pbmclapply(
                      bam.files.check,
                      ##c(fastq_files, fastq_files.neg), 
                        mc.cores = 10, 
                        function(f){
                          cnt_f = paste0(f, ".cnt")
                          ## if(TRUE){
                          if(!file.exists(cnt_f)){
                            cnt = system(paste0("gunzip -c ", f, "| wc -l"), intern = TRUE)
                            cnt_dt = data.table(f, as.numeric(cnt)/4)
                            fwrite(cnt_dt, cnt_f, sep = "\t", col.names = FALSE)
                            cnt_dt
                          }else{
                            fread(cnt_f, sep = "\t")
                          }
                        }
                      )
                    )
    
    setnames(fq_dt, c("fastq", "count"))
    fq_dt.tmp = make_dt(fq_dt$fastq, group_depth = 1)
    fq_dt.tmp$fastq_count = fq_dt$count
    fq_dt = fq_dt.tmp
    
    fq_dt$file
    fq_dt$name <- c("ctcf_mock_rep_1",
                    "input_mock_rep_1",
                    "ctcf_infected_rep_1",
                    "input_infected_rep_1",
                    "ctcf_mock_rep_2",
                    "ctcf_mock_rep_3",
                    "input_mock_rep_2",
                    "input_mock_rep_3",
                    "ctcf_infected_rep_2",
                    "ctcf_infected_rep_3",
                    "input_infected_rep_2",
                    "input_infected_rep_3")
    
    fq_dt$batch <- 1
    
    file_tag <- "infected"
    negative_tag <- "mock"
    
    p_dt.pos = make_dt(peak_files)
    p_dt.pos$treatment = file_tag
    p_dt.pos$batch <- "1"
    p_dt.pos$name <- c("ctcf_infected_rep_1",
                       "ctcf_infected_rep_2",
                       "ctcf_infected_rep_3")
    
    qdt.bam = make_dt(bam_files)
    qdt.bam$treatment = file_tag
    qdt.bam$name <- c("ctcf_mock_rep_1",
                      "ctcf_infected_rep_1",
                      "ctcf_mock_rep_2",
                      "ctcf_mock_rep_3",
                      "ctcf_infected_rep_2",
                      "ctcf_infected_rep_3")
    qdt.bam$batch <- 1
    qdt.bam$treatment <- c("mock",
                           "infected",
                           "mock",
                           "mock",
                           "infected",
                           "infected")
    
    qdt.bw = make_dt(bw_files)
    ##qdt.bw$treatment = file_tag
    qdt.bw$treatment = c("mock",
                         "infected",
                         "mock",
                         "mock",
                         "infected",
                         "infected")
    
    qdt.bw$name <- c("ctcf_mock_rep_1",
                     "ctcf_infected_rep_1",
                     "ctcf_mock_rep_2",
                     "ctcf_mock_rep_3",
                     "ctcf_infected_rep_2",
                     "ctcf_infected_rep_3")
    
    p_dt.neg = make_dt(peak_files.neg)
    p_dt.neg$treatment = negative_tag
    p_dt.neg$batch <- 1
    p_dt.neg$name <- c("ctcf_mock_rep_1",
                       "ctcf_mock_rep_2",
                       "ctcf_mock_rep_3")
    
    qdt.bw.neg = make_dt(bw_files.neg)
    qdt.bw.neg$treatment = negative_tag
    qdt.bw.neg$name <- c("input_mock_rep_1",
                         "input_infected_rep_1",
                         "input_mock_rep_2",
                         "input_mock_rep_3",
                         "input_infected_rep_2",
                         "input_infected_rep_3")
    
    qdt.bam.neg = make_dt(bam_files.neg)
    qdt.bam.neg$treatment = negative_tag
    qdt.bam.neg$name <- c("input_mock_rep_1",
                          "input_infected_rep_1",
                          "input_mock_rep_2",
                          "input_mock_rep_3",
                          "input_infected_rep_2",
                          "input_infected_rep_3")
    qdt.bam.neg$batch <- 1
    qdt.bam.neg$treatment <- c("mock",
                           "infected",
                           "mock",
                           "mock",
                           "infected",
                           "infected")
    
    stopifnot(!any(duplicated(p_dt.pos$name)))
    stopifnot(!any(duplicated(p_dt.neg$name)))
    
    message("load peaks...")
    
    ##peak_load_FUN = easyLoad_narrowPeak
    peak_grs = seqsetvis::easyLoad_narrowPeak(peak_files,
                                              file_names = c("ctcf_infected_rep_1",
                                                             "ctcf_infected_rep_2",
                                                             "ctcf_infected_rep_3"))
    ##peak_grs = mode_var$peak_load_FUN(p_dt.pos$file)
    names(peak_grs) = factor(p_dt.pos$name, levels=p_dt.pos$name)
    lengths(peak_grs)
    ##S4Vectors::mcols(peak_grs)
    
    
    peak_dt.pos = ssvFeatureBars(peak_grs, return_data = TRUE)
    setnames(peak_dt.pos, c("count", "group"), c("peak_count", "name"))
    peak_dt.pos$name <- c("ctcf_infected_rep_1",
                          "ctcf_infected_rep_2",
                          "ctcf_infected_rep_3")
    peak_dt.pos$name <- factor(peak_dt.pos$name, levels=peak_dt.pos$name)
    
    ##peak_grs.neg = mode_var$peak_load_FUN(p_dt.neg$file)
    peak_grs.neg = seqsetvis::easyLoad_narrowPeak(peak_files.neg,
                                                  file_names = c("ctcf_mock_rep_1",
                                                                 "ctcf_mock_rep_2",
                                                                 "ctcf_mock_rep_3"))
    names(peak_grs.neg) = factor(p_dt.neg$name, levels = p_dt.neg$name)
    lengths(peak_grs.neg)
    
    peak_dt.neg = ssvFeatureBars(peak_grs.neg, return_data = TRUE)
    setnames(peak_dt.neg, c("count", "group"), c("peak_count", "name"))
    peak_dt.neg$name <- c("ctcf_mock_rep_1",
                          "ctcf_mock_rep_2",
                          "ctcf_mock_rep_3")
    peak_dt.neg$name <- factor(peak_dt.neg$name, levels=peak_dt.neg$name)
    
    peak_o = rev(order(lengths(peak_grs)))
    
    message("overlap peaks...")
    olaps_all = ssvOverlapIntervalSets(peak_grs[peak_o])
    p_peak_counts_all = ssvFeatureBars(peak_grs[peak_o], 
                                       bar_colors = "black", 
                                       show_counts = FALSE, 
                                       counts_text_colors = "gray60") +
                                       guides(fill = "none") +
                                       scale_y_continuous(labels = function(x)x/1e3) +
                                       labs(y = "peak count (k)", title = "all peaks", subtitle = paste("counts:", formatC(max(lengths(peak_grs)), big.mark = ",", format = "d"), "max"), x = "")
    
    p_peak_overlaps_all = ssvFeatureBinaryHeatmap(olaps_all, 
                                                  raster_approximation = TRUE) + 
                                       labs(title = "all peaks", 
                                            subtitle = paste("overlaps:", 
                                                             formatC(length(olaps_all), 
                                                                     big.mark = ",", 
                                                                     format = "d"), 
                                                             "regions"))
    
    ##Joe's recommendations for parameters:
    ##Joe says 0.5 is the same as min_number=3
    min_fraction = c(0.5)
    min_number = 3
    view_size = c(2e3)
    
    message("consensus peaks...")
    
    olaps = ssvConsensusIntervalSets(peak_grs[peak_o], 
                                     min_fraction = min_fraction, 
                                     min_number = min(min_number, 
                                                      length(peak_o)))
    
    p_peak_counts_consenus = ssvFeatureBars(olaps, 
                                            bar_colors = "black", 
                                            show_counts = FALSE, 
                                            counts_text_colors = "gray60" ) +
                                            guides(fill = "none") +
                                            scale_y_continuous(labels = function(x)x/1e3) +
                                            labs(y = "peak count (k)", 
                                                 title = "consensus peaks only", 
                                                 subtitle = paste("counts:", 
                                                                  formatC(max(colSums(as.data.frame(mcols(olaps)))), 
                                                                          big.mark = ",", 
                                                                          format = "d"), "max"))
    
    p_peak_overlaps_consensus = ssvFeatureBinaryHeatmap(olaps, 
                                                        raster_approximation = TRUE) +
                                                        labs(title = "consensus peaks only", 
                                                             subtitle = paste("overlaps:", 
                                                                              formatC(length(olaps), 
                                                                                      big.mark = ",", 
                                                                                      format = "d"), 
                                                                              "regions")
                                                             )
    
    pg_peaks = cowplot::plot_grid(p_peak_counts_all + 
                                  theme(axis.text.x = element_blank()), 
                                  p_peak_overlaps_all + 
                                  theme(axis.text.x = element_blank()), 
                                  p_peak_counts_consenus, 
                                  p_peak_overlaps_consensus, 
                                  rel_heights = c(2, 3))
    
    output.dir <- "/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip"
    
    ggsave(paste0(output.dir,
                  "/peak_numbers.", 
                  file_tag, 
                  ".pdf"), 
           pg_peaks, 
           width = 8, 
           height = 9)
    
    olaps = ssvRecipes::sampleCap(olaps, 5e3)
    #filter for canonical and no chrM
    olaps = subset(olaps, grepl("chr[0-9XY]+$", seqnames ))
    options(mc.cores = 20)
    
    qdt.bw$batch <- 1
    qdt.bw[, facet := gsub("_", "\n", name)]
    message("center query regions...")
    query_gr.uncentered = resize(olaps, 5e3, fix = "center")
    
    bw_dt.uncentered = ssvFetchBigwig(qdt.bw, 
                                      query_gr.uncentered, 
                                      win_size = 100, 
                                      win_method = "summary",
                                      return_data.table = TRUE, 
                                      n_region_splits = 20)
    
    query_gr = unique(resize(centerGRangesAtMax(bw_dt.uncentered, 
                                                query_gr.uncentered), 
                                                view_size, 
                                                fix = 'center')
                      )
    
    message("fetch query regions...")
    bw_dt = ssvFetchBigwig(bw_files,
                           ##qdt.bw, 
                           query_gr, 
                           win_size = 100, 
                           win_method = "summary",
                           return_data.table = TRUE, 
                           n_region_splits = 20)
    bw_dt = append_ynorm(bw_dt)
    ##unique(bw_dt$sample)
    
    bw_dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_Control_HFL_AntiCTCF_S31_L005_R1_001.trimmed.R1.sorted_norm_unstranded.bw",name:="ctcf_mock_rep_1",]
    bw_dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_VZV_AntiCTCF_S33_L005_R1_001.trimmed.R1.sorted_norm_unstranded.bw",name:="ctcf_infected_rep_1",]        
    bw_dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/HFL_con_CTCF_set2_S2_L001_R1_001.trimmed.R1.sorted_norm_unstranded.bw",name:="ctcf_mock_rep_2",]                                                      
    bw_dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/HFL_con_CTCF_set3_S6_L001_R1_001.trimmed.R1.sorted_norm_unstranded.bw",name:="ctcf_mock_rep_3",]                                                      
    bw_dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/HFL_VZV_CTCF_set2_S4_L001_R1_001.trimmed.R1.sorted_norm_unstranded.bw",name:="ctcf_infected_rep_2",]                                                      
    bw_dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/HFL_VZV_CTCF_set3_S8_L001_R1_001.trimmed.R1.sorted_norm_unstranded.bw",name:="ctcf_infected_rep_3",]                                                      

    ## set.seed(0)
    ## clust_dt = ssvSignalClustering(bw_dt, fill_ = "y_norm", max_cols = Inf, facet_ = "facet")
    ##negative
    
    qdt.bw.neg$batch = 1
    qdt.bw.neg[, facet := gsub("_", "\n", name)]
    
    bw_dt.neg = ssvFetchBigwig(bw_files.neg,
                               ##qdt.bw.neg, 
                               query_gr, 
                               win_size = 100, 
                               win_method = "summary",
                               return_data.table = TRUE, 
                               n_region_splits = 20)
    bw_dt.neg = append_ynorm(bw_dt.neg)
    ##unique(bw_dt.neg$sample)
    bw_dt.neg[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_Control_HFL_Input_S30_L005_R1_001.trimmed.R1.sorted_norm_unstranded.bw", name:="input_mock_rep_1",]
    bw_dt.neg[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_VZV_Input_S32_L005_R1_001.trimmed.R1.sorted_norm_unstranded.bw", name:="input_infected_rep_1",]        
    bw_dt.neg[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/HFL_con_in_set2_S1_L001_R1_001.trimmed.R1.sorted_norm_unstranded.bw", name:="input_mock_rep_2",]                                                     
    bw_dt.neg[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/HFL_con_in_set3_S5_L001_R1_001.trimmed.R1.sorted_norm_unstranded.bw", name:="input_mock_rep_3",]                                                     
    bw_dt.neg[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/HFL_VZV_in_set2_S3_L001_R1_001.trimmed.R1.sorted_norm_unstranded.bw", name:="input_infected_rep_2",]                                                     
    bw_dt.neg[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_bowtie2_alignment/hg38/joe_boyd_norm_unstranded_bws/HFL_VZV_in_set3_S7_L001_R1_001.trimmed.R1.sorted_norm_unstranded.bw", name:="input_infected_rep_3",]                                                     
    
    bw_dt.all = rbind(bw_dt.neg, bw_dt)
    
    bw_dt.all[, y_relative := y / max(y), .(id)]
    #count reads
    # Rsamtools::countBam(bam_files[1], which = olaps[1:10])
    qdt.bam.all = rbind(qdt.bam, qdt.bam.neg)
    qdt.bam.all$batch <- 1
    
    ##qdt.bam.all$file
    
    ##qdt.bam.all$name <- c("ctcf_mock_rep_1",
    ##                      "ctcf_infected_rep_1",
    ##                      "ctcf_mock_rep_2",
    ##                      "ctcf_mock_rep_3",
    ##                      "ctcf_infected_rep_2",
    ##                      "ctcf_infected_rep_3",
    ##                      "input_mock_rep_1",
    ##                      "input_infected_rep_1",
    ##                      "input_mock_rep_2",
    ##                      "input_mock_rep_3",
    ##                      "input_infected_rep_2",
    ##                      "input_infected_rep_3")
    
    
    ##Check fq_dt before merging
    intersect(fq_dt$name, qdt.bam.all$name)
    
    fq_dt.check = merge(fq_dt, qdt.bam.all[, .(name, treatment)], by = "name")
    
    fq_dt = merge(fq_dt, qdt.bam.all[, .(name, treatment)], by = "name")
    
    ##fq_dt$file

    peak_dt.pos <- rbind(peak_dt.pos, peak_dt.neg)
    
    peak_dt.neg <- peak_dt.pos
    
    peak_dt.pos$name <- c("ctcf_infected_rep_1",
                          "ctcf_infected_rep_2",
                          "ctcf_infected_rep_3",
                          "ctcf_mock_rep_1",
                          "ctcf_mock_rep_2",
                          "ctcf_mock_rep_3")

    peak_dt.neg$name <- c("input_infected_rep_1",
                          "input_infected_rep_2",
                          "input_infected_rep_3",
                          "input_mock_rep_1",
                          "input_mock_rep_2",
                          "input_mock_rep_3")
    
    peak_dt = merge(rbind(peak_dt.pos,peak_dt.neg), 
                    qdt.bam.all[, .(name, treatment)], 
                    by = "name")
    
    message("fetch read counts...")
    ##dt = ssvFetchBamPE(qdt.bam.all, 
    ##                   query_gr, 
    ##                   return_unprocessed = TRUE, 
    ##                   n_region_splits = 50)
    
    dt.in = ssvFetchBam(c(bam_files,bam_files.neg),
                     query_gr, 
                     return_unprocessed = TRUE, 
                     n_region_splits = 50)
    
    dt <- dt.in
    ##I need to manually add name and treatment fields to dt
    basename(as.character(unique(dt.in$sample)))
    
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_Control_HFL_AntiCTCF_S31_L005_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "ctcf_mock_rep_1",]
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_VZV_AntiCTCF_S33_L005_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "ctcf_infected_rep_1",]      
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_con_CTCF_set2_S2_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "ctcf_mock_rep_2",]                                                      
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_con_CTCF_set3_S6_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "ctcf_mock_rep_3",]                                                      
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_VZV_CTCF_set2_S4_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "ctcf_infected_rep_2",]                                                      
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_VZV_CTCF_set3_S8_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "ctcf_infected_rep_3",]                                                      
    
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_Control_HFL_AntiCTCF_S31_L005_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "mock",]
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_VZV_AntiCTCF_S33_L005_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "infected",]      
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_con_CTCF_set2_S2_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "mock",]                                                      
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_con_CTCF_set3_S6_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "mock",]                                                      
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_VZV_CTCF_set2_S4_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "infected",]                                                      
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_VZV_CTCF_set3_S8_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "infected",] 
    
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_Control_HFL_Input_S30_L005_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "input_mock_rep_1",]  
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_VZV_Input_S32_L005_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "input_infected_rep_1",]           
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_con_in_set2_S1_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "input_mock_rep_2",]                                                     
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_con_in_set3_S5_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "input_mock_rep_3",]                                                        
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_VZV_in_set2_S3_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "input_infected_rep_2",]                                                        
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_VZV_in_set3_S7_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", name := "input_infected_rep_3",] 
    
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_Control_HFL_Input_S30_L005_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "mock",]  
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/Frietze_Cohrs_11497_180409_SNL128_0174_ACC3JAACXX_VZV_Input_S32_L005_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "infected",]           
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_con_in_set2_S1_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "mock",]                                                     
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_con_in_set3_S5_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "mock",]                                                        
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_VZV_in_set2_S3_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "infected",]                                                        
    dt[sample=="/slipstream/home/mmariani/projects/vzv_interactions/vzv_cohrs_ctcf_chip/vzv_cohrs_ctcf_chip_marked_duplicates/hg38/HFL_VZV_in_set3_S7_L001_R1_001.trimmed.R1.sorted.marked_duplicates.bam", treatment := "infected",] 
    
    
    ##dt.box = dt[, .(N = length(unique(qname))), .(id, name, treatment, sample)]
    dt.box = dt[, .(N = length(unique(read_id))), .(id, name, treatment, sample)]
    
    unique(dt.box[, .(name, treatment, sample)])
    
    unique(dt.box$name)
    
    dt.box.casted <- dcast(dt.box, 
                           id~name, 
                           value.var = "N", 
                           fill = 0)
    
    dt.box_filled = melt(dt.box.casted,
                         id.vars = "id", 
                         value.name = "N", 
                         variable.name = "name")
    
    dt.box = merge(dt.box_filled, unique(dt.box[, .(name, treatment, sample)]), by = "name")
    
    unique(dt.box$sample)
    unique(dt.box$treatment)
    unique(dt.box$name)
    
    ##mapped_counts = sapply(rbind(qdt.bam, qdt.bam.neg)$file, function(f){
    mapped_counts = sapply(rbind(bam_files, bam_files.neg), function(f){
        stats = Rsamtools::idxstatsBam(f)
        stats = subset(stats, grepl("chr[0-9XY]+$", seqnames))
        sum(stats[,3])
    })
    
    dt.box$mapped_reads = mapped_counts[dt.box$sample]
    
    dt.box[, frip := N/mapped_reads]
    
    ##name_lev = dt.box[, median(N) , .(name)][rev(order(V1))]$name
    
    name_lev <- c("ctcf_mock_rep_1",
                  "ctcf_mock_rep_2",     
                  "ctcf_mock_rep_3",
                  "ctcf_infected_rep_1",
                  "ctcf_infected_rep_2",  
                  "ctcf_infected_rep_3",  
                  "input_mock_rep_1",     
                  "input_mock_rep_2",     
                  "input_mock_rep_3",     
                  "input_infected_rep_1",
                  "input_infected_rep_2",
                  "input_infected_rep_3")
    
    dt.box$name = factor(dt.box$name, levels = name_lev)
    
    peak_dt$name = factor(peak_dt$name, levels = name_lev)
    
    fq_dt$name = factor(fq_dt$name, levels = name_lev)
    
    message("make plots...")
    
    p_fq1 = ggplot(fq_dt, 
                   aes(x = name, y = fastq_count, fill = treatment)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      scale_y_continuous(labels = function(x)x/1e6) +
      labs(y = "FASTQ reads (M)", x = "")
    
    p_input1 = ggplot(unique(dt.box[, .(name, mapped_reads, treatment)]), 
                      aes(x = name, y = mapped_reads, fill = treatment)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      scale_y_continuous(labels = function(x)x/1e6) +
      labs(y = "BAM mapped reads (M)", x = "")
    
    p_peaks1 = ggplot(peak_dt, aes(x = name, y = peak_count, fill = treatment)) +
      geom_bar(stat = 'identity') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      scale_y_continuous(labels = function(x)x/1e3) +
      labs(y = "peaks (k)", x = "")
    
    p_reads1 = ggplot(dt.box, 
                      aes(x = name, y = N, color = treatment)) +
      geom_boxplot(outlier.shape = NA) +
      coord_cartesian(ylim = quantile(dt.box$N, c(0.1, 0.96))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      labs(y = "reads per peak", x = "")
    
    p_frip1 = ggplot(dt.box, 
                     aes(x = name, y = frip, color = treatment)) +
      geom_boxplot(outlier.shape = NA) +
      coord_cartesian(ylim = quantile(dt.box$frip, c(0.1, 0.96))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      labs(y = "FRIP per peak", x = "")
    
    tmp = dt.box[, .(N = sum(N), mapped_reads = unique(mapped_reads)), .(sample, treatment, name)]
    tmp[, frip := N / mapped_reads]
    p_fripSum1 = ggplot(tmp, 
                        aes(x = name, y = frip, fill = treatment)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      labs(subtitle = paste(sum(width((query_gr)))/3.2e9, 
                            "of genome covered by peaks"), 
                            x = "", 
                            y = "FRIP")
    
    ##name_lev_frip = as.character(dt.box[, median(frip) , .(name)][rev(order(V1))]$name)
    
    name_lev_frip <- c("ctcf_mock_rep_1",
                       "ctcf_mock_rep_2",
                       "ctcf_mock_rep_3",
                       "ctcf_infected_rep_1",
                       "ctcf_infected_rep_2",
                       "ctcf_infected_rep_3",
                       "input_mock_rep_1",
                       "input_mock_rep_2",
                       "input_mock_rep_3",
                       "input_infected_rep_1",
                       "input_infected_rep_2",
                       "input_infected_rep_3")
    
    dt.box$name = factor(dt.box$name, levels = name_lev_frip)
    fq_dt$name = factor(fq_dt$name, levels = name_lev_frip)
    peak_dt$name = factor(peak_dt$name, levels = name_lev_frip)
    
    p_fq2 = ggplot(fq_dt, 
                   aes(x = name, y = fastq_count, fill = treatment)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      scale_y_continuous(labels = function(x)x/1e6) +
      labs(y = "FASTQ reads (M)", x = "")
    
    p_input2 = ggplot(unique(dt.box[, .(name, mapped_reads, treatment)]), 
                      aes(x = name, y = mapped_reads, fill = treatment)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            plot.margin = margin(.01, .01, .01, .1, unit = "npc"))  +
      scale_y_continuous(labels = function(x)x/1e6) +
      labs(y = "BAM mapped reads (M)", x = "")
    
    p_peaks2 = ggplot(peak_dt, aes(x = name, y = peak_count, fill = treatment)) +
      geom_bar(stat = 'identity') +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), 
            plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      scale_y_continuous(labels = function(x)x/1e3) +
      labs(y = "peaks (k)", x = "")
    
    p_reads2 = ggplot(dt.box, aes(x = name, y = N, color = treatment)) +
      geom_boxplot(outlier.shape = NA) +
      coord_cartesian(ylim = quantile(dt.box$N, c(0.1, 0.96))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      labs(y = "reads per peak", x = "")
    
    
    p_frip2 = ggplot(dt.box, aes(x = name, y = frip, color = treatment)) +
      geom_boxplot(outlier.shape = NA) +
      coord_cartesian(ylim = quantile(dt.box$frip, c(0.1, 0.96))) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      labs(y = "FRIP per peak", x = "")
    
    tmp = dt.box[, .(N = sum(N), mapped_reads = unique(mapped_reads)), .(sample, treatment, name)]
    tmp[, frip := N / mapped_reads]
    p_fripSum2 = ggplot(tmp, 
                        aes(x = name, y = frip, fill = treatment)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      labs(subtitle = paste(sum(width((query_gr)))/3.2e9, "of genome covered by peaks"), x = "", y = "FRIP")
    
    message("assemble plots...")
    col1 = list(p_fq1, p_input1, p_reads1, p_peaks1, p_frip1, p_fripSum1)
    col2 = list(p_fq2, p_input2, p_reads2, p_peaks2, p_frip2, p_fripSum2)
    
    col1 = ssvRecipes::sync_width(col1)
    col2 = ssvRecipes::sync_width(col2)
    
    pg = cowplot::plot_grid(ncol = 2,
                            cowplot::plot_grid(plotlist = col1, ncol = 1),
                            cowplot::plot_grid(plotlist = col2, ncol = 1)
    )
    ggsave(paste0(output.dir,"/frip_boxes.", file_tag, ".pdf"), 
           pg, 
           width = 13.3, 
           height = 18)
    
    unique(bw_dt$sample)

    bw_dt.all$name = factor(bw_dt.all$name, levels = name_lev_frip)
    bw_dt.all$facet = bw_dt.all$name
    levels(bw_dt.all$facet) = gsub("_", "\n", levels(bw_dt.all$facet))
    
    set.seed(0)
    
    clust_dt = ssvSignalClustering(bw_dt.all, 
                                   fill_ = "y_relative", 
                                   max_cols = Inf, 
                                   facet_ = "facet", 
                                   max_rows = Inf)
    
    toplot_id = ssvRecipes::sampleCap(unique(clust_dt$id), 500)
    
    p_heat = ssvSignalHeatmap(clust_dt[id %in% toplot_id],
                              fill_ = "y_relative", 
                              max_cols = Inf, 
                              facet_ = "facet", show_cluster_bars = FALSE)
    
    p_heat = p_heat + 
      labs(x = paste(view_size, "bp view size"), fill = "relative pileup") + 
      theme(axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            legend.position = "bottom")
    
    
    
    
    
    
    ref_gr = rtracklayer::import.gff("/slipstream/home/joeboyd/gencode.v35.annotation.gtf.gz", format = "gtf")
    head(ref_gr)
    table(ref_gr$type)
    gene_gr = reduce(subset(ref_gr, type == "gene"))
    exon_gr = reduce(subset(ref_gr, type == "exon" & tag == "basic"))
    intron_gr = setdiff(gene_gr, exon_gr)
    tx_gr = subset(ref_gr, type == "transcript" & tag == "basic")
    tss_gr = flank(tx_gr, 1e3, start = TRUE, both = TRUE)
    tts_gr = flank(tx_gr, 1e3, start = FALSE, both = TRUE)
    artifact_gr = rtracklayer::import.bed("/slipstream/home/joeboyd/R/qc_cutnrun/reference/blacklist_hg38.bed")
    anno_grs = rev(list(artifact = artifact_gr, tss = tss_gr, tts = tts_gr, exon = exon_gr, intron = intron_gr, genebody = gene_gr))
    
    
    
    anno_clust_dt = ssvFetchGRanges(anno_grs, resize(query_gr[toplot_id], view_size, fix = "center"), return_data.table = TRUE)
    anno_clust_dt$id = factor(anno_clust_dt$id, levels = levels(clust_dt$id))
    anno_clust_dt$sample = factor(anno_clust_dt$sample, levels = rev(names(anno_grs)))
    p_heat_anno = ggplot(anno_clust_dt, aes(x = x, y = id, fill = y>0)) +
      geom_raster() +
      facet_wrap(~sample, nrow = 1) +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
      scale_fill_manual(values = c("FALSE" = "gray80", "TRUE" = "gray20")) +
      scale_x_continuous(labels = function(x)x/1e3) +
      labs(fill = "feature overlap", x= "kbp", y = "peak region")
    
    clust_sizes = unique(clust_dt[, .(cluster_id, id)])[, .N, .(cluster_id)][rev(order(cluster_id))]
    clust_sizes[, xmin := 0]
    clust_sizes[, xmax := 1]
    clust_sizes[, ymin := c(0, cumsum(N)[-length(N)])]
    clust_sizes[, ymax := cumsum(N)]
    clust_sizes[, col := as.character(cluster_id%%2)]
    clust_sizes$facet = ""
    p_clust_boxes = ggplot(clust_sizes, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax)) + 
      geom_rect(aes(fill = col), color = "black") +
      geom_text(aes(label = cluster_id, x = (xmin+xmax)/2, y = (ymin+ymax)/2)) +
      scale_fill_manual(values = c("0" = "gray70", "1" = 'gray90')) +
      guides(fill = 'none') +
      coord_cartesian(expand = FALSE) +
      theme_void() + 
      facet_grid(~facet)
    
    
    p_heat_l = ssvRecipes::sync_height(list(p_clust_boxes, p_heat, p_heat_anno))
    pg_heat_top = cowplot::plot_grid(plotlist = p_heat_l, rel_widths = c(1, length(peak_grs)+2, length(anno_grs)+1), nrow = 1)
    
    dt.heat = copy(dt.box)
    assign_dt = unique(clust_dt[, .(id, cluster_id)])
    dt.heat = merge(dt.heat, assign_dt, by = "id")
    tmp = dt.heat[, .(N = sum(N), mapped_reads = unique(mapped_reads)), .(sample, treatment, name, cluster_id)]
    tmp[, frip := N / mapped_reads]
    
    w_frac = lapply(split(assign_dt$id, assign_dt$cluster_id), function(x){
      sum(width(query_gr[x]))  /3.2e9
    })
    w_frac = unlist(w_frac)
    tmp[, cluster_label := paste0("cluster", cluster_id, "\n", w_frac[cluster_id])]
    
    
    p_fripHeat = ggplot(tmp, 
                        aes(x = name, y = frip, fill = treatment)) +
      geom_bar(stat = "identity") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), plot.margin = margin(.01, .01, .01, .1, unit = "npc")) +
      labs(subtitle = paste(sum(width((query_gr)))/3.2e9, "of genome covered by peaks"), x = "", y = "FRIP") +
      facet_grid(cluster_label~., scales = "free_y", switch = "y") +
      theme(strip.text.y = element_text(angle = 0), strip.placement = "outside")
    # p_fripHeat
    # g_fripHeat = ssvRecipes::sync_width(list(p_fripHeat, p_heat))[[1]]
    
    peak_dt
    fq_dt
    
    mapped_dt = unique(dt.box[, .(name, mapped_reads)])
    
    meta_dt = merge(fq_dt, peak_dt[, .(name, peak_count)], by = "name")
    meta_dt = merge(meta_dt, mapped_dt[, .(name, mapped_reads)], by = "name")
    
    fwrite(meta_dt, paste0(output.dir,"/sample_info.", file_tag, ".csv"), sep = ",")
    
    region_frip_dt = dcast(dt.box[, .(id, name, N, treatment, frip)], id~name, value.var = "frip")
    setkey(region_frip_dt, id)
    assign_dt
    setkey(assign_dt, id)
    bed_frip_towrite = query_gr
    # bed_frip_towrite$id = NULL
    bed_peak_towrite = query_gr
    # bed_peak_towrite$id = NULL
    
    mcols(bed_frip_towrite) = region_frip_dt[.(names(bed_frip_towrite))]
    # bed_frip_towrite$id = NULL
    
    message("write files...")
    bed_frip_towrite$cluster_id = assign_dt[.(names(bed_frip_towrite))]$cluster_id
    # rtracklayer::export.bed(bed_frip_towrite, paste0("consensus_regions_with_FRIP.", file_tag, ".bed"))
    
    bed_frip_towrite_dt = as.data.table(bed_frip_towrite)
    bed_frip_towrite_dt$strand = "."
    bed_frip_towrite_dt$score = 0
    extra_cols = setdiff(colnames(mcols(bed_frip_towrite)), "id")
    extra_cols = gsub("-", ".", extra_cols)
    fwrite(bed_frip_towrite_dt[, c("seqnames", "start", "end", "id", "score", "strand", extra_cols), with = FALSE], 
           file = paste0(output.dir,"/consensus_regions_with_FRIP.", file_tag, ".txt"),
           sep = "\t")
    
    bed_peak_towrite$cluster_id = assign_dt[.(names(bed_peak_towrite))]$cluster_id
    bed_peak_towrite_dt = as.data.table(bed_peak_towrite)
    bed_peak_towrite_dt$strand = "."
    bed_peak_towrite_dt$score = 0
    extra_cols = setdiff(colnames(mcols(bed_peak_towrite)), "id")
    extra_cols = gsub("-", ".", extra_cols)
    fwrite(bed_peak_towrite_dt[, c("seqnames", "start", "end", "id", "score", "strand", extra_cols), with = FALSE], 
           file = paste0(output.dir,"/consensus_regions_with_peak_call.", file_tag, ".txt"),
           sep = "\t")
    # rtracklayer::export.bed(bed_peak_towrite, paste0("consensus_regions_with_peak_call.", file_tag, ".bed"))
    
    message("features overlap")
    
    to_overlap_grs = c(peak_grs, peak_grs.neg)
    peak_grs.anno = lapply(to_overlap_grs, function(x){
      x$anno = "intergenic"
      for(i in seq_along(anno_grs)){
        olaps = findOverlaps(x, anno_grs[[i]])
        x$anno[queryHits(olaps)] = names(anno_grs)[i]
        
      }
      x
      # table(x$anno)
    })
    peak_grs.anno_cnt = lapply(peak_grs.anno, function(x){
      tab = table(x$anno)
      data.table(feature = names(tab), count = as.numeric(tab))
    })
    
    anno_cnt = rbindlist(peak_grs.anno_cnt, idcol = "sample")
    anno_cnt[, sample_cnt := paste0(sample, "\n", sum(count)), .(sample)]
    anno_cnt[, fraction := count / sum(count), .(sample)]
    anno_cnt$sample = factor(anno_cnt$sample, levels = name_lev_frip)
    anno_cnt = anno_cnt[order(sample)]
    anno_cnt$sample_cnt = factor(anno_cnt$sample_cnt, 
                                 levels = unique(anno_cnt$sample_cnt))
    
    p_features = ggplot(anno_cnt, 
                        aes(x = sample_cnt, 
                            y = fraction, 
                            fill = feature)) +
      geom_bar(stat = "identity") +
      labs(x = "sample\npeak count") +
      theme(panel.background = element_blank(), 
            panel.grid.major.y = element_line(color = "black"), 
            panel.grid.minor.y = element_line(color = "black"),
            panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank()) +
      scale_fill_brewer(palette = "Set1") +
      labs(title = "Feature overlaps for peaks") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    
    
    dt.heat_summary = dt.heat[, .(mean_frip = mean(frip), median_frip = median(frip)), .(name, cluster_id)]
    dt.heat_summary[, txt_mean := round(mean_frip*1e6, 1)]
    dt.heat_summary[, txt_median := round(median_frip*1e6, 1)]
    dt.heat_summary[, txt := paste0(round(mean_frip*1e6, 1), "\n", round(median_frip*1e6, 1))]
    
    p_clust_text = ggplot(dt.heat_summary, aes(x = name, y = factor(cluster_id), label = txt)) +
      geom_text(size = 5, color = NA) +
      geom_text(data = dt.heat_summary, aes(x = name, cluster_id+.1, label = txt_mean), color = "red", vjust = 0) +
      geom_text(data = dt.heat_summary, aes(x = name, cluster_id-.1, label = txt_median), color = "blue", vjust = 1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1), panel.background = element_blank(), panel.grid = element_blank()) +
      labs(y = "cluster", x= "", title = "mean FRIP e6", subtitle = "median FRIP e6") +
      theme(plot.title = element_text(size = 14, color = "red"), plot.subtitle = element_text(size = 14, color = "blue"))
    
    pg_heat = cowplot::plot_grid(ncol = 1,
                                 pg_heat_top,
                                 cowplot::plot_grid(nrow = 1, p_fripHeat, p_clust_text)) 
    ggsave(paste0(output.dir,"/frip_heatmap.", file_tag, ".pdf"), pg_heat, width = 8+length(peak_grs), height = 14)
    ggsave(paste0(output.dir,"/feature_overlap.", file_tag, ".pdf"), p_features, width = 2+length(peak_grs)*.6, height = 6)
    
    fwrite(anno_cnt, 
           file = paste0(output.dir,"/feature_counts.", file_tag, ".csv"),
           sep = ",")
    
    # utr_gr = subset(ref_gr, gene_name == "IKZF1" & type == "UTR")
    # split(utr_gr, utr_gr$transcript_name)
    message("done")
    
    ##setwd(base_dir)

    