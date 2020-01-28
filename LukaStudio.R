# LukaStudio v0.9.2
# This script exploits multiple threads for parallel computing:
# This script has been tested under ubuntu 19.10, to set up the environment with terminal:
  # sudo apt update && sudo apt install r-base r-cran-rjava -y
  # sudo apt install libxml2-dev libssl-dev libcurl4-openssl-dev libxml2-dev -y
  # sudo -i R
  # install.packages(c("filesstrings", "R.utils", "vcfR", "rvest", "dplyr", "tidyr", "doParallel", "foreach", "writexl"))
# If you are using it under other systems, make sure that all packages are properly installed
# Put ".vcf" or ".vcf.gz" files in the same folder with this script before starting


# Load R packages
lapply(c("vcfR", "parallel", "foreach", "doParallel", "R.utils", 
         "dplyr", "tidyr", "filesstrings", "rlang", "writexl"), 
       require, character.only = TRUE)
print("Welcome to LukaStudio, the system is loading...")


# Set the folder containing current R script as working directory
setwd(".")
print(paste("Current working dir:", getwd()))
if (dir.exists("original_files")==FALSE){
  dir.create("original_files")
}
if (dir.exists("cache")==FALSE){
  dir.create("cache")
}


# Set CPU cores for parallel computing
numCores <- detectCores(all.tests = FALSE, logical = TRUE)
print(paste("Parallel computing:", numCores, "cores will be used for data processing"))
Sys.sleep(1)


# Convert all VCF files to csv files, decompress the .vcf.gz files if necessary
gzfiles <- list.files(pattern=".vcf.gz")
if (is_empty(gzfiles)==FALSE){
  registerDoParallel(numCores)
  foreach (gzfile = gzfiles) %dopar% {
    gunzip(gzfile, remove=FALSE, overwrite=TRUE)
    file.move(gzfile, "./original_files/c_gz_files", overwrite=TRUE)
  }
  print("Decompressed all vcf.gz files")
}
rm("gzfiles")
vcffiles <- list.files(pattern=".vcf")
hcvcfs <- list.files(pattern=".HC.vcf")
mdlvcs <- list.files(pattern=".MDLVC.vcf")
if (is_empty(vcffiles)==FALSE){
  print("Scanning for all .vcf files...")
  print(vcffiles)
  print("Extracting variants info from VCF files...(this might take a while)")
  if (is_empty(hcvcfs)==FALSE){
    registerDoParallel(numCores)
    foreach (hcvcf = hcvcfs) %dopar% {
      allvariants <- read.vcfR(hcvcf, verbose=FALSE)
      filename <- paste0(hcvcf, ".hc.csv")
      write.table(cbind(allvariants@fix, allvariants@gt), 
                  file=filename, sep=",", row.name=FALSE)
      file.move(hcvcf, "./original_files/a_hcvcf_files", overwrite=TRUE)
    }
  }
  rm("hcvcfs")
  if (is_empty(mdlvcs)==FALSE){
    registerDoParallel(numCores)
    foreach (mdlvc = mdlvcs) %dopar% {
      allvariants <- readLines(mdlvc)
      allvariants <- as.data.frame(allvariants)
      allvariants <-allvariants[-grep('#',allvariants$allvariants),]
      filename <- paste0(mdlvc, ".txt")
      write.table(allvariants, file = filename, sep = "\t", row.names = FALSE)
      file.move(mdlvc, "./original_files/b_mdlvc_files", overwrite=TRUE)
    }
    mdlvcs <- list.files(pattern=".MDLVC.vcf.txt")
    registerDoParallel(numCores)
    foreach (mdlvc = mdlvcs) %dopar% {
      allvariants <- read.table(mdlvc,header=T,sep="\t");
      is.data.frame(allvariants)
      allvariants <- separate(allvariants, "x",
                        into=c("CHROM","POS","ID", "REF", "ALT", "QUAL", 
                               "FILTER", "INFO", "FORMAT", "Code"),sep="\t")
      file.move(mdlvc, "./cache/a_txt_files", overwrite=TRUE)
      haplotype <- filter(allvariants, INFO == "VCcaller=HaplotypeCaller;GT:AD:GQ:PL:SAC")
      mdlvccall <- filter(allvariants, INFO != "VCcaller=HaplotypeCaller;GT:AD:GQ:PL:SAC")
      filename1 <- paste0(mdlvc, ".haplotype.csv")
      filename2 <- paste0(mdlvc, ".mdlvccall.csv")
      write.table(haplotype, file=filename1, sep=",", row.names = FALSE)
      write.table(mdlvccall, file=filename2, sep=",", row.names = FALSE)
    }
    rm("mdlvcs")
  }
  print("All variants info is extracetd")
}
rm("vcffiles")
Sys.sleep(2)


# Translate "Genotype" part
rawcsvs <- list.files(pattern=".csv")
rawcsvas <- list.files(pattern=".hc.csv")
rawcsvbs <- list.files(pattern=".haplotype.csv")
rawcsvcs <- list.files(pattern=".mdlvccall.csv")
if (is_empty(rawcsvs)==FALSE){
  print("Translating variants info...(part 1 of 2)")
  if (is_empty(rawcsvas)==FALSE){
    registerDoParallel(numCores)
    foreach (rawcsva = rawcsvas) %dopar% {
      csv1 <- read.csv(rawcsva, header = TRUE)
      if(nrow(csv1)>1){
        colnames(csv1) <- c("CHROM","POS","ID", "REF", "ALT", "QUAL", 
                            "FILTER", "INFO", "FORMAT", "Code")
        csv1 <- separate(csv1, "Code",
                         into=c("Genotype", 
                                "Allelic_depths", 
                                "Genotype_Quality", 
                                "Phred-scaled_likelihoods", 
                                "Number_of_supporting_strand"), sep=":")
      }
      filename <- paste0(rawcsva, ".hctrans.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
      file.move(rawcsva, "./cache/b_raw_hc_csv", overwrite=TRUE)
    }
  }
  rm("rawcsvas")
  if (is_empty(rawcsvbs)==FALSE){
    registerDoParallel(numCores)
    foreach (rawcsvb = rawcsvbs) %dopar% {
      csv1 <- read.csv(rawcsvb, header = TRUE)
      if(nrow(csv1)>1){
        csv1 <- separate(csv1, "FORMAT",
                         into=c("Genotype", 
                                "Allelic_depths", 
                                "Genotype_Quality", 
                                "Phred_scaled_likelihoods", 
                                "Supporting_allele_counts"), sep=":")
        csv1 <- separate(csv1, "INFO", into=c("INFO", "FORMAT"), sep=";")
      }
      filename <- paste0(rawcsvb, ".haptrans.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
      file.move(rawcsvb, "./cache/c_raw_hap_csv", overwrite=TRUE)
    }
  }
  rm("rawcsvbs")
  if (is_empty(rawcsvcs)==FALSE){
    registerDoParallel(numCores)
    foreach (rawcsvc = rawcsvcs) %dopar% {
      csv1 <- read.csv(rawcsvc, header = TRUE)
      if(nrow(csv1)>1){
        csv1 <- separate(csv1, "Code",
                         into=c("Genotype", "Ref_Alt_depth"), sep=":")
      }
      filename <- paste0(rawcsvc, ".mdltrans.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
      file.move(rawcsvc, "./cache/d_raw_mdl_csv", overwrite=TRUE)
    }
  }
  rm("rawcsvcs")
  print("All variants info is translated (part 1 of 2)")
}
rm("rawcsvs")
Sys.sleep(2)


# Translate "INFO" part
medcsvs <- list.files(pattern=".csv")
medcsvas <- list.files(pattern=".hctrans.csv")
medcsvbs <- list.files(pattern=".haptrans.csv")
medcsvcs <- list.files(pattern=".mdltrans.csv")
if (is_empty(medcsvs)==FALSE){
  print("Translating variants info...(part 2 of 2)")
  if (is_empty(medcsvas)==FALSE){
    registerDoParallel(numCores)
    foreach (medcsva = medcsvas) %dopar% {
      csv1 <- read.csv(medcsva, header = TRUE)
      if(nrow(csv1)>1){
        csv1$INFO <- gsub(';BaseQRankSum=', ':BaseQRankSum=', csv1$INFO)
        csv1$INFO <- gsub(';MQRankSum=', ':MQRankSum=', csv1$INFO)
        csv1$INFO <- gsub(';ReadPosRankSum=', ':ReadPosRankSum=', csv1$INFO)
        csv1$INFO <- gsub('AN=2;DB', 'AN=2::DB', csv1$INFO)
        csv1$INFO <- gsub(';DB', ':DB', csv1$INFO)
        csv1$Sample <- rep(medcsva,nrow(csv1))
        csv1 <- separate(csv1, "INFO",
                         into=c("Allele_Count",
                                "Allele_Frequency",
                                "Allele_Number", 
                                "Read_Depth",
                                "Phred-scaled_p-value",
                                "MLEAC",
                                "MLEAF",
                                "Mapping_Quality",
                                "Mapping_Quality0",
                                "Number_of_alternate_alleles_discovered",
                                "Quality_by_Depth",
                                "Symmetric Odds Ratio"), sep=";")
        csv1 <- separate(csv1, "Allele_Number",
                         into=c("Allele_Number",
                                "BaseQRankSum", 
                                "dbSNP_Membership"), sep=":")
        csv1 <- separate(csv1, "Mapping_Quality0",
                         into=c("Mapping_Quality0", 
                                "MQRankSum"), sep=":")
        csv1 <- separate(csv1, "Quality_by_Depth",
                         into=c("Quality_by_Depth", 
                                "ReadPosRankSum"), sep=":")
      }
      filename <- paste0(medcsva, ".hctrim.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
      file.move(medcsva, "./cache/e_med_hc_csv", overwrite=TRUE)
    }
  }
  rm("medcsvas")
  Sys.sleep(2)
  if (is_empty(medcsvbs)==FALSE){
    registerDoParallel(numCores)
    foreach (medcsvb = medcsvbs) %dopar% {
      csv1 <- read.csv(medcsvb, header = TRUE)
      if(nrow(csv1)>1){
        csv1$ID <- gsub('.;', '', csv1$ID)
        csv1$ID <- gsub('-', '', csv1$ID)
        csv1$Sample <- rep(medcsvb,nrow(csv1))
        colnames(csv1) <- c("CHROM","POS","ID","REF","ALT","QUAL","FILTER","VCcaller","FORMAT","Genotype","Allelic_depths","Genotype_Quality",
                            "PL: Normalized Phred-scaled likelihoods for genotypes as defined in the VCF specification",
                            "SAC: Number of reads on the forward and reverse strand supporting each allele",
                            "Code","Sample")
      }
      filename <- paste0(medcsvb, ".haptrim.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
      file.move(medcsvb, "./cache/f_med_hap_csv", overwrite=TRUE)
    }
  }
  rm("medcsvbs")
  Sys.sleep(2)
  if (is_empty(medcsvcs)==FALSE){
    registerDoParallel(numCores)
    foreach (medcsvc = medcsvcs) %dopar% {
      csv1 <- read.csv(medcsvc, header = TRUE)
      if(nrow(csv1)>1){
        csv1$ID <- gsub('.', '', csv1$ID)
        csv1$Sample <- rep(medcsvc,nrow(csv1))
        csv1 <- separate(csv1, "INFO",
                         into=c("TYPE", 
                                "VCcaller", 
                                "FSAF", 
                                "FSAR", 
                                "FSRF", 
                                "FSRR", 
                                "Read_depth", 
                                "Allelic_Frequency", 
                                "Alternate_allele_observed", 
                                "PoN", 
                                "control_VAF", 
                                "Strand_Bias_Test"), sep=";")
        csv1 <- separate(csv1, "control_VAF",
                         into=c("control_VAF", 
                                "controlVAF_sd"), sep=",")
      }
      filename <- paste0(medcsvc, ".mdltrim.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
      file.move(medcsvc, "./cache/g_med_mdl_csv", overwrite=TRUE)
    }
  }
  rm("medcsvcs")
  print("All variants info is translated (part 2 of 2)")
  Sys.sleep(2)
}
rm("medcsvs")


# Trim data
medwelcsvs <- list.files(pattern=".csv")
medwelcsvas <- list.files(pattern=".hctrim.csv")
medwelcsvbs <- list.files(pattern=".haptrim.csv")
medwelcsvcs <- list.files(pattern=".mdltrim.csv")
if (is_empty(medwelcsvs)==FALSE){
  print("Triming data...(this might take a while)")
  if (is_empty(medwelcsvas)==FALSE){
    registerDoParallel(numCores)
    foreach (medwelcsva = medwelcsvas) %dopar% {
      csv1 <- read.csv(medwelcsva, header = TRUE, na.strings=c("","NA"))
      if(nrow(csv1)>1){
        csv1$Allele_Count <- gsub('AC=', '', csv1$Allele_Count)
        csv1$Allele_Frequency <- gsub('AF=', '', csv1$Allele_Frequency)
        csv1$Allele_Number <- gsub('AN=', '', csv1$Allele_Number)
        csv1$BaseQRankSum <- gsub('BaseQRankSum=', '', csv1$BaseQRankSum)
        csv1$dbSNP_Membership <- gsub('DB', 'YES', csv1$dbSNP_Membership)
        csv1$Read_Depth <- gsub('DP=', '', csv1$Read_Depth)
        csv1$MQRankSum <- gsub('MQRankSum=', '', csv1$MQRankSum)
        csv1$ReadPosRankSum <- gsub('ReadPosRankSum=', '', csv1$ReadPosRankSum)
        csv1$Sample <- gsub('.HC.vcf.hc.csv.hctrans.csv', '', csv1$Sample)
        csv1 <- csv1[c(30, 1, 2, 4, 5, 6, 7, 3, 25, 26, 27, 28, 29, 8, 9, 10, 
                       11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23)]
      }
      filename <- paste0(medwelcsva, ".hctrim.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
      file.move(medwelcsva, "./cache/h_hc_trim_csv", overwrite=TRUE)
    }
  }
  rm("medwelcsvas")
  if (is_empty(medwelcsvbs)==FALSE){
    registerDoParallel(numCores)
    foreach (medwelcsvb = medwelcsvbs) %dopar% {
      csv1 <- read.csv(medwelcsvb, header = TRUE, na.strings=c("","NA"))
      if(nrow(csv1)>1){
        csv1$VCcaller <- gsub('VCcaller=', '', csv1$VCcaller)
        csv1$Sample <- gsub('.MDLVC.vcf.txt.haplotype.csv.haptrans.csv', '', csv1$Sample)
        csv1 <- csv1[c(16, 1, 2, 4, 5, 6, 7, 3, 10, 11, 12, 13, 14, 8)]
      }
      filename <- paste0(medwelcsvb, ".hatrim.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
      file.move(medwelcsvb, "./cache/i_hap_trim_csv", overwrite=TRUE)
    }
  }
  rm("medwelcsvbs")
  if (is_empty(medwelcsvcs)==FALSE){
    registerDoParallel(numCores)
    foreach (medwelcsvc = medwelcsvcs) %dopar% {
      csv1 <- read.csv(medwelcsvc, header = TRUE, na.strings=c("","NA"))
      if(nrow(csv1)>1){
        csv1$VCcaller <- gsub('VCcaller=', '', csv1$VCcaller)
        csv1$Read_depth <- gsub('DP=', '', csv1$Read_depth)
        csv1$Allelic_Frequency <- gsub('AF=', '', csv1$Allelic_Frequency)
        csv1$control_VAF <- gsub('PNAF=', '', csv1$control_VAF)
        csv1$control_VAF <- gsub('-', '', csv1$control_VAF)
        csv1$Sample <- gsub('.MDLVC.vcf.txt.mdlvccall.csv.mdltrans.csv', '', csv1$Sample)
        csv1 <- csv1[c(24, 1, 2, 4, 5, 6, 7, 3, 22, 14, 23, 8, 15, 18, 19, 
                       16, 17, 10, 11, 12, 13, 20, 9)]
      }
      filename <- paste0(medwelcsvc, ".mdtrim.csv")
      write.table(csv1, file=filename, sep=",", row.names = FALSE)
      file.move(medwelcsvc, "./cache/j_mdl_trim_csv", overwrite=TRUE)
    }
  }
  rm("medwelcsvcs")
  print("All variants info is reorganized")
}
rm("medwelcsvs")
Sys.sleep(2)


# Filter data
welcsvcs <- list.files(pattern=".mdtrim.csv")
if (is_empty(welcsvcs)==FALSE){
  print("Filtering data...(this might take a while)")
  registerDoParallel(numCores)
  foreach (welcsvc = welcsvcs) %dopar% {
    csv1 <- read.csv(welcsvc, header = TRUE)
    if(nrow(csv1)>1){
      csv1 <- mutate(csv1, VAFtoCVAF=Allelic_Frequency/control_VAF)
      csv1$VAFtoCVAF <- gsub('Inf', 'NA', csv1$VAFtoCVAF)
      csv1 <- filter(csv1 ,Allelic_Frequency >2 & Allelic_Frequency <45)
      csv1 <- filter(csv1 ,VAFtoCVAF >2 | VAFtoCVAF == "NA")
      csv1 <- csv1[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 24, 
                     15, 16, 17, 18, 19, 20, 21, 22, 23)]
    }
    filename <- paste0(welcsvc, ".mdlfilter.csv")
    write.table(csv1, file=filename, sep=",", row.names = FALSE)
    file.move(welcsvc, "./cache/k_mdl_filter_csv", overwrite=TRUE) 
  }
  print("All variants info filtered")
}
rm("welcsvcs")
Sys.sleep(2)


# Merge and export
results <- list.files(pattern=".csv")
resultas <- list.files(pattern=".hctrim.csv")
resultbs <- list.files(pattern=".hatrim.csv")
resultcs <- list.files(pattern=".mdlfilter.csv")
if (is_empty(results)==FALSE){
  print("Exporting final results...")
  if (is_empty(resultas)==FALSE){
    csv1 <- lapply(resultas, read.csv)
    csv1 <- bind_rows(csv1)
  }
  rm("resultas")
  if (is_empty(resultbs)==FALSE){
    csv2 <- lapply(resultbs, read.csv)
    csv2 <- bind_rows(csv2)
    csv2 <- filter(csv2, Sample != "NA")
    csv2 <- csv2[c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14)]
  }
  rm("resultbs")
  if (is_empty(resultcs)==FALSE){
    csv3 <- lapply(resultcs, read.csv)
    csv3 <- bind_rows(csv3)
  }
  rm("resultcs")
  registerDoParallel(numCores)
  foreach (result = results) %dopar% {
    file.move(result, "./cache/l_merge_csv", overwrite=TRUE)
  }
  now <- Sys.time()
  filename <- paste0(format(now, "%Y-%m-%d_%H:%M:%S_"), "data_export.xlsx")
  sheets <- list("HC(germline)" = csv1, "Haplotype" = csv2, "MDLVC(filtered)" = csv3)
  write_xlsx(sheets, filename)
  rm("now", "filename", "csv1", "csv2", "csv3", "sheets")
  print("All done! Good luck analysing!")
}
rm("results", "numCores")

