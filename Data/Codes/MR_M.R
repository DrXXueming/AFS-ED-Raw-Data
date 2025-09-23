library(TwoSampleMR)
library(dplyr)
library(ggplot2)
library(htmltools)
library(htmlwidgets)
library(httpuv)
library(ieugwasr)
library(knitr)
library(MRInstruments)
library(phenoscanner)
library(FastTraitR)
library(plinkbinr)

DIR<-setwd("/Users/emmetthui/Desktop/0/Jiani/AFS_ED/(male)3_1")
#file.create(file.path(DIR, '#confounder_SNPs.txt'))

#exp<-'ukb-b-6591'
#exp<-'ebi-a-GCST90000047'
exp<-'ebi-a-GCST90000046'
#exp<-'ebi-a-GCST90000045'

out<-'ebi-a-GCST006956'
#out<-'finn-b-ERECTILE_DYSFUNCTION'

#data(gwas_catalog)
#GC_dat <-subset(gwas_catalog,grepl("Peyrot WJ", Author) &Phenotype == "Attention deficit hyperactivity disorder")
#exp_dat <-format_data(GC_dat)

exp_dat <- extract_instruments(outcomes = exp,p1=5e-07,clump=FALSE)

unclump_snps <- ieugwasr::ld_clump(dat = dplyr::tibble(rsid = exp_dat$SNP, pval = exp_dat$pval.exposure, id = exp_dat$id.exposure),
                                   clump_kb = 5000,
                                   clump_r2 = 0.01,
                                   plink_bin = get_plink_exe(),
                                   bfile = file.path('/Users/emmetthui/Desktop/0/Jiani/LD/1kg.v3', 'EUR'))
exp_dat <- exp_dat %>%
  dplyr::inner_join(unclump_snps, by = c("SNP" = "rsid")) %>%
  dplyr::select(names(.))

exp_dat<-exp_dat[exp_dat$beta.exposure^2/exp_dat$se.exposure^2>10,]

snp_with_trait <- FastTraitR::look_trait(rsids = exp_dat$SNP, out_file = 'check_SNPs_trait.csv')
snp_with_trait_save <- snp_with_trait %>%
  arrange(trait) %>%
  select(trait) %>%
  distinct()  
writeLines(snp_with_trait_save$trait, 'check_SNPs_trait.txt')  # 保存到文件


confounders <- readLines("#confounder_SNPs.txt")
snp_with_trait$trait <- tolower(snp_with_trait$trait)  # 确保 trait 列文本均为小写
for (confounder in confounders) {
  snp_with_trait <- snp_with_trait[!grepl(tolower(confounder), snp_with_trait$trait),]
}
snp_with_trait <- dplyr::distinct(snp_with_trait, rsid, .keep_all = FALSE)  # 去重
exp_dat <- exp_dat %>%
  dplyr::inner_join(snp_with_trait, by = c("SNP" = "rsid")) %>%
  dplyr::select(names(exp_dat))

out_dat <- extract_outcome_data(snps = exp_dat$SNP,outcomes = out)

dat <- harmonise_data(exposure_dat = exp_dat,outcome_dat = out_dat)

mr_presso_result <- run_mr_presso(dat)
write.csv(mr_presso_result[[1]]$
            `MR-PRESSO results`$
            `Outlier Test`, file = "outlier_SNPs.csv")

mr_result<-mr(dat)

write.csv(generate_odds_ratios(mr_result), file = "MR-Result.csv", row.names = FALSE)
pFilter <- 1  

draw_forest_map <- function(inputFile = null, forestFile = null, forestCol = null) {
  
  rt <- read.csv(inputFile, header = T, sep = ",", check.names = F)
  row.names(rt) <- rt$method
  rt <- rt[rt$pval < pFilter,]
  method <- rownames(rt)
  or <- sprintf("%.3f", rt$"or")
  orLow <- sprintf("%.3f", rt$"or_lci95")
  orHigh <- sprintf("%.3f", rt$"or_uci95")
  OR <- paste0(or, "(", orLow, "-", orHigh, ")")
  pVal <- ifelse(rt$pval < 0.001, "<0.001", sprintf("%.3f", rt$pval))
  
  
  pdf(file = forestFile, width = 7, height = 4.6)
  n <- nrow(rt)
  nRow <- n + 1
  ylim <- c(1, nRow)
  layout(matrix(c(1, 2), nc = 2), width = c(3.5, 2))
  
 
  xlim <- c(0, 3)
  par(mar = c(4, 2.5, 2, 1))
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = F, xlab = "", ylab = "")
  text.cex <- 0.8
  text(0, n:1, method, adj = 0, cex = text.cex)
  text(1.9, n:1, pVal, adj = 1, cex = text.cex); text(1.9, n + 1, 'pvalue', cex = 1, font = 2, adj = 1)
  text(3.1, n:1, OR, adj = 1, cex = text.cex); text(2.7, n + 1, 'OR', cex = 1, font = 2, adj = 1)
  
  par(mar = c(4, 1, 2, 1), mgp = c(2, 0.5, 0))
  xlim <- c(min(as.numeric(orLow) * 0.975, as.numeric(orHigh) * 0.975, 0.9), max(as.numeric(orLow), as.numeric(orHigh)) * 1.025)
  plot(1, xlim = xlim, ylim = ylim, type = "n", axes = F, ylab = "", xaxs = "i", xlab = "OR")
  arrows(as.numeric(orLow), n:1, as.numeric(orHigh), n:1, angle = 90, code = 3, length = 0.05, col = "darkblue", lwd = 3)
  abline(v = 1, col = "black", lty = 2, lwd = 2)
  boxcolor <- ifelse(as.numeric(or) > 1, forestCol, forestCol)
  points(as.numeric(or), n:1, pch = 15, col = boxcolor, cex = 2)
  axis(1)
  dev.off()
}

draw_forest_map(inputFile = file.path(DIR,"MR-Result.csv"), forestFile = "forest_map.pdf", forestCol = "red")

write.csv(mr_heterogeneity(dat), file = "heterogeneity.csv", row.names = FALSE)

write.csv(mr_pleiotropy_test(dat), file = "pleiotropy.csv", row.names = FALSE)

write.csv(directionality_test(dat), file = "directionality.csv", row.names = FALSE)

pdf(file = "pic.scatter_plot.pdf", width = 7.5, height = 7); mr_scatter_plot(mr_result, dat); dev.off()  
res_single <- mr_singlesnp(dat)
pdf(file = "pic.forest.pdf", width = 7, height = 7); mr_forest_plot(res_single); dev.off()  
pdf(file = "pic.funnel_plot.pdf", width = 7, height = 6.5); mr_funnel_plot(singlesnp_results = res_single); dev.off() 
pdf(file = "pic.leaveoneout.pdf", width = 10, height = 7); mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat)); dev.off()  

#mr_report(dat)
