## plot_BD_Beta_together?

## fecet version?

ALL50K_BD_Beta_posterior <- read.csv2("output/ALL50K/ALL50K_BD_Beta.log", sep = "") %>% mutate(DataType="ALL50K")
intron_BD_Beta_posterior <- read.csv2("output/intron/intron_BD_Beta.log", sep = "") %>% mutate(DataType="intron")
CDS_BD_Beta_posterior <- read.csv2("output/CDS/CDS_BD_Beta.log", sep = "") %>% mutate(DataType="CDS")
lnc_RNA_BD_Beta_posterior <- read.csv2("output/lnc_RNA/lnc_RNA_BD_Beta.log", sep = "") %>% mutate(DataType="lnc_RNA")
other_BD_Beta_posterior <- read.csv2("output/other/other_BD_Beta.log", sep = "") %>% mutate(DataType="other")
pseudogene_BD_Beta_posterior <- read.csv2("output/pseudogene/pseudogene_BD_Beta.log", sep = "") %>% mutate(DataType="pseudogene")
UTR_BD_Beta_posterior <- read.csv2("output/UTR/UTR_BD_Beta.log", sep = "") %>% mutate(DataType="UTR")
unannotated_BD_Beta_posterior <- read.csv2("output/unannotated/unannotated_BD_Beta.log", sep = "") %>% mutate(DataType="unannotated")


#combine BD_Beta_data

combined_BD_Beta <- rbind(ALL50K_BD_Beta_posterior, intron_BD_Beta_posterior, CDS_BD_Beta_posterior,
                       lnc_RNA_BD_Beta_posterior, other_BD_Beta_posterior,  pseudogene_BD_Beta_posterior,
                       unannotated_BD_Beta_posterior)

class(combined_BD_Beta)

subset_combined_BD_Beta <- combined_BD_Beta %>% select(birth_rate, death_rate, DataType) 

subset_combined_BD_Beta$birth_rate <- as.numeric(as.character(subset_combined_BD_Beta$birth_rate))

subset_combined_BD_Beta$death_rate <- as.numeric(as.character(subset_combined_BD_Beta$death_rate))

melted_BD_Beta <- melt(subset_combined_BD_Beta)

ggplot(melted_BD_Beta, aes(x=value, fill=variable)) + geom_density(adjust=1.5, alpha=.4) + facet_wrap(.~DataType) + theme(legend.position = c(0.80,0.80))
