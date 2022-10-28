# Install RevGadget if you haven't done so already
#library(devtools)
#install_github("revbayes/RevGadgets")

library(RevGadgets)
library(ggplot2)
library(dplyr)
library(data.table)

# specify the output files
ALL50K_speciation_time_file <- "output/ALL50K/ALL50K_EBD_speciation_times.log"
ALL50K_speciation_rate_file <- "output/ALL50K/ALL50K_EBD_speciation_rates.log"
ALL50K_extinction_time_file <- "output/ALL50K/ALL50K_EBD_extinction_times.log"
ALL50K_extinction_rate_file <- "output/ALL50K/ALL50K_EBD_extinction_rates.log"


intron_speciation_time_file <- "output/intron/intron_EBD_speciation_times.log"
intron_speciation_rate_file <- "output/intron/intron_EBD_speciation_rates.log"
intron_extinction_time_file <- "output/intron/intron_EBD_extinction_times.log"
intron_extinction_rate_file <- "output/intron/intron_EBD_extinction_rates.log"


CDS_speciation_time_file <- "output/CDS/CDS_EBD_speciation_times.log"
CDS_speciation_rate_file <- "output/CDS/CDS_EBD_speciation_rates.log"
CDS_extinction_time_file <- "output/CDS/CDS_EBD_extinction_times.log"
CDS_extinction_rate_file <- "output/CDS/CDS_EBD_extinction_rates.log"


lnc_RNA_speciation_time_file <- "output/lnc_RNA/lnc_RNA_EBD_speciation_times.log"
lnc_RNA_speciation_rate_file <- "output/lnc_RNA/lnc_RNA_EBD_speciation_rates.log"
lnc_RNA_extinction_time_file <- "output/lnc_RNA/lnc_RNA_EBD_extinction_times.log"
lnc_RNA_extinction_rate_file <- "output/lnc_RNA/lnc_RNA_EBD_extinction_rates.log"


other_speciation_time_file <- "output/other/other_EBD_speciation_times.log"
other_speciation_rate_file <- "output/other/other_EBD_speciation_rates.log"
other_extinction_time_file <- "output/other/other_EBD_extinction_times.log"
other_extinction_rate_file <- "output/other/other_EBD_extinction_rates.log"


UTR_speciation_time_file <- "output/UTR/UTR_EBD_speciation_times.log"
UTR_speciation_rate_file <- "output/UTR/UTR_EBD_speciation_rates.log"
UTR_extinction_time_file <- "output/UTR/UTR_EBD_extinction_times.log"
UTR_extinction_rate_file <- "output/UTR/UTR_EBD_extinction_rates.log"


pseudogene_speciation_time_file <- "output/pseudogene/pseudogene_EBD_speciation_times.log"
pseudogene_speciation_rate_file <- "output/pseudogene/pseudogene_EBD_speciation_rates.log"
pseudogene_extinction_time_file <- "output/pseudogene/pseudogene_EBD_extinction_times.log"
pseudogene_extinction_rate_file <- "output/pseudogene/pseudogene_EBD_extinction_rates.log"




# read in and process rates
ALL50K_rates <- processDivRates(speciation_time_log = ALL50K_speciation_time_file,
                         speciation_rate_log = ALL50K_speciation_rate_file,
                         extinction_time_log = ALL50K_extinction_time_file,
                         extinction_rate_log = ALL50K_extinction_rate_file,
                         burnin = 0.25,
                         summary = "median")  %>% mutate(DataType="ALL50K")

# read in and process rates
intron_rates <- processDivRates(speciation_time_log = intron_speciation_time_file,
                         speciation_rate_log = intron_speciation_rate_file,
                         extinction_time_log = intron_extinction_time_file,
                         extinction_rate_log = intron_extinction_rate_file,
                         burnin = 0.25,
                         summary = "median") %>% mutate(DataType="intron")

# read in and process rates
CDS_rates <- processDivRates(speciation_time_log = CDS_speciation_time_file,
                         speciation_rate_log = CDS_speciation_rate_file,
                         extinction_time_log = CDS_extinction_time_file,
                         extinction_rate_log = CDS_extinction_rate_file,
                         burnin = 0.25,
                         summary = "median") %>% mutate(DataType="CDS")

# read in and process rates
pseudogene_rates <- processDivRates(speciation_time_log = pseudogene_speciation_time_file,
                         speciation_rate_log = pseudogene_speciation_rate_file,
                         extinction_time_log = pseudogene_extinction_time_file,
                         extinction_rate_log = pseudogene_extinction_rate_file,
                         burnin = 0.25,
                         summary = "median") %>% mutate(DataType="pseudogene")

# read in and process rates
lnc_RNA_rates <- processDivRates(speciation_time_log = lnc_RNA_speciation_time_file,
                         speciation_rate_log = lnc_RNA_speciation_rate_file,
                         extinction_time_log = lnc_RNA_extinction_time_file,
                         extinction_rate_log = lnc_RNA_extinction_rate_file,
                         burnin = 0.25,
                         summary = "median") %>% mutate(DataType="lnc_RNA")

# read in and process rates
other_rates <- processDivRates(speciation_time_log = other_speciation_time_file,
                         speciation_rate_log = other_speciation_rate_file,
                         extinction_time_log = other_extinction_time_file,
                         extinction_rate_log = other_extinction_rate_file,
                         burnin = 0.25,
                         summary = "median") %>% mutate(DataType="other")

# read in and process rates
UTR_rates <- processDivRates(speciation_time_log = UTR_speciation_time_file,
                               speciation_rate_log = UTR_speciation_rate_file,
                               extinction_time_log = UTR_extinction_time_file,
                               extinction_rate_log = UTR_extinction_rate_file,
                               burnin = 0.25,
                               summary = "median") %>% mutate(DataType="UTR")


# plot rates through time
p <- plotDivRates(ALL50K_rates) +
        xlab("Millions of years ago") +
        ylab("Rate per million years")

p

ggsave("pseudogene_EBD.png", p)

mean(extinction_rate_file)
