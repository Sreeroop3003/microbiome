# CRAN libraries
library(tidyverse) # for data cleaning/wrangling
library(ggplot2) # for visualization

# Bioconductor libraries
library(dada2)
library(phyloseq)
library(Biostrings)
library(DECIPHER)
library(phangorn)
library(gridExtra)

# list amplicons into forward/reverse reads
amplicons_F <- sort(list.files(path = "data/bpd_microbiome", # location of your amplicons
                               pattern = "_1", # forward reads
                               full.names = TRUE))

amplicons_R <- sort(list.files(path = "data/bpd_microbiome",
                               pattern = "_2", # reverse reads
                               full.names = TRUE))
sample_names <- basename(amplicons_F) %>% # get file names
  word(start = 1L, sep = "_") # get the first word as the sample names

# plotting 1st observation
f_plot <- plotQualityProfile(amplicons_F[1:3]) + labs(x = "Sequence Position")
r_plot <- plotQualityProfile(amplicons_R[1:3]) + labs(x = "Sequence Position")

grid.arrange(f_plot,r_plot,ncol = 1)

# creating directory for bpd_filtered reads

if(!file_test("-d", "data/bpd_filtered")) #"if there is no directory `data/bpd_filtered`,
  dir.create("data/bpd_filtered") # "create one"

# creating file path
# "data/bpd_filtered/sample_names_*_bpd_filtered.fastq.gz"

bpd_filtered_F <- file.path("data", "bpd_filtered", paste0(sample_names,"_F_bpd_filtered.fastq.gz"))
bpd_filtered_R <- file.path("data", "bpd_filtered", paste0(sample_names,"_R_bpd_filtered.fastq.gz"))

tnf_summary <- filterAndTrim(amplicons_F, bpd_filtered_F, # input and output
                             amplicons_R, bpd_filtered_R,
                             
                             # trimming
                             trimLeft=10, # trim the first n observation from each reads
                             truncLen=c(240,160), # truncate reads after this position; c(Forward/Reverse)
                             
                             # filtering standard
                             maxN=0, maxEE=c(2,2), # max expected error (maxEE) = 2
                             truncQ=2, rm.phix=TRUE,
                             
                             # additional setting
                             compress=TRUE, # whether outputs should be compressed
                             multithread=FALSE) # default for Windows, Mac can use `multithread=TRUE`
# error model for forward reads
error_F <- learnErrors(bpd_filtered_F) #input: file path for bpd_filtered reads
# error model for reverse reads
error_R <- learnErrors(bpd_filtered_R)
# infer sequence variants
dada_F <- dada(bpd_filtered_F, err = error_F, verbose = FALSE) 
dada_R <- dada(bpd_filtered_R, err = error_R, verbose = FALSE)

merged <- mergePairs(dadaF = dada_F, # dada result
                     derepF = bpd_filtered_F, # path of bpd_filtered reads
                     dadaR = dada_R,
                     derepR = bpd_filtered_R)


seqtab <- makeSequenceTable(merged)


seqtab[1:3,1:3]


rownames(seqtab) <- sample_names
dim(seqtab)


seqtab_nochim <- removeBimeraDenovo(seqtab,
                                    method = "consensus")
dim(seqtab_nochim)

getN <- function(x) {
  sum(getUniques(x))
}

# generate the summary
track_reads <- data.frame(row.names = sample_names,
                          raw_reads = tnf_summary[,1],
                          bpd_filtered = tnf_summary[,2],
                          ASVs_F = sapply(dada_F, getN), # sapply to apply the function to all rows (sample)
                          ASVs_R = sapply(dada_F, getN),
                          joined = sapply(merged, getN),
                          no_chimera = rowSums(seqtab_nochim)) # form row-column sums

taxa <- assignTaxonomy(seqtab_nochim, # sequence table
                       refFasta = "data/bpd_filtered/silva_nr_v138_train_set.fa.gz", # reference sequence
                       minBoot = 50) # minimal bootstrap confidence; default to 50

# add taxa until species level
# note that it may result NA if there is no exact match in the reference
taxa <- addSpecies(taxa, # result from assignTaxonomy
                   refFasta = "data/bpd_filtered/silva_species_assignment_v138.fa.gz")

# check taxonomic assignment
taxa_print <- taxa
rownames(taxa_print) <- NULL # for visualization purpose

head(taxa_print)  

# obtain sequence
seqs <- getSequences(seqtab_nochim)
names(seqs) <- seqs # This propagates to the tip labels of the tree

# do msa 
alignment <- AlignSeqs(DNAStringSet(seqs), anchor=NA)

phang.align <- phyDat(as(alignment, "matrix"), type="DNA") # converting into phyDat format
dm <- dist.ml(phang.align) # compute pairwaise distances
treeNJ <- NJ(dm) # construct neighbor-joining tree
fit <- pml(treeNJ, data = phang.align) # computes likelihood of phylogenetic tree

fitGTR <- update(fit, k=4, inv=0.2) # re-fit a model
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))
fitGTR
gender <- substr(sample_names, start = 1, stop = 1) # get the first character
gender

subject <- substr(sample_names, 2, 2) # get the second character
subject

day <- sapply(strsplit(sample_names, "D"), `[`, 2) %>% # separate by 'D'; get the 2nd value [2]
  as.integer() # convert to integer
day

host_age <- c("62","53","28","32","41","58","65","61","24","28","48","60","52","54","73","69","65","51","38","38","59","60","35","55","56","58","38","61","31","49","70","54","60","62","49","58","81","87","41","25","36","34","44","30","50","68","43","40","26","30","65","60","30","31","34","40","41","50","62","53","53","56","26","55","50","44","53","39","56","40","57","59","32","49","59","38","37","28","55","62","36","32","58","53","59","47","71","57","56","67","72","63","68","54","69","61","73","63","60","27","62","57","32","65","55","54","53","58","68","69","59","56","50","36","72","74","73","23","27","54","53","37","27","70","29","53","55","29","53","48","56","29","31","52","42","29","58")
host_age
sex <- c("Male","Female","Female","Female","Female","Female","Male","Female","Female","Female","Female","Male","Female","Female","Female","Female","Female","Female","Male","Male","Male","Female","Female","Female","Female","Female","Female","Female","Male","Male","Female","Male","Female","Female","Female","Female","Female","Male","Female","Female","Female","Female","Female","Female","Female","Male","Female","Female","Female","Female","Female","Male","Female","Male","Male","Female","Female","Male","Female","Female","Female","Female","Female","Male","Male","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Female","Male","Male","Male","Male","Female","Male","Female","Female","Male","Female","Female","Female","Female","Male","Female","Male","Female","Female","Male","Male","Male","Female","Female","Male","Male","Female","Female","Female","Female","Female","Female","Male","Male","Female","Female","Female","Female","Female","Female","Male","Female","Male","Female","Female","Female","Male","Female","Female","Female","Female","Female","Female","Female","Female","Male","Male","Male","Female","Male","Female","Female")
sex
bmi <- c("24.69","24.96","19.84","21.91","29.83","38.89","27.19","20.67","15.95","23.4","20.37","27.95","21.3","36.9","26.04","23.38","25.37","22.59","28.24","25.06","25.42","28.71","33.52","25.79","20.63","31.32","28.34","20.35","22.24","28.3","34.75","34.82","30.78","32","23.03","19.72","27.61","22.64","26.57","31.64","21.92","19.52","34.33","35.71","22.24","30.27","32.43","21.95","25.01","17.18","28.19","25.84","21.14","31.47","24.9","34.01","28.34","44.72","36.27","24.8","42.57","29.29","24.69","32.28","31.32","21.6","23.92","40.68","29.08","39.45","42.28","31.47","20.59","20.51","25.01","31.24","24.39","23.49","27.45","24.65","42.91","38.03","29.53","18.51","31.01","25.79","29.05","26.58","19.87","28.82","35.87","22.14","37.49","54.08","25.11","36.41","24.13","21.25","22.13","21.22","28.5","22.59","29.05","38.69","30.66","31.64","31.55","25.8","28.75","29.05","23.66","43.85","33.14","34.77","25.39","32.48","22.65","25.8","24.03","23.43","38.41","25.84","23.4","25.12","46.6","24.21","39.1","25.1","27.89","36.61","31.26","20.34","22.04","29.23","24","22.14","27.44")
bmi
class <- c("Healthy","BPD","Healthy","Healthy","BPD","BPD","Healthy","BPD","BPD","Healthy","BPD","BPD","BPD","BPD","Healthy","Healthy","BPD","Healthy","BPD","BPD","BPD","Healthy","BPD","BPD","BPD","Healthy","Healthy","Healthy","Healthy","BPD","Healthy","Healthy","BPD","Healthy","BPD","BPD","Healthy","BPD","BPD","Healthy","BPD","BPD","BPD","BPD","BPD","BPD","BPD","BPD","Healthy","Healthy","Healthy","BPD","BPD","BPD","Healthy","BPD","BPD","BPD","Healthy","BPD","BPD","BPD","Healthy","BPD","BPD","Healthy","Healthy","BPD","BPD","BPD","BPD","BPD","Healthy","BPD","BPD","BPD","Healthy","Healthy","BPD","BPD","BPD","BPD","BPD","Healthy","Healthy","Healthy","BPD","BPD","BPD","Healthy","BPD","Healthy","BPD","BPD","Healthy","BPD","Healthy","BPD","BPD","BPD","BPD","BPD","Healthy","BPD","Healthy","BPD","BPD","BPD","Healthy","BPD","BPD","BPD","BPD","BPD","BPD","BPD","Healthy","BPD","Healthy","BPD","BPD","Healthy","Healthy","BPD","BPD","Healthy","BPD","BPD","BPD","BPD","BPD","Healthy","Healthy","BPD","BPD","BPD","BPD")
class
length(bmi)
# combine metadata
seq_data <- data.frame(Age = host_age, 
                       Gender = sex, 
                       BMI = bmi,
                       Class = class,
                       row.names = sample_names)
head(seq_data)
sample_names
# create phylosec object
ps <- phyloseq(otu_table(seqtab_nochim, taxa_are_rows=FALSE),
               sample_data(seq_data),
               tax_table(taxa),
               phy_tree(fitGTR$tree))

# remove mock sample
ps <- prune_samples(sample_names(ps) != "Mock", ps)

ps

# get sequences
dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps) # connect dna data to taxa names

dna


ps <- merge_phyloseq(ps, dna)

# change taxa names into shorter id (ASVn)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))

otu_table(ps)[,1:3]

refseq(ps)[1:3]


ps
getSequences(otu_table(ps))[1:5]

#############################################################################################################

rank_names(ps)
## [1] "Kingdom" "Phylum"  "Class"   "Order"   "Family"  "Genus"   "Species"

# rows - column: 234 ASVs - 7 Taxonomic Rank 
dim(tax_table(ps))
## [1] 234   7

table(tax_table(ps)[,"Phylum"])
## 
##  Actinobacteriota      Bacteroidota  Campilobacterota     Cyanobacteria 
##                 7                20                 1                 3 
##      Deinococcota        Firmicutes   Patescibacteria    Proteobacteria 
##                 1               193                 1                 7 
## Verrucomicrobiota 
##                 1

# available genus: 49 (not including NA) 
get_taxa_unique(ps, taxonomic.rank = "Genus")
genus <- tax_table(ps)[,'Genus']
length(genus)
head <- ps_agg <- tax_glom(ps, "Genus", NArm = TRUE) 

# original data
ps

ps_tree <-  plot_tree(ps, method = "treeonly",
                      ladderize = "left",
                      title = "Before Agglomeration")

ps_agg_tree <-  plot_tree(ps_agg, method = "treeonly",
                          ladderize = "left",
                          title = "After Agglomeration")
ps_tree
ps_agg_tree

#############################################################################################################


ps_log <- transform_sample_counts(ps_agg, function(x) log(1 + x))
# making ordinate
out_wuf_log <- ordinate(ps_log, 
                        method = "MDS", # for PCoA 
                        distance = "wunifrac") # weighted Unifrac distance
# prepare eigen values to adjust axis
evals <- out_wuf_log$values$Eigenvalues

plot_pcoa<- plot_ordination(ps_log, out_wuf_log, color = "Gender") +
  geom_text(aes(label = sample_names(ps_log)), size = 3, nudge_y = 0.02) +
  labs(col = "Gender") +
  # to adjust axis length based on eigen values (variance it contains)
  coord_fixed(sqrt(evals[2] / evals[1])) 
plot_pcoa

dat <- otu_table(ps)
typeof(dat)
dat
OTUdf <- as.data.frame(otu_table(ps))
#OTUdf_tax <- setNames(OTUdf, c(genus))
write.csv(OTUdf,"OTU_Table_Unnamed.csv")

#plot_richness(ps.rarefied, x="body.site", measures=c("Observed", "Shannon")) +
#  geom_boxplot() +
#  theme_classic() +
#  theme(strip.background = element_blank(), axis.text.x.bottom = element_text(angle = -90))

#ps.rarefied = rarefy_even_depth(ps, rngseed=1, sample.size=1103, replace=F)
ps_cut <- prune_samples(sample_names(ps_agg) != "F3D142", ps_agg)


ps_relav <- transform_sample_counts(ps_cut, 
                                    function(x){x / sum(x)})

# before transformation
ps_cut@otu_table[1:6, 1:3]

ps_relav@otu_table[1:6, 1:3]
ps_relav@otu_table
# convert to data frame for easier access
tax_table <- as.data.frame(ps_relav@tax_table@.Data)
tax_table

tax_table_df <- as.data.frame(tax_table)
#OTUdf_tax <- setNames(OTUdf, c(genus))
write.csv(tax_table_df,"tax_table.csv")

# phylum
unique(tax_table$Phylum)

ps_relav@sam_data

# DIY function
plot_abundance <- function(x = physeq, # phyloseq data
                           title = "",
                           Facet = "Phylum", # taxa rank for facets
                           Category = "class", # categorical features for x axis
                           Color = "Phylum",
                           legend = "none"
) {
  
  mphyseq <- psmelt(x)
  mphyseq <- subset(mphyseq, Abundance > 0)
  
  ggplot(data = mphyseq, 
         mapping = aes_string(x = Category,
                              y = "Abundance",
                              color = Color, fill = Color)
  ) +
    geom_violin(fill = NA) +
    geom_point(size = 1, alpha = 0.3, 
               position = position_jitter(width = 0.3)) +
    facet_wrap(facets = Facet, ncol = 3) + 
    scale_y_log10() +
    labs(title = title) +
    theme(legend.position = legend)
}

plot
#save(ps, file = "C:/Users/george-pc/Documents/ps_george.RData")

# plotting abundance
plot <- plot_abundance(ps_relav, 
                       "Microbial Abundance on All Phylum")
plot

# subset taxa for Phylum "Firmicutes"
ps_firm <-  subset_taxa(ps_relav, Phylum == "Firmicutes")

plot_ordo <- plot_abundance(ps_firm,
                            title = "Microbial Abundance on Firmicutes",
                            Facet = "Order")
plot_ordo

plot2 <- plot_bar(ps_relav, # using all data (relative abundance)
                  fill = "Phylum", # fill colour by Phylum.
                  x = "reorder(Sample, Day)") + 
  labs(title = "Microbial Abundance of Murine Gut",
       subtitle = "Phylum Distribution",
       x = NULL) # to remove x-axis title
plot2

plot3 <- plot_bar(ps_firm, # using only firmicutes data
                  fill = "Order", # fill colour by order.
                  x = "reorder(Sample, Day)") + 
  labs(title = "Microbial Abundance of Murine Gut",
       subtitle = "Phylum Firmicutes",
       x = NULL) # to remove x-axis title

plot3

#############################################################################################################


# load library
library(caret)
library(randomForest)

# additional library for data tidying later
library(tidyr)
ps_cut@otu_table[1:3,1:3]
ps_cut_log <- transform_sample_counts(ps_cut, function(x) log(1 + x))

ps_cut_log@otu_table[1:3,1:3]  
# get age label & microbiome data
data_ml <- data.frame(age = sample_data(ps_cut_log)$When, otu_table(ps_cut_log))  
head(data_ml)


set.seed(100)
idx_train <- sample(nrow(sample_data(ps_cut_log)), size = nrow(sample_data(ps_cut_log))*0.8)

# splitting train-test
training <- data_ml[idx_train,]
testing <- data_ml[-idx_train,]

rfFit <- train(age ~ ., 
               data = training,
               method = "rf")

rfFit

rfClasses <- predict(rfFit, 
                     newdata = testing)

# other options use: `confusionMatrix()`
table(pred = rfClasses, 
      actual = testing$age)
importance(rfFit$finalModel)

# obtain taxonomy of all microbial sample
tax_rf <- tax_table(ps_cut_log)
# filter microbes with the highest importance
tax_rf[which.max(importance(rfFit$finalModel)),]
## Taxonomy Table:     [1 taxa by 7 taxonomic ranks]:
##       Kingdom    Phylum       Class        Order             Family            
## ASV27 "Bacteria" "Firmicutes" "Clostridia" "Oscillospirales" "Oscillospiraceae"
##       Genus           Species
## ASV27 "Oscillibacter" NA

# get microbial abundance
imp_abd <- as.vector(otu_table(ps_cut_log)[,"ASV27"])

# combine with sample data
imp_df <- data.frame(sample_data(ps_cut_log),
                     abund = imp_abd)

# plotting
imp_plot<- ggplot(imp_df, aes(x = abund)) + 
  geom_density(aes(fill = When),
               alpha = 0.5) +
  labs(title = "Abundance of Discriminative Species",
       subtitle = "Oscillibacter sp.",
       x = "Abundance",
       y = "Density of samples",
       fill = "Age") +
  # below is for aesthetics
  theme_minimal()
imp_plot


# obtain the importance of microbes on rf model
rf_imp <- data.frame(importance(rfFit$finalModel))
# obtain 5 microbes with highest importance
rf_imp %>% 
  mutate(ASV = rownames(.)) %>% 
  arrange(desc(MeanDecreaseGini)) %>% 
  head(5)
##   MeanDecreaseGini   ASV
## 1        0.4050516 ASV27
## 2        0.3336478 ASV16
## 3        0.3068830 ASV74
## 4        0.3012801 ASV68
## 5        0.2809128 ASV66
# saving ASV to object
asv_imp <- rf_imp %>% 
  mutate(ASV = rownames(.)) %>% 
  arrange(desc(MeanDecreaseGini)) %>% 
  head(5) %>% 
  pull(ASV)
# obtain taxonomy for those ASVs
data.frame(tax_rf) %>% 
  filter(rownames(.) %in% asv_imp)
##        Kingdom           Phylum          Class             Order
## ASV16 Bacteria       Firmicutes     Clostridia    Lachnospirales
## ASV27 Bacteria       Firmicutes     Clostridia   Oscillospirales
## ASV66 Bacteria Actinobacteriota Actinobacteria Bifidobacteriales
## ASV68 Bacteria       Firmicutes     Clostridia    Lachnospirales
## ASV74 Bacteria       Firmicutes     Clostridia     Clostridiales
##                   Family                         Genus Species
## ASV16    Lachnospiraceae Lachnospiraceae_NK4A136_group    <NA>
## ASV27   Oscillospiraceae                 Oscillibacter    <NA>
## ASV66 Bifidobacteriaceae               Bifidobacterium    <NA>
## ASV68    Lachnospiraceae                 Acetatifactor    <NA>
## ASV74     Clostridiaceae   Clostridium_sensu_stricto_1    <NA>

# get microbial abundance
imp_abd_asv <- data.frame(otu_table(ps_cut_log)[,asv_imp])

# combine with sample data
imp_df_asv <- data.frame(sample_data(ps_cut_log)) %>% 
  cbind(imp_abd_asv) %>% 
  pivot_longer(cols = asv_imp,
               names_to = "ASV", values_to = "Abundance")
head(imp_df_asv, 10)
## # A tibble: 10 x 6
##    Subject Gender   Day When  ASV   Abundance
##    <fct>   <fct>  <int> <fct> <chr>     <dbl>
##  1 3       F          0 Early ASV27      4.76
##  2 3       F          0 Early ASV16      6.80
##  3 3       F          0 Early ASV74      4.86
##  4 3       F          0 Early ASV68      3.40
##  5 3       F          0 Early ASV66      3.22
##  6 3       F          1 Early ASV27      5.08
##  7 3       F          1 Early ASV16      6.86
##  8 3       F          1 Early ASV74      0   
##  9 3       F          1 Early ASV68      4.50
## 10 3       F          1 Early ASV66      0

# plotting
imp_plot_asv <- ggplot(imp_df_asv, aes(x = Abundance)) + 
  geom_density(aes(fill = When),
               alpha = 0.5) +
  facet_wrap(~ASV) +
  labs(title = "Abundance of Discriminative Species",
       # subtitle = "Oscillibacter sp.",
       x = "Abundance",
       y = "Density of samples",
       fill = "Age") +
  # below is for aesthetics
  theme_minimal()
imp_plot_asv  



#----------------------------------------------------------------------------------------


ps
head(sample_data(ps), 25)
library("phyloseq")
head(sample_data(ps)$Gender, 25)
library("DESeq2")
packageVersion("DESeq2")
coldata <- read.table("C:/Users/karth/Documents/colData.txt", header = T, sep= '\t')
head(coldata) 
diagdds = phyloseq_to_deseq2(ps, colData(coldata), ~ class)
diaGgdds = DESeq(diagdds, test="Wald", fitype="parametric")



seq_data

