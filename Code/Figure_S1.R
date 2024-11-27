#Code to replicate Figure S1 results
#updated Friday, November 1st, 2024

#setup----

library(tidyr)            #1.3.0
library(dplyr)            #1.1.4
library(ggplot2)          #3.4.4
library(ggpubr)           #0.6.0
library(vegan)            #2.6.4
library(phyloseq)         #1.46.0
library(ape)              #5.7.1
library(ggExtra)          #0.10.1
library(Maaslin2)         #1.8.0
library(NatParksPalettes) #0.2.0
library(reshape2)         #1.4.4
library(rstatix)          #0.7.2
library(scales)           #1.3.0

##set working directory
setwd('/Users/annadeveaux/Box Sync/Functional consequences of microbiome variation in CDI susceptibility/')

##set seed
set.seed(42)

##define %ni%
`%ni%` <- Negate(`%in%`)

##define plot function
theme_plot <- function() {
  theme_classic() +
    theme(legend.position = "right",
          axis.title.x = element_blank(),
          axis.title.y = element_text(size=15),
          axis.text.y = element_text(size=15),
          axis.text.x = element_text(size=12, angle=45, hjust=1),
          strip.text.x = element_text(size=15),
          panel.border = element_rect(fill=NA, colour="black", size=0.5))
}

##load metadata
h.meta <- read.delim('Metadata/241031_Human_Metadata.csv',
                     sep = ',')
row.names(h.meta) <- h.meta$Sample

##load antibiotic exposure data
abx <- read.csv('Metadata/240122_Antibiotic_Exposures.csv')
rownames(abx) <- abx$Sample
abx$Sample <- NULL
abx[is.na(abx)] <- 0

##load phyloseq object with human taxonomic data (created in Figure 1 code)
h.ps <- readRDS('ManuscriptAnalyses/Human/MetaPhlAn4/ps_human_m4_final.rds')

##define plot colors
pal <- natparks.pals('Triglav')[c(3, 5, 4, 1)]



#taxonomic alpha diversity (Fig S1A)----

##prep data----
##create wide species data frame
h.ps.spp <- tax_glom(h.ps, taxrank = 'Species')
h.ps.spp.df <- psmelt(h.ps.spp)
h.ps.spp.df.wide <- dcast(h.ps.spp.df,
                          Sample ~ Species,
                          value.var = "Abundance")
row.names(h.ps.spp.df.wide) <- h.ps.spp.df.wide$Sample
h.ps.spp.df.wide$Sample <- NULL

##compute metrics----
##convert relative abundance measurements for FILTERED data to integers, preserving precision to 0.00001%
##ignore warning about singletons - this applies to 16S
h.ps.int <- transform_sample_counts(h.ps, function(x) trunc(x * 100000))
adiv <- estimate_richness(h.ps.int, measures = c('Observed', 'Shannon'))

##append alpha diversity metrics to phyloseq object
sample_data(h.ps)$Richness <- adiv$Observed
sample_data(h.ps)$Shannon <- adiv$Shannon

##gather metadata to long format
h.ps.df <- data.frame(sample_data(h.ps))
h.ps.df.a <- gather(h.ps.df, Alpha_Measure, Value, Richness:Shannon, factor_key = TRUE)
h.ps.df.a$Group <- factor(h.ps.df.a$Group,
                          levels = c("CDI", "Carrier", "IC", "HV"))

##plot shannon diversity----
##check whether variables are normally distributed
shapiro.test(filter(h.ps.df.a, Alpha_Measure == "Richness")$Value) #p = 0.0015, nonnormal
shapiro.test(filter(h.ps.df.a, Alpha_Measure == "Shannon")$Value) #p = 1.1e-05, nonnormal

##compute adjusted p-values between groups to add to plots
stat.test2<- h.ps.df.a %>%
  filter(Primary_Sample == 1, Alpha_Measure == "Shannon") %>% 
  wilcox_test(Value ~ Group) %>% 
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance() %>%
  add_y_position() %>% 
  mutate(Alpha_Measure = 'Shannon')
stat.test2$y.position <- c(4.5,4.5,4.5,4.5,4.5,4.5)

p.spec.shannon <- h.ps.df.a %>% 
  filter(Primary_Sample == 1, Alpha_Measure == 'Shannon') %>% 
  ggplot(aes(x = Group, y = Value)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Group), alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.7, color = 'black') +
  scale_fill_manual(values = pal) +
  theme_plot() +
  theme(legend.position = 'none',
        aspect.ratio=1.5) +
  ylab('Shannon diversity') +
  stat_pvalue_manual(stat.test2, label = 'p.adj', hide.ns = T, step.increase = 0.1) +
  scale_y_continuous(limits = c(0, NA)) +
  scale_x_discrete(labels = c('CDI', 'Carrier', 'IC', 'HV'))



#taxonomic beta diversity (Fig S1B)----
#this will plot bray-curtis dissimilarity at the species level. primary samples are used, except
#for the case where a single patient was used for two donor samples in the engraftment screen, in
#which case both are plotted.

##prep species data frame
h.ps.spp.df.wide2 <- h.ps.spp.df.wide

##compute dissimilarity metrics
h.ps.spp.bray.all <- vegdist(filter(h.ps.spp.df.wide2,
                                    row.names(h.ps.spp.df.wide2) %in% (filter(h.meta,
                                                                              Primary_Sample != 0)$Sample)),
                             method = 'bray')

##compute pcoa
h.ps.spp.pcoa.all <- pcoa(h.ps.spp.bray.all)

##extract PCoA1 and PCoA2
h.spp.all.pcoa1 <- 100 * (h.ps.spp.pcoa.all$values$Eigenvalues[1]/sum(h.ps.spp.pcoa.all$values$Eigenvalues))
h.spp.all.pcoa2 <- 100 * (h.ps.spp.pcoa.all$values$Eigenvalues[2]/sum(h.ps.spp.pcoa.all$values$Eigenvalues))

##prep plot----
h.ps.spp.pcoa.plot.all <- data.frame(h.ps.spp.pcoa.all$vectors[,])
h.ps.spp.meta.all <- filter(h.ps.spp.df.wide2,
                            row.names(h.ps.spp.df.wide2) %in% (filter(h.meta,
                                                                      Primary_Sample != 0)$Sample))
h.ps.spp.meta.all$Sample <- row.names(h.ps.spp.meta.all)
h.ps.spp.meta.all <- left_join(h.ps.spp.meta.all, h.meta, by = 'Sample')
h.ps.spp.meta.all$Axis.1 <- h.ps.spp.pcoa.plot.all$Axis.1
h.ps.spp.meta.all$Axis.2 <- h.ps.spp.pcoa.plot.all$Axis.2
row.names(h.ps.spp.meta.all) <- h.ps.spp.meta.all$Sample
h.ps.spp.meta.all$Group <- factor(h.ps.spp.meta.all$Group,
                                  levels = c('CDI', 'Carrier', 'IC', 'HV'))

##compute adjusted p-values between groups to add to plots
stat.test.marginal1 <- h.ps.spp.meta.all %>%
  wilcox_test(Axis.1 ~ Group) %>% 
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance() %>%
  add_y_position()

stat.test.marginal2<- h.ps.spp.meta.all %>%
  wilcox_test(Axis.2 ~ Group) %>% 
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance() %>%
  add_y_position()

stat.test.marginal3 <- rbind(stat.test.marginal1, stat.test.marginal2)

##PERMANOVA----
##between all 4 groups
adonis2(h.ps.spp.bray.all ~ Group, data = h.ps.spp.meta.all,
        permutations=9999, method="bray") # p = 1e-04

##between HV and Hospitalized groups
adonis2(h.ps.spp.bray.all ~ Hospitalized, data = h.ps.spp.meta.all,
        permutations = 9999, method = 'bray') # p = 1e-04

##plot----
h.p.pcoa.spp <- ggMarginal(h.ps.spp.meta.all %>% 
                             ggplot(aes(x = Axis.1, y = Axis.2, color = Group)) +
                             stat_ellipse(geom = 'polygon', alpha = 0.02, lwd = 1.2) +
                             geom_point(size = 1.5, alpha = 0.8) +
                             scale_color_manual(values = pal) +
                             theme_plot() +
                             theme(legend.position = 'none',
                                   axis.title.x = element_text(size=15),
                                   axis.text.x = element_text(size=15, angle=45, hjust=1),
                                   aspect.ratio = 0.8) +
                             labs(x = paste("PCoA 1 (", round(h.spp.all.pcoa1, 1), '%)', sep = ''),
                                  y = paste("PCoA 2 (", round(h.spp.all.pcoa2, 1), '%)', sep = '')),
                           type = 'boxplot', groupFill = T, alpha = 1)



#antimicrobial resistance gene (ARG) analysis (Fig S1C-D)----

#look at overall ARG burden between groups
##load ShortBRED output----
arg.full <- read.delim('ManuscriptAnalyses/Human/shortBRED/241101_shortBRED_out_final.txt')
arg.wide <- dcast(arg.full, Sample ~ Family, value.var = 'Count')
row.names(arg.wide) <- arg.wide$Sample
arg.wide <- arg.wide[, c(-1)]

##sum counts to get RPKM per sample
arg.wide$RPKM <- rowSums(arg.wide)

##get ARG richness
arg.wide$Richness <- specnumber(arg.wide[, 1:ncol(arg.wide)-1])

##add columns to metadata frame
arg.wide.meta <- select(arg.wide, RPKM, Richness)
arg.wide.meta$Sample <- row.names(arg.wide.meta)
arg.wide.meta <- left_join(filter(h.meta,
                                  Primary_Sample == 1), arg.wide.meta, by = 'Sample')
arg.wide.meta$Group <- factor(arg.wide.meta$Group,
                              levels = c('CDI', 'Carrier', 'IC', 'HV'))

##check whether variables are normally distributed
shapiro.test(arg.wide.meta$RPKM) #p = <2.2-16, nonnormal
shapiro.test(arg.wide.meta$Richness) #p = 4.42e-08, nonnormal

##compute adjusted p-values between groups to add to plots
stat.test4 <- arg.wide.meta %>%
  wilcox_test(RPKM ~ Group) %>% 
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance() %>%
  add_y_position()
stat.test4$y.position <- 5

stat.test5<- arg.wide.meta %>%
  wilcox_test(Richness ~ Group) %>% 
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance() %>%
  add_y_position()
stat.test5$y.position <- 275

##plot richness (S1C)----
p.arg.richness <- arg.wide.meta %>% 
  ggplot(aes(x = Group, y = Richness)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Group), alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.7, color = 'black') +
  scale_fill_manual(values = pal) +
  theme_plot() +
  theme(legend.position = 'none',
        aspect.ratio = 1.5) +
  ylab('Unique ARGs') +
  stat_pvalue_manual(stat.test5, label = 'p.adj', hide.ns = T, step.increase = 0.12) +
  scale_x_discrete(labels = c('CDI', 'Carrier', 'IC', 'HV'))

##plot RPKM (S1D)----
p.arg.rpkm <- arg.wide.meta %>% 
  ggplot(aes(x = Group, y = RPKM)) +
  geom_boxplot(outlier.shape = NA, aes(fill = Group), alpha = 0.8) +
  geom_jitter(width = 0.15, size = 0.7, color = 'black') +
  scale_fill_manual(values = pal) +
  theme_plot() +
  theme(legend.position = 'none',
        aspect.ratio = 1.5) +
  ylab('ARGs (RPKM)') +
  stat_pvalue_manual(stat.test4, label = 'p.adj', hide.ns = T, step.increase = 0.075) +
  scale_y_continuous(trans = 'log10', labels = label_comma()) +
  scale_x_discrete(labels = c('CDI', 'carrier', 'IC', 'HV'))



#cdtB positivity in isolates (Fig S1E)

##plot----
p.cdtB <- h.meta %>% 
  filter(Colonized == 1 & Primary_Sample == 1) %>%
  filter(is.na(CdtB) == FALSE) %>% 
  mutate(CdtB_presence = factor(if_else(CdtB == 1, 'cdtB+', 'cdtB-'),
                                levels = c('cdtB+', 'cdtB-'))) %>% 
  ggplot(aes(x = Group, fill = CdtB_presence)) +
  geom_bar(position = "stack") +
  labs(x = NULL, fill = "cdtB Presence") +
  theme_plot() +
  theme(legend.position = 'right',
        aspect.ratio = 0.4) +
  coord_flip()



#tcdB type distribution (Fig S1F-G)----

tcdb.map <- read.delim('Scripts/Anna/files_for_github/data/241125_tcdB_subtype_map.txt')

meta.tcdb <- tcdb.map %>% 
  left_join(h.meta) %>% 
  filter(is.na(tcdB_subtype_grouped) == FALSE)

##prep data
h.ps.spp.df.wide3 <- h.ps.spp.df.wide
h.ps.spp.df.wide3$Sample <- row.names(h.ps.spp.df.wide3)
h.ps.spp.df.wide3 <- h.ps.spp.df.wide3 %>%
  filter(Sample %in% meta.tcdb$Sample) %>% 
  filter(Sample %in% (filter(h.meta, Primary_Sample == TRUE))$Sample)
h.ps.spp.df.wide3$Sample <- NULL

##compute dissimilarity metrics
h.ps.spp.bray.tcdb <- vegdist(h.ps.spp.df.wide3,
                              method = 'bray')

##compute pcoa
h.ps.spp.pcoa.tcdb <- pcoa(h.ps.spp.bray.tcdb)

##extract PCoA1 and PCoA2
h.spp.tcdb.pcoa1 <- 100 * (h.ps.spp.pcoa.tcdb$values$Eigenvalues[1]/sum(h.ps.spp.pcoa.tcdb$values$Eigenvalues))
h.spp.tcdb.pcoa2 <- 100 * (h.ps.spp.pcoa.tcdb$values$Eigenvalues[2]/sum(h.ps.spp.pcoa.tcdb$values$Eigenvalues))

##prep plot----
h.ps.spp.pcoa.plot.tcdb <- data.frame(h.ps.spp.pcoa.tcdb$vectors[,])
h.ps.spp.meta.tcdb <- h.ps.spp.df.wide3
h.ps.spp.meta.tcdb$Sample <- row.names(h.ps.spp.meta.tcdb)
h.ps.spp.meta.tcdb <- h.ps.spp.meta.tcdb %>% 
  left_join(meta.tcdb)
h.ps.spp.meta.tcdb$Axis.1 <- h.ps.spp.pcoa.plot.tcdb$Axis.1
h.ps.spp.meta.tcdb$Axis.2 <- h.ps.spp.pcoa.plot.tcdb$Axis.2
row.names(h.ps.spp.meta.tcdb) <- h.ps.spp.meta.tcdb$Sample
h.ps.spp.meta.tcdb$Group <- factor(h.ps.spp.meta.tcdb$Group,
                                levels = c('Carrier', 'CDI'))

# h.ps.spp.meta.tcdb$group <- as.character(h.ps.spp.meta.tcdb$group)
h.ps.spp.meta.tcdb$tcdB_subtype_grouped <- factor(h.ps.spp.meta.tcdb$tcdB_subtype_grouped,
                                                  levels = c('B1', 'B2', 'B3', 'B5'))
cols <- c('#C77CFF', '#7CAE00', '#F8766D', '#00BFC4')

##run PERMANOVA
adonis2(h.ps.spp.bray.tcdb ~ tcdB_subtype_grouped, data = h.ps.spp.meta.tcdb,
        permutations = 9999, method = 'bray') # p = 0.5002

##plot Fig S1F
p.tcdb.bar <- h.ps.spp.meta.tcdb %>% 
  mutate(Cultured_Cd_Clade = factor(Cultured_Cd_Clade,
                                    levels = c('None', '5', '4', '2', '1'))) %>%
  ggplot(aes(x = Cultured_Cd_Clade)) +
  scale_fill_manual(values = cols) +
  geom_bar(aes(fill = tcdB_subtype_grouped),
           position = "stack") +
  # scale_fill_manual(values = c('#1C1C1C', '#4B8DA5', '#A0A0A0', '#8B5B7E', '#D9D9D9')) +
  labs(x = "Clade", y = "Genome count", fill = "tcdB type") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_plot() +
  theme(legend.position = 'right',
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        aspect.ratio = 0.8) +
  coord_flip()

##plot Fig S1G----
p.tcdb.pcoa <- h.ps.spp.meta.tcdb %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = tcdB_subtype_grouped)) +
  scale_color_manual(values = cols) +
  geom_point(size = 1.5, alpha = 0.8, aes(shape = Group)) +
  stat_ellipse(geom = 'polygon', alpha = 0.02, lwd = 1.2) +
  theme_plot() +
  theme(legend.position = 'right',
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        aspect.ratio = 1)



#Run generalized linear mixed models on CDI vs. carriers (species level) (Fig S1H)----

#Test taxa associations between CDI patients and Cd carriers
Maaslin2(input_data = h.ps.spp.df.wide2,
         input_metadata = filter(h.meta, Colonized == 1),
         output = 'ManuscriptAnalyses/Human/MetaPhlAn4/MaAsLin2/CDI_vs_carrier/Species_default',
         fixed_effects = c('Group'),
         random_effects = c('Participant_ID'),
         normalization = 'NONE')

##load results
cdi.vs.col.sp <- read.delim('Scripts/Anna/files_for_github/data/241125_CDI_v_carrier_sp_significant_results.tsv') %>% 
  dplyr::rename(Species = feature) %>% 
  mutate(log_qval = -log10(qval), direction = coef > 0, comparison = 'cdi.vs.col')

##plot significant associations (q <= 0.05)----
p.spp.cdi.vs.col <- cdi.vs.col.sp %>% 
  filter(qval <= 0.05) %>% 
  ggplot(aes(x = coef, y = reorder(Species, coef))) +
  geom_errorbarh(aes(xmax = coef + stderr, xmin = coef - stderr, height = 0), size = 0.7) +
  geom_point(color = '#18678d', shape = 16, alpha = 0.8, aes(size = log_qval)) +
  theme_minimal() +
  theme(legend.position = 'bottom',
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10),
        panel.border = element_rect(fill = NA, color = 'black', size = 1),
        aspect.ratio = 1.5) +
  geom_vline(xintercept=0, colour='black',size=0.5, linetype="dotted") +
  labs(x = 'Model coefficient (CDI relative to Carriers)')



#Random forest model results (Fig S1I)----
##refer to figure 1 code for performance of the best model, which was then formatted into a matrix
##in adobe illustrator.



#Run generalized linear mixed models on ceph exposure (species level) (Fig S1H)----

abx.meta <- merge(filter(h.meta, Cohort != "HH"),
                  select(abx, contains(c("all", "vanc", "ceph", "carb", "fluoro"))),
                  by = 'row.names')
rownames(abx.meta) <- abx.meta$Row.names
abx.meta$Row.names <- NULL

#Test taxa associations with recent cephalosporin exposure in colonized patients
Maaslin2(input_data = h.ps.spp.df.wide2,
         input_metadata = filter(abx.meta, Colonized == 1),
         output = 'ManuscriptAnalyses/Human/MetaPhlAn4/MaAsLin2/Antibiotic_exposures/Species_ceph_fixed',
         fixed_effects = c('ceph_30d'),
         random_effects = c('Participant_ID'),
         normalization = 'NONE')

##load results
ceph.sp <- read.delim('Scripts/Anna/files_for_github/data/241125_ceph_sp_significant_results.tsv') %>% 
  dplyr::rename(Species = feature) %>% 
  mutate(log_qval = -log10(qval), direction = coef > 0, comparison = 'ceph.vs.not')

##plot significant associations (q <= 0.05)
p.ceph.sp <- ceph.sp %>% 
  filter(qval <= 0.05) %>% 
  ggplot(aes(x = coef, y = reorder(Species, coef))) +
  geom_errorbarh(aes(xmax = coef + stderr, xmin = coef - stderr, height = 0), size = 0.7) +
  geom_point(color = '#18678d', shape = 16, alpha = 0.8, aes(size = log_qval)) +
  theme_minimal() +
  theme(legend.position = 'bottom',
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10, face = 'italic'),
        axis.text.x = element_text(size = 10),
        panel.border = element_rect(fill = NA, color = 'black', size = 1),
        aspect.ratio = 1.5) +
  geom_vline(xintercept=0, colour='black',size=0.5, linetype="dotted") +
  labs(x = 'Model coefficient (Cephalosporin exposure)')

