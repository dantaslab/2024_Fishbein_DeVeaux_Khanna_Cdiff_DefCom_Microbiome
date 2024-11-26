#Code to replicate Figure S3 R plots
#updated Tuesday, November 26th, 2024

#setup----

library(plyr)             #1.8.9
library(dplyr)            #1.1.4
library(phyloseq)         #1.46.0
library(reshape2)         #1.4.4
library(ggplot2)          #3.4.4
library(ggrepel)          #0.9.5   
library(scales)           #1.3.0
library(vegan)            #2.6-4
library(ape)              #5.7-1
library(NatParksPalettes) #0.2.0
library(broom)            #1.0.5
library(rstatix)          #0.7.2
library(ggpubr)           #0.6.0

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

##load human metadata
h.meta <- read.delim('Scripts/Anna/files_for_github/data/241031_Human_Metadata.csv',
                     sep = ',')
row.names(h.meta) <- h.meta$Sample

##load defcom metadata

dc.meta <- read.delim('Scripts/Anna/files_for_github/data/241126_DefCom_Metadata_fltd_clean_v3_AD.txt')
row.names(dc.meta) <- dc.meta$SampleName

##load phyloseq object with human taxonomic data (created in Figure 1 code)
h.ps <- readRDS('Scripts/Anna/files_for_github/data/ps_human_m4_final.rds')

##load phyloseq object with defcom taxonomic data (created in Figure 4 code)
dc.ps <- readRDS('Scripts/Anna/files_for_github/data/ps_defcoms_m4_final.rds')



#Colonized patient genus prev/relab plot (S3A)----

##genus
##create wide genus data frame and filter by Colonized patients (n = 149, including repeat samples)
h.ps.genus <- tax_glom(h.ps, taxrank = 'Genus')
h.ps.genus.df <- psmelt(h.ps.genus)
h.ps.genus.df.wide <- dcast(h.ps.genus.df, Sample ~ Genus, value.var = "Abundance")
h.ps.genus.df.wide.col <- h.ps.genus.df.wide %>% 
  filter(Sample %in% (filter(h.meta, Colonized == 1))$Sample)
row.names(h.ps.genus.df.wide.col) <- h.ps.genus.df.wide.col$Sample
h.ps.genus.df.wide.col$Sample <- NULL

##compute prevalence by species
h.ps.genus.col.prev <- h.ps.genus.df.wide.col
h.ps.genus.col.prev[h.ps.genus.col.prev != 0] = 1
h.ps.genus.col.prev.sums <- as.data.frame(colSums(h.ps.genus.col.prev)) %>% 
  transmute(Prevalence = `colSums(h.ps.genus.col.prev)`)
h.ps.genus.col.prev.sums$Genus <- row.names(h.ps.genus.col.prev.sums)

##determine average relative abundance of genera among colonized patients----

h.ps.genus.col.relab <- as.data.frame(colMeans(h.ps.genus.df.wide.col)) %>% 
  transmute(AvgRelAb = `colMeans(h.ps.genus.df.wide.col)`)
h.ps.genus.col.relab$Genus <- row.names(h.ps.genus.col.relab)

##plot defcom2 taxa prevalence amongst Colonized patients----
##load defcom1/2 metaphlan4 species map
dc.map <- read.delim('Scripts/Anna/files_for_github/data/241105_defcom_species_map.txt')
dc.map$Species <- gsub('s__', '', dc.map$Species)

h.spp <- readRDS('Scripts/Anna/files_for_github/data/spp_human_final.RDS')

dc.map <- left_join(dc.map, select(as.data.frame(h.spp), c('Species', 'Genus')),
                    by = 'Species')

##manually assign genera to species without representation in metaphlan database
dc.map$Genus[is.na(dc.map$Genus) & dc.map$Species == "Hungatella_effluvii"] <- 'Hungatella'
dc.map$Genus[is.na(dc.map$Genus) & dc.map$Species == "Streptococcus_caecimuris"] <- 'Streptococcus'
dc.map$Genus[is.na(dc.map$Genus) & dc.map$Species == "Blautia_coccoides"] <- 'Blautia'

##select unique genera with their assigned defcom(s)
dc.map.genus <- dc.map %>% 
  select(c('Genus', 'Community')) %>% 
  group_by(Genus) %>% 
  distinct()

##if multiple species assigned to the same genera are present, this will only select the "both" designation,
##if it exists for that genus
dc.map.genus2 <- dc.map.genus[!duplicated(dc.map.genus$Genus), ]

##combine genus-level data frames
h.ps.genus.prev.relab <- left_join(h.ps.genus.col.prev.sums,
                                   h.ps.genus.col.relab, by = 'Genus')
dc.map.genus.prev.relab <- left_join(h.ps.genus.prev.relab,
                                     dc.map.genus2, by = 'Genus')

##assign a new label for genera represented in defcom(s)
dc.map.genus.prev.relab <- dc.map.genus.prev.relab %>%
  arrange(desc(Prevalence)) %>% 
  mutate(present = ifelse(!is.na(Community),
                          'present', 'not present')) %>% 
  mutate(label = ifelse(Genus %in% dplyr::slice(dc.map.genus.prev.relab, 1:15)$Genus &
                          present == 'present',
                        Genus, "")) %>% 
  mutate(label2 = ifelse(present == 'present',
                         Genus, ""))

##plot (filtering for >0.01% relab for visualization)
p.prev.relab <- dc.map.genus.prev.relab %>%
  filter(Prevalence > 0 & AvgRelAb > 0.01) %>%
  ggplot(aes(x = 100 * (Prevalence/nrow(h.ps.genus.col.prev)), y = AvgRelAb, color = Community)) +
  geom_point(aes(alpha = present, size = present)) +
  geom_text_repel(aes(label = label2), color = 'black',
                  arrow = arrow(length = unit(0.02, "npc")),
                  max.overlaps = 30,
                  size = 2,
                  force = 5) +
  scale_alpha_discrete(range = c(0.4, 0.8)) +
  scale_size_discrete(range = c(1.5, 3)) +
  scale_color_manual(values = c('black', '#7F9A4E', '#77429A')) +
  scale_y_continuous(trans = 'log10', labels = label_number(accuracy = 0.01)) +
  theme_plot() +
  # theme(legend.position = c(.8, .3),
  #       legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'),
  #       axis.title.x = element_text(size=15)) +
  labs(x = 'Prevalence (%)', y = 'Average Relative Abundance (%)') +
  xlim(0, 100)



#DefCom day 28 distances overlaid on human data (S3B)----

d28.spp <- subset_samples(dc.ps, (Day == 28))
d28.genus <- tax_glom(d28.spp, taxrank = 'Genus')
g.d28.stacked <- psmelt(d28.genus)

g.means.dc1 <- ddply(filter(g.d28.stacked,
                            MicrobiomeType == 'DefCom1'), ~ Genus,
                     function(x) c(mean = mean(x$Abundance)))

g.means.dc2 <- ddply(filter(g.d28.stacked,
                            MicrobiomeType == 'DefCom2'), ~ Genus,
                     function(x) c(mean = mean(x$Abundance)))

g.means_average <- g.means.dc1 %>%
  full_join(g.means.dc2, by = "Genus", suffix = c("_community1", "_community2")) %>%
  mutate(average_relab = (mean_community1 + mean_community2)/2)

g.dc.means <- as.data.frame(g.means_average) %>%
  select(!average_relab) %>% 
  melt() %>% 
  dplyr::mutate(Community = case_when((variable == 'mean_community1') ~ 'DefCom1',
                                      (variable == 'mean_community2') ~ 'DefCom2')) %>% 
  dplyr::rename(Sample = Community,
                Abundance = value) %>% 
  select(!variable)

full.ps.genus.df <- full_join(h.ps.genus.df,
                              g.dc.means)

full.ps.genus.df.wide <- dcast(full.ps.genus.df,
                               Sample ~ Genus,
                               value.var = "Abundance")

full.ps.genus.df.wide[is.na(full.ps.genus.df.wide)] <- 0
row.names(full.ps.genus.df.wide) <- full.ps.genus.df.wide$Sample
full.ps.genus.df.wide$Sample <- NULL

full.ps.genus.df.wide.filt <- filter(full.ps.genus.df.wide,
                                     row.names(full.ps.genus.df.wide) %in% (filter(h.meta,
                                                                                   Colonized == 1 & Primary_Sample != 0)$Sample) |
                                       row.names(full.ps.genus.df.wide) %in% c('DefCom1', 'DefCom2'))

full.ps.genus.bray <- vegdist(full.ps.genus.df.wide.filt)

##compute pcoa
full.ps.genus.pcoa <- pcoa(full.ps.genus.bray)

##extract PCoA1 and PCoA2
full.genus.pc1 <- 100 * (full.ps.genus.pcoa$values$Eigenvalues[1]/sum(full.ps.genus.pcoa$values$Eigenvalues))
full.genus.pc2 <- 100 * (full.ps.genus.pcoa$values$Eigenvalues[2]/sum(full.ps.genus.pcoa$values$Eigenvalues))

##prep plot
full.ps.genus.pcoa.plot <- data.frame(full.ps.genus.pcoa$vectors[,])
full.genus.plot.meta <- as.data.frame(row.names(full.ps.genus.df.wide.filt)) %>% 
  dplyr::rename(Sample = `row.names(full.ps.genus.df.wide.filt)`) %>%
  left_join(.,
            select(h.meta, c('Sample', 'Group'))) %>% 
  mutate(Group = case_when((Sample == 'DefCom1') ~ 'DefCom1',
                           (Sample == 'DefCom2') ~ 'DefCom2',
                           (Group == 'CDI') ~ 'CDI',
                           (Group == 'Carrier') ~ 'Carrier'))
full.genus.plot.meta$Axis.1 <- full.ps.genus.pcoa.plot$Axis.1
full.genus.plot.meta$Axis.2 <- full.ps.genus.pcoa.plot$Axis.2
row.names(full.genus.plot.meta) <- full.genus.plot.meta$Sample
full.genus.plot.meta$Group <- factor(full.genus.plot.meta$Group,
                                     levels = c('Carrier', 'CDI', 'DefCom1', 'DefCom2'))
full.genus.plot.meta <- full.genus.plot.meta %>% 
  mutate(Type = case_when((Sample == 'DefCom1') ~ 'DefCom',
                          (Sample == 'DefCom2') ~ 'DefCom',
                          (Group == 'CDI') ~ 'patient',
                          (Group == 'Carrier') ~ 'patient'))

##PERMANOVA
##between colonized and defcoms
adonis2(full.ps.genus.bray ~ Type, data = full.genus.plot.meta,
        permutations=9999, method="bray") # p = 0.245

full.genus.plot.meta$Type <- factor(full.genus.plot.meta$Type,
                                    levels = c('patient', 'DefCom'))

##plot
##not colored by FMT donor
p.full.genus.plot.meta <- full.genus.plot.meta %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = Group)) +
  geom_point(aes(size = Type),
             alpha = 0.8) +
  # scale_size_continuous(range = c(2, 4)) +
  # scale_alpha_continuous(range = c(0.3, 0.95)) +
  scale_color_manual(values = c(natparks.pals('Triglav')[c(3, 5)], 'darkolivegreen4', '#77439A')) +
  theme_plot() +
  theme(axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        aspect.ratio = 1.5) +
  labs(x = paste("PCoA 1 (", round(full.genus.pc1, 1), '%)', sep = ''),
       y = paste("PCoA 2 (", round(full.genus.pc2, 1), '%)', sep = ''))



#PCoA of Bray-Curtis distances stratified by DefCom (S3C)----

colsnames <- c("DefCom1", "DefCom2")
colscolz <- c("darkolivegreen4", "#77439A")
names(colscolz) <- colsnames

d28.ord <- ordinate(d28.spp, 'PCoA')

##extract PCoA1 and PCoA2
d28.spp.pc1 <- 100 * (d28.ord$values$Eigenvalues[1]/sum(d28.ord$values$Eigenvalues))
d28.spp.pc2 <- 100 * (d28.ord$values$Eigenvalues[2]/sum(d28.ord$values$Eigenvalues))

p2 <- plot_ordination(d28.spp,
                      d28.ord,
                      type = 'samples',
                      color = 'MicrobiomeType',
                      shape = 'Experiment',
                      justDF = TRUE)

p.defcom.pcoa <- p2 %>%
  ggplot(aes(x = Axis.1, y = Axis.2,
             color = MicrobiomeType,
             shape = Experiment,
             label = MouseIDExp)) +
  geom_point(aes(size = 5, alpha = 0.7)) +
  theme_linedraw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        aspect.ratio = 1.5) +
  labs(x = paste("PCoA 1 (", round(d28.spp.pc1, 1), '%)', sep = ''),
       y = paste("PCoA 2 (", round(d28.spp.pc2, 1), '%)', sep = '')) +
  scale_color_manual(values = colscolz)



#Bray-Curtis distance comparison (S3D)----

d28.dist <- phyloseq::distance(d28.spp,
                               method = 'bray')
d28.dist.df <- tidy(d28.dist)

d28.dist.labels <- data.frame(sample_data(d28.spp)[, c('SampleName', 'MicrobiomeType', 'Experiment')])
colnames(d28.dist.labels)[1] <- 'item1'
d28.dist.df1 <- merge(d28.dist.df, d28.dist.labels,
                      by = 'item1')
colnames(d28.dist.df1) <- c('item1', 'item2', 'distance', 'Cohort1', 'Experiment1')
colnames(d28.dist.labels)[1] <- 'item2'
d28.dist.df2 <- merge(d28.dist.df1, d28.dist.labels,
                      by = 'item2')
colnames(d28.dist.df2) <- c('item1', 'item2', 'distance', 'Cohort1', 'Experiment1', 'Cohort2', 'Experiment2')
d28.dist.df2$SameCohort <- d28.dist.df2$Cohort1 == d28.dist.df2$Cohort2
d28.dist.df2$DifferentReplicate <- d28.dist.df2$Experiment1 != d28.dist.df2$Experiment2
d28.dist.between <- d28.dist.df2$distance[d28.dist.df2$SameCohort == FALSE]
d28.dist.within <- d28.dist.df2$distance[d28.dist.df2$SameCohort == TRUE]
d28.dist.bw.rep <- d28.dist.df2$distance[(d28.dist.df2$SameCohort == TRUE) &
                                           (d28.dist.df2$DifferentReplicate == TRUE)]
d28.dists <- c(d28.dist.between,
               d28.dist.within,
               d28.dist.bw.rep)
d28.dist.annot <- c(rep('Between DefCom', length(d28.dist.between)),
                    rep('Within DefCom', length(d28.dist.within)),
                    rep('Between Replicate', length(d28.dist.bw.rep)))

plot.d28.dist <- data.frame(d28.dists,
                            d28.dist.annot)

stat.test.d28.dist <- plot.d28.dist %>% 
  wilcox_test(d28.dists ~ d28.dist.annot) %>% 
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance() %>%
  add_y_position()

p.d28.dist <- plot.d28.dist %>% 
  ggplot(aes(x = d28.dist.annot,
             y = d28.dists)) +
  # geom_jitter(alpha = 0.3) +
  geom_violin() +
  # geom_boxplot(outlier.shape = NA) +
  # stat_compare_means(ref.group = 'Between DefCom') +
  theme_linedraw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 1) +
  xlab("Group Comparisons") +
  ylab("Beta diversity (Bray-Curtis Dissimilarity)") +
  stat_pvalue_manual(stat.test.d28.dist,
                     label = 'p.adj', hide.ns = T, step.increase = 0.1)



#DefCom avg genus relab across experiments (Fig S3E)----

##compute relative abundance at genus level for each experiment----

dc.ps.genus <- tax_glom(dc.ps, taxrank = 'Genus')
dc.ps.genus.df <- psmelt(dc.ps.genus)
dc.ps.genus.df.wide <- dcast(dc.ps.genus.df,
                          Sample ~ Genus,
                          value.var = "Abundance")
row.names(dc.ps.genus.df.wide) <- dc.ps.genus.df.wide$Sample
dc.ps.genus.df.wide$Sample <- NULL

##defcom1, experiment A
ps.a.dc1.relab.genus <- dc.ps.genus.df.wide %>% 
  filter(row.names(.) %in% (filter(dc.meta, (MicrobiomeType == 'DefCom1' & Day == 28) & Experiment == 'A'))$SampleName) %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  dplyr::rename(AvgRelAb = '.') %>% 
  mutate(Experiment = 'A', Community = 'DefCom1')
ps.a.dc1.relab.genus$Genus <- rownames(ps.a.dc1.relab.genus)
rownames(ps.a.dc1.relab.genus) <- NULL

##defcom1, experiment B
ps.b.dc1.relab.genus <- dc.ps.genus.df.wide %>% 
  filter(row.names(.) %in% (filter(dc.meta, (MicrobiomeType == 'DefCom1' & Day == 28) & Experiment == 'B'))$SampleName) %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  dplyr::rename(AvgRelAb = '.') %>% 
  mutate(Experiment = 'B', Community = 'DefCom1')
ps.b.dc1.relab.genus$Genus <- rownames(ps.b.dc1.relab.genus)
rownames(ps.b.dc1.relab.genus) <- NULL

##defcom1, experiment C
ps.c.dc1.relab.genus <- dc.ps.genus.df.wide %>% 
  filter(row.names(.) %in% (filter(dc.meta, (MicrobiomeType == 'DefCom1' & Day == 28) & Experiment == 'C'))$SampleName) %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  dplyr::rename(AvgRelAb = '.') %>% 
  mutate(Experiment = 'C', Community = 'DefCom1')
ps.c.dc1.relab.genus$Genus <- rownames(ps.c.dc1.relab.genus)
rownames(ps.c.dc1.relab.genus) <- NULL

##defcom2, experiment A
ps.a.dc2.relab.genus <- dc.ps.genus.df.wide %>% 
  filter(row.names(.) %in% (filter(dc.meta, (MicrobiomeType == 'DefCom2' & Day == 28) & Experiment == 'A'))$SampleName) %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  dplyr::rename(AvgRelAb = '.') %>% 
  mutate(Experiment = 'A', Community = 'DefCom2')
ps.a.dc2.relab.genus$Genus <- rownames(ps.a.dc2.relab.genus)
rownames(ps.a.dc2.relab.genus) <- NULL

##defcom2, experiment B
ps.b.dc2.relab.genus <- dc.ps.genus.df.wide %>% 
  filter(row.names(.) %in% (filter(dc.meta, (MicrobiomeType == 'DefCom2' & Day == 28) & Experiment == 'B'))$SampleName) %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  dplyr::rename(AvgRelAb = '.') %>% 
  mutate(Experiment = 'B', Community = 'DefCom2')
ps.b.dc2.relab.genus$Genus <- rownames(ps.b.dc2.relab.genus)
rownames(ps.b.dc2.relab.genus) <- NULL

##defcom2, experiment C
ps.c.dc2.relab.genus <- dc.ps.genus.df.wide %>% 
  filter(row.names(.) %in% (filter(dc.meta, (MicrobiomeType == 'DefCom2' & Day == 28) & Experiment == 'C'))$SampleName) %>% 
  colMeans() %>% 
  as.data.frame() %>% 
  dplyr::rename(AvgRelAb = '.') %>% 
  mutate(Experiment = 'C', Community = 'DefCom2')
ps.c.dc2.relab.genus$Genus <- rownames(ps.c.dc2.relab.genus)
rownames(ps.c.dc2.relab.genus) <- NULL

relab.genus.all <- rbind(ps.a.dc1.relab.genus, ps.b.dc1.relab.genus, ps.c.dc1.relab.genus,
                         ps.a.dc2.relab.genus, ps.b.dc2.relab.genus, ps.c.dc2.relab.genus)
relab.genus.all$Experiment <- factor(relab.genus.all$Experiment,
                                     levels = c('C', 'B', 'A'))

##manually add back genera that did not show up in either community in vivo but were present in input material
##for at least one community
genera.zero <- data.frame(
  AvgRelAb = rep(c(0, 0, 0, 0, 0, 0), 4),
  Experiment = rep(c('A', 'B', 'C', 'A', 'B', 'C'), 4),
  Community = rep(c(rep('DefCom1', 3), rep('DefCom2', 3)), 4),
  Genus = c(rep('Klebsiella', 6),
            rep('Lacticaseibacillus', 6),
            rep('Limosilactobacillus', 6),
            rep('Streptococcus', 6))
)

relab.genus.all.v2 <- rbind(relab.genus.all,
                            genera.zero)


##filter out genera in the metaphlan table with 0% abundance and were not in either input gavage material
p.relab.genus.all <- relab.genus.all.v2 %>% 
  filter(Genus %ni% c('Achromobacter', 'Bacillus', 'Clostridioides', 'GGB33432')) %>% 
  ggplot(aes(x = Genus, y = Experiment, size = log10(AvgRelAb), color = Experiment)) +
  geom_point() +
  facet_wrap(~Community, ncol = 1) +
  theme_linedraw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size=15),
        axis.text.x = element_text(size=12, angle=45, hjust=1),
        strip.text.x = element_text(size=15),
        aspect.ratio = 0.2) +
  scale_color_manual(values = c('darkblue', 'darkgreen', 'darkred'))



