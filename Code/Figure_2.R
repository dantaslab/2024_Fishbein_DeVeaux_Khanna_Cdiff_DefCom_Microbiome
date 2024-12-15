#Code to replicate Figure 2C & D results
#updated Thursday, October 31st, 2024

##load libraries
library(dplyr)            #v1.1.4
library(vegan)            #v2.6-4
library(phyloseq)         #v1.46.0
library(ape)              #v5.7-1
library(spatstat.geom)    #v3.1-2
library(ggrepel)          #v0.9.5   
library(reshape2)         #v1.4.4
library(NatParksPalettes) #0.2.0
library(ggpubr)           #0.6.0

##set working directory
setwd('/Users/annadeveaux/Box Sync/Functional consequences of microbiome variation in CDI susceptibility/')

##load rds
ps <- readRDS('Scripts/Anna/files_for_github/data/ps_human_m4_final.rds')

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

##define plot colors
pal <- natparks.pals('Triglav')[c(3, 5, 4, 1)]

##load metadata
h.meta <- read.delim('Scripts/Anna/files_for_github/data/241031_Human_Metadata.csv',
                     sep = ',')
row.names(h.meta) <- h.meta$Sample

##create wide (filtered) species data frame
ps.spp <- tax_glom(ps, taxrank = 'Species')
ps.spp.df <- psmelt(ps.spp)
ps.spp.df.wide <- dcast(ps.spp.df,
                        Sample ~ Species,
                        value.var = "Abundance")

##prep species data frame
ps.spp.df.wide2 <- ps.spp.df.wide
row.names(ps.spp.df.wide2) <- ps.spp.df.wide2$Sample
ps.spp.df.wide2$Sample <- NULL

##compute dissimilarity metrics
ps.spp.bray.all <- vegdist(filter(ps.spp.df.wide2,
                                  row.names(ps.spp.df.wide2) %in% (filter(h.meta,
                                                                          Primary_Sample != 0)$Sample)),
                           method = 'bray')

##compute pcoa
ps.spp.pcoa.all <- pcoa(ps.spp.bray.all)

##extract PC1 and PC2
spp.all.pc1 <- 100 * (ps.spp.pcoa.all$values$Eigenvalues[1]/sum(ps.spp.pcoa.all$values$Eigenvalues))
spp.all.pc2 <- 100 * (ps.spp.pcoa.all$values$Eigenvalues[2]/sum(ps.spp.pcoa.all$values$Eigenvalues))

##prep plot
ps.spp.pcoa.plot.all <- data.frame(ps.spp.pcoa.all$vectors[,])
ps.spp.meta.all <- filter(ps.spp.df.wide2,
                          row.names(ps.spp.df.wide2) %in% (filter(h.meta,
                                                                  Primary_Sample != 0)$Sample))
ps.spp.meta.all$Sample <- row.names(ps.spp.meta.all)
ps.spp.meta.all <- left_join(ps.spp.meta.all, h.meta, by = 'Sample')
ps.spp.meta.all$Axis.1 <- ps.spp.pcoa.plot.all$Axis.1
ps.spp.meta.all$Axis.2 <- ps.spp.pcoa.plot.all$Axis.2
row.names(ps.spp.meta.all) <- ps.spp.meta.all$Sample
ps.spp.meta.all$Group <- factor(ps.spp.meta.all$Group,
                                levels = c('CDI', 'Carrier', 'IC', 'HV'))

##plot PCoA colored by FMT donor
p.pcoa.spp.donor <- ps.spp.meta.all %>% 
  ggplot(aes(x = Axis.1, y = Axis.2, color = Group)) +
  stat_ellipse(geom = 'polygon', alpha = 0.02, lwd = 1.2) +
  geom_point(aes(alpha = Engraftment_Donor, size = Engraftment_Donor, fill = Group,
                 shape = Lethality_Outcome)) +
  scale_size_continuous(range = c(2, 4)) +
  scale_alpha_continuous(range = c(0.3, 0.95)) +
  scale_color_manual(values = pal) +
  theme_plot() +
  theme(legend.position = 'right',
        axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        aspect.ratio = 0.8) +
  labs(x = paste("PCoA 1 (", round(spp.all.pc1, 1), '%)', sep = ''),
       y = paste("PCoA 2 (", round(spp.all.pc2, 1), '%)', sep = ''))

#compute centroid of healthy humans cluster
hh.x.spp <- mean(filter(ps.spp.meta.all, Cohort == "HH")$Axis.1)
hh.y.spp <- mean(filter(ps.spp.meta.all, Cohort == "HH")$Axis.2)

#compute distances between FMT donors and centroid of HH points
hh.fmt.dist.spp <- ps.spp.meta.all %>% 
  filter(Engraftment_Donor == 1) %>% 
  mutate(dist = crossdist(Axis.1, Axis.2, hh.x.spp, hh.y.spp))

##plot distances of donors from HVs by engraftment CFU/mL
p.hh.dist.cfu <- ggscatter(hh.fmt.dist.spp,
                           x = 'CFU_per_mL',
                           y = 'dist',
                           color = 'Lethality_Outcome',
                           palette = c('red', 'black'),
                           ggtheme = theme_plot(),
                           add = "reg.line",
                           add.params = list(color = "grey", fill = "lightgray"),
                           xlab = 'CFU per mL',
                           ylab = 'Distance from HV cluster') +
  xscale("log10", .format = TRUE) +
  stat_cor() +
  geom_text_repel(aes(label = Donor_ID)) +
  theme(axis.title.x = element_text(size=15),
        axis.text.x = element_text(size=15, angle=45, hjust=1),
        aspect.ratio = 0.8,
        legend.position = c(.85, .2),
        legend.background = element_rect(colour = 'black', fill = 'white', linetype='solid'))
