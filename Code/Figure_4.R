#Code to replicate Figure 4 R plots
#updated Tuesday, November 26th, 2024

#setup----

library(plyr)             #1.8.9
library(dplyr)            #1.1.4
library(phyloseq)         #1.46.0
library(reshape2)         #1.4.4
library(ggplot2)          #3.4.4
library(NatParksPalettes) #0.2.0
library(Maaslin2)         #1.8.0

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
dc.meta <- read.delim('Scripts/Anna/files_for_github/data/241126_DefCom_Metadata_fltd_clean_v3_AD.txt')
row.names(dc.meta) <- dc.meta$SampleName

###prep MetaPhlAn4 output----
##format raw output for phyloseq
dc.m4.raw <- read.csv('Scripts/Anna/files_for_github/data/241126_defcom_m4_filtered.csv',
                        check.names = F, row.names = 1)

##filter out rows containing strain-level annotations that are new to metaphlan4
dc.m4.raw <- dc.m4.raw %>% filter(!grepl('t__',rownames(.)))

##take only rows that include s__ (species level annotation)
dc.m4.sp <- dc.m4.raw[grepl('s__', rownames(dc.m4.raw)), ]

##filter out species calls below < 0.1% abundance (within sample)
dc.m4.sp.filt <- dc.m4.sp
dc.m4.sp.filt[dc.m4.sp.filt < 0.1] <- 0

##remove rows that sum to zero (lowly abundant taxa that were filtered out in the previous step)
dc.m4.sp.filt <- dc.m4.sp.filt[rowSums(dc.m4.sp.filt[ ]) > 0,]

##renormalize abundance to 100%
dc.m4.sp.filt.renorm <- t((t(dc.m4.sp.filt) / rowSums(t(dc.m4.sp.filt))) * 100)

##create phyloseq taxa table
dc.spp <- data.frame(Names = rownames(dc.m4.sp.filt.renorm))
dc.spp <- data.frame(do.call('rbind', strsplit(as.character(dc.spp$Names),'|',fixed=TRUE)))
rownames(dc.spp) <- rownames(dc.m4.sp.filt.renorm)
colnames(dc.spp) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

##remove prefixes from taxa
dc.spp$Domain <- gsub('k__', '', dc.spp$Domain)
dc.spp$Phylum <- gsub('p__', '', dc.spp$Phylum)
dc.spp$Class <- gsub('c__', '', dc.spp$Class)
dc.spp$Order <- gsub('o__', '', dc.spp$Order)
dc.spp$Family <- gsub('f__', '', dc.spp$Family)
dc.spp$Genus <- gsub('g__', '', dc.spp$Genus)
dc.spp$Species <- gsub('s__', '', dc.spp$Species)

##convert to matrix
dc.spp <- as.matrix(dc.spp)

##make phyloseq object
dc.ps <- phyloseq(otu_table(dc.m4.sp.filt.renorm, taxa_are_rows = TRUE),
                  sample_data(dc.meta),
                  tax_table(dc.spp))

##species
##create wide species data frame
dc.ps.spp <- tax_glom(dc.ps, taxrank = 'Species')
dc.ps.spp.df <- psmelt(dc.ps.spp)
dc.ps.spp.df.wide <- dcast(dc.ps.spp.df,
                           Sample ~ Species,
                           value.var = "Abundance")
row.names(dc.ps.spp.df.wide) <- dc.ps.spp.df.wide$Sample
dc.ps.spp.df.wide$Sample <- NULL

##genus
##create wide genus data frame
dc.ps.genus <- tax_glom(dc.ps, taxrank = 'Genus')
dc.ps.genus.df <- psmelt(dc.ps.genus)
dc.ps.genus.df.wide <- dcast(dc.ps.genus.df,
                             Sample ~ Genus,
                             value.var = "Abundance")
row.names(dc.ps.genus.df.wide) <- dc.ps.genus.df.wide$Sample
dc.ps.genus.df.wide$Sample <- NULL



#Stacked barplots comparing community composition (Figure 4A)----

##get representative DefCom stacked barplots (4B)
d28.spp <- subset_samples(dc.ps, (Day == 28 & MicrobiomeType != 'DefCom1_Blautiaobeum'))
d28.stacked <- psmelt(d28.spp)

medians.dc1 <- ddply(filter(d28.stacked,
                            MicrobiomeType == 'DefCom1'), ~ Species,
                     function(x) c(median = median(x$Abundance)))

medians.dc2 <- ddply(filter(d28.stacked,
                            MicrobiomeType == 'DefCom2'), ~ Species,
                     function(x) c(median = median(x$Abundance)))

medians_average <- medians.dc1 %>%
  full_join(medians.dc2, by = "Species", suffix = c("_community1", "_community2")) %>%
  mutate(average_relab = (median_community1 + median_community2)/2)

remainder.dc1 <- medians.dc1[medians.dc1$median <= 1, ]$Species
remainder.dc2 <- medians.dc2[medians.dc2$median <= 1, ]$Species

remainder.shared <- intersect(remainder.dc1,
                              remainder.dc2)

d28.stacked[d28.stacked$Species %in% remainder.shared, ]$Species <- 'Other'

d28.stacked$Species <- factor(d28.stacked$Species,
                              levels = c((arrange(medians_average, desc(average_relab)))$Species, 'Other'))

p.d28.stacked <- d28.stacked %>% 
  arrange(Abundance) %>% 
  ggplot(aes(x = MicrobiomeType, y = Abundance, fill = Species)) +
  geom_bar(stat = 'identity', position = 'fill', alpha = 0.8) +
  theme_linedraw() +
  theme(panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.line = element_line(),
        panel.grid.minor.y = element_blank(),
        panel.border = element_blank(),
        aspect.ratio = 2.5) +
  xlab('') +
  scale_fill_manual(values = c(natparks.pals('Torres'),
                               natparks.pals('CapitolReef'), natparks.pals('CraterLake'), 'grey'))



#Differential pathway expression (Figure 4F)----

##load cecal sample metadata
cec.meta <- read.delim('Scripts/Anna/files_for_github/data/240515_DefCom_repeat_sample_metadata.txt')
row.names(cec.meta) <- cec.meta$Sample

##load pathway output
cec.path.rna.raw <- read.delim('Scripts/Anna/files_for_github/data/240520_RNA_humann3_relab_pathabundance_merged.tsv', sep = '\t',
                               check.names = F, row.names = 1)

names(cec.path.rna.raw) <- gsub("_Abundance$", "", names(cec.path.rna.raw))

##omit rows that break down the species associations (contain '|' delimiter)
cec.path.rna <- cec.path.rna.raw[!grepl('\\|', rownames(cec.path.rna.raw)), ]

##remove unintegrated and unmapped pathway calls
cec.path.rna.pruned <- cec.path.rna %>% 
  filter(row.names(.) != 'UNMAPPED' & row.names(.) != 'UNINTEGRATED')

##renormalize pruned abundance calls to 100
cec.path.rna.pruned.renorm <- t((t(cec.path.rna.pruned) / rowSums(t(cec.path.rna.pruned))) * 100)

Maaslin2(cec.path.rna.pruned.renorm,
         cec.meta,
         'Metatranscriptomics/HUMAnN3/MaAsLin2_RNA_Pathways_community_strain_cage_sex_random_FILTERED',
         fixed_effects = c('Community'),
         random_effects = c('Strain', 'Cage', 'Sex'),
         normalization = 'NONE')

##load significant RNA pathway differences
cec.path.diff <- read.delim('Scripts/Anna/files_for_github/data/241125_RNAseq_significant_results.tsv') %>% 
  dplyr::rename(Pathway = feature) %>%
  mutate(log_qval = -log10(qval), direction = coef > 0)

##plot----
##manually trimmed pathway names in illustrator
p.cec.path.diff <- cec.path.diff %>% 
  filter(qval <= 0.05) %>% 
  ggplot(aes(x = coef, y = reorder(Pathway, coef))) +
  geom_errorbarh(aes(xmax = coef + stderr, xmin = coef - stderr, height = 0), size = 0.7) +
  geom_point(aes(size = log_qval, color = direction)) +
  scale_color_manual(values = c("darkolivegreen4","#77439A")) +
  theme_linedraw() +
  theme(legend.position = 'bottom',
        axis.title.x = element_text(size = 12),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size=10, angle=45, hjust=1),
        panel.border = element_rect(fill = NA, color = 'black', size = 1)) +
  geom_vline(xintercept=0, colour='black',size=0.5, linetype="dotted") +
  labs(x = 'DefCom2 expression relative to DefCom1') +
  coord_flip()
