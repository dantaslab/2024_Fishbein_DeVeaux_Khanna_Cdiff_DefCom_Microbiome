#Code to replicate Figure S4 R plots
#updated Monday, November 25th, 2024

#setup----

library(ggplot2)
library(ggpubr)
library(ggthemes)
library(plyr)
library(dplyr)
library(phyloseq)
library(Maaslin2)
library(rstatix)

##set working directory
setwd('/Users/annadeveaux/Box Sync/Functional consequences of microbiome variation in CDI susceptibility/')

##set seed
set.seed(42)

##load metadata
dc.meta <- read.delim('Scripts/Anna/files_for_github/data/241126_DefCom_Metadata_fltd_clean_v3_AD.txt')
row.names(dc.meta) <- dc.meta$SampleName



#Calprotectin cecal concentration plot (Fig S4D)----

dred <- c("darkolivegreen4", "#77439A", "#424242")
dred_name <- c('DefCom1', 'DefCom2', 'DefCom1 + Bo')
names(dred) <- dred_name
cal.df <- read.delim('Scripts/Anna/files_for_github/data/241022_Calprotectin_DC1DC2Bo_R_input.txt')

stat.test.cal <- cal.df %>% 
  wilcox_test(Calprotectin ~ Arm) %>% 
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance() %>%
  add_y_position()

p.cal <- cal.df %>% 
  ggplot(aes(x = Arm, y = Calprotectin)) +
  geom_boxplot(aes(fill = Arm)) +
  geom_point(position = position_dodge(width = 0.75)) +
  ylab("Calprotectin Cecal Concentration (pg/mg)") +
  stat_pvalue_manual(stat.test.cal,
                     label = 'p.adj', hide.ns = T, step.increase = 0.1) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x =element_text(size = 10, angle = 35, hjust = 1),
        strip.text.x = element_text(size = 12),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
        aspect.ratio = 1.5) +
  scale_fill_manual(values = dred)

# ggsave('ManuscriptDocs/Figures/FigureS4/241125_calprotectin_adjust.pdf',
#        plot = p.cal,
#        width = 3, height = 3.5,
#        units = 'in', device = "pdf")



#TcdB cecal concentration plot (Fig S4E)----

tcdb.df <- read.delim('Scripts/Anna/files_for_github/data/241022_TcdB_DC1DC2Bo_R_input.txt')

stat.test.tcdb <- tcdb.df %>% 
  wilcox_test(TcdB ~ Arm) %>% 
  adjust_pvalue(method = 'bonferroni') %>%
  add_significance() %>%
  add_y_position()

p.tcdb <- tcdb.df %>% 
  ggplot(aes(x = Arm, y = TcdB)) +
  geom_boxplot(aes(fill = Arm)) +
  geom_point(position = position_dodge(width=0.75))+
  ylab("TcdB Cecal Concentration (ng/mg)") +
  stat_pvalue_manual(stat.test.tcdb,
                     label = 'p.adj', hide.ns = T, step.increase = 0.1) +
  theme_classic() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        axis.text.x =element_text(size = 10, angle = 35, hjust = 1),
        strip.text.x = element_text(size = 12),
        panel.border = element_rect(fill = NA, colour = "black", size = 0.5),
        aspect.ratio = 1.5) +
  scale_fill_manual(values = dred)

# ggsave('ManuscriptDocs/Figures/FigureS4/241125_tcdb_adjust.pdf',
#        plot = p.tcdb,
#        width = 3, height = 3.5,
#        units = 'in', device = "pdf")



#Pathway differences between communities (Fig S4G)----

dc.path.raw <- read.delim('Scripts/Anna/files_for_github/data/241127_DefComs_humann3_ra_v2mod.txt',
                      row.names = 1) %>% 
  select(dc.meta$SampleName)

metacyc.map <- read.delim('Scripts/Anna/files_for_github/data/231205_Pwys_metacyc_directparent.txt',
                          header = F, row.names = 1) %>% 
  mutate(PWY_name = row.names(.)) %>% 
  rename(Agg_pwy = V2) %>% 
  select(c(PWY_name, Agg_pwy))

##omit rows that break down the species associations (contain '|' delimiter)
dc.path <- dc.path.raw[!grepl('\\|', rownames(dc.path.raw)), ]

##prep data frame of aggregated pathways, filter out unmapped and unintegrated calls
dc.path.agg <- dc.path %>% 
  filter(row.names(.) != 'UNMAPPED' & row.names(.) != 'UNINTEGRATED')

##renormalize pruned pathway abundance calls to 100
dc.path.agg.renorm <- t(t(dc.path.agg)/rowSums(t(dc.path.agg)))*100

##filter out pathway calls <0.01% abundance (within sample) after initial renormalization
dc.path.agg.filt <- dc.path.agg.renorm
dc.path.agg.filt[dc.path.agg.filt < 0.01] <- 0

##renormalize abundance back to 100%
dc.path.agg.filt.renorm <- t(t(dc.path.agg.filt)/rowSums(t(dc.path.agg.filt)))*100
dc.path.agg.filt.renorm.2 <- dc.path.agg.filt.renorm
rownames(dc.path.agg.filt.renorm.2) <- sub("^(.*?):.*", "\\1", rownames(dc.path.agg.filt.renorm.2))

##convert pathway calls to descriptive labels
dc.path.agg.rename <- merge(metacyc.map,
                            dc.path.agg.filt.renorm.2,
                            by = 'row.names') %>% 
  select(!c('Row.names', 'PWY_name'))

##sum total abundance for a given pathway into one value
dc.path.agg.hdf <- ddply(dc.path.agg.rename,
                         .(Agg_pwy),
                         numcolwise(sum))
row.names(dc.path.agg.hdf) <- dc.path.agg.hdf$Agg_pwy
dc.path.agg.hdf$Agg_pwy <- NULL

ps.dc.path <- phyloseq(otu_table(dc.path.agg.hdf,
                                 taxa_are_rows = T),
                       sample_data(dc.meta))
ps.dc.path.df <- psmelt(ps.dc.path)

fit <- Maaslin2(dc.path.agg.hdf,
         dc.meta,
         'ManuscriptAnalyses/DefCom/DefCom_code_FigS3_4_S4/241106_defcom_humann_test3',
         fixed_effects = c("MicrobiomeType"),
         transform = "LOG",
         random_effects = c("Cage", "Sex", "Strain"),
         normalization = "NONE")

humann <- fit$results[fit$results$qval<1e-15,]
humann$feature <- factor(humann$feature,
                             levels = unique(humann$feature[order(humann$coef)]))

p.humann <- humann %>% 
  ggplot(aes(x = feature, y = coef)) +
  geom_point(aes(size=2*-log10(qval))) +
  theme_linedraw()+
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 15)) +
  xlab("Pathways") +
  ylab("Coefficient") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 8,
      colour = "black"
    ),
    axis.title.x = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.title.y = element_text(
      size = 14,
      face = "bold",
      colour = "black"
    ),
    axis.text.y = element_text(size = 12, colour = "black"),
    axis.line.x = element_line(color="black", size = 0.3),
    axis.line.y = element_line(color="black", size = 0.3),
    panel.border = element_rect(colour = "black", fill = NA, size = 0.3),
    legend.position = "bottom",
    plot.margin = margin(20, 20, 20, 80, "pt"))

# ggsave('ManuscriptAnalyses/DefCom/DefCom_code_FigS3_4_S4/241127_pwy_defcoms_dna.pdf',
#        plot = p.humann,
#        width = 7, height = 6,
#        units = 'in', device = "pdf")




