# Random forest analyses for Fishbein & DeVeaux et al., 'Commensal support of Clostridioides
# difficile colonization suppresses virulence'. (2024)    
# Last updated Tuesday, Jan 23rd, 2024
#                                                                              
# Bioconductor version 3.18 (BiocManager 1.30.22), R 4.3.1 (2023-06-16)
#
# Adapted from Ferreiro et al., 'Gut microbiome composition may be an indicator of preclinical
# Alzheimer's disease.' (2023)

#setup----

##load packages
library(dplyr)     #v1.1.4
library(tidyr)     #v1.3.0
library(ggplot2)   #v3.4.4
library(phyloseq)  #v1.46.0
library(Boruta)    #v8.0.0
library(caret)     #v6.0-94
library(VIM)       #v6.2.2
library(stringr)   #v1.5.1
library(egg)       #v0.4.5
library(scales)    #v1.3.0
library(reshape2)  #1.4.4
library(ggpubr)    #0.6.0
library(rstatix)   #0.7.2



##general----
##set working directory
# setwd('your directory')

##set seed
set.seed(42)

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



##load metadata----
meta <- read.delim("path/240124_Human_Metadata.csv",
                   sep = ",",
                   row.names = 1)
samples <- row.names(meta)
meta$Sample <- row.names(meta)



##load antibiotic exposure data----
abx <- read.csv('path/240122_Antibiotic_Exposures.csv') %>% 
  filter(Sample %in% samples)
rownames(abx) <- abx$Sample
abx$Sample <- NULL
abx[is.na(abx)] <- 0



##load MetaPhlAn4 output----
##format raw output for phyloseq
m4.raw <- read.delim('path/MetaPhlAn4_raw_output.csv',
                     sep = ",",
                     check.names = F,
                     row.names = 1)

##filter out rows containing strain-level annotations that are new to metaphlan4
m4.raw <- m4.raw %>% filter(!grepl('t__',rownames(.)))

##take only rows that include s__ (species level annotation)
m4.sp <- m4.raw[grepl('s__', rownames(m4.raw)), ]

##filter out species calls below < 0.1% abundance (within sample)
m4.sp.filt <- m4.sp
m4.sp.filt[m4.sp.filt < 0.1] <- 0

##remove rows that sum to zero (lowly abundant taxa that were filtered out in the previous step)
m4.sp.filt <- m4.sp.filt[rowSums(m4.sp.filt[ ]) > 0,]

##renormalize abundance to 100%
m4.sp.filt.renorm <- t((t(m4.sp.filt) / rowSums(t(m4.sp.filt))) * 100)

##create phyloseq taxa table
spp <- data.frame(Names = rownames(m4.sp.filt.renorm))
spp <- data.frame(do.call('rbind', strsplit(as.character(spp$Names),'|',fixed=TRUE)))
rownames(spp) <- rownames(m4.sp.filt.renorm)
colnames(spp) <- c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

##remove prefixes from taxa
spp$Domain <- gsub('k__', '', spp$Domain)
spp$Phylum <- gsub('p__', '', spp$Phylum)
spp$Class <- gsub('c__', '', spp$Class)
spp$Order <- gsub('o__', '', spp$Order)
spp$Family <- gsub('f__', '', spp$Family)
spp$Genus <- gsub('g__', '', spp$Genus)
spp$Species <- gsub('s__', '', spp$Species)

##convert to matrix
spp <- as.matrix(spp)

##make phyloseq object
ps <- phyloseq(otu_table(m4.sp.filt.renorm, taxa_are_rows = TRUE), #1201 taxa
               sample_data(meta),
               tax_table(spp))

##convert relative abundance metrics to integers
ps.int <- transform_sample_counts(ps, function(x) trunc(x*100000))


#feature selection----
##prepare feature data for training----
##metadata - just colonized, primary samples with corresponding toxigenic isolates ( n=52 CDI, n=78 carriage)
meta.rf <- meta %>% 
  filter(Colonized == 1 & Primary_Sample == 1) %>% 
  filter(is.na(Cultured_Cd_Clade) != 'TRUE' & Race != 'unknown') %>% 
  arrange()



##taxa
##access taxonomic abundance data (note transform of matrix)
ps.int.spp <- tax_glom(ps.int, taxrank = 'Species')
ps.int.spp.df <- psmelt(ps.int.spp)

taxa.caret.df <- dcast(ps.int.spp.df, Sample ~ Species, value.var = "Abundance")
row.names(taxa.caret.df) <- taxa.caret.df$Sample
taxa.caret.df$Sample <- NULL

##filter by only colonized samples (and primary = 1 per patient)
taxa.caret.df.col <- taxa.caret.df %>% 
  filter(row.names(.) %in% meta.rf$Sample)
taxa.caret.col <- as.matrix(taxa.caret.df.col)
taxa.caret.col <- taxa.caret.col[order(row.names(taxa.caret.col)), ]

##remove columns that sum to zero
taxa.caret.col <- taxa.caret.col[, colSums(taxa.caret.col) > 0] #removes 160 taxa



##prepare metadata for training----
##append antibiotic exposure data (of prevalence >=10%) to existing metadata
abx2 <- abx
abx2$Sample <- row.names(abx2)
meta.caret <- merge(meta.rf,
                    select(abx2, c('Sample', 'all_abx_30d_binary', 'vancomycin_30d', 'ceph_30d',
                                   'fluoroquinolone_30d', 'carbapenem_30d')))

##define metadata variables to include for training
RF.vars <- c('Group', 'Age', 'Sex', 'Race', 'CdtB', 'all_abx_30d_binary', 'vancomycin_30d', 'ceph_30d',
             'fluoroquinolone_30d', 'carbapenem_30d','Sample')
meta.caret <- meta.caret[, RF.vars]

##update factor levels for ease of interpretation
meta.caret2 <- meta.caret %>% 
  mutate(all_abx_30d_binary = ifelse(all_abx_30d_binary == 1, 'Abx past 30d+', 'Abx past 30d-'),
         vancomycin_30d = ifelse(vancomycin_30d == 1, 'Vancomycin past 30d+', 'Vancomycin past 30d-'),
         ceph_30d = ifelse(ceph_30d == 1, 'Cephalosporin past 30d+', 'Cephalosporin past 30d-'),
         fluoroquinolone_30d = ifelse(fluoroquinolone_30d == 1, 'Fluoroquinolone past 30d+',
                                      'Fluoroquinolone past 30d-'),
         carbapenem_30d = ifelse(carbapenem_30d == 1, 'Carbapenem past 30d+', 'Carbapenem past 30d-'),
         CdtB = ifelse(CdtB == 1, 'cdtB+', 'cdtB-'),
         Group = ifelse(Group == 'CDI', 'CDI', 'Carrier'))
row.names(meta.caret2) <- meta.caret2$Sample
meta.caret2$Sample <- NULL

##update feature names for ease of interpretation
meta.caret2 <- meta.caret2 %>% dplyr::rename(All_Abx_30d = all_abx_30d_binary,
                                             IV_Vanc_30d = vancomycin_30d,
                                             Ceph_30d = ceph_30d,
                                             FQN_30d = fluoroquinolone_30d,
                                             Carb_30d = carbapenem_30d,
                                             cdtB = CdtB)

##convert categorical variables to factors
factor.vars <- c('Group', 'Sex', 'Race', 'cdtB', 'All_Abx_30d', 'IV_Vanc_30d', 'Ceph_30d', 'FQN_30d', 'Carb_30d')
meta.caret2[, factor.vars] <- lapply(meta.caret2[, factor.vars], factor)



#feature selection----

##partition out training cohort for feature selection
set.seed(42)
train_idx <- createDataPartition(meta.caret2$Group, p = 0.6, list=FALSE)

##access metadata and taxonomic abundances for training cohort
meta.train <- meta.caret2[train_idx, ]
taxa.train <- taxa.caret.col[train_idx, ]

##confirm samples match in each data frame for the training cohort
gplots::venn(list(taxa.train = rownames(taxa.train), meta.train = rownames(meta.train)))

##merge sample class identity with feature data
class <- subset(meta.train, select = c('Group'))
taxa.train.wclass <- merge(class, taxa.train, by = 'row.names')
rownames(taxa.train.wclass) <- taxa.train.wclass$Row.names
taxa.train.wclass$Row.names <- NULL



# FUNCTION: callBoruta helper function (called by runFeatureSelection() below)
# Required packages: Boruta, stringr
# Arguments:
#  taxa.data (dataframe)  = taxa.train.wclass (tax abundances for training cohort)
#  seed.Boruta (int)      = will be passed iteratively in defined range
#
# Return:
#  list of names of feature-selected taxa for current iteration

callBoruta <- function(taxa.data, seed.Boruta){
  set.seed(seed.Boruta)
  
  taxa.boruta <- Boruta(Group~., data = taxa.data, maxRuns = 500, doTrace = 0)
  taxa.boruta.fix <- TentativeRoughFix(taxa.boruta)
  
  taxa.boruta.df <- data.frame('boruta' = taxa.boruta.fix$finalDecision)
  
  taxanames <- str_replace_all(rownames(subset(taxa.boruta.df, 
                                               boruta == 'Confirmed')), "`", "")
  
}



# FUNCTION: runFeatureSelection - iterative feature selection function
# Required packages: Boruta, stringr
# Arguments:
#  taxa.train.data (dataframe)   = taxa.train.wclass (passed to callBoruta())
#  seed.range (num list)         = range for iteration, e.g. 1:100.
#
# Return:
#  dataframe summarizing frequency at which unique taxa were feature-selected
#  across all iterations (random seeds).

runFeatureSelection <- function(taxa.train.data, seed.range){
  
  fs.taxa.train <- vector('list', length(seed.range))
  for (i in seed.range){
    fs.taxa.train[[i]] <- callBoruta(taxa.train.data, i)
  }
  
  fs.taxa.train.summ <- data.frame(table(unlist(fs.taxa.train)))
}



##carry out iterative feature selection
fs.taxa.no.cd <- runFeatureSelection(select(taxa.train.wclass, -Clostridioides_difficile), 1:100)

##filter for taxa selected in >=25% of iterations
fs.taxa.no.cd.top25 <- subset(fs.taxa.no.cd, Freq >= 25, select = 'Var1')

##Subset taxa abundance data to these selected taxa (for entire cohort including
##training and validation sets -- will be re-partitioned later using same index)
taxa.caret.boruta.no.cd.top25 <- taxa.caret.col[, colnames(taxa.caret.col) %in% fs.taxa.no.cd.top25$Var1]



#train random forest classifiers----

##define the test harness. Within the training cohort, will train using 10-fold cross-validation. 
cv10 <- trainControl(method = 'cv', number = 10, classProbs = T, savePredictions = T)

##define categorical variables that should not be normalized / scaled.
varsnot2norm <- c('Group', 'Sex', 'Race', 'cdtB', 'All_Abx_30d', 'IV_Vanc_30d', 'Ceph_30d', 'FQN_30d', 'Carb_30d')

##prep data frames
##without c diff >=25%
meta.wtax.no.cd <- merge(meta.caret2, taxa.caret.boruta.no.cd.top25, all = TRUE, by = 'row.names')
rownames(meta.wtax.no.cd) <- meta.wtax.no.cd$Row.names
meta.wtax.no.cd$Row.names <- NULL
meta.wtax.no.cd$Group <- factor(meta.wtax.no.cd$Group,
                                levels = c('CDI', 'Carrier'))



# FUNCTION: train_rf_models trains a classifier on the provided training data  
# and test harness. It does this iteratively (100x) on random 80:20 partitions of
# the training cohort, testing on the entire validation cohort at each iteration. 
# Predictive results are collated.
# Required packages: caret
# Arguments:
#  data (dataframe)         = data subset generated by createSubsets()
#  control.harness          = cv10. Control harness generated by caret::trainControl()
#  data.name.string (str)   = Identifier for the data subset provided
#  varsNot2Norm (chr list)  = List of categorical variables that shouldn't be
#                             normalized / scaled.
#  shuffle.class (logical)  = TRUE if class labels should be shuffled during
#                             model training to generate null performance 
#                             parameter distributions.
#
# Return: 
#  list(prediction results [on training cohort], prediction results [on 
#            validation cohort], variable importances)

train_rf_models <- function(data, control.harness, data.name.string, 
                            varsNot2Norm, shuffle.class) {
  #Cross Validation (within training cohort) 
  out <- list()
  varimportance <- list() 
  
  #Validation Set 
  out.val <- list()
  
  #Separate into train/test and validation subsets, using same index as before.
  set.seed(42)
  data1_idx <- createDataPartition(data$Group, p = 0.6, list = FALSE)
  data1 <- data[data1_idx, ]
  data.val <- data[-data1_idx, ]
  
  for (i in 1:100) {
    #Create random partition of training cohort (80:20). 
    set.seed(i)
    
    train_idx <- createDataPartition(data1$Group, p = 0.8, list = FALSE)
    
    data.train <- data1[train_idx, ]
    data.test <- data1[-train_idx, ]
    
    # Optionally shuffle class labels in training data (data.train).
    if (shuffle.class == TRUE) {
      data.train$Group <- sample(data.train$Group)
    }
    
    #Pre-process data (center and scale)
    preprocessrule <- preProcess(data.train[, !(colnames(data.train) %in% varsNot2Norm)], 
                                 method = c('center', 'scale'))
    data.train.p <- predict(preprocessrule, data.train)
    data.test.p <- predict(preprocessrule, data.test)
    
    data.val.p <- predict(preprocessrule, data.val)
    
    
    #train model
    set.seed(42)
    fit.rf <- train(Group~., data = data.train.p, method = 'rf',
                    metric = 'Accuracy', trControl = control.harness)
    
    #Make predictions for this iteration's test set and store performance measures.
    predictions.rf <- predict(fit.rf, data.test.p)
    cM <- confusionMatrix(predictions.rf, data.test.p$Group,
                          positive = 'CDI')
    pred.results <- c('CDI-CDI' = cM$table[1,1], 
                      'CDI-Carrier' = cM$table[2,1], 
                      'Carrier-CDI' = cM$table[1,2], 
                      'Carrier-Carrier' = cM$table[2,2], 
                      cM$overall, 
                      cM$byClass,
                      'Data' = data.name.string, 
                      'Seed' = i)
    
    out[[i]] <- pred.results
    
    #Make predictions for retained VALIDATION set and store performance measures
    predictions.val.rf <- predict(fit.rf, data.val.p)
    cM.val <- confusionMatrix(predictions.val.rf, data.val.p$Group,
                              positive = 'CDI')
    pred.results.val <- c('CDI-CDI' = cM.val$table[1,1], 
                          'CDI-Carrier' = cM.val$table[2,1], 
                          'Carrier-CDI' = cM.val$table[1,2], 
                          'Carrier-Carrier' = cM.val$table[2,2], 
                          cM.val$overall, 
                          cM.val$byClass,
                          'Data' = data.name.string, 
                          'Seed' = i)
    
    out.val[[i]] <- pred.results.val
    
    
    #Find important vars
    fit.rf.importance <- varImp(fit.rf, scale=FALSE)
    fit.rf.importance.df <- fit.rf.importance$importance
    var.importance <- fit.rf.importance.df$Overall
    names(var.importance) <- row.names(fit.rf.importance.df)
    
    varimportance[[i]] <- var.importance
  }
  
  #Return 
  out.df <- data.frame(do.call('rbind', out))
  out.val.df <- data.frame(do.call('rbind', out.val))
  varimportance.df <- data.frame(do.call('rbind', varimportance))
  
  allout <- list('Pred.Results.CV'=out.df, 
                 'Pred.Results.Val'=out.val.df, 
                 'Var.Importance'=varimportance.df)
}



##train models (without shuffling class labels)
##without c diff
rf.no.cd <- train_rf_models(meta.wtax.no.cd, cv10, 'No C diff',
                            varsnot2norm, FALSE)



##train model for null distributions (shuffle class labels)
##without c diff
rfN.no.cd <- train_rf_models(meta.wtax.no.cd, cv10, 'No C diff',
                             varsnot2norm, TRUE)



##collate results----
##test
all.pred.CV <- rf.no.cd$Pred.Results.CV

all.pred.Val <- rf.no.cd$Pred.Results.Val

no.cd.VarImp <- list(rf.no.cd$Var.Importance) %>% 
  Reduce(function(d1, d2) full_join(d1, d2), .)


##null
NULL.all.pred.CV <- rfN.no.cd$Pred.Results.CV

NULL.all.pred.Val <- rfN.no.cd$Pred.Results.Val

NULL.no.cd.VarImp <- list(rfN.no.cd$Var.Importance) %>% 
  Reduce(function(d1, d2) full_join(d1, d2), .)



##evaluate predictive performance metrics----
##gather performance metrics for plotting
all.pred.Val2 <- all.pred.Val %>%
  gather('Performance_Measure', 
         'Performance_Measure_Value', c('Accuracy', 'Sensitivity', 'Specificity'))

##determine the order of mean accuracy amongst models
accuracy_order <- all.pred.Val2 %>%
  group_by(Data) %>%
  filter(Performance_Measure == 'Accuracy') %>%
  summarise(mean(as.numeric(Performance_Measure_Value))) %>%
  arrange(desc(`mean(as.numeric(Performance_Measure_Value))`)) %>% 
  select(Data)

##order factor by mean accuracy (decreasing)
all.pred.Val2$Data <- factor(all.pred.Val2$Data,
                             levels = accuracy_order$Data)

##plot
all.pred.Val2 %>% 
  ggplot(aes(x = Data, y = as.numeric(Performance_Measure_Value))) +
  geom_boxplot() +
  facet_wrap(~ Performance_Measure) +
  coord_flip() +
  theme_plot() #+
  # stat_compare_means(label = 'p.signif')



##assess feature importance----
##summarise variable importance
##test
no.cd.VarImp.mean <- no.cd.VarImp %>% 
  summarise_all(mean)

no.cd.VarImp.sd <- no.cd.VarImp %>% 
  summarise_all(sd)

##null
NULL.no.cd.VarImp.mean <- NULL.no.cd.VarImp %>% 
  summarise_all(mean)

NULL.no.cd.VarImp.sd <- NULL.no.cd.VarImp %>% 
  summarise_all(sd)

##convert to long and merge into a single data frame
##test
no.cd.VarImp.mean.long <- gather(no.cd.VarImp.mean, Feature, Feature.MeanImp,
                                 Age:Veillonella_parvula,
                                 factor_key = T)

no.cd.VarImp.sd.long <- gather(no.cd.VarImp.sd, Feature, Feature.SDImp,
                               Age:Veillonella_parvula,
                               factor_key = T)

no.cd.VarImp.summ <- left_join(no.cd.VarImp.mean.long,
                               no.cd.VarImp.sd.long,
                               by = 'Feature') %>% 
  dplyr::mutate(distribution = 'TEST')

##null
NULL.no.cd.VarImp.mean.long <- gather(NULL.no.cd.VarImp.mean, Feature, Feature.NULLMeanImp,
                                      Age:Veillonella_parvula,
                                      factor_key = T)

NULL.no.cd.VarImp.sd.long <- gather(NULL.no.cd.VarImp.sd, Feature, Feature.NULLSDImp,
                                    Age:Veillonella_parvula,
                                    factor_key = T)

NULL.no.cd.VarImp.summ <- left_join(NULL.no.cd.VarImp.mean.long,
                                    NULL.no.cd.VarImp.sd.long,
                                    by = 'Feature') %>% 
  dplyr::mutate(distribution = 'NULL')

##combine
no.cd.VarImp.summ <- merge(no.cd.VarImp.mean.long, no.cd.VarImp.sd.long)
no.cd.VarImp.summ <- merge(no.cd.VarImp.summ, NULL.no.cd.VarImp.mean.long)
no.cd.VarImp.summ <- merge(no.cd.VarImp.summ, NULL.no.cd.VarImp.sd.long)

##define feature labels
varimp.labels <- c('Ceph_30dCephalosporin.past.30d.' = 'Cephalosporin (30d)',
                   'Clostridium_innocuum' = 'Clostridium innocuum',
                   'Blautia_wexlerae' = 'Blautia wexlerae',
                   'Veillonella_parvula' = 'Veillonella parvula',
                   'Fusobacterium_nucleatum' = 'Fusobacterium nucleatum',
                   'Enterocloster_clostridioformis' = 'Enterocloster clostridioformis',
                   'Sexmale' = 'Sex',
                   'All_Abx_30dAbx.past.30d.' = 'Any antibiotics (30d)',
                   'cdtBcdtB.' = 'cdtB',
                   'Racewhite' = 'Race',
                   'Carb_30dCarbapenem.past.30d.' = 'Carbapenem (30d)',
                   'IV_Vanc_30dVancomycin.past.30d.' = 'IV Vancomycin (30d)',
                   'FQN_30dFluoroquinolone.past.30d.' = 'Fluoroquinolone (30d)')

##plot - filter out Race with Pacific Islander as a level because it has zero variance (one participant)
p.varimp <- no.cd.VarImp.summ %>%
  filter(Feature != 'Racepacific.islander') %>% 
  ggplot(aes(x=reorder(Feature, Feature.MeanImp), y=Feature.MeanImp))+
  geom_col(fill = 'grey')+
  geom_point(color='black')+
  geom_errorbar(aes(ymin = Feature.MeanImp - Feature.SDImp, 
                    ymax = Feature.MeanImp + Feature.SDImp), width=0.7, size=0.8)+
  geom_point(aes(x=reorder(Feature, Feature.MeanImp), y=Feature.NULLMeanImp), 
             color='maroon', alpha=0.80)+
  geom_errorbar(aes(ymin = Feature.NULLMeanImp - Feature.NULLSDImp, 
                    ymax = Feature.NULLMeanImp + Feature.NULLSDImp), width=0.5, size=0.6, 
                color='maroon', alpha=0.80)+
  theme_classic()+
  theme(axis.text.y = element_text(size=10),
        axis.text.x = element_text(size=10),
        axis.title.x = element_text(size=10),
        axis.title.y = element_blank(),
        plot.title = element_blank(),
        panel.border = element_rect(fill=NA, colour="black", size=0.5))+
  #facet_wrap(~DataSet, nrow=2)+
  scale_x_discrete(labels = varimp.labels)+
  coord_flip()



##check for significance----
imp.features <- as.character(filter(no.cd.VarImp.summ, Feature.MeanImp > Feature.NULLMeanImp)$Feature)

no.cd.VarImp.comb <- rbind(mutate(no.cd.VarImp, run = 'empirical'),
                           mutate(NULL.no.cd.VarImp, run = 'null'))

no.cd.VarImp.comb.long <- melt(no.cd.VarImp.comb) %>% 
  filter(variable %in% imp.features) %>% 
  dplyr::rename(Feature.Importance = value,
                Feature = variable)

stat.test.imp <- no.cd.VarImp.comb.long %>% 
  mutate(group = str_c(run, '_', Feature)) %>% 
  wilcox_test(Feature.Importance ~ group) %>% 
  adjust_pvalue(method = 'bonferroni') %>% 
  add_significance()

stat.test.imp.filt <- filter(stat.test.imp,
                             group1 %in% (str_c('empirical', '_', imp.features)) &
                               group2 %in% (str_c('null', '_', imp.features)))
