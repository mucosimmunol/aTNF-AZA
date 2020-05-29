### CEDTherapy complete script

####---- FIG1: Responder Statistics + Generate Metadata ----####

library(tidyverse)
library(readxl)
library(lubridate)
library(data.table)
library(ggpubr)
library(forcats)
library(cowplot)
library(Maaslin)
library(phyloseq)
library(vegan)
library(DESeq2)
library(lme4)
library(Maaslin)

####---- Fig.1 A Response to Therapy ----####
responderstat <- readRDS("metadata/responder_stat.rds")
fig1A <- ggbarplot(responderstat, "disther", "value",fill="Responder", 
                   color = "Responder",label = T,lab.col="white",lab.nb.digits = 2,
                   lab.pos = "in",palette="Paired", xlab="", ylab="Number of patients")+
  theme_pubr(base_size=8)


####---- Fig.1 B/C Clinical Baselinecharestics ----####


#read data
all <- readRDS("metadata/sample_metadata.rds")
all<- all[,2:28]

x <- dplyr::filter(all) %>%
  select(DIS, Case, CDAI,MayoScore, CRP, Leuko, Ferritin, Calprotectin, TP, Therapy, Responder=Responder_predes) %>%
  filter(TP<=2) %>%
  mutate(Responder=factor(Responder, labels = c("Non-Responder", "Deep Remission")))%>%
  gather(key = variable, value=value, CDAI:Calprotectin) %>%
  filter(!(Case %in% c(49,56))) %>%
  arrange(Case)
x$TP <- factor(x$TP, labels = c("baseline", "timepoint 1", "timepoint 2"))
x$Therapy <- fct_recode(x$Therapy, AZA="aza", 'anti-TNF'="biol")

# Crohn's Disease
p_calpro <- ggplot(filter(x, variable=="Calprotectin", DIS=="mc"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 0.8)+
  theme_pubr(base_size=9, legend="none")+
  labs(x = "", y = "Calprotectin [µg/g]")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_sqrt(expand = c(0, 0), limits = c(0, NA))


p_crp <- ggplot(filter(x, variable=="CRP", DIS=="mc"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 0.8)+
  theme_pubr(base_size=9, legend="none")+
  labs(x = "", y = "CRP [mg/dL]")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_ferritin <- ggplot(filter(x, variable=="Ferritin", DIS=="mc"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 0.8)+
  theme_pubr(base_size=9, legend="none")+
  labs(x = "", y = "Ferritin [�g/L]")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_leuco <- ggplot(filter(x, variable=="Leuko", DIS=="mc"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 0.8)+
  theme_pubr(base_size=9, legend="none")+
  labs(x = "", y = "Leucocytes [G/L]")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

p_cdai <- ggplot(filter(x, variable=="CDAI", DIS=="mc"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 0.8)+
  theme_pubr(base_size=9, legend="none")+
  labs(x = "", y = "CDAI")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))

# ulcerative colitis
p_calpro2 <- ggplot(filter(x, variable=="Calprotectin", DIS=="cu"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 1)+
  theme_pubr(base_size=9, legend="none")+
  labs(x = "", y = "Calprotectin [µg/g]")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_sqrt(expand = c(0, 0), limits = c(0, NA))


p_crp2 <- ggplot(filter(x, variable=="CRP", DIS=="cu"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 0.8)+
  theme_pubr(base_size=9, legend="none")+
  labs(x = "", y = "CRP [mg/dL]")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_ferritin2 <- ggplot(filter(x, variable=="Ferritin", DIS=="cu"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 0.8)+
  theme_pubr(base_size=9, legend="none")+
  labs(x = "", y = "Ferritin [µg/L]")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

p_leuco2 <- ggplot(filter(x, variable=="Leuko", DIS=="cu"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 0.8)+
  theme_pubr(base_size=9, legend="none")+
  labs(x = "", y = "Leucocytes [G/L]")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))


p_mayo <- ggplot(filter(x, variable=="MayoScore", DIS=="cu"), aes(x = factor(TP), y = value))+
  geom_point(aes(shape = Therapy),  size = 2)+
  geom_line(aes(group = Case, color = Responder), size = 0.8)+
  theme_pubr(base_size=9)+
  labs(x = "", y = "Mayo Score")+
  scale_color_brewer(palette="Pastel1")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))


####---- FIG2: Baseline-charakteristics diversity ---####
meta <- all
meta$IMGM <- sapply(str_split(meta$IMGM, "-"), `[`, 2)
rownames(meta) <- meta$IMGM
meta$Responder <- factor(meta$Responder_predes, labels = c("Non-Remission","Remission"))
meta$TP <- fct_recode(factor(meta$TP), baseline = "0", "timepoint 1" = "1", "timepoint 2" = "2", "exclude" = "3")
otu <- readRDS("16S/seqtab_nochim.rds")
tax <- readRDS("16S/species.rds")

phylo <- phyloseq(otu_table(otu, taxa_are_rows = F),
                  tax_table(tax),
                  sample_data(meta))
sample_names(phylo) <- paste0("sample", sample_names(phylo))
phylo <- prune_samples(sample_data(phylo)$TP!=99, phylo)

#Filter phyloseq
phylo <- prune_samples(sample_sums(phylo) > 1000, phylo)
phylo.ns <- prune_taxa(taxa_sums(phylo) > 1, phylo) 
phylo.filt <- filter_taxa(phylo.ns, function(x) sum(x>8) > (0.1 * length(x)), TRUE)

# filter for baseline
phylo.bl <- subset_samples(phylo, TP =="baseline")
phylo.filt.bl <- subset_samples(phylo.filt, TP=="baseline")
phylo.rel <- transform_sample_counts(phylo.filt.bl, function(x) x/sum(x)*100)

# a diversity @ baseline, compare Disease groups
adat <- plot_richness(phylo.bl, measures=c("Chao1", "Shannon", "InvSimpson"), x = "DIS")
adat$data$DIS <- factor(adat$data$DIS, labels = c("UC", "CD"))
p <- ggplot(filter(adat$data, variable=="Shannon"),  aes(x = DIS, y = value))+geom_point()+
  geom_boxplot(alpha=0.4, fill = "darkgrey")+theme_pubr(base_size=9)+
  labs(x = "", y = "Shannon index")+
  lims(y = c(0,5))

# Test for significance using t-tests
t.test(value~DIS, data = filter(adat$data, variable=="Shannon"))


# a-diversity @ baseline between future responders/non-responders (within diseases)
adat2 <- adat$data
adat2$Disr <- paste0(adat2$DIS, "_", adat2$Responder)

p2 <- ggplot(filter(adat2, variable=="Shannon"), aes(x = DIS, y = value))+
  geom_point(aes(color=Responder), position=position_dodge(width = 0.7))+
  geom_boxplot(alpha=0.4, aes(fill = Responder))+theme_pubr(legend="right", base_size=9)+
  labs(x = "", y = "Shannon index")+
  lims(y = c(0,5))

# Test for significance using ANOVA
summary(aov(value~DIS+Responder+Therapy, data = filter(adat2, variable=="Shannon")))


# b-diversity pcoa plots of Bray Curtis by Disease and prospective response to Therapy
sample_data(phylo.rel)$DIS <- factor(sample_data(phylo.rel)$DIS, labels = c("UC", "CD"))

prel2 <- tax_glom(phylo.rel, taxrank = "Genus")
DistBC <- phyloseq::distance(prel2, method="bray")
ord <- ordinate(prel2, method = "PCoA", distance = DistBC)
BCplot <- plot_ordination(prel2, ord)
BCplot$data$Therapy <- factor(BCplot$data$Therapy, labels = c("AZA", "anti-TNF"))

pcoa_bray <- ggplot(BCplot$data, aes(x = Axis.1, y = Axis.2))+
  geom_point(aes(color=Responder), size = 2)+
  ggtitle("PCoA: Bray-Curtis (Genus)")+
  facet_grid(DIS~Therapy)+
  labs(x = BCplot$labels$x, y = BCplot$labels$y)+
  theme_bw(base_size=9)+
  theme(aspect.ratio = 15.4/42.6)

# PERMANOVA
df = as(sample_data(prel2), "data.frame")
df$Therapy <- as.factor(df$Therapy)
d = phyloseq::distance(prel2, "bray")
adonis_genus = adonis(d ~ Responder + DIS + Therapy, df)
adonis_genus

DistBCASV <- phyloseq::distance(phylo.rel, method="bray")
ordASV <- ordinate(phylo.rel, method = "PCoA", distance = DistBCASV)

BCplotASV <- plot_ordination(phylo.rel, ordASV)
BCplotASV$data$Therapy <- factor(BCplotASV$data$Therapy, labels = c("AZA", "anti-TNF"))

pcoa_bray_ASV <- ggplot(BCplotASV$data, aes(x = Axis.1, y = Axis.2))+
  geom_point(aes(color=Responder), size = 2)+
  ggtitle("PCoA: Bray-Curtis (ASVs)")+
  facet_grid(DIS~Therapy)+
  labs(x = BCplot$labels$x, y = BCplot$labels$y)+
  theme_bw(base_size=9)+
  theme(aspect.ratio = 15.4/42.6)

# PERMANOVA
df = as(sample_data(phylo.rel), "data.frame")
df$Therapy <- as.factor(df$Therapy)
d = phyloseq::distance(phylo.rel, "bray")
adonis_asv = adonis(d ~ Responder + DIS + Therapy, df)
adonis_asv



####---- FIG3: Longitudinal analysis ---####
# taxplots on phylum level
psp <- tax_glom(phylo.ns, "Phylum")
pspr <- transform_sample_counts(psp, function(x) x/sum(x)*100)
pspr <- filter_taxa(pspr, function(x) sum(x>3) > (0.1 * length(x)), TRUE)
psp2 <- psmelt(pspr)

psp3 <- reshape2::dcast(psp2, Case+DIS+Therapy+Responder+TP~Phylum, value.var="Abundance", fun.aggregate = median)
phyla <- names(psp3)[6:9]

mc <- filter(psp3, DIS=="mc")
mc_aza <- filter(mc, Therapy=="aza")
mc_biol <- filter(mc, Therapy=="biol")

cu <- filter(psp3, DIS=="cu")
cu_aza <- filter(cu, Therapy=="aza")
cu_biol <- filter(cu, Therapy=="biol")

lmer.comp <- function(df, dis, ther) {
  r <- list()
  for (pn in phyla) {
    print(pn)
    mixed.lmer <- lmer(get(pn)~TP+Responder+(1|Case), data=df, REML = F)
    p<-plot(fitted(mixed.lmer), residuals(mixed.lmer, type = "pearson"))
    sing <- isSingular(mixed.lmer)
    null.lmer <- lmer(get(pn)~TP+(1|Case), data=df, REML=F)
    t<-anova(mixed.lmer, null.lmer)
    x<-summary(mixed.lmer)
    r[[pn]] <- c(paste0(dis, "_", ther), t$`Pr(>Chisq)`[2], x$coefficients[3,1], sing)
  }
  r2 <- as.data.frame(r)
  colnames(r2) <- phyla
  rownames(r2) <- c("dataset", "pval_chi", "Coefficient_Responder", "Singular")
  r2 <- t(r2) %>% as.data.frame()
  r2 <- rownames_to_column(r2, var = "Phylum")
  return(r2)
}

res_mcaza <- lmer.comp(mc_aza, "mc", "aza")
res_mcaza
res_mctnf <- lmer.comp(mc_biol, "mc", "biol") 
res_mctnf
res_cuaza <- lmer.comp(cu_aza, "cu", "aza") 
res_cuaza
res_cutnf <- lmer.comp(cu_biol, "cu", "biol") 
res_cutnf
res <- rbind(res_mcaza, res_mctnf, res_cuaza, res_cutnf)
res$pval_chi <- as.character(res$pval_chi) %>% as.numeric()
res$Coefficient_Responder <- as.character(res$Coefficient_Responder) %>% as.numeric()

res <- res[res$pval_chi < 0.05,]
res <- res[res$Singular==FALSE,]
res

####---- Phylum composition ----####
#mit Differenzierung nach Erkrankung
ps.phylum.dis <- pspr
ps.phylum.dis <- subset_samples(ps.phylum.dis, TP !="exclude")
sample_data(ps.phylum.dis)$DIS <- factor(sample_data(ps.phylum.dis)$DIS)
sample_data(ps.phylum.dis)$Therapy <- factor(sample_data(ps.phylum.dis)$Therapy)
sample_data(ps.phylum.dis)$TP <- factor(sample_data(ps.phylum.dis)$TP)
sample_data(ps.phylum.dis)$Responder <- factor(sample_data(ps.phylum.dis)$Responder)
sample_data(ps.phylum.dis)$Group1 <-factor(paste0(sample_data(ps.phylum.dis)$DIS, "_", 
                                                  sample_data(ps.phylum.dis)$Therapy, "_", 
                                                  sample_data(ps.phylum.dis)$TP, "_",
                                                  sample_data(ps.phylum.dis)$Responder))
ps.phylum.dis <- merge_samples(ps.phylum.dis, sample_data(ps.phylum.dis)$Group1)
ps.phylum.dis <- transform_sample_counts(ps.phylum.dis, function(x) x/sum(x)*100)

topgenus = names(sort(taxa_sums(ps.phylum.dis), decreasing = TRUE)[1:4])
ps.phylum.dis <- prune_taxa(topgenus, ps.phylum.dis)
p.phyla2 <- plot_bar(ps.phylum.dis, x="Sample", y="Abundance")
p.phyla2$data$DIS <- factor(p.phyla2$data$DIS, labels = c("cu", "mc"))
p.phyla2$data$TP <- factor(p.phyla2$data$TP, labels = c("baseline", "timepoint 1", "timepoint 2")) 
              
pdat2 <- p.phyla2$data %>%
  dplyr::select(Sample,Phylum, Abundance,DIS, TP, Therapy, Responder) %>%
  dplyr::mutate(Therapy = factor(Therapy, labels = c("Azathioprine", "anti-TNFa")),
                Responder = factor(Responder, labels = c("Non-Responder", "Deep Remission")))

p.phyla2cu <- ggplot(filter(pdat2, DIS=="cu"), aes(x=TP, y=Abundance))+geom_col(aes(fill = Phylum))+facet_grid(Therapy~Responder)+
  theme_pubr(base_size=9.6, legend = "right")+
  labs(x ="", y = "Relative Abundance (%)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p.phyla2cu

p.phyla2mc <- ggplot(filter(pdat2, DIS=="mc"), aes(x=TP, y=Abundance))+geom_col(aes(fill = Phylum))+facet_grid(Therapy~Responder)+
  theme_pubr(base_size=9.6, legend = "right")+
  labs(x ="", y = "Relative Abundance (%)")+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p.phyla2mc

####---- qiime2 - adaptor for longitudinal analysis----####
# Export data for qiime2
# Preparations
phylo <- phyloseq(otu_table(otu, taxa_are_rows = F),
                  tax_table(tax),
                  sample_data(meta))

tax<-as(tax_table(phylo),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL


# Export feature/OTU table
# As a biom file
sample_names(phylo) <- paste0("sam", sample_names(phylo))
library(biomformat);packageVersion("biomformat")

otu<-t(as(otu_table(phylo),"matrix"))
otu_biom<-make_biom(data=otu)

tax<-as(tax_table(phylo),"matrix")
tax_cols <- colnames(tax)
tax<-as.data.frame(tax)
tax$taxonomy<-do.call(paste, c(tax[tax_cols], sep=";"))
for(co in tax_cols) tax[co]<-NULL

smd <- sample_data(phylo)
sampleid <- rownames(smd)
smd2 <- cbind(sampleid, smd)
smd2 <- dplyr::select(smd2)
smd2$TherResp <- paste0(smd2$Therapy, smd2$Responder_predes, sep = "_")

# read from qiime2 output
mc_shannon <- readxl::read_excel("figures/fig3_qiime2-data.xlsx", sheet="MC_alpha") %>%
  mutate(Responder=factor(Responder, labels=c("Non-Remission","Remission")),
         Therapy = fct_rev(factor(Therapy, labels = c("anti-TNFa", "AZA"))))


mc_BC <- readxl::read_excel("figures/fig3_qiime2-data.xlsx", sheet="MC_BC") %>%
  mutate(Responder=factor(Responder, labels=c("Non-Remission","Remission")),
         Therapy = fct_rev(factor(Therapy, labels = c("anti-TNFa", "AZA"))))

cu_shannon <- readxl::read_excel("figures/fig3_qiime2-data.xlsx", sheet="CU_alpha") %>%
  mutate(Responder=factor(Responder, labels=c("Non-Remission","Remission")),
         Therapy = fct_rev(factor(Therapy, labels = c("anti-TNFa", "AZA"))))


cu_BC <- readxl::read_excel("figures/fig3_qiime2-data.xlsx", sheet="CU_BC") %>%
  mutate(Responder=factor(Responder, labels=c("Non-Remission","Remission")),
         Therapy = fct_rev(factor(Therapy, labels = c("anti-TNFa", "AZA"))))



# generate figures and combine into Fig. 3
mc_shannon_plot<-ggplot(mc_shannon, aes(x = Therapy, y=Difference,fill=Responder))+
  geom_boxplot()+
  labs(y = "Pairwise difference in Shannon index",x="")+
  scale_fill_brewer(type = "qual", palette = 3)+
  theme_pubr(base_size=9.6, legend = "right")
cu_shannon_plot<-ggplot(cu_shannon, aes(x = Therapy, y=Difference,fill=Responder))+
  geom_boxplot()+
  labs(y = "Pairwise difference in Shannon index",x="")+
  scale_fill_brewer(type = "qual", palette = 3)+
  theme_pubr(base_size=9.6, legend = "right")

mc_bc_plot<-ggplot(mc_BC, aes(x = Therapy, y=Distance,fill=Responder))+
  geom_boxplot()+
  labs(y = "Pairwise Bray-Curtis distance",x="")+
  scale_fill_brewer(type = "qual", palette = 3)+
  theme_pubr(base_size=9.6, legend = "right")
cu_bc_plot<-ggplot(cu_BC, aes(x = Therapy, y=Distance,fill=Responder))+
  geom_boxplot()+
  labs(y = "Pairwise Bray-Curtis distance",x="")+
  scale_fill_brewer(type = "qual", palette = 3)+
  theme_pubr(base_size=9.6, legend = "right")

fig3_1 <- cowplot::plot_grid(mc_shannon_plot, mc_bc_plot, align= "hv", ncol = 2, labels =c("A","B"))
fig3_2<- cowplot::plot_grid(p.phyla2mc,NULL, align= "hv", ncol = 2, labels =c("C",NULL))

fig3_3 <- cowplot::plot_grid(cu_shannon_plot, cu_bc_plot, align= "hv", ncol = 2, labels =c("D","E"))
fig3_4<- cowplot::plot_grid(p.phyla2cu,NULL, align= "hv", ncol = 2, labels =c("F",NULL))

fig3<-cowplot::plot_grid(fig3_1, fig3_2,fig3_3,fig3_4, ncol = 1, rel_heights = c(1,1.3,1,1.3))


####---- Figure 4: Maaslin for Genus Abundance ----####

ps_genus <- tax_glom(phylo, taxrank = "Genus")
ps_genus <- transform_sample_counts(ps_genus, function(x) x/sum(x))

formatMaaslin <- function(phyloseqobj, taxlevel) {
  cat("Transform Sample counts \n")
  # clean up phyloseqobject
  phyloseqobj <- transform_sample_counts(phyloseqobj, function(x) x / sum(x)*100)
  sample_data(phyloseqobj)$TP <- as.numeric(factor(sample_data(phyloseqobj)$TP, labels = c(1,2,3,4,99)))
  phyloseqobj <- prune_samples(sample_data(phyloseqobj)$TP <4, phyloseqobj)
  sample_data(phyloseqobj)$Responder <- factor(sample_data(phyloseqobj)$Responder, labels = c("NR", "R"))
  sample_data(phyloseqobj)$TPR <- paste0(sample_data(phyloseqobj)$TP, sample_data(phyloseqobj)$Responder)
  
  
  if (taxlevel!="Genus") {
    cat(paste("summarizing at:", taxlevel, "\n"))
    phyloseqobj <- tax_glom(phyloseqobj, taxrank = taxlevel, NArm=T)
  }
  
  phyloseqobjCD <- subset_samples(phyloseqobj, DIS =="mc")
  phyloseqobjUC <- subset_samples(phyloseqobj, DIS =="cu")
  
  #extract otu_table from phyloseq object and shift rownames = OTUs into separate column for merging
  dat.sharedCD <- as.data.frame(t(otu_table(phyloseqobjCD)))
  otusCD <- rownames(dat.sharedCD)
  dat.sharedCD$OTU <- otusCD
  
  #extract otu_table from phyloseq object and shift rownames = OTUs into separate column for merging
  dat.sharedUC <- as.data.frame(t(otu_table(phyloseqobjUC)))
  otusUC <- rownames(dat.sharedUC)
  dat.sharedUC$OTU <- otusUC
  
  
  #extract and format metadata and taxonomy
  cat("Extract Metadata and Taxonomy \n")
  
  metadatCD <- sample_data(phyloseqobjCD)
  metadatCD <- select(metadatCD,
                      Case, TP, Responder_predes)
  mCD<-t(metadatCD[1:nrow(metadatCD),])
  sample <- colnames(mCD)
  mCD<-rbind(sample, mCD)
  mCD <- as.data.frame(mCD)
  mCD <- rownames_to_column(mCD,var ="metavar")
  taxonomyCD <-as.data.frame(phyloseqobjCD@tax_table@.Data)
  
  metadatUC <- sample_data(phyloseqobjUC)
  metadatUC <- select(metadatUC,
                      Case, TP, Responder_predes)
  mUC<-t(metadatUC[1:nrow(metadatUC),])
  sample <- colnames(mUC)
  mUC<-rbind(sample, mUC)
  mUC <- as.data.frame(mUC)
  mUC <- rownames_to_column(mUC,var ="metavar")
  taxonomyUC <-as.data.frame(phyloseqobjUC@tax_table@.Data)
  
  
  # Merge abundance and taxonomy data
  cat("Merge abundance and taxonomy data \n")
  if (taxlevel == "Genus") {
    otutax <- paste(taxonomyCD$Phylum, taxonomyCD$Class, taxonomyCD$Order, taxonomyCD$Family, taxonomyCD$Genus, sep = "|")
  }
  if (taxlevel == "Family") {
    otutax <- paste(taxonomyCD$Phylum, taxonomyCD$Class, taxonomyCD$Order, taxonomyCD$Family, sep = "|")
  }
  if (taxlevel == "Order") {
    otutax <- paste(taxonomyCD$Phylum, taxonomyCD$Class, taxonomyCD$Order, sep = "|")
  }
  if (taxlevel == "Class") {
    otutax <- paste(taxonomyCD$Phylum, taxonomyCD$Class, sep = "|")
  }
  if (taxlevel == "Phylum") {
    otutax <- paste(taxonomyCD$Phylum, sep = "|")
  }
  
  otutax <- as.data.frame(otutax)
  otutax$OTU <- rownames(taxonomyCD)
  otutax$taxonomy <- gsub("-", "_", otutax$otutax)
  otutax <- select(otutax, -otutax)
  
  
  abundancesCD <- base::merge(otutax, dat.sharedCD, by = "OTU") %>%
    dplyr::select(-OTU)
  abundancesUC <- base::merge(otutax, dat.sharedUC, by = "OTU") %>%
    dplyr::select(-OTU)
  
  colnames(abundancesCD) <- NULL
  colnames(abundancesUC) <- NULL
  
  # output pcl files
  # Crohns
  cat("Write Results \n")
  write.table(mCD, paste0("maaslin-datasets/", taxlevel, "MC.prepare.pcl"), sep="\t", row.names = F, col.names = F,quote = F)
  write.table(abundancesCD, paste0("maaslin-datasets/", taxlevel, "MC.prepare.pcl"), sep="\t", append = T, row.names = F, col.names = F, quote=F)
  
  x <- read.table(file =  paste0("maaslin-datasets/", taxlevel, "MC.prepare.pcl"), sep = "\t", header = F)
  write.table(x, file =  paste0("maaslin-datasets/", taxlevel, "_MC.pcl"), sep="\t", row.names = F, col.names = F, quote = F)
  
  #Ulcerative Colitis
  write.table(mUC, paste0("maaslin-datasets/", taxlevel, "CU.prepare.pcl"), sep="\t", row.names = F, col.names = F,quote = F)
  write.table(abundancesUC, paste0("maaslin-datasets/", taxlevel, "CU.prepare.pcl"), sep="\t", append = T, row.names = F, col.names = F, quote=F)
  
  x <- read.table(file =  paste0("maaslin-datasets/", taxlevel, "CU.prepare.pcl"), sep = "\t", header = F)
  write.table(x, file =  paste0("maaslin-datasets/", taxlevel, "_CU.pcl"), sep="\t", row.names = F, col.names = F, quote = F)
  
  #Generate read config
  ## CD
  x <- last(mCD$metavar)
  z <- abundancesCD[1,1]
  write(paste0("Matrix: Metadata\nRead_PCL_Rows: -",x,"\n","\nMatrix: Abundance\nRead_PCL_Rows: ", z, "-"),
        file = paste0("maaslin-datasets/", taxlevel, "_MC.read.config"))
  
  ## Ulcerative colitis
  x <- last(mUC$metavar)
  z <- abundancesUC[1,1]
  write(paste0("Matrix: Metadata\nRead_PCL_Rows: -",x,"\n","\nMatrix: Abundance\nRead_PCL_Rows: ", z, "-"),
        file = paste0("maaslin-datasets/", taxlevel, "_UC.read.config"))
}


setwd("maaslin_results")

formatMaaslin(ps_genus, "Genus")
formatMaaslin(ps_genus, "Phylum")

#remove intermediary files
preparefiles <- list.files("maaslin-datasets")
preparefiles <- preparefiles[grepl("prepare", preparefiles)]
preparefiles <- paste0("maaslin-datasets/", preparefiles)
file.remove(preparefiles)
rm(preparefiles)

Maaslin("maaslin-datasets/Genus_MC.pcl", "maaslin_Genus_MC",strInputConfig="maaslin-datasets/Genus_MC.read.config", dSignificanceLevel = 0.05,
        dMinAbd=0.02, dMinSamp=0.1, strRandomCovariates=c("Case"))

Maaslin("maaslin-datasets/Genus_CU.pcl", "maaslin_Genus_CU",strInputConfig="maaslin-datasets/Genus_UC.read.config", dSignificanceLevel = 0.05,
        dMinAbd=0.02, dMinSamp=0.1, strRandomCovariates=c("Case"))

Maaslin("maaslin-datasets/Phylum_MC.pcl", "maaslin_Phylum_MC",strInputConfig="maaslin-datasets/Phylum_MC.read.config", dSignificanceLevel = 0.05,
        dMinAbd=0.02, dMinSamp=0.1, strRandomCovariates=c("Case"))

Maaslin("maaslin-datasets/Phylum_CU.pcl", "maaslin_Phylum_CU",strInputConfig="maaslin-datasets/Phylum_UC.read.config", dSignificanceLevel = 0.05,
        dMinAbd=0.02, dMinSamp=0.1, strRandomCovariates=c("Case"))



### FIG4 Maaslin: Results from Maaslin
# read in taxa list from maaslin
setwd("maaslin-datasets/") # use this line when analysing no_calpro data
cu_responders <- read.delim("maaslin_Genus_CU/Genus_CU-Responder.txt") %>%
  filter(P.value < 0.05)
mc_responders <- read.delim("maaslin_Genus_MC/Genus_MC-Responder.txt") %>%
  filter(P.value < 0.05)

# read in phyloseq object used in maaslin
ps_rel <- transform_sample_counts(ps_genus, function(x) x/sum(x)*100)
keeptax_cu <- paste(cu_responders$Feature)
keeptax_mc <- paste(mc_responders$Feature)

getkeep <- function(keepvector) {
  keepvector <- unique(keepvector) %>% as.data.frame()
  names(keepvector) <- "keep"
  keepvector <- tidyr::separate(data = keepvector,
                                col="keep", 
                                into = c("Phylum", "Order", "Class", "Family", "Genus"),
                                sep = "\\|")
  keepvector
}
keeptax_cu <- getkeep(keeptax_cu)
keeptax_mc <- getkeep(keeptax_mc)
keeptax <- (unique(c(keeptax_cu$Genus, keeptax_mc$Genus)))

pdat <- psmelt(ps_rel)
pdat2 <- filter(pdat, Genus %in% keeptax)

pdat2 <- filter(pdat2, TP !="exclude")
pdat2$Therapy <- factor(pdat2$Therapy, labels = c("Azathioprine", "anti-TNFa"))
plot_cu1 <- ggplot(filter(pdat2, DIS=="cu", Genus==keeptax_cu$Genus[1]), aes(x = TP, y = Abundance, fill = Responder))+
  geom_point(position = position_dodge(width =0.75), aes(color=Responder))+
  facet_wrap(~Therapy)+
  geom_boxplot(position = position_dodge(preserve = "single"), alpha = 0.5)+
  theme_bw()+
  labs(x = "", title = keeptax_cu$Genus[1])+
  scale_y_sqrt()+
  theme(axis.text.x=element_text(angle=90))

plot_cu2 <- ggplot(filter(pdat2, DIS=="cu", Genus==keeptax_cu$Genus[2]), aes(x = TP, y = Abundance, fill = Responder))+
  geom_point(position = position_dodge(width =0.75), aes(color=Responder))+
  facet_wrap(~Therapy)+
  geom_boxplot(position = position_dodge(preserve = "single"), alpha = 0.5)+
  theme_bw()+
  labs(x = "", title = keeptax_cu$Genus[2])+
  scale_y_sqrt()+
  theme(axis.text.x=element_text(angle=90))

plot_mc1 <- ggplot(filter(pdat2, DIS=="mc", Genus==keeptax_mc$Genus[1]), aes(x = TP, y = Abundance, fill = Responder))+
  geom_point(position = position_dodge(width =0.75), aes(color=Responder))+
  facet_wrap(~Therapy)+
  geom_boxplot(position = position_dodge(preserve = "single"), alpha = 0.5)+
  theme_bw()+
  labs(x = "", title = keeptax_mc$Genus[1])+
  scale_y_sqrt()+
  theme(axis.text.x=element_text(angle=90))

plot_mc2 <- ggplot(filter(pdat2, DIS=="mc", Genus==keeptax_mc$Genus[2]), aes(x = TP, y = Abundance, fill = Responder))+
  geom_point(position = position_dodge(width =0.75), aes(color=Responder))+
  facet_wrap(~Therapy)+
  geom_boxplot(position = position_dodge(preserve = "single"), alpha = 0.5)+
  theme_bw()+
  labs(x = "", title = keeptax_mc$Genus[2])+
  scale_y_sqrt()+
  theme(axis.text.x=element_text(angle=90))


plot_dis <- ggpubr::ggarrange(plot_mc1, plot_mc2, plot_cu1, plot_cu2,
                              ncol = 2, nrow=2,legend = "right", 
                              common.legend = TRUE, labels = c("A", "", "B", ""))

####---- Fig4B - Deseq2 ----####

####---- preparations ----####
library(DESeq2)
library(phyloseq)

# Function for geometric mean, because of zero counts
geo_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

alph <- 0.1 # alpha for filtering of results table


# get dataset without singletons
cedt <- phylo.ns

# filter very low abundance taxa that occur only in a few samples
ps <- tax_glom(cedt, taxrank = "Genus") 
ps <- filter_taxa(ps, function(x) sum(x>3)>(0.1*length(x)), TRUE) 
ps <- prune_samples(sample_data(ps)$TP!="exclude", ps)

# format metadata
sample_data(ps)$TP <- factor(sample_data(ps)$TP)
sample_data(ps)$Case <- factor(sample_data(ps)$Case)
sample_data(ps)$Responder <- factor(sample_data(ps)$Responder, labels=c(0,1))
sample_data(ps)$TPR <- factor(paste0(sample_data(ps)$TP, 
                                     "_", 
                                     sample_data(ps)$Responder))

####---- baseline vs first FU: Dis, Therapy corrected----####
dat5 <- subset_samples(ps, TP=="baseline"| TP=="timepoint 1")
ced = phyloseq_to_deseq2(dat5,~DIS+Therapy+Responder)
geoMeans = apply(counts(ced), 1, geo_mean)
ced_dds = estimateSizeFactors(ced, geoMeans = geoMeans)
ced_dds = DESeq(ced_dds, test="Wald", fitType="local")

# Getting Results: 
res = results(ced_dds, contrast = c("Responder", "1", "0"), cooksCutoff = FALSE)
if(min(res$padj)> alph) {
  warning("No significant change!")
} else{
  res = res[order(res$padj, na.last=NA), ]
  sigtab = res[(res$padj < alph), ]
  sigtab = cbind(as(sigtab, "data.frame"), as(tax_table(ps)[rownames(sigtab), ], "matrix"))
  sigtab_nr <- tibble::rownames_to_column(sigtab, var = "OTU")
  sigtab_nr
}

#klebsiella is detected
sigtab5 <- sigtab_nr
####---- 8) baseline vs second FU: Dis corrected----####


####---- Plot results from analysis - Klebsiella ----####
psr <- transform_sample_counts(ps, function(x) x/sum(x)*100)
ps7 <- subset_taxa(physeq = psr,Genus %in% sigtab5$Genus)
sample_names(ps7) <- paste0("s", sample_names(ps7)) # avoiding problems with "numbers-only" sample names
ps_unmerged5 <- psmelt(ps7)
ps_unmerged5$TPR <- factor(paste0(ps_unmerged5$TP, "_", ps_unmerged5$Responder))
ps_unmerged5$Responder <- factor(ps_unmerged5$Responder, labels = c("Non-Remission", "Remission"))
ps_unmerged5$Therapy <- factor(ps_unmerged5$Therapy, labels = c("AZA", "anti-TNF"))
ps_unmerged5$DIS <- factor(ps_unmerged5$DIS, labels = c("UC", "CD"))
ps_unmerged5$TP <- factor(ps_unmerged5$TP, labels = c("baseline",
                                                      "timepoint 1",
                                                      "timepoint 2"))


# only in CD, filter for CD patients
# Proteus is skipped, because not meaningful data
# especially in Azathioprine subgroup
plot_subset<-dplyr::filter(ps_unmerged5, DIS=="CD", Therapy =="AZA")
klebsiella <- ggplot(filter(ps_unmerged5, DIS=="CD", Therapy == "AZA"), aes(x = factor(TP), y = Abundance, fill = Responder))+
  geom_boxplot(na.rm = TRUE, 
               position = position_dodge(preserve = "single", width = 0.9),alpha=0.7)+
  geom_point(pch=21, position=position_jitterdodge())+
  facet_grid(~Therapy)+
  theme_bw(base_size=12)+
  #coord_cartesian(ylim=c(0,15))+
  scale_y_sqrt(expand=c(0,0),labels = scales::number_format(accuracy = 0.1))+
  labs(x = "", y = "Abundance (%)", fill="")+
  theme(axis.text.x = element_text(angle = 90, vjust=0.5, hjust=1))+
  expand_limits(y=c(0,75))

klebsiella
ggsave("figures/suppfig2_deseq2_klebsiella.pdf", plot = klebsiella,
       device=cairo_pdf, width = 3.5)

# statistic test for klebsiella abundance alone (anova with Tukey)
library(multcomp)
x2 <- aov(formula = Abundance~DIS+Therapy+TP+Responder,
          data= ps_unmerged5)
x <- aov(formula = Abundance~Responder,
         data= ps_unmerged5)

ps_unmerged5$groups <- mutate(ps_unmerged5, paste0(DIS, Therapy, TP, Responder))

x2res <- summary(glht(x2, linfct=mcp(Responder="Tukey")))
x2res
pval<-x2res$test$pvalues[1]