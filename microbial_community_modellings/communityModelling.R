library(MicrobiomeAGORA)
library(ggplot2)

source("dat/init.agora.R")
rm(agora)

mic <- new("Microbiome",
           uniq.table.file = "dat/16s/otutab.txt",
           agora.mapping.file = "dat/16s/agora.mapping.output.m8",
           sample.description.file = "dat/16s/sample.info.csv"
)

mic <- filter.mapping(mic,allowed.identity.deviation=0,method.resolve.multiple = "user") # always choose 1 if asked
mic <- create.AgoraTable(mic)
mic <- filter.samples(mic, min.seqs = 250, max.unclassified = 0.5)

mic.comm <- simulate_community_metabolism.FB2(mbo = mic,
                                              agora = readRDS("dat/Agora-1.02-Western_corrected_v4_new.RDS"),
                                              abunCutoff = 0.001, scale.boundaries = 1)
mic.comm.sum <- summary.comm.mods(mic.comm)
#saveRDS(mic.comm.sum, file = "dat/mic.comm.sum_jcc2.RDS")

#
# STATS
#
dt <- merge(mic.comm.sum$met.int, mic@sample.description, by = "sample")
dt <- merge(dt, mic.comm.sum$comm.gr, by = "sample")
dt[,flux.c := flux / c.gr]
dt[,o.flux.c := o.flux / c.gr]

library(ggplot2)
met <- "EX_but(e)"

dt.tmp <- rbind(copy(dt), copy(dt)[,Therapy := "overall"])
# Note: The Overall plot includes "cu" and "mc". Defining DIS = "cu" is only for the purpose to include all data points 
# in this panel. Layout of plots was later edited using Inkscape
dt.tmp[Therapy == "overall", DIS := "cu"] 
p <- ggplot(dt.tmp[TP==0 & rxn == met], aes(Responder, o.flux.c)) +
  geom_boxplot(fill = "grey", outlier.shape = NA) + geom_point(size = 0.5) + labs(y = "Butyrate production [mmol/gDW]") +
  theme_bw() +
  facet_grid(DIS~Therapy)
p # Fig ???

# Overall
wilcox.test(dt[TP==0 & Responder == "Responder" & rxn == met, o.flux.c], 
            dt[TP==0 & Responder == "Non-Responder" & rxn == met, o.flux.c])

# aza all DIS
wilcox.test(dt[TP==0 & Responder == "Responder"     & rxn == met & Therapy == "aza", o.flux.c], 
            dt[TP==0 & Responder == "Non-Responder" & rxn == met & Therapy == "aza", o.flux.c])

# biol all DIS
wilcox.test(dt[TP==0 & Responder == "Responder"     & rxn == met & Therapy == "biol", o.flux.c], 
            dt[TP==0 & Responder == "Non-Responder" & rxn == met & Therapy == "biol", o.flux.c])

# MC all treatments
wilcox.test(dt[TP==0 & Responder == "Responder"     & DIS=="mc" & rxn == met, o.flux.c], 
            dt[TP==0 & Responder == "Non-Responder" & DIS=="mc" & rxn == met, o.flux.c])

# MC & aza
wilcox.test(dt[TP==0 & Responder == "Responder"     & rxn == met & DIS=="mc" & Therapy == "aza", o.flux.c], 
            dt[TP==0 & Responder == "Non-Responder" & rxn == met & DIS=="mc" & Therapy == "aza", o.flux.c])

# MC & biol
wilcox.test(dt[TP==0 & Responder == "Responder"     & rxn == met & DIS=="mc" & Therapy == "biol", o.flux.c], 
            dt[TP==0 & Responder == "Non-Responder" & rxn == met & DIS=="mc" & Therapy == "biol", o.flux.c])

# cu all treatmens
wilcox.test(dt[TP==0 & Responder == "Responder"     & DIS=="cu" & rxn == met, o.flux.c], 
            dt[TP==0 & Responder == "Non-Responder" & DIS=="cu" & rxn == met, o.flux.c])

# cu & aza
wilcox.test(dt[TP==0 & Responder == "Responder"     & rxn == met & DIS=="cu" & Therapy == "aza", o.flux.c], 
            dt[TP==0 & Responder == "Non-Responder" & rxn == met & DIS=="cu" & Therapy == "aza", o.flux.c])

# cu & biol
wilcox.test(dt[TP==0 & Responder == "Responder"     & rxn == met & DIS=="cu" & Therapy == "biol", o.flux.c], 
            dt[TP==0 & Responder == "Non-Responder" & rxn == met & DIS=="cu" & Therapy == "biol", o.flux.c])



