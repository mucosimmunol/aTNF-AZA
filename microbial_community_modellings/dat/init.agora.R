library(sybil)
library(stringr)
library(data.table)
source("dat/correct_common_errors_agora4.R")
source("dat/repair.metRxn.names.R")
sybil::SYBIL_SETTINGS("SOLVER","cplexAPI"); ok <- 1

#~~~~~~~~~~~~~~#
# Western diet #
#~~~~~~~~~~~~~~#
#agora <- readRDS("dat/Agora-1.02-Western.RDS")
a1 <- readRDS("dat/Agora-1.02-Western_pt1.RDS")
a2 <- readRDS("dat/Agora-1.02-Western_pt2.RDS")
a3 <- readRDS("dat/Agora-1.02-Western_pt3.RDS")
a4 <- readRDS("dat/Agora-1.02-Western_pt4.RDS")
agora <- c(a1, a2, a3, a4)
agora <- lapply(agora, repair.metRxn.names)


agora <- lapply(agora, correct_common_errors_agora)

dtc <- fread("dat/exrxn.repair.csv")
dtc <- dtc[!duplicated(exrxn)]
agora <- lapply(agora, function(x) {
  dtc.tmp <- copy(dtc)[exrxn %in% x@react_id]
  
  if(nrow(dtc.tmp)>0)
  x <- changeBounds(x, react = dtc.tmp$exrxn, lb = dtc.tmp$lb)
  
  return(x)
})

saveRDS(agora, file = "dat/Agora-1.02-Western_corrected_v4_new.RDS")