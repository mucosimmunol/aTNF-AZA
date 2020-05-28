library(stringr)

repair.metRxn.names <- function(mod) {
  mod@met_id <- str_replace_all(mod@met_id,"__40__","(")
  mod@met_id <- str_replace_all(mod@met_id,"__41__",")")
  mod@met_id <- str_replace_all(mod@met_id,"__91__","[")
  mod@met_id <- str_replace_all(mod@met_id,"__93__","]")
  
  mod@react_id <- str_replace_all(mod@react_id,"__40__","(")
  mod@react_id <- str_replace_all(mod@react_id,"__41__",")")
  mod@react_id <- str_replace_all(mod@react_id,"__91__","[")
  mod@react_id <- str_replace_all(mod@react_id,"__93__","]")
  
  return(mod)
}