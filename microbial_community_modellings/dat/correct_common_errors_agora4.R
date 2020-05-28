correct_common_errors_agora <- function(mod,excl.rules=integer(0)) {
  #correct_common_errors_agora <- function(mod,excl.rules=96) {
  i <- 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 1 - Stochiometry in D-Lac+Na / L-mal+H Antiporter (ID: MALLACDt)
    rxn.id <- which(mod@react_id=="MALLACDt")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -2
      mod@S[h.e.id,rxn.id] <- 2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 2 - sodium proton antiporter (H:NA is 1:2) (ID: NAt3_1) Should not be 1:2!!!
    # OK
    rxn.id <- which(mod@react_id=="NAt3_1")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 3 (amino acid / H+ symport (1:2)!!)
    # OK
    rxn.id <- which(mod@react_id=="ALAt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 4 (carboxylic acid / H+ symport (1:2)!!)
    # OK
    rxn.id <- which(mod@react_id=="L_LACt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 5 (carboxylic acid / H+ symport (1:2)!!)
    # ?
    rxn.id <- which(mod@react_id=="FORt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 6 (carboxylic acid / H+ symport (1:2)!!)
    # OK
    rxn.id <- which(mod@react_id=="D_LACt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 7 (amino acid / H+ symport (1:2)!!)
    # ?
    rxn.id <- which(mod@react_id=="GLYt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 8 (Phosphate/Na+ Symport: Bacteria probably have the PMID: 23398154. 
    # Citation: "NaPi-IIc are insensitive to membrane potential, no net charge is translocated, and transport 
    # is mediated with a 2:1 stoichiometry")
    # ?
    rxn.id <- which(mod@react_id=="PIt7")
    h.c.id <- which(mod@met_id=="na1[c]")
    h.e.id <- which(mod@met_id=="na1[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 9 (amino acid / H+ symport (1:2)!!)
    # OK
    rxn.id <- which(mod@react_id=="THRt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 10 (carboxylic acid / H+ symport (1:2)!!)
    # OK
    rxn.id <- which(mod@react_id=="L_LACt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 11 (amino acid / H+ symport (1:2)!!)
    # ?
    rxn.id <- which(mod@react_id=="GLYt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 12 H+ is on wrong side of the reaction "FDHfdx" + Stoichiometries of fdx.. are wrong
    # new stoichiometry based on curated reaction description: http://www.rhea-db.org/reaction?id=46952
    rxn.id <- which(mod@react_id=="FDHfdx")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    nad.id <- which(mod@met_id=="nad[c]")
    nadh.id <- which(mod@met_id=="nadh[c]")
    for.id <- which(mod@met_id=="for[c]")
    
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
      mod@S[fox.id,rxn.id] <- 2
      mod@S[frd.id,rxn.id] <- -2
      mod@S[nad.id,rxn.id] <- 1
      mod@S[nadh.id,rxn.id] <- -1
      mod@S[for.id,rxn.id] <- 2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 13 (amino acid / H+ symport (1:2)!!)
    # ?
    rxn.id <- which(mod@react_id=="DSERt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 14 (amino acid / H+ symport (1:2)!!)
    # OK
    rxn.id <- which(mod@react_id=="SERt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 15 (amino acid / H+ symport (1:2)!!)
    # ?
    rxn.id <- which(mod@react_id=="ASPt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 16 (amino acid / H+ symport (1:2)!!)
    # ?
    rxn.id <- which(mod@react_id=="ASPt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 17 - sodium proton antiporter (H:NA is 1:1) (ID: NAt3) Should not be 1:1!!!
    # ?
    rxn.id <- which(mod@react_id=="NAt3")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 18 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="MALt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 19 (organic acid / H+ symport (1:2)!!) CHECK IF THIS IS TRUE / NECESSARY? NOT!
    #rxn.id <- which(mod@react_id=="GLYC3Pt")
    #h.c.id <- which(mod@met_id=="h[c]")
    #h.e.id <- which(mod@met_id=="h[e]")
    #if(length(rxn.id) == 1) {
    #  mod@S[h.c.id,rxn.id] <- 2
    #  mod@S[h.e.id,rxn.id] <- -2
    #}
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 20 (phosphate acid / H+ symport (1:4) - to make it electroneural) (only import)
    rxn.id <- which(mod@react_id=="PIt6b")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 4
      mod@S[h.e.id,rxn.id] <- -4
      mod@lowbnd[rxn.id] <- 0
      mod@react_rev[rxn.id] <- FALSE
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERRRO 21 Reversibility of Trp-Synthase (it is not reversible)
    rxn.id <- which(mod@react_id=="TRPS2r")
    if(length(rxn.id) == 1) {
      mod@lowbnd[rxn.id] <- 0
      mod@react_rev[rxn.id] <- FALSE
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERRRO 22 Reversibility of phosphate:H+ sympoter. No much known about kinetics and mechanisms. Only Import has been
    # reported so far. (https://doi.org/10.1016/j.bios.2011.06.024)
    rxn.id <- which(mod@react_id=="GLYC3Pt")
    if(length(rxn.id) == 1) {
      #mod@lowbnd[rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 23 H+ is on wrong side of the reaction "FDNADOX_H" AND Stoichiometriy of ferredoxin is wrong
    rxn.id <- which(mod@react_id=="FDNADOX_H")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
      mod@S[h.e.id,rxn.id] <- 0
      mod@S[fox.id,rxn.id] <- 2
      mod@S[frd.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 24 Diffusion of formate is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="FORt")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      #mod@lowbnd[rxn.id] <- 0
      #mod@uppbnd[rxn.id] <- 0
      
      #mod@S[h.c.id,rxn.id] <- 2
      #mod@S[h.e.id,rxn.id] <- -2
      mod <- addReact(mod,
                      id = "FORt2r",
                      met = c("h[e]","for[e]","h[c]","for[c]"),
                      Scoef = c(-2,-1,2,1),
                      reversible = T,
                      reactName = "formate transport in via proton symport")
      mod@react_attr[which(mod@react_id=="FORt2r"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="FORt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 25 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="GLUt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 26 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="MALt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 27 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="ABUTt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 28 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="GLNt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 29 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="FUMt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 30 Diffusion of fumarate is unlikely (at least in high amounts) and causes futile cycles
    # replace with fumarate symport (2H+)
    rxn.id <- which(mod@react_id=="FUMt")
    if(length(rxn.id) == 1) {
      
      if(length(which(mod@react_id=="FUMt2_2"))==0) {
        mod <- addReact(mod, id = "FUMt2_2",
                        met = c("h[e]","fum[e]","h[c]","fum[c]"),
                        Scoef = c(-2,-1,2,1),
                        reversible = F,
                        reactName = "Fumarate transport via proton symport (2 H)")
        mod@react_attr[which(mod@react_id=="FUMt2_2"),] <- c(annotation="",
                                                             notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      }
      
      mod <- rmReact(mod,rxn.id)
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 31 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="3MOPt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 32 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="ILEt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 33 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="AKGt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 34 H+ is on wrong side of the reaction "BTCOADH" (doi: 10.1128/JB.01417-07) AND 
    # http://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec=1.3.1.109
    # Correct stoichiometry from: http://www.rhea-db.org/reaction?id=46960
    rxn.id <- which(mod@react_id=="BTCOADH")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 0
      mod@S[fox.id,rxn.id] <- -2
      mod@S[frd.id,rxn.id] <- 2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 35 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="DALAt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 36 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="SUCCt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 37 Diffusion of succinate is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="SUCCt")
    if(length(rxn.id) == 1) {
      
      
      if(length(which(mod@react_id=="SUCCt2_2"))==0) {
        mod <- addReact(mod, id = "SUCCt2_2",
                        met = c("h[e]","succ[e]","h[c]","succ[c]"),
                        Scoef = c(-2,-1,2,1),
                        reversible = T,
                        reactName = "succinate transport via proton symport (2 H)")
        mod@react_attr[which(mod@react_id=="SUCCt2_2"),] <- c(annotation="",
                                                              notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
        mod@lowbnd[which(mod@react_id=="SUCCt2_2")] <- -1000
        
      }
      
      mod <- rmReact(mod,rxn.id)
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 38 H+ is on wrong side of the reaction "FDNADOX_H" AND Stoichiometriy of ferredoxin is wrong
    rxn.id <- which(mod@react_id=="OOR2r")
    h.c.id <- which(mod@met_id=="h[c]")
    
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
      mod@S[fox.id,rxn.id] <- 2
      mod@S[frd.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 39 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="TRPt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 40 H+ is on wrong side of the reaction "POR4" AND Stoichiometriy of ferredoxin is wrong#
    # Plus: Forward and backwards reaction seem to be catalysed by different enzymes, while the pyruvate lysase ist much more active.
    # http://www.jbc.org/content/252/8/2657.full.pdf
    rxn.id <- which(mod@react_id=="POR4")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
      mod@S[fox.id,rxn.id] <- 2
      mod@S[frd.id,rxn.id] <- -2

      mod@uppbnd[rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERRRO 41 Reversibility of phosphoenolpyruvate carboxykinase (it is reversible, however the kinetics are 
    # very unfavourable [doi:10.1016/j.ymben.2004.03.001])
    rxn.id <- which(mod@react_id=="PPCKr")
    if(length(rxn.id) == 1) {
      mod@lowbnd[rxn.id] <- 0
      mod@react_rev[rxn.id] <- FALSE
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 42 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="ASNt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 43 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="r1088")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 44 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="ABUTt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 45 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="BUTt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 46 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="BUTt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 47 H+ is missing in the reaction "GLFRDO" AND Stoichiometriy of ferredoxin is wrong
    rxn.id <- which(mod@react_id=="GLFRDO")
    h.c.id <- which(mod@met_id=="h[c]")
    
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[fox.id,rxn.id] <- -2
      mod@S[frd.id,rxn.id] <- 2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 48 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="CITt10i")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 49 H+ has wrong Stiochiometry in the reaction "FXXRDO" AND Stoichiometriy of ferredoxin is wrong
    rxn.id <- which(mod@react_id=="FXXRDO")
    h.c.id <- which(mod@met_id=="h[c]")
    
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 8
      mod@S[fox.id,rxn.id] <- -6
      mod@S[frd.id,rxn.id] <- 6
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 50 H+ has wrong Stiochiometry in the reaction "SULR" - MAybe not
    rxn.id <- which(mod@react_id=="SULR")
    h.c.id <- which(mod@met_id=="h[c]")
    if(length(rxn.id) == 1) {
      #mod@S[h.c.id,rxn.id] <- 3
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 51 H+ has wrong Stiochiometry in the reaction "FRDOr" AND Stoichiometriy of ferredoxin is wrong
    rxn.id <- which(mod@react_id=="FRDOr")
    h.c.id <- which(mod@met_id=="h[c]")
    
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
      mod@S[fox.id,rxn.id] <- 2
      mod@S[frd.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 52 H+ is on wrong side of the reaction "IOR3" AND Stoichiometriy of ferredoxin is wrong
    rxn.id <- which(mod@react_id=="IOR3")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 1
      mod@S[fox.id,rxn.id] <- -2
      mod@S[frd.id,rxn.id] <- 2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 53 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="TYRt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 54 H+ is on wrong side of the reaction "HYD4" AND Stoichiometriy of ferredoxin is wrong
    # also: It is reversible! doi.org/10.1016/0005-2728(90)90044-5
    rxn.id <- which(mod@react_id=="HYD4")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -2
      mod@S[fox.id,rxn.id] <- 2
      mod@S[frd.id,rxn.id] <- -2
      mod@react_rev[rxn.id] <- TRUE
      mod@lowbnd[rxn.id] <- -1000
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 55 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="ARGt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 56 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="ORNt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 57 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="SUCCt2_3r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 58 Diffusion of Guanine is unlikely (at least in high amounts) and causes futile cycles
    # replace with proton-symport
    rxn.id <- which(mod@react_id=="GUAt")
    rxn.id <- which(mod@react_id=="GUAt")
    if(length(rxn.id) == 1) {
      
      if(length(which(mod@react_id=="GUAt2"))==0) {
        mod <- addReact(mod, id = "GUAt2",
                        met = c("h[e]","gua[e]","h[c]","gua[c]"),
                        Scoef = c(-1,-1,1,1),
                        reversible = F,
                        reactName = "guanine transport in via proton symport")
        mod@react_attr[which(mod@react_id=="GUAt2"),] <- c(annotation="",
                                                           notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      }
      mod <- rmReact(mod,rxn.id)
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 59 Stoichiometry of Protons incorrect
    rxn.id <- which(mod@react_id=="NADH6")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
      mod@S[h.e.id,rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 60 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="GLCNt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 61 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="DDGLCNt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 62 Reaction DMPPS got reclassified and updated in 2016
    rxn.id <- which(mod@react_id=="DMPPS")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    nadh.id <- which(mod@met_id=="nadh[c]")
    nad.id <- which(mod@met_id=="nad[c]")
    if(length(rxn.id) == 1 & length(fox.id) == 1) {
      #mod@S[h.c.id,rxn.id] <- -2
      #mod@S[fox.id,rxn.id] <- 2
      #mod@S[frd.id,rxn.id] <- -2
      
      #mod@S[nadh.id,rxn.id] <- 0
      #mod@S[nad.id,rxn.id] <- 0
    }
    if(length(rxn.id) == 1 & length(fox.id) == 0) { 
      # This is the case if fdxox and fdrd are not part of the model, but the reaction DMPPS is. In this case: 
      # Block the reaction
      #mod@lowbnd[rxn.id] <- 0
      #mod@uppbnd[rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 63 Reaction DMPPS2r got reclassified and updated in 2016 and merged with DMPPS
    rxn.id <- which(mod@react_id=="DMPPS2r")
    
    if(length(rxn.id) == 1) {
      #mod@lowbnd[rxn.id] <- 0
      #mod@uppbnd[rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 64 carboxylic acid / H+ symport  should be ration 1:2 - Mg2+ is transported in this reaction (CITt10)
    # as well, however, MG2+ can diffuse across membrane barrier by reaction "MGt5"
    rxn.id <- which(mod@react_id=="CITt10")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 65 carboxylic acid / H+ symport  should be ration 1:2 - Mg2+ is transported in this reaction (ICITt10)
    # as well, however, MG2+ can diffuse across membrane barrier by reaction "MGt5"
    rxn.id <- which(mod@react_id=="ICITt10")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 66 
    # UPDATE: Check if that reaction is essential for some AGORAs
    rxn.id <- which(mod@react_id=="NADH6")
    if(length(rxn.id) == 1) {
      #mod@lowbnd[rxn.id] <- 0
      #mod@uppbnd[rxn.id] <- 0
      #mod <- rmReact(mod,rxn.id)
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 67 Oxaloacetate synthase (citrate as substrate) does not seem to exist in chemical reaction databases. However, 
    # the opposite direction exists, which is the "citrate synthase" and is irreversible.
    rxn.id <- which(mod@react_id=="OAASr")
    h2o.id <- which(mod@met_id=="h2o[c]")
    accoa.id <- which(mod@met_id=="accoa[c]")
    oaa.id <- which(mod@met_id=="oaa[c]")
    h.id <- which(mod@met_id=="h[c]")
    cit.id <- which(mod@met_id=="cit[c]")
    coa.id <- which(mod@met_id=="coa[c]")
    if(length(rxn.id) == 1) {
      mod@S[h2o.id,rxn.id] <- -1
      mod@S[accoa.id,rxn.id] <- -1
      mod@S[oaa.id,rxn.id] <- -1
      
      mod@S[h.id,rxn.id] <- 1
      mod@S[cit.id,rxn.id] <- 1
      mod@S[coa.id,rxn.id] <- 1
      
      mod@react_rev[rxn.id] <- FALSE
      mod@lowbnd[rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 68 Stoichiometries incorrect in "MTHFRfdx"
    rxn.id <- which(mod@react_id=="MTHFRfdx")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    nadh.id <- which(mod@met_id=="nadh[c]")
    nad.id <- which(mod@met_id=="nad[c]")
    if(length(rxn.id) == 1 & length(fox.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[fox.id,rxn.id] <- -2
      mod@S[frd.id,rxn.id] <- 2
      
      mod@S[nadh.id,rxn.id] <- 0
      mod@S[nad.id,rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 69 H+ stoichiometry wrong in "HYDFDNrfdx" (https://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec=1.12.1.4)
    rxn.id <- which(mod@react_id=="HYDFDNrfdx")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1 & length(fox.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -3
      mod@S[fox.id,rxn.id] <- 2
      mod@S[frd.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 70 H+ stoichiometry wrong in "HYDFDN2rfdx" (https://www.ebi.ac.uk/intenz/query?cmd=SearchEC&ec=1.12.1.4)
    rxn.id <- which(mod@react_id=="HYDFDN2rfdx")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1 & length(fox.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -3
      mod@S[fox.id,rxn.id] <- 2
      mod@S[frd.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 71 H+ stoichiometry wrong in "CODH_ACS"
    rxn.id <- which(mod@react_id=="CODH_ACS")
    h.c.id <- which(mod@met_id=="h[c]")
    fox.id <- which(mod@met_id=="fdxox[c]")
    frd.id <- which(mod@met_id=="fdxrd[c]")
    if(length(rxn.id) == 1 & length(fox.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
      mod@S[fox.id,rxn.id] <- 2
      mod@S[frd.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 72 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="CITt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 73 Diffusion of R)-3-hydroxybutanoate is unlikely (at least in high amounts) and causes futile cycles
    # replace with proton symport
    rxn.id <- which(mod@react_id=="BHBtp")
    if(length(rxn.id) == 1) {
      #mod@lowbnd[rxn.id] <- 0
      #mod@uppbnd[rxn.id] <- 0
      
      if(length(which(mod@react_id=="BHBt"))==0) {
        mod <- addReact(mod, id = "BHBt",
                        met = c("h[e]","bhb[e]","h[c]","bhb[c]"),
                        Scoef = c(-1,-1,1,1),
                        reversible = T,
                        reactName = "(R)-3-Hydroxybutanoate transport via H+ symport")
        mod@react_attr[which(mod@react_id=="BHBt"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
        mod@lowbnd[which(mod@react_id=="BHBt")] <- -1000
        
      }
      
      mod <- rmReact(mod,rxn.id)
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 74 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="ARGt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 75 ASPNH4L is most likely (nearly) irreversible.
    rxn.id <- which(mod@react_id=="ASPNH4L")
    if(length(rxn.id) == 1) {
      #mod@lowbnd[rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 76 "GLXS" is the reverse version of  (http://identifiers.org/ec-code/2.3.3.9). However, EC 2.3.3.9 is irreversable
    rxn.id <- which(mod@react_id=="GLXS")
    h2o.id <- which(mod@met_id=="h2o[c]")
    accoa.id <- which(mod@met_id=="accoa[c]")
    glx.id <- which(mod@met_id=="glx[c]")
    h.id <- which(mod@met_id=="h[c]")
    mal.id <- which(mod@met_id=="mal_L[c]")
    coa.id <- which(mod@met_id=="coa[c]")
    if(length(rxn.id) == 1) {
      mod@S[h2o.id,rxn.id] <- -1
      mod@S[accoa.id,rxn.id] <- -1
      mod@S[glx.id,rxn.id] <- -1
      
      mod@S[h.id,rxn.id] <- 1
      mod@S[mal.id,rxn.id] <- 1
      mod@S[coa.id,rxn.id] <- 1
      
      mod@react_rev[rxn.id] <- FALSE
      mod@lowbnd[rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 77 (carboxylic acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="ACACt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      #mod@S[h.c.id,rxn.id] <- 2
      #mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 78 Adjusting stoichiometries of FDH2 (Formate dyhydrogenase (quinine-8)) according to 
    rxn.id <- which(mod@react_id=="FDH2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    for.c.id <- which(mod@met_id=="for[c]")
    co2.c.id <- which(mod@met_id=="co2[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
      mod@S[h.e.id,rxn.id] <- 0
      mod@S[for.c.id,rxn.id] <- -1
      mod@S[co2.c.id,rxn.id] <- 1
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 79 "ASPNH4L" is irreversable
    rxn.id <- which(mod@react_id=="ASPNH4L")
    if(length(rxn.id) == 1) {
      mod@lowbnd[rxn.id] <- 0
      mod@react_rev[rxn.id] <- FALSE
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 80 correcting stoichiometries of COF420H according to
    # http://identifiers.org/ec-code/1.12.98.1
    rxn.id <- which(mod@react_id=="COF420H")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
      mod@S[h.e.id,rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 81 (Phosphate/Na+ Symport: Bacteria probably have the PMID: 23398154. 
    # Citation: "NaPi-IIc are insensitive to membrane potential, no net charge is translocated, and transport 
    # is mediated with a 2:1 stoichiometry")
    # OK
    rxn.id <- which(mod@react_id=="PIt7ir")
    h.c.id <- which(mod@met_id=="na1[c]")
    h.e.id <- which(mod@met_id=="na1[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 82 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="PROt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 83 Diffusion of proline is unlikely (at least in high amounts) and causes futile cycles. Replace by H+ symport
    rxn.id <- which(mod@react_id=="PROPAT4te")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      
      if(length(which(mod@react_id=="PROt2r"))==0) {
        mod <- addReact(mod, id = "PROt2r",
                        met = c("h[e]","pro_L[e]","h[c]","pro_L[c]"),
                        Scoef = c(-2,-1,2,1),
                        reversible = F,
                        reactName = "L proline transport in via proton symport")
        mod@react_attr[which(mod@react_id=="PROt2r"),] <- c(annotation="",
                                                            notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      }
      mod <- rmReact(mod,rxn.id)
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 84 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="PROt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 85 correcting stoichiometries of COF420_NADP_OX according to
    # http://identifiers.org/ec-code/1.5.1.40
    rxn.id <- which(mod@react_id=="COF420_NADP_OX")
    h.c.id <- which(mod@met_id=="h[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -2
    }
    
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 86 - massive inconsistencies in reactions r0318 & r0220
    rxn.id <- which(mod@react_id=="r0220")
    HC01668.id <- which(mod@met_id=="HC01668[c]")
    ppa.id <- which(mod@met_id=="ppa[c]")
    atp.id <- which(mod@met_id=="atp[c]")
    ppi.id <- which(mod@met_id=="ppi[c]")
    h.id <- which(mod@met_id=="h[c]")
    if(length(rxn.id) == 1) {
      mod@S[HC01668.id,rxn.id] <- 0
      mod@S[ppa.id,rxn.id] <- -1
      mod@S[atp.id,rxn.id] <- -1
      mod@S[ppi.id,rxn.id] <- 1
      mod@S[h.id,rxn.id] <- 0
    }
    rxn.id <- which(mod@react_id=="r0318")
    if(length(rxn.id) == 1) {
      #mod@lowbnd[rxn.id] <- 0
      #mod@uppbnd[rxn.id] <- 0
      
      mod <- rmReact(mod,rxn.id)
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 87 - Electron inbalance in FHL (formate-hydrogen lyase)
    rxn.id <- which(mod@react_id=="FHL")
    nad.id <- which(mod@met_id=="nad[c]")
    nadh.id <- which(mod@met_id=="nadh[c]")
    h2.id <- which(mod@met_id=="h2[c]")
    h.id <- which(mod@met_id=="h[c]")
    if(length(rxn.id) == 1) {
      mod@S[nad.id,rxn.id] <- -1
      mod@S[nadh.id,rxn.id] <- 1
      mod@S[h2.id,rxn.id] <- 0
      mod@S[h.id,rxn.id] <- 0
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 88 - Uncertainty in reversebility of OCBT
    rxn.id <- which(mod@react_id=="OCBT")
    if(length(rxn.id) == 1) {
      #mod@lowbnd[rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 89 - the reaction PPC (Phosphoenolpyruvate carboxylase) is irreversible and forms oaa. Replace "PPCr" with "PPC"
    # Ref: http://www.brenda-enzymes.info/enzyme.php?ecno=4.1.1.31
    # and : doi/10.1073/pnas.0913127107
    rxn.id <- which(mod@react_id=="PPCr")
    if(length(rxn.id) == 1) {
      
      if(length(which(mod@react_id=="PPC"))==0) {
        mod <- addReact(mod, id = "PPC",
                        met = c("h2o[c]","co2[c]","pep[c]","h[c]","pi[c]","oaa[c]"),
                        Scoef = c(-1,-1,-1,1,1,1),
                        reversible = F,
                        reactName = "(R)-3-Hydroxybutanoate transport via H+ symport")
        mod@react_attr[which(mod@react_id=="PPC"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      }
      
      mod <- rmReact(mod,rxn.id)
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 90 - the reaction CITRH (Citrullinase) is reported to be reversible in old 50s papers. But under physiological
    # conditions, this reaction is most likely (nearly) irreversible. The reverse direction also causes futile cycles in 
    # paired models. It's also ir according to MetaCyc.
    rxn.id <- which(mod@react_id=="CITRH")
    if(length(rxn.id) == 1) {
      mod@lowbnd[rxn.id] <- 0
      mod@react_rev[rxn.id] <- FALSE
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 91 ethanol can usually pass the cell membranes by diffusion - thus proton symport is unlikely and causes
    # futile cycles
    rxn.id <- which(mod@react_id=="ETOHt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 0
      mod@S[h.e.id,rxn.id] <- 0
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 92 Nucleotide deaminases are usually irreversible - so is DCMPDA
    rxn.id <- which(mod@react_id=="DCMPDA")
    if(length(rxn.id) == 1) {
      mod@lowbnd[rxn.id] <- 0
      mod@react_rev[rxn.id] <- FALSE
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 93 Stoichiometry of Protons incorrect
    rxn.id <- which(mod@react_id=="NADH8")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -5
      mod@S[h.e.id,rxn.id] <- 4
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 94 (indole / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="INDOLEt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 95 Phosphoketolase is most likely (nearly) irreversible
    # ref: 10.1128/JB.00713-12 
    rxn.id <- which(mod@react_id=="PKL")
    if(length(rxn.id) == 1) {
      mod@lowbnd[rxn.id] <- 0
      mod@react_rev[rxn.id] <- FALSE
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 96 - correcting stoichiometry of AHSERL3 according to
    # http://www.ebi.ac.uk/intenz/query?q=2.5.1.49&submit=Search
    #print(i)
    rxn.id <- which(mod@react_id=="AHSERL3")
    trdrd.id <- which(mod@met_id=="trdrd[c]")
    trdox.id <- which(mod@met_id=="trdox[c]")
    h2s.id <- which(mod@met_id=="h2s[c]")
    h.id <- which(mod@met_id=="h[c]")
    tsul.id <- which(mod@met_id=="tsul[c]")
    so3.id <- which(mod@met_id=="so3[c]")
    if(length(rxn.id) == 1) {
      #mod@S[trdrd.id,rxn.id] <- 0
      #mod@S[trdox.id,rxn.id] <- 0
      #mod@S[h2s.id,rxn.id] <- -1
      #mod@S[h.id,rxn.id] <- 1
      #mod@S[tsul.id,rxn.id] <- 0
      #mod@S[so3.id,rxn.id] <- 0
      mod <- addReact(mod, id = "AHSERL3",
                      met = c("achms[c]","tsul[c]","so3[c]","hcys_L[c]","ac[c]","h[c]"),
                      Scoef = c(-1,-1,1,1,1,1),
                      reversible = T,
                      reactName = "O-Acetyl-L-homoserine acetate-lyase")
      mod@react_attr[which(mod@react_id=="AHSERL3"),] <- c(annotation="",
                                                           notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="AHSERL3")] <- -1000
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 97 Stoichiometry of Protons incorrect
    rxn.id <- which(mod@react_id=="NADH5")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      #mod@S[h.c.id,rxn.id] <- -5
      #mod@S[h.e.id,rxn.id] <- 4
    }
  }
  i <- i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 98 (amino acid / H+ symport (1:2)!!)
    rxn.id <- which(mod@react_id=="LYSt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 99 Diffusion of tyrosine is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="TYRt")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "TYRt2r",
                      met = c("h[e]","tyr_L[e]","h[c]","tyr_L[c]"),
                      Scoef = c(-2,-1,2,1),
                      reversible = T,
                      reactName = "Tyrosine transport via proton symport")
      mod@react_attr[which(mod@react_id=="TYRt2r"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="TYRt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 100 There's no evidence (literature) for the existence of an lysine:proton antiporter. (to our best knowledge)
    # -> Replace with proton symporter
    rxn.id <- which(mod@react_id=="LYSt3r")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "LYSt2r",
                      met = c("h[e]","lys_L[e]","h[c]","lys_L[c]"),
                      Scoef = c(-2,-1,2,1),
                      reversible = T,
                      reactName = "L-lysine reversible transport via proton symport")
      mod@react_attr[which(mod@react_id=="LYSt2r"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="LYSt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 101 Diffusion of ornithine is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="ORNt")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "ORNt2r",
                      met = c("h[e]","orn[e]","h[c]","orn[c]"),
                      Scoef = c(-2,-1,2,1),
                      reversible = T,
                      reactName = "orntithine reversible transport in via proton symport")
      mod@react_attr[which(mod@react_id=="ORNt2r"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="ORNt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 102 H+ is on wrong side of the reaction "CYSS3r" 
    # http://www.genome.jp/dbget-bin/www_bget?rn:R04859
    rxn.id <- which(mod@react_id=="CYSS3r")
    h.c.id <- which(mod@met_id=="h[c]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- -1
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 103 Diffusion of Propionate is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="PPAtr")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "PPAt2r",
                      met = c("h[c]","ppa[c]","h[e]","ppa[e]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "Propionate transport via proton symport reversible")
      mod@react_attr[which(mod@react_id=="PPAt2r"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="PPAt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1

  if(!(i %in% excl.rules)) {
    # ERROR 104 Diffusion of Nicotinamide is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="NCAMUP")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "NCAMt2r",
                      met = c("h[c]","ncam[c]","h[e]","ncam[e]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "Nicotinamide transport by proton symport reversible")
      mod@react_attr[which(mod@react_id=="NCAMt2r"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="NCAMt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 105 Diffusion of hypoxanthine is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="HYXNti")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "HXANt2r",
                      met = c("h[e]","hxan[e]","h[c]","hxan[c]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "hypoxanthine reversible transport via proton symport")
      mod@react_attr[which(mod@react_id=="HXANt2r"),] <- c(annotation="",
                                                           notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="HXANt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 106 Diffusion of hypoxanthine is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="HYXNt")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "HXANt2r",
                      met = c("h[e]","hxan[e]","h[c]","hxan[c]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "hypoxanthine reversible transport via proton symport")
      mod@react_attr[which(mod@react_id=="HXANt2r"),] <- c(annotation="",
                                                           notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="HXANt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 107 Diffusion of thymidine is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="THYMDtr2")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "THMDt2r",
                      met = c("h[e]","thymd[e]","h[c]","thymd[c]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "thymidine transport in via proton symport reversible")
      mod@react_attr[which(mod@react_id=="THMDt2r"),] <- c(annotation="",
                                                           notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="THMDt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 108 Diffusion of D-Arabinose is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="ARABt")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "ARABDt2",
                      met = c("h[e]","arab_D[e]","h[c]","arab_D[c]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "D-Arabinose transport via proton symport")
      mod@react_attr[which(mod@react_id=="ARABDt2"),] <- c(annotation="",
                                                           notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="ARABDt2")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 109 Diffusion of thymidine is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="THYMDt")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "THMDt2r",
                      met = c("h[e]","thymd[e]","h[c]","thymd[c]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "thymidine transport in via proton symport reversible")
      mod@react_attr[which(mod@react_id=="THMDt2r"),] <- c(annotation="",
                                                           notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="THMDt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 110 Diffusion of inosine is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="INSt")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "INSt2",
                      met = c("h[e]","ins[e]","h[c]","ins[c]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "inosine transport in via proton symport, reversible")
      mod@react_attr[which(mod@react_id=="INSt2"),] <- c(annotation="",
                                                           notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="INSt2")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 111 Diffusion of guanosine is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="GSNt")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "GSNt2r",
                      met = c("h[e]","gsn[e]","h[c]","gsn[c]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "guanosine transport in via proton symport reversible")
      mod@react_attr[which(mod@react_id=="GSNt2r"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="GSNt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 112 Diffusion of Nicotinic acid is unlikely (at least in high amounts) and causes futile cycles
    rxn.id <- which(mod@react_id=="NACUP")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "NACt2r",
                      met = c("h[c]","nac[c]","h[e]","nac[e]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "Nicotinic acid transport by proton symport reversible")
      mod@react_attr[which(mod@react_id=="NACt2r"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="NACt2r")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 113 potassium proton antiporter is only activated in bacteria to get potassium out of the cell (under hugh salinity)
    # Thus -> not reversible
    # Plus add K difffusion if not already there
    # ref: DOI:10.1074/jbc.M600333200
    rxn.id <- which(mod@react_id=="Kt3r")
    if(length(rxn.id) == 1) {
      mod@lowbnd[rxn.id] <- 0
      mod@react_rev[rxn.id] <- FALSE
      
      mod <- addReact(mod,
                      id = "Kt1r",
                      met = c("k[e]","k[c]"),
                      Scoef = c(-1,1),
                      reversible = T,
                      reactName = "potassium transport via uniport reversible (facilitated diffusion)")
      mod@react_attr[which(mod@react_id=="Kt1r"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="Kt1r")] <- -1000
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 114 - Stochiometry in putrescine proton symport should be 1:2
    rxn.id <- which(mod@react_id=="PTRCt2")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 115 transport of rmn for 2 protons is wrong. Use the 1:1 variant
    rxn.id <- which(mod@react_id=="RMNt2_1")
    if(length(rxn.id) == 1) {
      mod <- addReact(mod,
                      id = "RMNt2",
                      met = c("h[e]","rmn[e]","h[c]","rmn[c]"),
                      Scoef = c(-1,-1,1,1),
                      reversible = T,
                      reactName = "L-rhamnose transport via proton symport")
      mod@react_attr[which(mod@react_id=="RMNt2"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
      mod@lowbnd[which(mod@react_id=="RMNt2")] <- -1000
      
      mod <- rmReact(mod,rxn.id)
      
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 116 - Stochiometry in putrescine proton symport should be 1:2
    rxn.id <- which(mod@react_id=="PTRCt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 117 Thioredoxin reductase (EC 1.8.1.9) should not be reversible
    rxn.id <- which(mod@react_id=="TRDRr")
    if(length(rxn.id) == 1) {
      mod@lowbnd[rxn.id] <- 0
      mod@react_rev[rxn.id] <- FALSE
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 118 - Stochiometry in cys proton symport should be 1:2
    rxn.id <- which(mod@react_id=="CYSt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  if(!(i %in% excl.rules)) {
    # ERROR 119 - Stochiometry in Meso-2,6-Diaminoheptanedioate proton symport should be 1:2
    rxn.id <- which(mod@react_id=="26DAPt2r")
    h.c.id <- which(mod@met_id=="h[c]")
    h.e.id <- which(mod@met_id=="h[e]")
    if(length(rxn.id) == 1) {
      mod@S[h.c.id,rxn.id] <- 2
      mod@S[h.e.id,rxn.id] <- -2
    }
  }
  i <-  i + 1
  
  
  
  
  
  # if(!(i %in% excl.rules)) {
  #   # NOT an error, but for FC-identification I need to remove K+ diffusion
  #   # remove this change for real model corrections
  #   rxn.id <- which(mod@react_id=="Kt1r")
  #   if(length(rxn.id) == 1) {
  #     mod <- addReact(mod,
  #                     id = "Kt3r",
  #                     met = c("h[e]","k[c]","h[c]","k[e]"),
  #                     Scoef = c(-1,-1,1,1),
  #                     reversible = F,
  #                     reactName = "potassium proton antiporter")
  #     mod@react_attr[which(mod@react_id=="Kt3r"),] <- c(annotation="",
  #                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
  #     mod@lowbnd[which(mod@react_id=="Kt3r")] <- 0
  #     
  #     mod <- rmReact(mod,rxn.id)
  #     
  #   }
  # }
  # i <-  i + 1
  
  
  
  # A few Gap-filling reactions are required to facilitate growth in a few corrected models
  
  if(mod@mod_desc %in% c("Clostridium_sporosphaeroides_DSM_1294",
                         "Stoquefichus_massiliensis_AP9",
                         "Streptococcus_pleomorphus_DSM_20574",
                         "Faecalibacterium_prausnitzii_A2_165")) {
    mod <- addReact(mod, id = "H2CO3D",
                    met = c("h2o[c]","co2[c]","h[c]","hco3[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "carboxylic acid dissociation")
    mod@react_attr[which(mod@react_id=="H2CO3D"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="H2CO3D")] <- -1000
    
  }
  
  if(mod@mod_desc == "Lactobacillus_animalis_KCTC_3501") {
    mod <- addReact(mod,
                    id = "OCBT",
                    met = c("cbp[c]","orn[c]","h[c]","citr_L[c]","pi[c]"),
                    Scoef = c(-1,-1,1,1,1),
                    reversible = T,
                    reactName = "ornithine carbamoyltransferase")
    mod@react_attr[which(mod@react_id=="OCBT"),] <- c(annotation="",
                                                      notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="OCBT")] <- -1000
    
  }
  
  if(mod@mod_desc == "Bacteroides_sp_1_1_30") {

    mod <- addReact(mod,
                    id = "NH4tb",
                    met = c("nh4[e]","nh4[c]"),
                    Scoef = c(-1,1),
                    reversible = T,
                    reactName = "NH4tb")
    mod@react_attr[which(mod@react_id=="NH4tb"),] <- c(annotation="",
                                                       notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="NH4tb")] <- -1000
    
    # HURRA!!! Diese Reaktion wars!!!
    mod <- addReact(mod,
                    id = "DTMPK",
                    met = c("atp[c]","dtmp[c]","adp[c]","dtdp[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = F,
                    reactName = "dTMP kinase")
    mod@react_attr[which(mod@react_id=="DTMPK"),] <- c(annotation="",
                                                       notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
  }
  
  if(mod@mod_desc == "Clostridium_sp_L2_50") {
    mod <- addReact(mod,
                    id = "CITL",
                    met = c("cit[c]","ac[c]","oaa[c]"),
                    Scoef = c(-1,1,1),
                    reversible = F,
                    reactName = "Citrate lyase")
    mod@react_attr[which(mod@react_id=="CITL"),] <- c(annotation="",
                                                      notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
    
  }
  
  if(mod@mod_desc %in% c("Lactobacillus_fermentum_ATCC_14931","Lactobacillus_fermentum_IFO_3956")) {
    mod <- addReact(mod,
                    id = "MALCOAPYRCT",
                    met = c("malcoa[c]","pyr[c]","accoa[c]","oaa[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Malonyl CoA pyruvate carboxytransferase")
    mod@react_attr[which(mod@react_id=="MALCOAPYRCT"),] <- c(annotation="",
                                                             notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="MALCOAPYRCT")] <- -1000
    
  }
  
  if(mod@mod_desc %in% c("Eggerthella_lenta_DSM_2243")) {
    mod <- addReact(mod,
                    id = "ILEt2r",
                    met = c("h[e]","ile_L[e]","h[c]","ile_L[c]"),
                    Scoef = c(-2,-1,2,1),
                    reversible = T,
                    reactName = "L-isoleucine reversible transport via proton symport")
    mod@react_attr[which(mod@react_id=="ILEt2r"),] <- c(annotation="",
                                                             notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="ILEt2r")] <- -1000
  }
  
  if(mod@mod_desc %in% c("Helicobacter_pullorum_MIT_98_5489")) {
    # 1
    mod <- addReact(mod,
                    id = "ACACT1r",
                    met = c("accoa[c]","aacoa[c]","coa[c]"),
                    Scoef = c(-2,1,1),
                    reversible = T,
                    reactName = "Acetyl Coenzyme A C-Acetyltransferase")
    mod@react_attr[which(mod@react_id=="ACACT1r"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Tryptophan metabolism</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="ACACT1r")] <- -1000
    
    # 2
    mod <- addReact(mod,
                    id = "ACACt2",
                    met = c("acac[e]","h[e]","acac[c]","h[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Acetoacetate Transport via Proton Symport")
    mod@react_attr[which(mod@react_id=="ACACt2"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="ACACt2")] <- -1000
    
    # 3
    mod <- addReact(mod,
                    id = "GHMT2r",
                    met = c("ser_L[c]","thf[c]","gly[c]","h2o[c]","mlthf[c]"),
                    Scoef = c(-1,-1,1,1,1),
                    reversible = T,
                    reactName = "Glycine Hydroxymethyltransferase, Reversible")
    mod@react_attr[which(mod@react_id=="GHMT2r"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Glycine, serine, alanine, and threonine metabolism</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="GHMT2r")] <- -1000
    
    # 4
    mod <- addReact(mod,
                    id = "OCOAT1r",
                    met = c("acac[c]","succoa[c]","aacoa[c]","succ[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Glycine Hydroxymethyltransferase, Reversible")
    mod@react_attr[which(mod@react_id=="OCOAT1r"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Valine, leucine, and isoleucine metabolism</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="OCOAT1r")] <- -1000
    
  }
  
  if(mod@mod_desc %in% c("Helicobacter_cinaedi_CCUG_18818")) {
    # 1
    mod <- addReact(mod,
                    id = "ACACT1r",
                    met = c("accoa[c]","aacoa[c]","coa[c]"),
                    Scoef = c(-2,1,1),
                    reversible = T,
                    reactName = "Acetyl Coenzyme A C-Acetyltransferase")
    mod@react_attr[which(mod@react_id=="ACACT1r"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Tryptophan metabolism</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="ACACT1r")] <- -1000
    
    # 2
    mod <- addReact(mod,
                    id = "ACACt2",
                    met = c("acac[e]","h[e]","acac[c]","h[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Acetoacetate Transport via Proton Symport")
    mod@react_attr[which(mod@react_id=="ACACt2"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="ACACt2")] <- -1000
    
    # 3
    mod <- addReact(mod,
                    id = "OCOAT1r",
                    met = c("acac[c]","succoa[c]","aacoa[c]","succ[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Glycine Hydroxymethyltransferase, Reversible")
    mod@react_attr[which(mod@react_id=="OCOAT1r"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Valine, leucine, and isoleucine metabolism</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="OCOAT1r")] <- -1000
    
  }
  
  if(mod@mod_desc %in% c("Helicobacter_winghamensis_ATCC_BAA_430")) {
    # 1
    mod <- addReact(mod,
                    id = "ACACT1r",
                    met = c("accoa[c]","aacoa[c]","coa[c]"),
                    Scoef = c(-2,1,1),
                    reversible = T,
                    reactName = "Acetyl Coenzyme A C-Acetyltransferase")
    mod@react_attr[which(mod@react_id=="ACACT1r"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Tryptophan metabolism</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="ACACT1r")] <- -1000
    
    # 2
    mod <- addReact(mod,
                    id = "ACACt2",
                    met = c("acac[e]","h[e]","acac[c]","h[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Acetoacetate Transport via Proton Symport")
    mod@react_attr[which(mod@react_id=="ACACt2"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="ACACt2")] <- -1000
    
    # 3
    mod <- addReact(mod,
                    id = "OCOAT1r",
                    met = c("acac[c]","succoa[c]","aacoa[c]","succ[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Glycine Hydroxymethyltransferase, Reversible")
    mod@react_attr[which(mod@react_id=="OCOAT1r"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Valine, leucine, and isoleucine metabolism</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="OCOAT1r")] <- -1000
    
    # 4 
    mod <- addReact(mod,
                    id = "SERt2r",
                    met = c("h[e]","ser_L[e]","h[c]","ser_L[c]"),
                    Scoef = c(-2,-1,2,1),
                    reversible = T,
                    reactName = "L-serine reversible transport via proton symport")
    mod@react_attr[which(mod@react_id=="SERt2r"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="SERt2r")] <- -1000
  }
  
  if(mod@mod_desc %in% c("Campylobacter_hominis_ATCC_BAA_381")) {
    # 1
    mod <- addReact(mod,
                    id = "26DAPt2r",
                    met = c("26dap_M[e]","h[e]","26dap_M[c]","h[c]"),
                    Scoef = c(-1,-2,1,2),
                    reversible = T,
                    reactName = "Meso-2,6-Diaminoheptanedioate transport via proton symport, reversible")
    mod@react_attr[which(mod@react_id=="26DAPt2r"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="26DAPt2r")] <- -1000
    
    # 2
    mod <- addReact(mod,
                    id = "CYSt2r",
                    met = c("cys_L[e]","h[e]","cys_L[c]","h[c]"),
                    Scoef = c(-1,-2,1,2),
                    reversible = T,
                    reactName = "L-cysteine reversible transport via proton symport")
    mod@react_attr[which(mod@react_id=="CYSt2r"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="CYSt2r")] <- -1000
    
    # 3
    mod <- addReact(mod,
                    id = "GLYt4r",
                    met = c("gly[e]","na1[e]","gly[c]","na1[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "L-cysteine reversible transport via proton symport")
    mod@react_attr[which(mod@react_id=="GLYt4r"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="GLYt4r")] <- -1000

  }
  
  if(mod@mod_desc %in% c("Methanobrevibacter_smithii_ATCC_35061")) {
    # 1
    mod <- addReact(mod,
                    id = "PEPCK",
                    met = c("gtp[c]","oaa[c]","co2[c]","gdp[c]","pep[c]"),
                    Scoef = c(-1,-1,1,1,1),
                    reversible = F,
                    reactName = "Phosphoenolpyruvate Carboxykinase (GTP)")
    mod@react_attr[which(mod@react_id=="PEPCK"),] <- c(annotation="",
                                                          notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Glycolysis/gluconeogenesis</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
    # 2
    mod <- addReact(mod,
                    id = "PFL",
                    met = c("coa[c]","pyr[c]","accoa[c]","for[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = F,
                    reactName = "Pyruvate formate lyase")
    mod@react_attr[which(mod@react_id=="PFL"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Pyruvate metabolism</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
    # 3
    mod <- addReact(mod,
                    id = "PYK",
                    met = c("adp[c]","h[c]","pep[c]","atp[c]","pyr[c]"),
                    Scoef = c(-1,-1,-1,1,1),
                    reversible = F,
                    reactName = "Pyruvate Kinase")
    mod@react_attr[which(mod@react_id=="PYK"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Glycolysis/gluconeogenesis</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
  }
  
  if(mod@mod_desc %in% c("Gordonibacter_pamelaeae_7_10_1_bT_DSM_19378")) {
    # 1
    mod <- addReact(mod,
                    id = "PPCK",
                    met = c("atp[c]","oaa[c]","adp[c]","co2[c]","pep[c]"),
                    Scoef = c(-1,-1,1,1,1),
                    reversible = F,
                    reactName = "Phosphoenolpyruvate carboxykinase")
    mod@react_attr[which(mod@react_id=="PPCK"),] <- c(annotation="",
                                                       notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Citric acid cycle</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
    # 2
    mod <- addReact(mod,
                    id = "XPPT",
                    met = c("prpp[c]","xan[c]","ppi[c]","xmp[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = F,
                    reactName = "Xanthine phosphoribosyltransferase")
    mod@react_attr[which(mod@react_id=="XPPT"),] <- c(annotation="",
                                                     notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Purine catabolism<</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
    # 3
    mod <- addReact(mod,
                    id = "PYK",
                    met = c("adp[c]","h[c]","pep[c]","atp[c]","pyr[c]"),
                    Scoef = c(-1,-1,-1,1,1),
                    reversible = F,
                    reactName = "Pyruvate Kinase")
    mod@react_attr[which(mod@react_id=="PYK"),] <- c(annotation="",
                                                     notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Glycolysis/gluconeogenesis</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
  }
  
  if(mod@mod_desc %in% c("Mogibacterium_timidum_ATCC_33093")) {
    # 1
    mod <- addReact(mod,
                    id = "PYDXabc",
                    met = c("atp[c]","h2o[c]","pydx[e]","adp[c]","h[c]","pydx[c]"),
                    Scoef = c(-1,-1,-1,1,1,1),
                    reversible = F,
                    reactName = "Pyridoxal transport via ECF transporter")
    mod@react_attr[which(mod@react_id=="PYDXabc"),] <- c(annotation="",
                                                      notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
    # 2
    mod <- addReact(mod,
                    id = "PYK",
                    met = c("adp[c]","h[c]","pep[c]","atp[c]","pyr[c]"),
                    Scoef = c(-1,-1,-1,1,1),
                    reversible = F,
                    reactName = "Pyruvate Kinase")
    mod@react_attr[which(mod@react_id=="PYK"),] <- c(annotation="",
                                                     notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Glycolysis/gluconeogenesis</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
  }
  
  if(mod@mod_desc %in% c("Bilophila_wadsworthia_3_1_6")) {
    # 1
    mod <- addReact(mod,
                    id = "TALA",
                    met = c("g3p[c]","s7p[c]","e4p[c]","f6p[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Transaldolase")
    mod@react_attr[which(mod@react_id=="TALA"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Pentose phosphate pathway</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="TALA")] <- -1000
    
    # 2
    mod <- addReact(mod,
                    id = "r0220",
                    met = c("coa[c]","HC01668[c]","amp[c]","ppcoa[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Propinol Adenylate:Coa Ligase (AMP-Forming) Ec:6.2.1.17 Ec:6.2.1.1")
    mod@react_attr[which(mod@react_id=="r0220"),] <- c(annotation="",
                                                     notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Miscellaneous</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="r0220")] <- -1000
    
  }
  
  if(mod@mod_desc %in% c("Brevibacterium_massiliense_5401308")) {
    # 1
    mod <- addReact(mod,
                    id = "RBK_Dr",
                    met = c("atp[c]","rbl_D[c]","adp[c]","h[c]","ru5p_D[c]"),
                    Scoef = c(-1,-1,1,1,1),
                    reversible = T,
                    reactName = "D-ribulokinase (reversible)")
    mod@react_attr[which(mod@react_id=="RBK_Dr"),] <- c(annotation="",
                                                      notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Pentose phosphate pathway</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="RBK_Dr")] <- -1000
    
    # 2
    mod <- addReact(mod,
                    id = "SHCHCS3",
                    met = c("2s5epchc[c]","2shchc[c]","pyr[c]"),
                    Scoef = c(-1,1,1),
                    reversible = F,
                    reactName = "2-succinyl-6-hydroxy-2,4-cyclohexadiene-1-carboxylate synthase")
    mod@react_attr[which(mod@react_id=="SHCHCS3"),] <- c(annotation="",
                                                       notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Ubiquinone and other terpenoid-quinone biosynthesis</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    
  }
  
  if(mod@mod_desc %in% c("Fusobacterium_gonidiaformans_3_1_5R","Fusobacterium_gonidiaformans_ATCC_25563")) {
    # 1
    mod <- addReact(mod,
                    id = "ADNCNT3tc",
                    met = c("adn[e]","h[e]","adn[c]","h[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Adenosine Transport Cnt3")
    mod@react_attr[which(mod@react_id=="ADNCNT3tc"),] <- c(annotation="",
                                                        notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="ADNCNT3tc")] <- -1000
    
    # 2
    mod <- addReact(mod,
                    id = "DADNt2r",
                    met = c("dad_2[e]","h[e]","dad_2[c]","h[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Deoxyadenosine transport in via proton symport reversible")
    mod@react_attr[which(mod@react_id=="DADNt2r"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="DADNt2r")] <- -1000
    
    # 3
    mod <- addReact(mod,
                    id = "DRPAr",
                    met = c("2dr5p[c]","acald[c]","g3p[c]"),
                    Scoef = c(-1,1,1),
                    reversible = T,
                    reactName = "2-Deoxy-D-ribose-5-phosphate acetaldehyde-lyase")
    mod@react_attr[which(mod@react_id=="DRPAr"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Pentose phosphate pathway</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="DRPAr")] <- -1000
    
    # 4
    mod <- addReact(mod,
                    id = "INSt2",
                    met = c("h[e]","ins[e]","h[c]","ins[c]"),
                    Scoef = c(-1,-1,1,1),
                    reversible = T,
                    reactName = "Inosine Transport in via Proton Symport, Reversible")
    mod@react_attr[which(mod@react_id=="INSt2"),] <- c(annotation="",
                                                         notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Transport, extracellular</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="INSt2")] <- -1000
    
    # 5
    mod <- addReact(mod,
                    id = "r0570",
                    met = c("2dr1p[c]","2dr5p[c]"),
                    Scoef = c(-1,1),
                    reversible = T,
                    reactName = "2-Deoxy-D-Ribose 1-Phosphate 1, 5-Phosphomutase Ec:5.4.2.7")
    mod@react_attr[which(mod@react_id=="r0570"),] <- c(annotation="",
                                                       notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: Pentose phosphate pathway</p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
    mod@lowbnd[which(mod@react_id=="r0570")] <- -1000
  }
  
  # Add missing exchange reactions
  ex.mets.ind  <- grep("[e]",mod@met_id,fixed = T)
  ex.mets.ids  <- mod@met_id[ex.mets.ind]
  ex.mets.name <- mod@met_name[ex.mets.ind]
  
  present.exchanges <- str_match(mod@react_id[grep("^EX_.*\\(e\\)",mod@react_id)],"EX_(.*?)\\(")
  if(ncol(present.exchanges)==2){
    ind.new <- !(gsub("\\[.*\\]","",ex.mets.ids) %in% present.exchanges[,2])
    ex.mets.ind  <- ex.mets.ind[ind.new]
    ex.mets.ids  <- ex.mets.ids[ind.new]
    ex.mets.name <- ex.mets.name[ind.new]
  }
  
  if(length(ex.mets.ids)>0) {
    for(i in 1:length(ex.mets.ids)) {
      #cat("Adding new exchange reaction for: ",ex.mets.ids[i],"\n")
      mod <- addReact(model = mod,
                      id = paste0("EX_",gsub("\\[.*\\]","",ex.mets.ids[i]),"(e)"), 
                      met = ex.mets.ids[i],
                      Scoef = -1,
                      reversible = T, 
                      metComp = 2,
                      ub = 1000,
                      lb = 0,
                      reactName = paste0(ex.mets.name[i], " exchange"), 
                      metName = ex.mets.name[i])
    }
    mod@react_attr[which(mod@react_id==paste0("EX_",gsub("\\[.*\\]","",ex.mets.ids[i]),"(e)")),] <- c(annotation="",
                                                             notes="<notes>\n  <body xmlns=\"http://www.w3.org/1999/xhtml\">\n    <p>SUBSYSTEM: </p>\n    <p>Confidence Level: </p>\n  </body>\n</notes>")
  }
  
  return(mod)
}

sync_constraints <- function(agora) {
  ex.list <- list()
  for(i in agora) {
    ex.rxn.id <- grep("^EX_",i@react_id)
    dt <- data.table(rxn = i@react_id[ex.rxn.id],
                     lb  = i@lowbnd[ex.rxn.id],
                     mod = i@mod_desc)
    ex.list[[i@mod_desc]] <- dt
  }
  ex.list <- rbindlist(ex.list)
  ex.list[,min.lb := min(lb), by = "rxn"]
  
  # those entries indicate required changes:
  dt.cg <- ex.list[lb!=min.lb]
  print(dt.cg)
  if(nrow(dt.cg)>0) {
    for(i in 1:nrow(dt.cg)) {
      ind <- which(agora[[dt.cg[i,mod]]]@react_id == dt.cg[i,rxn])
      agora[[dt.cg[i,mod]]]@lowbnd[ind] <- dt.cg[i,min.lb]
    }
  }
  return(agora)
}
