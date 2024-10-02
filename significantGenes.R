
deseq <- read_csv('mouseedgeRcomp_updated.csv', show_col_types = FALSE)
deseq <- deseq[complete.cases(deseq), ]

pvalColumns <- c("Control.F-Acarbose.F_PValue", "Control.M-Acarbose.M_PValue", "Control.F-CR.F_PValue", 
                 "Control.M-CR.M_PValue", "GHRKOControl.M-GHRKO.M_PValue", "Control.F-Rapamycin.F_PValue", 
                 "Control.M-Rapamycin.M_PValue", "SnellDwarfMiceControl.M-SnellDwarfMice.M_PValue", 
                 "Control.F-Protandim.F_PValue", "Control.M-Protandim.M_PValue", "Control.F-alphaestradiol.F_PValue", 
                 "Control.M-alphaestradiol.M_PValue", "MRControl.M-MethionineRestriction.M_PValue")
statColumns <- c("Control.F-Acarbose.F_logCPM", "Control.M-Acarbose.M_logCPM", "Control.F-CR.F_logCPM", 
                 "Control.M-CR.M_logCPM", "GHRKOControl.M-GHRKO.M_logCPM", "Control.F-Rapamycin.F_logCPM", 
                 "Control.M-Rapamycin.M_logCPM", "SnellDwarfMiceControl.M-SnellDwarfMice.M_logCPM", 
                 "Control.F-Protandim.F_logCPM", "Control.M-Protandim.M_logCPM", "Control.F-alphaestradiol.F_logCPM", 
                 "Control.M-alphaestradiol.M_logCPM", "MRControl.M-MethionineRestriction.M_logCPM")
tablePval <- deseq[,pvalColumns]

pvalCheck <- tablePval < 0.05
rowsToKeep <- rowSums(pvalCheck == FALSE) == 0
deseq <- deseq[rowsToKeep, ]
