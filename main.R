#### Import training data and script ####
BRT.brACA   <- readRDS("Import/Training/BRT_brACA.Rds")
BRT.brACA.karid   <- readRDS("Import/Training/BRT_brACA_karid.Rds")
BRT.brACA.kwet   <- readRDS("Import/Training/BRT_brACA_kwet.Rds")
M.br.GDGT <- readRDS("Import/Training/M_brGDGT_ACADB.Rds")
M.br.GDGT.karid <- readRDS("Import/Training/M_brGDGT_ACADB_karid.Rds")
M.br.GDGT.kwet <- readRDS("Import/Training/M_brGDGT_ACADB_kwet.Rds")
BRT.param <- data.frame(read.csv(file="Import/Training/BRT_param.csv", sep=",",dec=".",header=T,row.names=1))        # GDGT indexe Ayrag
Cluster.prediction.ACADB.brGDGT <- readRDS("Import/Training/Cluster.prediction.ACADB.brGDGT.Rds")
MC <- readRDS("Import/Training/M_clim_ACADB.Rds")
MC.karid <- readRDS("Import/Training/M_clim_ACADB_karid.Rds")
MC.kwet <- readRDS("Import/Training/M_clim_ACADB_kwet.Rds")
set.seed(123)
source("Import/Script/BRT_script.R")

#### Import paleo data ####
GDGT.XRD <- data.frame(read.csv(file="Import/Paleo/XRD.csv", sep=",",dec=".",header=T,row.names=1))        # GDGT indexe Ayrag
XRD.metadata  <- data.frame(read.csv(file="Import/Paleo/XRD_metadata.csv", sep=",",dec=".",header=T,row.names=1))        # GDGT indexe Ayrag

#### Verify if the Fractional abundances of the cores are OK ####
GDGT.histo.plot.surf.core(Mcore = GDGT.XRD, Remove.7Me = T, W = 1600, H = 700, Save.path = "Figures/Hist_brGDGT_XRD.pdf")

#### BRT predictions ####
GDGT.XRD.conv <- GDGT.XRD[which(names(GDGT.XRD) %in% names(M.br.GDGT))]
XRD.brACA <- FT.core(Model.BRT = BRT.brACA, MCore = GDGT.XRD.conv, MAge = GDGT.XRD$Age, Save.path = "Results/XRD.csv")
XRD.brACA.karid <- FT.core(Model.BRT = BRT.brACA.karid, MCore = GDGT.XRD.conv, MAge = GDGT.XRD$Age, Save.path = "Results/XRD.csv")
XRD.brACA.kwet <- FT.core(Model.BRT = BRT.brACA.kwet, MCore = GDGT.XRD.conv, MAge = GDGT.XRD$Age, Save.path = "Results/XRD.csv")


#### Combined model ####
XRD.brACA <- readRDS("Results/XRD_BRT_brACA.Rds")[[2]]
XRD.brACA.karid <- readRDS("Results/XRD_BRT_karid.Rds")[[2]]
XRD.brACA.kwet <- readRDS("Results/XRD_BRT_kwet.Rds")[[2]]

pBRT.XRD <- Combine.ML.cluster(
  List.models = list(M1 = XRD.brACA, M2 = XRD.brACA.karid, M3 = XRD.brACA.kwet),
  Model.lab = c("ACADB", "K-warm/arid", "K-cold/wet"),
  Cluster.prediction = Cluster.prediction.ACADB.brGDGT,
  GDGT.paleo = GDGT.XRD, Time.in.k = T, Time.res = 1000, Time.lim = c(0, 8000),
  Highlight.combined = F, Facet = T, Only.best = F,
  Save.plot = "Figures/XRD_combined.pdf", H = 700, W = 500,
  Save.path = "Results/XRD_brACA_combined.Rds"
)

#### Test randomTF() ####
test.randomTF = F
if(test.randomTF == T){
  XRD.rTF.ACADB <- Plot.randomTF(MPsurf = M.br.GDGT, MPpaleo = GDGT.XRD.conv,  Mclim = MC, Database = "ACADB", Lake = "XRD", Plot.MAT = F, Plot.WAPLS = F, 
                                 Save.path = "Results/XRD_randomFT_ACADB.Rds", return.plot = T, Plot.RF = T, Plot.BRT = T,
                                 BRT.max_trees = mean(BRT.param$n_.trees.), BRT.learning_rate = mean(BRT.param$Learning.rate), BRT.tree_complexity = mean(BRT.param$Tree.complexity), BRT.minobs = mean(BRT.param$Minimal.number.of.observations), BRT.bag_fraction = 0.75,
                                 H = 500, W = 2000, Save.plot = "Figures/XRD_randomTF_ACADB.pdf")
  
  XRD.rTF.ACADB.comb <- Plot.randomTF(MPsurf = M.br.GDGT.karid, MPpaleo = GDGT.XRD.conv,  Mclim = MC.karid, Database = "ACADB-combined", Lake = "XRD", Plot.MAT = F, Plot.WAPLS = F, 
                                      BRT.max_trees = mean(BRT.param$n_.trees.), BRT.learning_rate = mean(BRT.param$Learning.rate), BRT.tree_complexity = mean(BRT.param$Tree.complexity), BRT.minobs = mean(BRT.param$Minimal.number.of.observations), BRT.bag_fraction = 0.75,
                                      Save.path = "Results/XRD_randomFT_ACADB_comb.Rds", return.plot = T, Plot.RF = T, Plot.BRT = T,
                                      H = 500, W = 2000, Save.plot = "Figures/XRD_randomTF_ACADB_comb.pdf")
  }


