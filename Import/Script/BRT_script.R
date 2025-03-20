#### Library ####
library(reshape2)
library(ggplot2)
library(patchwork)
library(randomForest)
library(dplyr)
library(palaeoSig)
library(caret)

#### Functions  ####

# This function plot the boxplot for each br-GDGT fractional abundance
# If you add Mtype, it make the groups for moss, soils and sediment...
# Graph présent chez Ding et al. 2015, Hopmans et al. 2004
# Si ajout iso.GDGT, Hopmans et al. 2004 
GDGT.histo.plot.surf.core <- function(Mcore, Msurf, Mtype, Select.type, Show.Plotly, Reorder.group, Global.box = F, Annot.size = 6,
                                      Keep.br, Leg.nb.lines, Leg.iso = F, Box.linewidth = .5, Leg.size = 13,
                                      Remove.8Me, Remove.7Me, Color.choice, Iso.GDGT, Leg.pos, Return.plot = F, Boxplot.title = NULL, Leg.box = F,
                                      Zoom1.comp = NULL, Zoom2.comp = NULL, Zoom1.Ymax, Zoom2.Ymax, Insert1.loc = NULL, Insert2.loc = NULL,
                                      Name.untype, Ymax, Save.path, W, H, Dot.pop, Overlap.OK){
  #### Initialization values ####
  if(missing(Mtype)){Mtype = NULL}
  if(missing(Msurf)){Msurf = NULL}
  if(missing(Mcore)){Mcore = NULL}
  if(missing(Reorder.group)){Reorder.group = NULL}
  if(missing(Show.Plotly)){Show.Plotly = F}
  if(missing(Save.path)){Save.path = NULL}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Keep.br)){Keep.br = NULL}
  if(missing(Leg.nb.lines)){Leg.nb.lines = NULL}
  if(missing(Dot.pop)){Dot.pop = NULL}
  if(missing(Overlap.OK)){Overlap.OK = F}
  if(missing(Name.untype)){Name.untype = "Other"}
  if(missing(Select.type)){Select.type = "Sample.type"}
  if(missing(Color.choice)){Color.choice = NULL}
  if(missing(Remove.8Me)){Remove.8Me = F}
  if(missing(Remove.7Me)){Remove.7Me = F}
  if(missing(Iso.GDGT)){Iso.GDGT = F}
  if(missing(Leg.pos)){Leg.pos = c(0.31,.86)}
  
  #### Extract br-GDGT ####
  MBRcore <- Mcore[,grep("^f.I", colnames(Mcore))]
  MBRsurf <- Msurf[,grep("^f.I", colnames(Msurf))]
  
  #### Remove 8 Me ####
  if(Remove.8Me == T){
    MBRcore <- MBRcore[,!grepl("_8Me", colnames(MBRcore))]
    MBRsurf <- MBRsurf[,!grepl("_8Me", colnames(MBRsurf))]
  }
  
  #### Remove 7 Me ####
  if(Remove.7Me == T){
    MBRcore <- MBRcore[,!grepl("_7Me", colnames(MBRcore))]
    MBRsurf <- MBRsurf[,!grepl("_7Me", colnames(MBRsurf))]
  }
  
  if(is.null(MBRsurf) == F){MBRsurf <- MBRsurf[setdiff(names(MBRsurf), names(MBRsurf)[which(colSums(MBRsurf) == 0)])]}
  if(is.null(MBRcore) == F){MBRcore <- MBRcore[setdiff(names(MBRcore), names(MBRcore)[which(colSums(MBRcore) == 0)])]}
  
  #### Test if the names are egual ####
  if(is.null(MBRcore) == F & is.null(MBRsurf) == F){
    '%nin%' <- Negate('%in%')
    test1 <- length(names(MBRcore)[names(MBRcore) %nin% names(MBRsurf)])
    test2 <- length(names(MBRsurf)[names(MBRsurf) %nin% names(MBRsurf)])
    
    if(test1 > 0 & test2 > 0){
      print("Names are different.")
      break
    }
    if(test1 == 0 & test2 == 0){
      print("Let's go !")
      Mplot <- rbind(MBRcore, MBRsurf)   
    }
    else{print("Some molecules are missing !")
      print(setdiff(names(MBRsurf), names(MBRcore)))
      print(setdiff(names(MBRcore), names(MBRsurf)))
    }
  }
  else{
    if(is.null(MBRcore) == F){Mplot = MBRcore[,grep("^f.I", colnames(MBRcore))]}
    if(is.null(MBRsurf) == F){Mplot = MBRsurf[,grep("^f.I", colnames(MBRsurf))]}
  }
  
  #### Extract iso-GDGT ####
  if(Iso.GDGT == T){
    #### Recuperation des iso GDGT dans M tot ####
    if(is.null(Mcore) == F){
      Mcore.iso <- cbind(Mcore[,grep("^GDGT", colnames(Mcore))], Crenarch = Mcore$Crenarch, Crenarch.p = Mcore$Crenarch.p)
      if(is.null(Mcore.iso$GDGT0.Crenar) == F){Mcore.iso <- subset(Mcore.iso, select = -c(GDGT0.Crenar))} 
      Mcore.iso <- Mcore.iso/rowSums(Mcore.iso)*100
      
      #### Stats iso-GDGTs ####
    }
    
    if(is.null(Msurf) == F){
      Msurf.iso <- cbind(Msurf[,grep("^GDGT", colnames(Msurf))], Crenarch = Msurf$Crenarch, Crenarch.p = Msurf$Crenarch.p)
      if(is.null(Msurf.iso$GDGT0.Crenar) == F){Msurf.iso <- Msurf.iso[!grepl("GDGT0.Crenar", names(Msurf.iso))]}
      Msurf.iso <- Msurf.iso/rowSums(Msurf.iso)*100
    }
    
    if(is.null(Mcore) == F & is.null(Msurf) == F){
      Ymax.iso = max(max(Mcore.iso), max(Msurf.iso))
      L.iso = rbind(Msurf.iso, Mcore.iso)}
    if(is.null(Mcore) == T & is.null(Msurf) == F){
      Ymax.iso = max(Msurf.iso)
      L.iso = Msurf.iso}
    if(is.null(Mcore) == F & is.null(Msurf) == T){
      Ymax.iso = max(Mcore.iso)
      L.iso = Mcore.iso}
    #### Mise en forme des labels des axes ####
    if(is.null(Mtype) == F){
      Mtype <- Mtype[Select.type]
      names(Mtype) <- "Sample.type"
      A.iso <- merge(L.iso, Mtype, all.x = T, by = 0, sort = F)
      row.names(A.iso)<- A.iso$Row.names
      A.iso <- A.iso[,-1]
      levels(A.iso$Sample.type)[length(levels(A.iso$Sample.type)) + 1] <- Name.untype
      A.iso$Sample.type[is.na(A.iso$Sample.type)] <- Name.untype
    }
    else{ 
      A.iso <- L.iso 
      A.iso["Sample.type"]="Samples"}
    
    
    
    Model.lab <- gsub("GDGT", "GDGT-", names(A.iso))
    # Model.lab <- gsub("0", "0]", Model.lab)
    # Model.lab <- gsub("1", "1]", Model.lab)
    # Model.lab <- gsub("2", "2]", Model.lab)
    # Model.lab <- gsub("3", "3]", Model.lab)
    # Model.lab <- gsub("4", "4]", Model.lab)
    #Model.lab <- Model.lab[-1]                  # on enleve Row.names()
    Model.lab <- Model.lab[-length(Model.lab)]     # on enleve Samples.type
    Model.lab <- Model.lab[-length(Model.lab)]     # on enleve Crenarch.p et plus tard on rajoute Crenarch'
    
    #### If little pop to dot activated ####
    if(is.null(Dot.pop) == F){ # Echantillons choisis comme points
      Little.pop.iso <- A.iso[A.iso$Sample.type %in% Dot.pop,]
      for(i in 1:length(Dot.pop)){
        New.lab.iso <- unique(A.iso$Sample.type)[grepl(Dot.pop[i], unique(A.iso$Sample.type))]
        Little.pop.iso$Sample.type <- gsub(Dot.pop[i], New.lab.iso, Little.pop.iso$Sample.type)
        if(Overlap.OK == T){A.iso[A.iso$Sample.type %in% New.lab.iso, names(A.iso[,-length(names(A.iso))])] <- 100}}
      
      Not.in.dot.iso <- setdiff(unique(as.character(A.iso$Sample.type)), Dot.pop)
      New.row.0.iso <- data.frame(Sample.type = Not.in.dot.iso)
      Little.pop.iso <- rbind.fill(Little.pop.iso, New.row.0.iso)
      Little.pop.iso[is.na(Little.pop.iso)] <- 100
      A.little.iso <- melt(Little.pop.iso, id ='Sample.type')
      Add.points.size.lim.iso <- geom_dotplot(data = A.little.iso, aes(x=variable, y=value, fill=Sample.type),
                                              binaxis='y', stackdir='center', 
                                              position = position_dodge(.85),
                                              show.legend = FALSE,
                                              binwidth = 1)
      A.iso$Sample.type <- as.character(A.iso$Sample.type)
      A.iso <- melt(A.iso, id ='Sample.type')}
    
    else{                   # Pas échatillon comme points
      Add.points.size.lim.iso <- NULL
      A.iso$Sample.type <- as.character(A.iso$Sample.type)
      A.iso <- melt(A.iso, id ='Sample.type')
    }
  }
  else{A.iso = NULL}
  
  
  #### Label type ####
  Mplot <- 100*Mplot
  AZER <- gsub("f.", "", names(Mplot))
  AZER <- gsub("_5Me", "", AZER)
  AZER <- gsub("_6Me", "\\'", AZER)
  AZER <- gsub("_7Me", "\\''", AZER)
  AZER <- gsub("_8Me", "\\'''", AZER)
  names(Mplot)  <- AZER
  if(missing(Ymax)){Ymax = max(Mplot)}
  
  #### Merging datas ####
  if(is.null(Mtype) == F){
    if(Iso.GDGT == F){
      Mtype <- Mtype[Select.type]
      names(Mtype) <- "Sample.type"
    }
    A <- merge(Mplot, Mtype, all.x = T, by = 0, sort = F)
    row.names(A)<- A$Row.names
    A <- A[,-1]
    levels(A$Sample.type)[length(levels(A$Sample.type)) + 1] <- Name.untype
    A[is.na(A)] <- Name.untype
  }
  
  else{ 
    A = Mplot 
    A["Sample.type"]="Samples"}
  
  #### Counting n of each samples ####
  M.type <- as.character(A$Sample.type)
  N.type <- unique(M.type)
  NT <- sapply(N.type, function(x) length(M.type[M.type == x]))
  A.save <- A
  for(i in 1:length(NT)){A$Sample.type <- gsub(names(NT)[i], paste(names(NT)[i], ", n = ", NT[i], sep = ""), A$Sample.type)}
  A$Sample.type <- as.factor(A$Sample.type)
  
  if(Iso.GDGT == T & Leg.iso == T){
    for(i in 1:length(NT)){A.iso$Sample.type <- gsub(names(NT)[i], paste(names(NT)[i], ", n = ", NT[i], sep = ""), A.iso$Sample.type)}
    A.iso$Sample.type <- as.factor(A.iso$Sample.type)
  }
  
  #### Little sample size br-GDGT ####
  if(is.null(Dot.pop) == F){ # Echantillons choisis comme points
    Little.pop <- A.save[A.save$Sample.type %in% Dot.pop,]
    for(i in 1:length(Dot.pop)){
      New.lab <- unique(A$Sample.type)[grepl(Dot.pop[i], unique(A$Sample.type))]
      Little.pop$Sample.type <- gsub(Dot.pop[i], New.lab, Little.pop$Sample.type)
      if(Overlap.OK == T){A[A$Sample.type %in% New.lab, names(A[,-length(names(A))])] <- 100}}
    
    Not.in.dot <- setdiff(levels(A$Sample.type), Dot.pop)
    New.row.0 <- data.frame(Sample.type = Not.in.dot)
    Little.pop <- rbind.fill(Little.pop, New.row.0)
    Little.pop[is.na(Little.pop)] <- 100
    
    B.little <- melt(Little.pop, id ='Sample.type')
    
    if(is.null(Keep.br) == F){
      Keep.br.little <- levels(B.little$variable)[Keep.br] 
      B.little <- B.little[which(B.little$variable %in% Keep.br.little),]}
    
    Add.points.size.lim <- geom_dotplot(data = B.little, aes(x=variable, y=value, fill=Sample.type),
                                        binaxis='y', stackdir='center', 
                                        position = position_dodge(.85),
                                        show.legend = FALSE,
                                        binwidth = .5)
    B <- melt(A, id ='Sample.type')
  }
  
  else{                   # Pas échatillon comme points
    Add.points.size.lim <- NULL
    B <- melt(A, id ='Sample.type')
  }
  
  #### Color settings ####
  if(is.null(Color.choice) == T){
    if(length(NT) == 1){Color.vec <- c("grey")}
    if(length(NT) == 2){Color.vec <- c("#74D43B", "#1EAEAE")}
    if(length(NT) == 3){Color.vec <- c("#74D43B", "#1EAEAE", "#F3C643")}
    if(length(NT) == 4){Color.vec <- c("#74D43B", "#1EAEAE", "#F3C643", "grey")}
    if(length(NT) == 5){Color.vec <- c("#74D43B", "#1EAEAE", "grey", "#F3C643",  "darkorange")}
    if(length(NT) == 6){Color.vec <- c("#74D43B", "#1EAEAE", "grey", "#F3C643",  "darkorange", "indianred2")}
    if(length(NT) == 7){Color.vec <- c("#74D43B", "#1EAEAE", "grey", "#F3C643",  "darkorange", "indianred2", "darksalmon")}
    if(length(NT) == 8){Color.vec <- c("#74D43B", "#1EAEAE", "grey", "#F3C643",  "darkorange", "indianred2", "darksalmon", "darkred")}
    if(length(NT) > 8){Color.vec <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(NT)}}
  else{Color.vec <- Color.choice}
  
  #### Other graphical settings ####
  if(Global.box == T){My_border <- element_rect(fill = NA)}
  else{My_border <- element_blank()}
  
  if(is.null(Boxplot.title) == F){My_title <- ggtitle(Boxplot.title)}
  else{My_title <- NULL}
  
  if(Leg.box == T){My_legbox <- element_rect(fill = "white", color = "grey")}
  else{My_legbox <- NULL}
  
  #### Change order of boxplot ####
  if(is.null(Reorder.group) == F){
    if(Iso.GDGT == T){
      Old.order.name <- unique(A.iso$Sample.type)[order(unique(A.iso$Sample.type))]
      A.iso$Sample.type <- factor(A.iso$Sample.type, levels = Old.order.name[Reorder.group], ordered = T)}
    
    B$Sample.type <- factor(B$Sample.type, levels = levels(B$Sample.type)[Reorder.group], ordered = T)
    Color.vec <- Color.vec[Reorder.group]}
  
  #### brGDGTs last settings ####
  if(Iso.GDGT == F){A.lab = ylab("Fractional Abundance (%)")}
  if(Iso.GDGT == T){A.lab = ylab("")}
  
  if(is.null(Keep.br) == F){
    Keep.br <- levels(B$variable)[Keep.br] 
    B <- B[which(B$variable %in% Keep.br),]}
  
  if(is.null(Leg.nb.lines) == F){
    Leg.guide <- guides(fill = guide_legend(nrow = Leg.nb.lines))
  }
  else{Leg.guide <- NULL}  
  
  if(Iso.GDGT == T & Leg.iso == T){
    Leg.pos.iso <- Leg.pos
    Leg.pos <- "none"
  }
  else{Leg.pos.iso <- "none"}
  
  #### Add ligne en tirets ####
  Xligne <- c(
    max(grep("III", levels(B$variable)))+0.5,
    max(grep("II", levels(B$variable)))+0.5)
  
  Xtext <- c(Xligne[1]/2-0.25, (Xligne[1]+Xligne[2])/2-0.25, (Xligne[2]+nlevels(B$variable))/2)
  My_lab <- c("Hexa.", "Penta.", "Tetra.")
  
  #### Insert LR ####
  get_inset <- function(B, Zoom1.Ymax, My_annot){
    pinsert <- ggplot(B, aes(x=variable, y=value, fill = Sample.type)) + 
      coord_cartesian(ylim = c(0, Zoom1.Ymax)) +
      geom_boxplot(outlier.shape = NA, position = position_dodge(.85), linewidth = Box.linewidth) +
      Add.points.size.lim + Leg.guide +
      scale_fill_manual(name = "Sample type", values = Color.vec
                        , labels = c("0" = "Foo", "1" = "Bar"))+
      annotate("text", x = 0.8, y = Zoom1.Ymax - 0.1*Zoom1.Ymax, label = My_annot, size = 5, fontface = 2)+
      theme(legend.position = "none", panel.background = element_blank(),
            legend.key = element_blank(), legend.background = element_blank(),
            axis.title = element_blank(),
            axis.text = element_text(size = 12), panel.grid = element_blank(),
            plot.margin = unit(x = c(1, 2, 2, 0),units="mm"), # Bas / Droite / Haut / Gauche
            axis.line = element_blank(),
            legend.text = element_text(size = 13),
            panel.border = element_rect(fill = NA),
            legend.title = element_text(size = 14),
            plot.background = element_blank())+
      A.lab
    return(pinsert)
  }
  
  if(is.null(Zoom1.comp) == F){
    B.insert <- B[which(B$variable %in% levels(B$variable)[Zoom1.comp]),]
    if(missing(Zoom1.Ymax)){Zoom1.Ymax = max(B.insert)}
    if(missing(Insert1.loc)){Insert1.loc = c(3, Xligne[1]-0.25, 10, 40)}
    
    My_insert1 <- get_inset(B.insert, Zoom1.Ymax, "B1")
    My_insert1 <- annotation_custom(ggplotGrob(My_insert1), xmin = Insert1.loc[1], xmax = Insert1.loc[2], ymin = Insert1.loc[3], ymax = Insert1.loc[4])
    Box1 <- geom_rect(xmin = min(Zoom1.comp)-0.5, xmax = Insert1.loc[2], ymin = -1, ymax = Zoom1.Ymax*2, fill = NA,
                      alpha = 1, color = "grey30", linewidth = .3, linetype = 2)
  }
  else{My_insert1 = NULL; Box1 <- NULL}
  
  if(is.null(Zoom2.comp) == F){
    B.insert <- B[which(B$variable %in% levels(B$variable)[Zoom2.comp]),]
    if(missing(Zoom2.Ymax)){Zoom2.Ymax = max(B.insert)}
    if(missing(Insert2.loc)){Insert2.loc = c(11, Xligne[2]-0.25, 10, 40)}
    
    My_insert2 <- get_inset(B.insert, Zoom2.Ymax, "B2")
    My_insert2 <- annotation_custom(ggplotGrob(My_insert2), xmin = Insert2.loc[1], xmax = Insert2.loc[2], ymin = Insert2.loc[3], ymax = Insert2.loc[4])
    Box2 <- geom_rect(xmin = min(Zoom2.comp)-0.5, xmax = Insert2.loc[2], ymin = -1, ymax = Zoom2.Ymax*2, fill = NA,
                      alpha = 1, color = "grey30", linewidth = .3, linetype = 2)
  }
  else{My_insert2 = NULL; Box2 <- NULL}
  
  #### Plot brGDGT ####
  p2 <- ggplot(B, aes(x=variable, y=value, fill = Sample.type)) + 
    coord_cartesian(ylim = c(0, Ymax)) + 
    xlab("brGDGTs") +
    geom_boxplot(outlier.shape = NA, position = position_dodge(.85), linewidth = Box.linewidth) +
    Add.points.size.lim + Leg.guide + My_title +
    geom_vline(xintercept = c(Xligne), lty="dotted")+
    annotate("text", x = Xtext, y = Ymax, label = My_lab, size = Annot.size, hjust = 0.25)+
    scale_fill_manual(name = "Sample type", values = Color.vec
                      , labels = c("0" = "Foo", "1" = "Bar"))+
    theme(legend.position = Leg.pos, panel.background = element_blank(),
          legend.key = element_blank(), legend.background = My_legbox,
          axis.title = element_text(size = 14), panel.border = My_border,
          axis.text = element_text(size = 12), panel.grid = element_blank(),
          plot.margin = unit(x = c(1, 2, 2, 0),units="mm"), # Bas / Droite / Haut / Gauche
          axis.line = element_line(colour = "black"), legend.text = element_text(size = Leg.size),
          legend.title = element_text(size = Leg.size*1.3),
          plot.background = element_blank())+
    A.lab + My_insert1 + My_insert2 + Box1 + Box2
  
  plot_build <- ggplot_build(p2)
  boxplot_stats <- plot_build$data[[1]]
  boxplot_stats <- round(boxplot_stats[c(13,3:5)], digits = 0)
  
  #### Plot isoGDGT ####
  if(Iso.GDGT == T){
    p1 <- ggplot(A.iso, aes(x=variable, y=value, fill=Sample.type)) + 
      coord_cartesian(ylim = c(0, Ymax.iso)) + 
      ylab("Fractional Abundance (%)")+
      xlab("isoGDGTs") +
      geom_boxplot(outlier.shape = NA, position = position_dodge(.85), linewidth = Box.linewidth) +
      scale_x_discrete(labels = c(parse(text = Model.lab), "Crenach\'"))+
      Add.points.size.lim.iso +
      scale_fill_manual(name = "Sample type", values = Color.vec
                        , labels = c("0" = "Foo", "1" = "Bar"))+
      theme(axis.text.x = element_text(angle = 25, hjust = 1),
            legend.position = Leg.pos.iso, panel.grid = element_blank(),
            panel.background = element_blank(), legend.background = element_blank(),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            axis.line = element_line(colour = "black"),
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 14),
            plot.margin = unit(x = c(0, 0, 0, 0),units="mm"), 
            plot.background = element_blank()
      )
    
    # p3 <- plot_grid(p1, p2, ncol = 2, axis="tb", align = "h", rel_widths = c(3/10, 7/10), labels = c("A","B"), label_size = 18, label_x = -0.01) # 2 graphs
    p3 <- p1 + p2 + plot_layout(ncol = 2, widths = c(3/10, 7/10)) + plot_annotation(tag_levels = 'A') &
      theme(plot.tag = element_text(size = 18, face = "bold", hjust = 0.5, vjust = -3.5), strip.background = element_blank(),
            plot.background = element_blank(), plot.margin = unit(x = c(0, 0, 0, 0),units="mm"))
    
  }
  
  #### Save html ####
  if(Show.Plotly == T){
    library(plotly)
    library(htmlwidgets)
    Save.plot.html <- gsub("pdf", "html", Save.path)
    Keep.name <- gsub(".*\\/", "", Save.plot.html)
    Path.root <- paste(gsub(Keep.name, "", Save.plot.html), "HTML_files/", sep = "")
    if(file.exists(Path.root) == F){dir.create(Path.root)}
    Save.plot.html <- paste(Path.root, Keep.name, sep = "")
    p1_ly <- ggplotly(p2)
    p1_ly <- p1_ly %>% layout(boxmode = "group", boxpoints = F, legend = list(font = list(size = 12)))
    options(warn = - 1) 
    saveWidget(p1_ly, file = Save.plot.html)
  }
  #### Save plots ####
  if(is.null(Save.path) == F){
    if(is.null(W) == F & is.null(H) == F){
      ggsave(Save.path, width = W*0.026458333, height = H*0.026458333, units = "cm")}
    else{ggsave(Save.path)}}
  else(
    if(Iso.GDGT == T){return(p3)}
    else{return(p2)})
  # if(is.null(p2) == F){View(ggplot_build(p1)$data[[1]])}
  
  if(Return.plot == T){
    if(Iso.GDGT == T){return(p3)}
    else{return(p2)}
  }
}

FT.core <- function(MCore, MAge, Model.WAPLS, Model.MAT, Fit.val, Model.RF, Model.BRT,
                    Save.tab, Save.plot, Save.RDS, H, W, Only.fit, LakeName, Select.clim,
                    Ecartype.curve, Model.param.show, Displot, Verbose, GDGT = T,
                    Zone.Clim.span, Zone.Temp, Save.path){
  #### Init param ####
  if(missing(Verbose)){Verbose = T}
  if(missing(Save.tab)){Save.tab = T}
  if(missing(Zone.Clim.span)){Zone.Clim.span = NULL}
  if(missing(Zone.Temp)){Zone.Temp = rep("U", length(Zone.Clim.span)/2)}
  if(missing(Select.clim)){Select.clim = NULL}
  if(missing(Model.WAPLS)){Model.WAPLS = NULL}
  if(missing(Model.MAT)){Model.MAT = NULL}
  if(missing(Model.RF)){Model.RF = NULL}
  if(missing(Model.BRT)){Model.BRT = NULL}
  if(missing(Save.path)){Save.tab = F; Save.path = NULL}
  if(missing(Displot)){Displot = F}
  if(missing(Zone.Clim.span)){Zone.OK = F}
  if(missing(Model.param.show)){Model.param.show = F}
  if(missing(Save.plot)){Save.plot = NULL}
  if(missing(Save.RDS)){Save.RDS = T}
  if(missing(W)){W = NULL}
  if(missing(H)){H = NULL}
  if(missing(Only.fit)){Only.fit = F}
  if(missing(LakeName)){LakeName = "Lake"}
  if(missing(Fit.val)){Fit.val = 0}
  if(missing(Ecartype.curve)){Ecartype.curve = c(F,F,F,F)}
  if(missing(MAge)){MAge = paste(LakeName, seq(1:nrow(MCore)), sep = "_")}
  
  #### Select le model type ####
  Keep.WAPLS <- deparse(substitute(Model.WAPLS))
  Keep.MAT <- deparse(substitute(Model.MAT))
  Keep.RF <- deparse(substitute(Model.RF))
  Keep.BRT <- deparse(substitute(Model.BRT))
  
  if(is.null(Model.WAPLS)== F){Model.type <- Model.WAPLS}
  if(is.null(Model.MAT)== F){Model.type <- Model.MAT}
  if(is.null(Model.RF)== F){Model.type <- Model.RF}
  if(is.null(Model.BRT)== F){Model.type <- Model.BRT}
  
  #### Save plots ####
  if(is.null(Save.plot) == F & Displot == T){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = FALSE)
    if(is.null(W) == F & is.null(H) == F){
      pdf(file = Save.plot, width = W*0.01041666666667, height = H*0.01041666666667)}
    else{pdf(file = Save.plot)}}
  
  #### Select param clim ####
  if(is.null(Select.clim) == F){
    Select.clim <- c(Select.clim, "Best.Param")
    Model.type <- Model.type[names(Model.type) %in% Select.clim]
    Model.type$Best.Param <- Model.type$Best.Param[row.names(Model.type$Best.Param) %in% Select.clim,]
    
    if(is.null(Model.BRT) == F){
      Model.BRT <- Model.BRT[names(Model.BRT) %in% Select.clim]
      Model.BRT$Best.Param <- Model.BRT$Best.Param[row.names(Model.BRT$Best.Param) %in% Select.clim,]
    }
    
    if(is.null(Model.MAT) == F){
      Model.MAT <- Model.MAT[names(Model.MAT) %in% Select.clim]
      Model.MAT$Best.Param <- Model.MAT$Best.Param[row.names(Model.MAT$Best.Param) %in% Select.clim,]
    }
    
    if(is.null(Model.WAPLS) == F){
      Model.WAPLS <- Model.WAPLS[names(Model.WAPLS) %in% Select.clim]
      Model.WAPLS$Best.Param <- Model.WAPLS$Best.Param[row.names(Model.WAPLS$Best.Param) %in% Select.clim,]
    }
    
    if(is.null(Model.RF) == F){
      Model.RF <- Model.RF[names(Model.RF) %in% Select.clim]
      Model.RF$Best.Param <- Model.RF$Best.Param[row.names(Model.RF$Best.Param) %in% Select.clim,]
    }
  }
  
  #### Graphical settings ####
  Zone.Temp <- rep(Zone.Temp, each=2)
  
  if(is.null(Zone.Clim.span) == T | length(Zone.Clim.span) != length(Zone.Temp)){
    Zone.OK = F; print("Something is wrong with the climate zones. Please check.")}
  else{Zone.OK = T}
  
  Tailleplot <-length(Model.type)-1
  if (Tailleplot <= 3){par(mfrow = c(1,Tailleplot))}
  if (Tailleplot == 4){par(mfrow = c(2,2))}
  if (Tailleplot >= 5 & Tailleplot <= 6){par(mfrow = c(2,3))}
  if (Tailleplot >= 7 & Tailleplot <= 9){par(mfrow = c(3,3))}
  Mmodel <- data.frame(Age = MAge)
  MModel.MAT <- data.frame(Age = MAge)
  MModel.RF <- data.frame(Age = MAge)
  MModel.BRT <- data.frame(Age = MAge)
  M.errors.WAPLS <- data.frame(Age = MAge)
  M.errors.MAT <- data.frame(Age = MAge) 
  M.errors.RF <- data.frame(Age = MAge)
  M.errors.BRT <- data.frame(Age = MAge) 
  Full.MAT.RDS <- list()
  Full.WAPLS.RDS <- list()
  Full.RF.RDS <- list()
  Full.BRT.RDS <- list()
  
  #### Loop on the climat param ####
  if(Verbose == F){
    library(lubridate)
    pb = txtProgressBar(min = 1, max = (length(Model.type)-1), width = 40, initial = 0,  style = 3) 
    init <- numeric((length(Model.type)-1))
    end <- numeric((length(Model.type)-1))
  }
  
  print(paste("Prediction for ", LakeName, " with the following models :", Keep.WAPLS, ", ", Keep.MAT, ", ", Keep.BRT, ", ",  Keep.RF, ".", sep = ""))
  for (i in 1:(length(Model.type)-1)){
    if(Verbose == F){init[i] <- Sys.time()}
    LabParamlim <- as.character(gsub("\\p{P}","", deparse(names(Model.type)[i]), perl = TRUE ))
    #### WAPLS ####
    if(is.null(Model.WAPLS) == F){
      NComp.WAPLS = Model.WAPLS$Best.Param[[3]][i]
      if(NComp.WAPLS < 2){NComp.WAPLS <- 2}
      if(Verbose == T){print(paste(round(i/(length(Model.WAPLS)), digits = 2)*100, "% done. The ", LabParamlim, " is modelling with the WAPLS and ", NComp.WAPLS, " parameters.", sep = ""))}
      Cor.WAPLS = predict(Model.WAPLS[[i]], MCore, npls = NComp.WAPLS, sse = T, nboot = 1000, verbose = F)
      Mmodel[i+1] <- cbind(Cor.WAPLS$fit[,NComp.WAPLS])
      colnames(Mmodel)[i+1] <- LabParamlim
      M.errors.WAPLS[i+1] <- cbind(Cor.WAPLS$SEP.boot[,NComp.WAPLS])
      colnames(M.errors.WAPLS)[i+1] <- LabParamlim
      Full.WAPLS.RDS[[i]] <- Cor.WAPLS
      names(Full.WAPLS.RDS)[[i]] <- LabParamlim}
    
    #### MAT ####
    if(is.null(Model.MAT) == F){
      NComp.MAT = Model.MAT$Best.Param[[3]][i]
      if(NComp.MAT < 4){NComp.MAT <- 4}
      if(Verbose == T){print(paste("The ", LabParamlim, " is modelling with the MAT and ", NComp.MAT, " analogues.", sep = ""))}
      Cor.MAT=predict(Model.MAT[[i]], MCore, k = NComp.MAT, sse = T, nboot = 1000, verbose = F)
      MModel.MAT[i+1] <- cbind(Cor.MAT$fit[,2])     # 2 = value-wm (weighted mean), 1 = normal-value
      M.errors.MAT[i+1] <- cbind(Cor.MAT$SEP.boot[,2]) # 2 = value-wm (weighted mean)
      colnames(MModel.MAT)[i+1] <- LabParamlim
      colnames(M.errors.MAT)[i+1] <- LabParamlim
      Full.MAT.RDS[[i]] <- Cor.MAT
      names(Full.MAT.RDS)[[i]] <- LabParamlim}
    
    #### Random forest ####
    if(is.null(Model.RF) == F){
      if(Verbose == T){print(paste("The ", LabParamlim, " is modelling with the RF.", sep = ""))}
      Cor.RF = predict(Model.RF[[i]], MCore, see = T)
      MSE <- Model.RF[[i]]$mse             # Pb dans le mse. Manquant pou Mongolie ?
      RMSE.RF <- sqrt(MSE[length(MSE)])
      MModel.RF <- cbind(MModel.RF, Cor.RF)
      colnames(MModel.RF)[i+1] <- LabParamlim
      Full.RF.RDS[[i]] <- Cor.RF
      names(Full.RF.RDS)[[i]] <- LabParamlim}
    
    #### Boosted Regression Tree (BRT) ####
    if(is.null(Model.BRT) == F){
      if(Verbose == T){print(paste("The ", LabParamlim, " is modelling with the BRT.", sep = ""))}
      MCore.i <- MCore[, names(MCore) %in% colnames(Model.BRT[[i]]$data$x.order)] # ATTENTION il manque les taxon < 0.1 pour les BRT
      Cor.BRT <- gbm::predict.gbm(Model.BRT[[i]], MCore.i, n.trees = Model.BRT[[i]]$gbm.call$best.trees, type="response")
      RMSE.BRT <- Model.BRT$Best.Param[i,2]
      
      MModel.BRT<- cbind(MModel.BRT, Cor.BRT)
      colnames(MModel.BRT)[i+1] <- LabParamlim
      Full.BRT.RDS[[i]] <- Cor.BRT
      names(Full.BRT.RDS)[[i]] <- LabParamlim
    }
    
    #### Plot résultats graphiques ####
    if(Displot == T){
      #### Val Min / Max ####
      ymin = 0
      ymax = 0
      if(Ecartype.curve[1] == F){
        
        if(is.null(Model.WAPLS) == F){
          ymin = min(Cor.WAPLS$fit[,NComp.WAPLS], na.rm = T)
          ymax = max(Cor.WAPLS$fit[,NComp.WAPLS], na.rm = T)}
        
        if(is.null(Model.MAT) == F & Ecartype.curve[2] == F){
          ymin = min(ymin, Cor.MAT$fit[,2], na.rm = T)
          ymax = max(ymax, Cor.MAT$fit[,2], na.rm = T)}
        
        if(is.null(Model.RF) == F & Ecartype.curve[3] == F){
          ymin = min(ymin, min(Cor.RF, na.rm = T), na.rm = T)
          ymax = max(ymax, max(Cor.RF, na.rm = T), na.rm = T)}
        
        if(is.null(Model.BRT) == F & Ecartype.curve[4] == F){
          ymin = min(ymin, min(Cor.BRT, na.rm = T), na.rm = T)
          ymax = max(ymax, max(Cor.BRT, na.rm = T), na.rm = T)}
        
        fullY <- abs(ymax) - abs(ymin)
        ymax = ymax + 0.05*fullY
      }
      else{
        if(is.null(Model.WAPLS) == F){
          ymin = min(Cor.WAPLS$fit[,NComp.WAPLS] - Cor.WAPLS$SEP.boot[,1], na.rm = T)
          ymax = max(Cor.WAPLS$fit[,NComp.WAPLS] + Cor.WAPLS$SEP.boot[,1], na.rm = T)}
        
        if(is.null(Model.MAT) == F){
          ymin = min(ymin, na.omit(Cor.MAT$fit[,2] - Cor.MAT$SEP.boot[,1]), na.rm = T)
          ymax = max(ymax, na.omit(Cor.MAT$fit[,2] + Cor.MAT$SEP.boot[,1]), na.rm = T)}
        
        if(is.null(Model.RF) == F){
          ymin = min(ymin, min(Cor.RF - RMSE.RF, na.rm = T), na.rm = T)
          ymax = max(ymax, max(Cor.RF + RMSE.RF, na.rm = T), na.rm = T)}
        
        if(is.null(Model.BRT) == F){
          ymin = min(ymin, c(Cor.BRT - RMSE.BRT), na.rm = T)
          ymax = max(ymax, max(Cor.BRT + RMSE.BRT, na.rm = T), na.rm = T)}
        
        fullY <- abs(ymax) - abs(ymin)
        ymax = ymax + 0.05*fullY
      }
      
      #### Plot Mean value ####
      if(is.null(Model.WAPLS) == F & Only.fit == F){Y <- Cor.WAPLS$fit[,NComp.WAPLS]}
      else{Y <- rep(NA, length(MAge))}
      
      plot(MAge, Y, 
           xlim = c(min(MAge),max(MAge)), 
           ylim = c(ymin, ymax), 
           type = "l", 
           ylab = LabParamlim, 
           xlab = "Time (yr cal BP)", 
           col = 1, 
           las = 0, lwd=1, bty="n")
      
      #### Legend ####
      if(Model.param.show == T){
        fullT <- abs(min(MAge)) + abs(max(MAge))
        x1 = 0.3*fullT
        x2 = 0.4*fullT
        x4 = 0.55*fullT
        x5 = 0.65*fullT
        
        ymax.lab1 <- ymax 
        ymax.lab2 <- ymax - 0.05*ymax
        
        #### Legend WAPLS ####
        if(is.null(Model.WAPLS) == F){
          mylabel1 = Keep.WAPLS
          mylabel2 = bquote(npls == .(Model.WAPLS[[1]][["npls"]]) ~ "," ~
                              k == .(Model.WAPLS[[length(Model.WAPLS)]][i,3]) ~ "," ~
                              italic(R)^2 == .(format(Model.WAPLS[[length(Model.WAPLS)]][i,4], digits = 2)) ~ "," ~
                              RMSE == .(format(Model.WAPLS[[length(Model.WAPLS)]][i,5], digits = 2)))
          
          text(x = x1, y = ymax.lab1, labels = mylabel1, font = 2, cex = 0.8, pos = 1)
          text(x = x2, y = ymax.lab1, labels = mylabel2, cex = 0.8, pos = 1)
        }
        
        #### Legend MAT ####
        if(is.null(Model.MAT) == F){
          mylabel1b = Keep.MAT
          mylabel3b = bquote(k == .(Model.MAT[[length(Model.MAT)]][i,3]) ~ "," ~
                               italic(R)^2 == .(format(Model.MAT[[length(Model.MAT)]][i,4], digits = 2)) ~ "," ~
                               RMSE == .(format(Model.MAT[[length(Model.MAT)]][i,5], digits = 2)))
          
          
          text(x = x1, y = ymax.lab2, labels = mylabel1b, col = "royalblue", font = 2, cex = 0.8)
          text(x = x2, y = ymax.lab2, labels = mylabel3b, cex = 0.8)
        }
        
        #### Legend RF ####
        if(is.null(Model.RF) == F){
          mylabel1b = Keep.RF
          mylabel2b = bquote(#Nb.tree == .(Model.RF[[i]][["call"]][["ntree"]]) ~ "," ~
            italic(R)^2 == .(format(Model.RF[[length(Model.RF)]][i,1], digits = 2)) ~ "," ~
              RMSE == .(format(Model.RF[[length(Model.RF)]][i,2], digits = 2)))
          
          
          text(x = x4, y = ymax.lab2, labels = mylabel1b, col = "darkorange", font = 2, cex = 0.8, pos = 1)
          text(x = x5, y = ymax.lab2, labels = mylabel2b, cex = 0.8, pos = 1)
        }
        
        
        #### Legend BRT ####
        if(is.null(Model.BRT) == F){
          mylabel1b = Keep.BRT
          mylabel2b = bquote(Nb.tree == .(Model.BRT[[i]][["call"]][["ntree"]]) ~ "," ~
                               italic(R)^2 == .(format(Model.BRT[[length(Model.BRT)]][i,1], digits = 2)) ~ "," ~
                               RMSE == .(format(Model.BRT[[length(Model.BRT)]][i,2], digits = 2)))
          
          
          text(x = x4, y = ymax.lab1, labels = mylabel1b, col = "darkgreen", font = 2, cex = 0.8, pos = 1)
          text(x = x5, y = ymax.lab1, labels = mylabel2b, cex = 0.8, pos = 1)
        }
      }
      
      #### Plot Model MAT add ####
      if (is.null(Model.MAT) == F & Only.fit == F){
        lines(MAge, Cor.MAT$fit[,2], col = "royalblue", lwd=1, las=0)}
      
      #### Plot Model RF add ####
      if(is.null(Model.RF) == F & Only.fit == F){
        lines(MModel.RF[[1]], Cor.RF, col = "darkorange", lwd=1, las=0)}
      
      #### Plot Model BRT add ####
      if(is.null(Model.BRT) == F & Only.fit == F){
        lines(MModel.BRT[[1]], Cor.BRT, col = "darkgreen", lwd=1, las=0)}
      
      #### Plot fitting ####
      if(Fit.val > 0){
        if (is.null(Model.WAPLS) == F){
          Curve.fit = lowess(Cor.WAPLS$fit[,NComp.WAPLS], f = Fit.val)
          lines(MAge, Curve.fit$y, col=1, lwd=2)}
        if (is.null(Model.MAT) == F){
          Curve.fit.MAT = lowess(Cor.MAT$fit[,2], f = Fit.val)
          lines(MAge, Curve.fit.MAT$y, col= "royalblue", lwd=2)}
        if (is.null(Model.RF) == F){
          Curve.fit.MAT = lowess(Cor.RF, f = Fit.val)
          lines(MAge, Curve.fit.MAT$y, col = "darkorange", lwd=2)}
        if (is.null(Model.BRT) == F){
          Curve.fit.MAT = lowess(Cor.BRT, f = Fit.val)
          lines(MAge, Curve.fit.MAT$y, col = "darkgreen", lwd=2)}
      }
      
      #### Interval WAPLS ####
      if (Ecartype.curve[1] == T & is.null(Model.WAPLS) == F){
        lines(MAge, Cor.WAPLS$fit[,NComp.WAPLS] + Cor.WAPLS$SEP.boot[,1], lwd=.4, lty = "dashed")
        lines(MAge, Cor.WAPLS$fit[,NComp.WAPLS] - Cor.WAPLS$SEP.boot[,1], lwd=.4, lty = "dashed")
      }
      
      #### Interval MAT ####
      if (Ecartype.curve[2] == T & is.null(Model.MAT) == F){
        lines(MAge, Cor.MAT$fit[,2] + Cor.MAT$SEP.boot[,1], lwd=.4, col = "royalblue", lty = "dashed")
        lines(MAge, Cor.MAT$fit[,2] - Cor.MAT$SEP.boot[,1], lwd=.4, col = "royalblue", lty = "dashed")
      }
      
      #### Interval RF ####
      if (Ecartype.curve[3] == T & is.null(Model.RF) == F){
        lines(MAge, Cor.RF + RMSE.RF, lwd=.4, col = "darkorange", lty = "dashed")
        lines(MAge, Cor.RF - RMSE.RF, lwd=.4, col = "darkorange", lty = "dashed")
      }
      
      #### Interval BRT ####
      if (Ecartype.curve[4] == T & is.null(Model.BRT) == F){
        lines(MAge, Cor.BRT + RMSE.BRT, lwd=.4, col = "darkgreen", lty = "dashed")
        lines(MAge, Cor.BRT - RMSE.BRT, lwd=.4, col = "darkgreen", lty = "dashed")
      }
      
      #### Plot climate zones ####
      if (Zone.OK == T){
        for(j in 1:length(Zone.Clim.span)){
          if(as.logical(j%%2) == T){        # seulement les j impairs
            a <- which(MAge>Zone.Clim.span[j])
            b <- which(MAge>Zone.Clim.span[j+1])
            #Ymax <- max(Cor.WAPLS$fit[a[1]:b[1],NComp] + Cor.WAPLS$SEP.boot[a[1]:b[1],NComp])
            #Ymin <- min(Cor.WAPLS$fit[a[1]:b[1],NComp] - Cor.WAPLS$SEP.boot[a[1]:b[1],NComp])
            chaud = rgb(1, 0, 0, 0.08)
            froid = rgb(0, 0, 1, 0.08)
            unknow = rgb(0.3, 0.3, 0.3, 0.08)
            if(Zone.Temp[j]=="C"){colTemp <- froid}
            if(Zone.Temp[j]=="W"){colTemp <- chaud}
            if(Zone.Temp[j]=="U"){colTemp <- unknow}
            polygon(c(Zone.Clim.span[j],Zone.Clim.span[j+1], Zone.Clim.span[j+1],Zone.Clim.span[j]), c(ymin, ymin, ymax, ymax), border = NA, col = colTemp)
          }}}
    }
    if(Verbose == F){
      end[i] <- Sys.time()
      setTxtProgressBar(pb, i)
      time <- round(seconds_to_period(sum(end - init)), 0)
      est <- (length(Model.type)-1) * (mean(end[end != 0] - init[init != 0])) - time
      remainining <- round(seconds_to_period(est), 0)
      cat(paste(" // Execution time:", time,
                " // Estimated time remaining:", remainining), "")}
    
  }
  
  if(Verbose == F){close(pb);library(beepr)} 
  
  #### Save / export data ####
  row.names(Mmodel) <- row.names(MCore)
  par(mfrow = c(1,1))
  if(Save.tab == T){
    Path.to.create.csv <- gsub("(.*/).*\\.csv.*","\\1", Save.path)
    dir.create(file.path(Path.to.create.csv), showWarnings = FALSE)
    
    if(is.null(Model.WAPLS) == F){
      #### Save WAPLS ####
      Save.Model.WAPLS <- gsub("\\.", "_", Keep.WAPLS)
      add.to.path <- paste("_", Save.Model.WAPLS, ".csv", sep = "")
      add.to.path.error <- paste("_SEP_", Save.Model.WAPLS, ".csv", sep = "")
      Save.path1 <- gsub("\\.csv", add.to.path, Save.path)
      Save.path.error <- gsub("\\.csv", add.to.path.error, Save.path)
      write.table(Mmodel, file = Save.path1, row.names=T, col.names=NA, sep=",", dec = ".")
      write.table(M.errors.WAPLS, file = Save.path.error, row.names=T, col.names=NA, sep=",", dec = ".")
    }}
  
  if(is.null(Model.MAT) == F){
    #### Save MAT ####
    if(Save.tab == T){
      Save.Model.WAPLS <- gsub("\\.", "_", Keep.MAT)
      add.to.path <- paste("_", Save.Model.WAPLS, ".csv", sep = "")
      add.to.path.error <- paste("_SEP_", Save.Model.WAPLS, ".csv", sep = "")
      Save.path2 <- gsub("\\.csv", add.to.path, Save.path)
      Save.path2.error <- gsub("\\.csv", add.to.path.error, Save.path)
      write.table(MModel.MAT, file = Save.path2, row.names=T, col.names=NA, sep=",", dec = ".")
      write.table(M.errors.MAT, file = Save.path2.error, row.names=T, col.names=NA, sep=",", dec = ".")
    }
    
    #### Add MAT to function return ####
    Mmodel = list(Mmodel, MModel.MAT, M.errors.WAPLS, M.errors.MAT)
    Save.DB.name <- gsub(".*\\.","", Keep.WAPLS)
    labtot <- c(Keep.WAPLS, Keep.MAT, paste("SEP.WAPLS", Save.DB.name, sep = "."), paste("SEP.MAT", Save.DB.name, sep = "."))
    names(Mmodel) <- labtot}
  
  if(is.null(Model.RF) == F){
    #### Save RF ####
    if(Save.tab == T){
      Save.Model.WAPLS <- gsub("\\.", "_", Keep.RF)
      add.to.path <- paste("_", Save.Model.WAPLS, ".csv", sep = "")
      Save.path2 <- gsub("\\.csv", add.to.path, Save.path)
      write.table(MModel.RF, file = Save.path2, row.names=T, col.names=NA, sep=",", dec = ".")
    }
    
    #### Add RF to function return ####
    Mmodel <- append(Mmodel, list(MModel.RF))
    names(Mmodel)[length(Mmodel)] <- Keep.RF
  }
  
  if(is.null(Model.BRT) == F){
    #### Save BRT ####
    if(Save.tab == T){
      Save.Model.BRT <- gsub("\\.", "_", Keep.BRT)
      add.to.path <- paste("_", Save.Model.BRT, ".csv", sep = "")
      Save.path2 <- gsub("\\.csv", add.to.path, Save.path)
      write.table(MModel.BRT, file = Save.path2, row.names=T, col.names = NA, sep=",", dec = ".")
    }
    
    #### Add BRT to function return ####
    Mmodel <- append(Mmodel, list(MModel.BRT))
    names(Mmodel)[length(Mmodel)] <- Keep.BRT
  }
  
  Total.model <- list(WAPLS = Full.WAPLS.RDS, MAT = Full.MAT.RDS, BRT = Full.BRT.RDS, RF = Full.RF.RDS)
  
  #### Save format RDS ####
  if(Save.RDS == T & is.null(Save.path) == F){
    if(GDGT == F & exists("Save.DB.name") == F){Save.DB.name <- gsub(".*\\.","", deparse(substitute(Model.WAPLS)))}
    if(GDGT == T){
      if(is.null(Model.BRT) == F){Save.DB.name <- gsub(".*\\.","", Keep.BRT)}
      if(is.null(Model.WAPLS) == F){Save.DB.name <- gsub(".*\\.","", Keep.WAPLS)}
      if(is.null(Model.MAT) == F){Save.DB.name <- gsub(".*\\.","", Keep.MAT)}
      if(is.null(Model.RF) == F){Save.DB.name <- gsub(".*\\.","", Keep.RF)}
    }
    
    Save.path.RDS <- paste(gsub("\\.csv", paste("_", Save.DB.name, sep = ""), Save.path), ".Rds", sep = "")
    Save.path.RDS.full <- paste(gsub("\\.csv", paste("_", Save.DB.name, "_full", sep = ""), Save.path), ".Rds", sep = "")
    Path.to.create2 <- gsub("(.*/).*\\.Rds.*","\\1", Save.path.RDS)
    dir.create(file.path(Path.to.create2), showWarnings = F)
    
    saveRDS(Total.model, Save.path.RDS.full)
    saveRDS(Mmodel, Save.path.RDS)}
  
  #### End ####
  if(is.null(Save.plot) == F){dev.off()}
  #return(Mmodel)
  return(Total.model)
}

Combine.ML.cluster <- function(Cluster.prediction, List.models, Model.lab, GDGT.paleo, Plot.y = "Age", Param.clim = "MAAT", 
                               Highlight.combined = F, Save.path = NULL,
                               Compare.curve = NULL, Cluster.prob = "K-warm/arid", Time.lim = NULL, Surf.val = NULL, 
                               Core.name = NULL, Plot.y.lab = Plot.y, Show.proba = T, Facet = F, Only.best = T, H = 900, W = 500, Save.plot = NULL){
  #### Cluster predictions ####
  Br.GDGT.paleo <- GDGT.paleo[grepl("f.I", names(GDGT.paleo)) & !grepl("_7Me", names(GDGT.paleo))]
  Br.GDGT.paleo <- data.frame(t(Br.GDGT.paleo))
  Br.GDGT.paleo <- apply(Br.GDGT.paleo, 2, MESS::round_percent)
  Br.GDGT.paleo <- data.frame(t(Br.GDGT.paleo/100))
  RF_class <- predict(Cluster.prediction, Br.GDGT.paleo)
  RF_class_prob <- predict(Cluster.prediction, Br.GDGT.paleo, type = "prob")
  Br.GDGT.paleo$Pred.cluster <- RF_class
  Br.GDGT.paleo$Plot.y <- GDGT.paleo[[Plot.y]]
  
  #### Matrix full ####
  List.models <- Map(function(df, Model){df$Model <- Model; return(df)}, List.models, Model.lab)
  M <- do.call(rbind, List.models)
  All.param <- setdiff(names(M), c(Plot.y, "Model"))
  M$Pred.cluster <- Br.GDGT.paleo$Pred.cluster[match(M[[Plot.y]], Br.GDGT.paleo$Plot.y)]
  M$Model <- factor(M$Model, ordered = T, levels = unique(M$Model))
  
  #### Matrix full (weighted) ####
  List.models[[1]] <- List.models[[1]][c(Plot.y, Param.clim)]
  Mw <- cbind(List.models[[1]], MAAT.karid = List.models[[2]][[Param.clim]], MAAT.kwet = List.models[[3]][[Param.clim]], RF_class_prob)
  names(Mw)[c(1:2)] <- c("Plot.y", "Param.clim")
  Mw$Combine.weighted <- Mw$MAAT.kwet*Mw$`K-cold/wet` + Mw$MAAT.karid*Mw$`K-warm/arid`
  Mw$Pred.cluster <- ifelse(Mw$`K-warm/arid` >= 0.5, "K-warm/arid", "K-cold/wet")
  Mw$Model <- "Combine-Weighted"
  Mw$Param.clim <- Mw$Combine.weighted
  Keep <- c(Plot.y, Param.clim, "Model", "Pred.cluster")
  M <- M[Keep]
  names(M)[c(1:2)] <- c("Plot.y", "Param.clim")
  M <- full_join(M, Mw[c(1,2,9,8)], by = join_by(Plot.y, Param.clim, Model, Pred.cluster))
  
  if(Only.best == T){
    M <- M[M$Model  %in% c("ACADB", "Combine-Weighted"),]  
  }
  M$Model <- factor(M$Model, ordered = T, levels = unique(M$Model))
  
  RF_class_prob <- data.frame(RF_class_prob)
  RF_class_prob$Plot.y <- GDGT.paleo[[Plot.y]]
  
  #### Graphical settings ####
  My_colors <- c("Combined" = "grey20", "Combine-Weighted" = "darkred", 
                 "ACADB" = "bisque3", "K-cold/wet" = "royalblue", "K-warm/arid" = "darkorange",
                 "MAAT_soil_Naaf" = "darkgreen", "MAAT_LSun" = "darkblue",
                 "MAAT_mr_DJ" = "darkolivegreen3", "MAAT_DJ_5Me" = "darkolivegreen3", "MAAT_NMSDB_mr5" = "aquamarine2")
  
  My_labs <- c("Combined" = "Combined", "Combine-Weighted" = "BRT (combined)", 
               "MAAT_soil_Naaf" = "MBT'5Me (Naafs et al. 2017)", "MAAT_LSun" = "MBT/CBT lake (Sun et al. 2010)",
               "ACADB" = "BRT (ACADB)", "K-cold/wet" = "BRT (K-cold/wet)", "K-warm/arid" = "BRT (K-warm/arid)",
               "MAAT_mr_DJ" = "MR (De Jonge)", "MAAT_DJ_5Me" = "MBT'5Me (De Jonge)", "MAAT_NMSDB_mr5" = "MR mong.")
  
  if(Facet == T){My_facet <- tidypaleo::facet_geochem_grid(vars(Model))}
  else{My_facet <- NULL; H <- H/2}
  
  if(is.null(Plot.y.lab) == F){
    if(is.null(Compare.curve) == F){
      Age.lab.2 <- Plot.y.lab; Age.lab.1 <- NULL}
    else{Age.lab.1 <- Plot.y.lab}
    
  }
  else{
    if(is.null(Compare.curve) == F){
      Age.lab.2 <- Plot.y; Age.lab.1 <- NULL}
    else{Age.lab.1 <- Plot.y}
  }
  
  if(Param.clim == "MAAT"){
    Clim.lab.1 <- "MAAT\nBRT (°C)"
    Clim.lab.2 <- "MAAT\nclassic (°C)"
  }
  else{Clim.lab.1 <- Param.clim; Clim.lab.2 <- Param.clim}
  
  if(is.null(Surf.val) == F){
    Surf.line <- geom_hline(yintercept = Surf.val, linetype = "dotdash", linewidth = 0.55, color = "black", alpha = 0.55)
  }
  else{Surf.line <- NULL}
  
  if(is.null(Core.name) == F){
    My_title <- ggtitle(Core.name)
  }
  else{My_title <- NULL}
  
  if(Highlight.combined == T){
    Double.line <- geom_line(data = Mw, aes(x = Plot.y, y = Combine.weighted, group = Model, color = Model), linewidth = .8)
  }
  else{Double.line <- NULL}
  
  if(is.null(Time.lim) == F){Xlim <- xlim(Time.lim)}
  else{Xlim <- NULL}
  
  if(is.null(Plot.y.lab) == T){
    X.title <- element_blank()
  }
  else{X.title <- element_text()}
  
  #### Comparative curves ####
  if(is.null(Compare.curve) == F){
    Mc <- GDGT.paleo[c(Plot.y, Compare.curve)]
    Mc$'Combine-Weighted' <- Mw$Combine.weighted
    Mc <- melt(Mc, id = Plot.y)
    Ticks.1 <- element_blank(); Text.1 <- element_blank()
    
    Mc$Pred.cluster <- Br.GDGT.paleo$Pred.cluster[match(Mc[[Plot.y]], Br.GDGT.paleo$Plot.y)]
    names(Mc)[c(1:2)] <- c("Plot.y", "Model")
    Mc$Model <- factor(Mc$Model, ordered = T, levels = unique(Mc$Model))
    p.comp <- ggplot(Mc, aes(x = Plot.y, y = value, group = Model))+
      My_facet+ Surf.line+ Xlim+
      geom_point(aes(color = Pred.cluster))+
      geom_line(aes(color = Model))+ xlab(Age.lab.2)+ ylab(Clim.lab.2)+
      Double.line +
      scale_color_manual(values = My_colors, label = My_labs, name = NULL)
  }
  else{Ticks.1 <- element_line(); Text.1 <- element_text()}
  
  #### Plot (main) ####
  p.main <- ggplot(M, aes(x = Plot.y, y = Param.clim, group = Model))+
    My_facet+ Surf.line+ Xlim+
    geom_point(aes(color = Pred.cluster))+
    geom_line(aes(color = Model))+ xlab(Age.lab.1)+ ylab(Clim.lab.1)+
    Double.line+
    scale_color_manual(values = My_colors, label = My_labs, name = NULL)+
    theme(axis.text.x = Text.1, axis.ticks.x = Ticks.1)
  
  #### Plot (proba. cluster) ####
  if(Show.proba == T){
    if(Cluster.prob == "K-warm/arid"){
      RF_class_prob$Prob.to.show <- RF_class_prob$K.warm.arid
      Estim.col <- "darkorange"}
    if(Cluster.prob == "K-cold/wet"){
      RF_class_prob$Prob.to.show <- RF_class_prob$K.cold.wet
      Estim.col <- "royalblue"}
    
    RF_class_prob$Prob.to.show <- round(RF_class_prob$Prob.to.show*100, digits = 0)
    p.prob <- ggplot(RF_class_prob, aes(x = Plot.y, y = Prob.to.show))+
      geom_area(fill = Estim.col, position = "stack")+ My_title+
      Xlim+ scale_y_continuous(limits = c(0,100), breaks = c(0,50,100))+
      xlab(NULL)+ylab("Cluster\nprob.")+
      theme(
        axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.title.x = element_blank())
    
    p.main <- p.prob / p.main# + plot_layout(heights = c(1, nlevels(M$Model)))
  }
  
  #### Merge plots ####
  if(is.null(Compare.curve) == F){
    p.main <- p.main / p.comp + plot_layout(guides = "collect")
  }
  
  p.main <- p.main&
    theme(panel.background = element_rect(fill = NA, color = "black"), legend.key = element_rect(fill = NA, color = NA),
          panel.grid = element_blank(), axis.title.x = X.title, 
          plot.background = element_blank(),
          plot.margin = unit(c(0,0,0,0), 'pt'),
          legend.position = "bottom")&
    guides(color = guide_legend(nrow = 2))
  
  #### Generalized combined plot ####
  if(is.null(Save.path) == F){
    Mexp <- data.frame(List.models[[1]][Plot.y])
    for(i in 1:length(All.param)){
      Clim.i <- All.param[i]
      Mw <- cbind(List.models[[1]], M2 = List.models[[2]][[Clim.i]], M1 = List.models[[3]][[Clim.i]], RF_class_prob)
      
      Mexp[i+1] <- Mw$M1*Mw$K.cold.wet + Mw$M2*Mw$K.warm.arid
      names(Mexp)[i+1] <- Clim.i
    }
    Mexp$Pred.cluster <- ifelse(Mw$K.warm.arid >= 0.5, "K-warm/arid", "K-cold/wet")
    saveRDS(Mexp, Save.path)
  }
  #### Export plot ####
  if(is.null(Save.plot) == F){ggsave(filename = Save.plot, p.main, width = W*0.026458333, height = H*0.026458333, units = "cm")}
  return(p.main)
}

Plot.randomTF <- function(MPsurf, MPpaleo, Mclim, Plot.MAT = F, Plot.WAPLS = F, Plot.RF = F, Plot.BRT = T, Bold = T, Database = NULL, Lake = NULL,
                          Nb.simul = 99, H = 300, W = 1500, Save.path = NULL, return.plot = F, Save.plot = NULL) {
  #### Settings ####
  DB.result <- NULL
  if(is.null(Lake) == T & is.null(Database) == F){print(paste("randomFT() made for", Database))}
  if(is.null(Lake) == F & is.null(Database) == F){print(paste("randomFT() made for", Database, "on lake", Lake))}
  
  #### MAT settings ####
  if(Plot.MAT == T){
    rlghr <- randomTF(spp = sqrt(MPsurf), 
                      env = Mclim,
                      fos = sqrt(MPpaleo), 
                      n = Nb.simul, fun = MAT, col = 1, k = 10, lean = F)
    
    PP1 <- autoplot(rlghr) + theme_light() + ggtitle(paste("MAT (", Database, ")", sep = ""))
    
    P95 <- sort(rlghr$sim.ex)[floor(length(rlghr$sim.ex)*0.95)]
    DB.MAT <- data.frame(round(cbind(t(rlghr$EX), P95 = P95, Var.max = rlghr$MAX), 2))
    if(Bold == T){DB.MAT[1,which(rlghr$EX >= round(P95, 2))] <- paste("\\textbf{", DB.MAT[1,which(rlghr$EX >= round(P95, 2))], "}", sep = "")}
    DB.MAT[1,which(rlghr$sig <= 0.01)] <- paste(DB.MAT[1,which(rlghr$sig <= 0.01)], "*", sep = "") # Rehfeld et al., 2016
    DB.MAT[1,which(rlghr$sig <= 0.1)] <- paste(DB.MAT[1,which(rlghr$sig <= 0.1)], "*", sep = "") # Rehfeld et al., 2016
    DB.MAT <- cbind(Database = Database, Model = "MAT", DB.MAT)
  }
  else{PP1 <- NULL; DB.MAT <- NULL}
  
  #### WAPLS settings ####
  if(Plot.WAPLS == T){
    rlghr <- randomTF(spp = sqrt(MPsurf), 
                      env = Mclim,
                      fos = sqrt(MPpaleo), 
                      n = Nb.simul, fun = WAPLS, col = 1, npls = 5, lean = F)
    
    PP2 <- autoplot(rlghr) + theme_light() + ggtitle(paste("WAPLS (", Database, ")", sep = ""))
    
    P95 <- sort(rlghr$sim.ex)[floor(length(rlghr$sim.ex)*0.95)]
    DB.WAPLS <- data.frame(round(cbind(t(rlghr$EX), P95 = P95, Var.max = rlghr$MAX), 2))
    if(Bold == T){DB.WAPLS[1,which(rlghr$EX >= round(P95, 2))] <- paste("\\textbf{", DB.WAPLS[1,which(rlghr$EX >= round(P95, 2))], "}", sep = "")}
    DB.WAPLS[1,which(rlghr$sig <= 0.01)] <- paste(DB.WAPLS[1,which(rlghr$sig <= 0.01)], "*", sep = "")
    DB.WAPLS[1,which(rlghr$sig <= 0.1)] <- paste(DB.WAPLS[1,which(rlghr$sig <= 0.1)], "*", sep = "")
    DB.WAPLS <- cbind(Database = Database, Model = "WAPLS", DB.WAPLS)
  }
  else{PP2 <- NULL; DB.WAPLS <- NULL}
  
  #### RF settings ####
  if(Plot.RF == T){
    rlghr <- randomTF(spp = sqrt(MPsurf), 
                      env = Mclim,
                      fos = sqrt(MPpaleo), 
                      n = Nb.simul, fun = randomForest, col = 1, ntree = 100, mtry = 2, na.action = na.roughfix)
    
    PP3 <- autoplot(rlghr) + theme_light() + ggtitle(paste("RF (", Database, ")", sep = ""))
    P95 <- sort(rlghr$sim.ex)[floor(length(rlghr$sim.ex)*0.95)]
    DB.RF <- data.frame(round(cbind(t(rlghr$EX), P95 = P95, Var.max = rlghr$MAX), 2))
    if(Bold == T){DB.RF[1,which(rlghr$EX >= round(P95, 2))] <- paste("\\textbf{", DB.RF[1,which(rlghr$EX >= round(P95, 2))], "}", sep = "")}
    DB.RF[1,which(rlghr$sig <= 0.01)] <- paste(DB.RF[1,which(rlghr$sig <= 0.01)], "*", sep = "")
    DB.RF[1,which(rlghr$sig <= 0.1)] <- paste(DB.RF[1,which(rlghr$sig <= 0.1)], "*", sep = "")
    DB.RF <- cbind(Database = Database, Model = "RF", DB.RF)
  }
  else{PP3 <- NULL; DB.RF <- NULL}
  
  #### BRT settings ####
  if(Plot.BRT == T){
    # rlghr <- randomTF(spp = sqrt(MPsurf), env = Mclim, fos = sqrt(MPpaleo), 
    rlghr <- randomTF(spp = MPsurf, env = Mclim, fos = MPpaleo,
                      n = Nb.simul, fun = train, col = 1, method = "gbm", 
                      verbose = F)
    
    PP4 <- autoplot(rlghr) + theme_light() + ggtitle(paste("BRT (", Database, ")", sep = ""))
    
    P95 <- sort(rlghr$sim.ex)[floor(length(rlghr$sim.ex)*0.95)]
    DB.BRT <- data.frame(round(cbind(t(rlghr$EX), P95 = P95, Var.max = rlghr$MAX), 2))
    if(Bold == T){DB.BRT[1,which(rlghr$EX >= round(P95, 2))] <- paste("\\textbf{", DB.BRT[1,which(rlghr$EX >= round(P95, 2))], "}", sep = "")}
    DB.BRT[1,which(rlghr$sig <= 0.01)] <- paste(DB.BRT[1,which(rlghr$sig <= 0.01)], "*", sep = "")
    DB.BRT[1,which(rlghr$sig <= 0.1)] <- paste(DB.BRT[1,which(rlghr$sig <= 0.1)], "*", sep = "")
    DB.BRT <- cbind(Database = Database, Model = "BRT", DB.BRT)
  }
  else{PP4 <- NULL; DB.BRT <- NULL}
  
  #### Export results ####
  DB.result <- rbind(DB.MAT, DB.WAPLS, DB.RF, DB.BRT)
  PP <- PP1 + PP2 + PP3 + PP4
  
  if(is.null(Save.path) == F){
    Path.to.create <- gsub("(.*/).*\\.Rds.*","\\1", Save.path)
    dir.create(file.path(Path.to.create), showWarnings = F)
    saveRDS(DB.result, Save.path)
  }
  if(is.null(Save.plot) == F){
    Path.to.create <- gsub("(.*/).*\\.pdf.*","\\1", Save.plot)
    dir.create(file.path(Path.to.create), showWarnings = F)
    ggsave(file = Save.plot, PP, width = W*0.01041666666667, height = H*0.01041666666667)}
  
  if(return.plot == F){return(DB.result)}
  else{return(PP)}
}