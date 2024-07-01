#install.packages("C:/Users/simon/OneDrive/Desktop/Simone/R/APP/NanoStringNorm_1.2.1.tar.gz")
library(later)
library(shinyWidgets)
library(shiny)
library(shinyFiles)
library(tidyverse)
library(pheatmap)
library(DT)
library(forestmodel) ##
library(ggsci)
library(NanoStringNorm)
#library(NanoStringDiff) 
library(DESeq2) ##
library(tidyverse)
library(apeglm)
library(openxlsx)
library(randomcoloR)
library(survival)
library(ggplot2)
library(edgeR) ##
library(naniar)
library(plotly)
library(RColorBrewer)
library(pals)
library(knitr) ##
library(ggfortify)
library(dbplyr)

# Imposta la cartella di lavoro
#setwd("~/Dropbox/RCC TESI/")
#se qualocsa non va riguardare i save passo passo
ui <- fluidPage( 
  
  headerPanel("Nanostring QC"),
  helpText("Analisi QC di file RCC prodotti da Nanostring nCounter®"),
  helpText("Le boxplot colorate di rosso indicano che il campione non ha superato gli sntandard richiesti dal filtraggio"),
  mainPanel(p("1 • Seleziona una cartella contenente i file RCC che vuoi analizzare"
              , style = "font-family: 'times'; font-si16pt"),
            p("2 • Scarica il report html contenente i risultati dei QC"
              , style = "font-family: 'times'; font-si16pt"),
            p("3 • Scarica le matrici relative ai valori di espressione normalizzati con o senza outliers"
              , style = "font-family: 'times'; font-si16pt")),
  
  # shinyDirButton("Btn_GetFolder", "Choose a folder" ,
  #                title = "Please select a folder:", multiple = FALSE,
  #                buttonType = "default", class = NULL),
  
  #fluidRow
  sidebarPanel(width = 12,
               
               #),
               #sidebarPanel( width = 12, 
               column(3,          
                      selectInput("background", ("Background"),
                                  choices = list("none" = "none", "mean" = "mean", "mean.2sd" = "mean.2sd"), 
                                  selected = "none")
               ),
               
               column(3,
                      selectInput("samplecontent",("SampleContent"),
                                  choices = list("none"="none", "housekeeping.sum"="housekeeping.sum", "housekeeping.geo.mean"="housekeeping.geo.mean",
                                                 "total.sum"="total.sum", "low.cv.geo.mean"="low.cv.geo.mean", "top.mean"="top.mean", "top.geo.mean."="top.geo.mean."),
                                  selected = "none")
               ),
               
               column(3,
                      selectInput("codecount",("CodeCount"),
                                  choices = list("none" = "none", "sum"="sum","geo.mean"="geo.mean"), 
                                  selected = "none")
               ),
               
               column(3,  
                      selectInput("othernorm",("OtherNorm"),
                                  choices = list("none" = "none", "quantile"="quantile","zscore"="zscore"), 
                                  selected = "none")
               ),
               fluidRow(    
                 shinyDirButton("Btn_GetFolder", "Scegli una Cartella" ,
                                title = "Please select a folder:", multiple = FALSE,
                                buttonType = "default", class = NULL, align = "right")
               ),
               fluidRow(
                 downloadLink("downloadReport", "Download Report"),
               ),
               fluidRow(
                 downloadLink("tablepre", "Download matrice pre-filtraggio"),
               ),
               fluidRow(
                 downloadLink("tablepost", "Download matrice post-filtraggio")
               )
  ),
  
  # downloadLink("downloadReport", "Download Report"),
  # downloadLink("tablepre", "Download matrice pre filtaggio"),
  # downloadLink("tablepost", "Download matrice post filtaggio"),  
  fluidRow( align = "center",
            textOutput("outputT"),
            width=12, solidheader=T
  ),
  
  fluidRow(
    tabsetPanel(id = "main_tabs",
    
    # Tab per la scheda "QC"
    tabPanel("NSolver's QC",
             plotOutput("MatriceNonFilt", width = "100%", height = "10px"),
             plotOutput("inputAnalQC", width = "100%", height = "400px"),
             plotOutput("inputAnalBD", width = "100%", height = "400px")
             
    ),
    
    # Tab per la scheda "Altri"
    tabPanel("Unfiltered Samples",
             plotOutput("inputAnalpre1", width = "100%", height = "400px"),
             plotOutput("inputAnalpre2", width = "100%", height = "400px"),
             plotOutput("inputAnalpre3", width = "100%", height = "400px"),
             plotOutput("inputAnalpre4", width = "100%", height = "400px"),
             plotOutput("inputAnalpre5", width = "100%", height = "400px"),
             plotOutput("inputAnalpre6", width = "100%", height = "400px"),
             plotOutput("inputAnalpre7", width = "100%", height = "400px"),
             plotOutput("Heat5", width = "100%", height = "400px"),
    ),
    #si dovrebbe creare un observe nel server che si aggionra ion modo dinamico andando a cancellare nel casio di
 #mancanza di outlier il tabpanel
 
            tabPanel("Filtered Samples",
             plotOutput("inputAnalpost1", width = "100%", height = "400px"),
             plotOutput("inputAnalpost2", width = "100%", height = "400px"),
             plotOutput("inputAnalpost3", width = "100%", height = "400px"),
             plotOutput("inputAnalpost4", width = "100%", height = "400px"),
             plotOutput("inputAnalpost5", width = "100%", height = "400px"),
             plotOutput("inputAnalpost6", width = "100%", height = "400px"),
             plotOutput("inputAnalpost7", width = "100%", height = "400px"),
             
             plotOutput("Heat", width = "100%", height = "400px"), # Heatmap anche qui
             
             #aggiunger i post filt
             )
      )
    
  )
 
)
#scaricare matrice in modo automatico così da averla sempre in input senza dover fare dei renderplot uno dentro l'altro
server <- function(input, output, session) {
  #infile=reactive({c(parseDirPath(roots = volumes, input$Btn_GetFolder))})
  
  # Variabile reattiva per tenere traccia degli outlier
  #outlier <- reactiveVal(NULL)
  
  volumes = c(home = "C:/Users/simon/")
  observe({  
    shinyDirChoose(input, "Btn_GetFolder", roots = volumes)
    
    
  }) #scelta dei file+
  
  output$MatriceNonFilt <- renderTable({
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    
    # Variabile reattiva per tenere traccia degli outlier
    #outlier<- NULL
    df <- NULL
    data <- NULL
    Normdata <- NULL
    Rawdata <- NULL
    RawdataDGE <- NULL
    NormdataDGE <- NULL
    NumberOfPatients <- NULL
    nin <- NULL
    idx <- NULL
    columnz <- NULL
    CodeClass <- NULL
    ctrls <- NULL
    outliers <- NULL
    mediana <- NULL
    treSUP <- NULL
    treMIN <- NULL
    out <- NULL
    Mypalette16 <- NULL
    count <- NULL
    colorz <- NULL
    matrix.raw.HK <- NULL
    nopts <- NULL
    outlier <- NULL
    ImagingQCFILT <- NULL
    
    ###scelgo file###
    df <- read.markup.RCC(
      rcc.path = parseDirPath(roots = volumes, input$Btn_GetFolder),
      rcc.pattern = "*.RCC|*.rcc",
      include = NULL, nprobes = -1,
      exclude = NULL
    )
    
    # c(data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
    #      nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
    #      Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,ImagingQCFILT)<-NULL
    
    data <- NanoStringNorm(
      x = df,
      anno = NA,
      CodeCount = input$codecount,
      Background = input$background,
      SampleContent = input$samplecontent,
      OtherNorm = input$othernorm,
      round.values = FALSE,
      is.log = FALSE,
      take.log = TRUE,
      return.matrix.of.endogenous.probes = FALSE
    )
    
    rownames(data$raw.data) <- rownames(data$normalized.data)
    Normdata<- data$normalized.data[,4:ncol(data$normalized.data)]
    Rawdata<- data$raw.data[,4:ncol(data$raw.data)]
    RawdataDGE <- Rawdata
    NormdataDGE <- Normdata
    
    NumberOfPatients = ncol(data$normalized.data[4:ncol(data$normalized.data)])
    nin<- df$header[15,]
    nin <- names(nin)
    idx <- order(nin)
    
    columnz = 3 + 1:NumberOfPatients
    CodeClass <- as.factor(df$x$CodeClass)
    ctrls<- which(CodeClass == "Negative" | CodeClass == "Positive")
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,
         file = "saved_objects.RData")
  }) 
  output$inputAnalQC<- renderPlot({ ###ImagingQC
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    
    print("QC")
    
    load("saved_objects.RData")
    
    
    ###imagingqc###
    ImagingQC<-(as.numeric(df$header["FovCounted",])/as.numeric(df$header["FovCount",]))*100
    names(ImagingQC)=colnames(df$header)
    ImagingQC<-ImagingQC[order(ImagingQC, decreasing = FALSE)]
    barplot(ImagingQC,las=2, cex.names = ifelse(NumberOfPatients > 95, 0.7,
                                                ifelse(NumberOfPatients > 80, 0.8,
                                                       ifelse(NumberOfPatients > 70, 0.9, 1))),
            main = ifelse(any(ImagingQC<75),paste(sum(ImagingQC<75), "Sample/s have less than 75% FOV"),"All samples have FOV higher than 75%"),
            col = ifelse(ImagingQC < 75, "red", "green3")
            ,lwd=2,
            xlab = "",
            ylab = "FOV %")
    abline(h= 75, lwd = 3, lty = 2, col ="black")
    
    
    #a<-recordPlot()
    ImagingQCFILT<-ImagingQC>75
    outlier<-ImagingQCFILT[ImagingQCFILT==F] #vanno agg insieme nomi e non 
    outliers<-paste0(names(ImagingQCFILT[ImagingQCFILT == FALSE]), ".RCC") #vanno agg insieme solo nomi
    #lista_creata[["a"]] <-a
    #outliers = c(paste0(names(outlier[outlier== FALSE]), ".RCC")
    controlli<- paste0(colnames(Normdata)[str_detect(colnames(Normdata), "ctrl|CTRL")],".RCC")
    
    
    #salvare la matrice ha senso nel momento che quello che salverai sarà direttamente
    #la matrice senza outliers trivati (forse), dovremo aggioranre tutti insieme ogni volta che si tratta di 
    #trobare degli outliers siccome dobbiamo aggiornare i parametri
    
    ###salvo file###
    df <- read.markup.RCC(
      rcc.path = parseDirPath(roots = volumes, input$Btn_GetFolder),
      rcc.pattern = "*.RCC|*.rcc",
      include = NULL, nprobes = -1,
      exclude = outliers
    )
    
    data <- NanoStringNorm(
      x = df,
      anno = NA,
      CodeCount = input$codecount,
      Background = input$background,
      SampleContent = input$samplecontent,
      OtherNorm = input$othernorm,
      round.values = FALSE,
      is.log = FALSE,
      take.log = TRUE,
      return.matrix.of.endogenous.probes = FALSE
    )
    
    rownames(data$raw.data) <- rownames(data$normalized.data)
    Normdata<- data$normalized.data[,4:ncol(data$normalized.data)]
    Rawdata<- data$raw.data[,4:ncol(data$raw.data)]
    RawdataDGE <- Rawdata
    NormdataDGE <- Normdata
    
    NumberOfPatients = ncol(data$normalized.data[4:ncol(data$normalized.data)])
    nin<- df$header[15,]
    nin <- names(nin)
    idx <- order(nin)
    
    columnz = 3 + 1:NumberOfPatients
    CodeClass <- as.factor(df$x$CodeClass)
    ctrls<- which(CodeClass == "Negative" | CodeClass == "Positive")
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,ImagingQCFILT,
         file = "saved_objects.RData") # metodo per salvarla in locale in modo da usare come input 
    
    print(outlier)
  })
  output$inputAnalBD<-renderPlot({
    
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    ###carico###
    print("BD")
    
    load("saved_objects.RData")
    
    ###bindingdensity###
    BindingDensity<-as.numeric(df$header["BindingDensity",])
    binding<-BindingDensity >0.1 & BindingDensity <2.25
    names(binding)=colnames(df$header) #associo nomi a chi rispetto o no parmaetri 
    names(BindingDensity)=colnames(df$header) #asoccio nomi alla binding density 
    BindingDensity<-BindingDensity[order(BindingDensity, decreasing = FALSE)]
    
    barplot(BindingDensity, las=2, cex.names = ifelse(NumberOfPatients > 95, 0.9,
                                                      ifelse(NumberOfPatients > 80, 1,
                                                             ifelse(NumberOfPatients > 70, 1.1, 1.2)))
            , ylim=c(0,3), 
            main = ifelse(all(binding<-BindingDensity >0.1 & BindingDensity <2.25),
                          "All samples meet the required Bindinding density","Some samples do not meet the Binding Density requirements"),
            col= ifelse(BindingDensity <0.1 | BindingDensity >2.25,"red","green3")
            ,lwd=2,
            xlab = "",
            ylab = "Binding Density")
    
    abline(h=0.1 , lwd = 3, lty = 2, col ="black")
    abline(h=2.25, lwd = 3, lty = 2, col ="black")
    
    outlier<-binding[binding==F] #vanno agg insieme nomi e non 
    outliers<-paste0(names(binding[binding == FALSE]), ".RCC") #vanno agg insieme solo nomi
    
    df <- read.markup.RCC(
      rcc.path = parseDirPath(roots = volumes, input$Btn_GetFolder),
      rcc.pattern = "*.RCC|*.rcc",
      include = NULL, nprobes = -1,
      exclude = outliers
    )
    
    data <- NanoStringNorm(x = df,#data,
                           anno = NA,
                           CodeCount = 'geo.mean',
                           Background = 'mean.2sd',
                           SampleContent = 'total.sum',
                           round.values = FALSE,
                           is.log = FALSE,
                           take.log = TRUE,
                           return.matrix.of.endogenous.probes = FALSE )
    
    
    rownames(data$raw.data) <- rownames(data$normalized.data)
    Normdata<- data$normalized.data[,4:ncol(data$normalized.data)]
    Rawdata<- data$raw.data[,4:ncol(data$raw.data)]
    RawdataDGE <- Rawdata
    NormdataDGE <- Normdata
    
    NumberOfPatients = ncol(data$normalized.data[4:ncol(data$normalized.data)])
    nin<- df$header[15,]
    nin <- names(nin)
    idx <- order(nin)
    
    columnz = 3 + 1:NumberOfPatients
    CodeClass <- as.factor(df$x$CodeClass)
    ctrls<- which(CodeClass == "Negative" | CodeClass == "Positive")
    
    
    mediana<-apply(as.matrix(Normdata),2, median)
    treSUP<-median(as.matrix(Normdata))+ mad(as.matrix(Normdata))
    treMIN<-median(as.matrix(Normdata))- mad(as.matrix(Normdata))
    out<-mediana>treMIN & mediana<treSUP
    #print(out)
    #print(outliers) #mi da dei valori fantasma .RCC .RCC
    outliers<-c(paste0(names(binding[binding== FALSE]), ".RCC"),paste0(names(ImagingQCFILT[ImagingQCFILT == FALSE]), ".RCC"))
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,ImagingQCFILT,
         file = "saved_objects.RData")
  })
  output$inputAnalpre1<- renderPlot({ ###
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    
    print("LineePlot")
    
    load("saved_objects.RData")
    
    Mypalette16 <- c("dodgerblue2", "#E31A1C", "yellow2", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
    
    colorz <- Mypalette16[1:length(ctrls)]
    for ( i in ctrls ) {
      plot (log2(as.numeric(data$raw.data [i, columnz ]
                            [idx]
      )),
      ylim = c(#log2(min(data$raw.data [ctrls, columnz])+0.1),
        ifelse(length(ctrls)>21,-3,
               ifelse(length(ctrls)>14,-2,-1)), #this line automatically gives you the correct amount of space for placing the legend
        log2(max(data$raw.data [ctrls, columnz]))),
      col = colorz[i-(ctrls[1]-1)],
      lwd=2,
      xlab = "",
      ylab = "Reads (log2 scaled)",
      cex.axis= 1,
      xaxt = "n",
      type = "l",
      las = 2,
      main="Positive and Negative Controls")
      par (new=TRUE)
    }
    axis(side=1, at=1:NumberOfPatients, labels = colnames(data$raw.data)[columnz][idx], las = 2, cex = 1, cex.axis= 0.7)
    legend("bottom", legend = data$raw.data$Name[ctrls], bty = "n", cex= 0.6, lwd = 3, col = colorz,
           ncol= 6)#modify here if you want to add space to the legend
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,
         file = "saved_objects.RData")
  })
  
  output$inputAnalpre2<- renderPlot({ ###
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    
    print("Boxplot")
    
    load("saved_objects.RData")
    
    par (mfrow = c(1,1), mar =c(5, 4, 2, 2) + 0.1, las=2)
    boxplot(log2(data$raw.data[,3+idx]), 
            ylim= c(0,round(max(log2(data$raw.data[,3+idx])), digits = -1)),
            las =2, cex.axis = 0.7, col= ifelse(mediana>treMIN & mediana<treSUP,"green3","red")
    )
    
    title ("Signal distributions - Raw data")
  })
  
  output$inputAnalpre3<- renderPlot({ ###
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    print("Signal distributions - Cartridge-corrected HK normalized data")
    load("saved_objects.RData")
    
    boxplot(data$normalized.data[,3+idx],ylim= c(0,round(max(log2(data$raw.data[,3+idx])), digits = -1)),
            las =2, cex.axis = 0.7, col= ifelse(mediana>treMIN & mediana<treSUP,"green3","red"))
    title ("Signal distributions - Cartridge-corrected HK normalized data")
    abline(h = treSUP)
    abline(h = treMIN)
    text(x = -1, y=treSUP+0.4, label= round(treSUP, 3), col= "darkblue", font= 2, xpd = TRUE)
    text(x = -1, y=treMIN+0.4, label= round(treMIN, 3), col= "darkblue", font= 2, xpd = TRUE)
    
    nopts = NumberOfPatients
    
    count = grep ("Housekeeping", data$normalized.data$Code.Class )  ### nr housekeepings
    colorz<- Mypalette16[1:length(count)]
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,
         file = "saved_objects.RData")
    
  })
  output$inputAnalpre4<- renderPlot({ ###
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    print("Housekeeping")
    load("saved_objects.RData")
    
    matrix.raw.HK <- data$raw.data [ count, ]
    
    boxplot(matrix.raw.HK [,columnz][idx],
            ylim = c(0,max(matrix.raw.HK[,columnz][idx])),
            names = colnames(data$raw.data)[columnz][idx] ,
            las=2,
            ylab = "Raw counts",
            cex = 0.8,
            cex.axis = 0.7,
            main= "Housekeeping genes distributions",
            col= ifelse(mediana>treMIN & mediana<treSUP,"green3","red")
    )
    
    write.table(apply(matrix.raw.HK [,columnz],2,quantile), file="HK.matrix.txt", sep="\t", quote = F, )
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,
         file = "saved_objects.RData")
  })
  output$inputAnalpre5<- renderPlot({ ###
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    print("matplot")
    load("saved_objects.RData")
    
    matplot (t(apply(matrix.raw.HK [,columnz][idx],2,quantile)),
             type = "l",
             pch =".",
             lwd=2,
             col = Mypalette16[1:5], #we are plotting the quantiles
             las = 2,
             xaxt = "n" ,
             ylab = "Raw counts",
             ylim= c(0,max(matrix.raw.HK[,columnz][idx])),
             cex.axis = 0.7,
             main= "Housekeeping genes distributions quantiles")
    
    axis(side=1, at=1:nopts, labels = colnames(data$raw.data)[columnz][idx], las = 2, cex = 0.8, cex.axis = 0.7)
    legend("top", legend = colnames(t(apply(matrix.raw.HK [,columnz][idx],2,quantile))), bty = "n", cex= 0.6, lwd = 3, col = Mypalette16[1:5], ncol= 5)
    columnz = 3 + 1:nopts ## nr experiments
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,
         file = "saved_objects.RData")
  })
  output$inputAnalpre6<- renderPlot({ ###
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    print("for2")
    load("saved_objects.RData")
    
    for ( i in count) {
      plot (as.numeric(data$raw.data [i, columnz ][idx]),
            ylim = c(min(data$raw.data [count, columnz]),max(data$raw.data [count, columnz])+5),
            col = colorz[which(count %in% i)],
            lwd=2,
            xlab = "",
            ylab = "Raw counts",
            xaxt = "n",
            type = "l",
            las = 2,
            cex.axis = 0.7,
            main= "Housekeeping genes raw expression values")
      par (new=TRUE)
    }
    axis(side=1, at=1:nopts, labels = colnames(data$raw.data)[columnz][idx], las = 2, cex = 0.8,cex.axis = 0.7)
    legend("top", legend = data$raw.data$Name [count], bty = "n", cex= 0.6, lwd = 3, col = colorz, ncol= 5)
    
    columnz = 3 + 1:nopts ## nr experiments
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,
         file = "saved_objects.RData")
  })
  
  output$inputAnalpre7<- renderPlot({ ###
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    print("bro")
    load("saved_objects.RData")
    for ( i in count) {
      plot (as.numeric(data$normalized.data [i, columnz ][idx]),
            ylim = c(min(data$normalized.data [count, columnz]),max(data$normalized.data [count, columnz])+5),
            col = colorz[which(count %in% i)],
            lwd=2,
            xlab = "",
            ylab = "Normalized counts (log2 scaled)",
            xaxt = "n",
            type = "l",
            las = 2,
            cex.axis = 0.7,
            main= "Housekeeping genes normalized expression values"
      )
      par (new=TRUE)
    }
    axis(side=1, at=1:nopts, labels = colnames(data$raw.data)[columnz][idx], las = 2, cex = 0.8, cex.axis = 0.7)
    legend("top", legend = data$raw.data$Name [count], bty = "n", cex= 0.6, lwd = 3, col = colorz, ncol= 5)
    
    print(outlier)
    out<-NULL
    mediana<-0
    treMIN<-0
    treSUP<-0
    mediana<-apply(as.matrix(Normdata),2, median)
    treSUP<-median(as.matrix(Normdata))+ mad(as.matrix(Normdata))
    treMIN<-median(as.matrix(Normdata))- mad(as.matrix(Normdata))
    out<-mediana>treMIN & mediana<treSUP
    controlli<- paste0(colnames(Normdata)[str_detect(colnames(Normdata), "ctrl|CTRL|STANDARD|standard|Standard|Control|CONTROL|control")],".RCC")
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,controlli,
         file = "saved_objects.RData")
    
  })
  output$Heat5<- renderPlot({ ###
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    print("Heat")
    load("saved_objects.RData")
    
    while(all(out)==FALSE){
      #tutti questi sono ariabili da dichiarare sia fuori che poi dentro al ciclo così da essere aggiornate sempre durante la sua esecuzione 
      mediana<-apply(as.matrix(Normdata),2, median)
      treSUP<-median(as.matrix(Normdata))+ mad(as.matrix(Normdata))
      treMIN<-median(as.matrix(Normdata))- mad(as.matrix(Normdata))
      out<-mediana>treMIN & mediana<treSUP
      print(ifelse(all(out==TRUE),"sono tutti veri", "rimane ancora qualche outlier")) #il controllo non funziona
      #newNorm<-Normdata[,out] #questo subset non ha senso vrearlo siccome dobbiamo fare altri sub set?
      outlier<-c(outlier,out[out==F])
      
      if (all(out)==FALSE) { 
        outliers = c(paste0(names(out[out == FALSE]), ".RCC"),outliers)
      }
      
      df <- read.markup.RCC(
        rcc.path = parseDirPath(roots = volumes, input$Btn_GetFolder),
        rcc.pattern = "*.RCC|*.rcc",
        include = NULL, nprobes = -1,
        exclude = outliers
      )
      
      data <- NanoStringNorm(x = df,#data,
                             anno = NA,
                             CodeCount = 'geo.mean',
                             Background = 'mean.2sd',
                             SampleContent = 'total.sum',
                             round.values = FALSE,
                             is.log = FALSE,
                             take.log = TRUE,
                             return.matrix.of.endogenous.probes = FALSE )
      
      
      rownames(data$raw.data) <- rownames(data$normalized.data)
      Normdata<- data$normalized.data[,4:ncol(data$normalized.data)]
      Rawdata<- data$raw.data[,4:ncol(data$raw.data)]
      RawdataDGE <- Rawdata
      NormdataDGE <- Normdata
    }
    
    CartridgeID<- as.data.frame(t(df$header["CartridgeID",]))
    
    Cartridges<- levels(as.factor(as.numeric(str_extract(CartridgeID$CartridgeID,"\\d+"))))
    
    for (n in Cartridges) {
      CartridgeID[grep(x = CartridgeID$CartridgeID, pattern = n),] <- paste("Batch", n)
    }
    CartridgeID$CartridgeID <- as.factor(CartridgeID$CartridgeID)
    
    ##ora mergi
    ptsinfo<- #cbind(# comment this line if no ptsinfo (xlsx) is available
      CartridgeID
    #,pheno) # comment this line if no ptsinfo (xlsx) is available
    NAMES <- rownames(ptsinfo)
    ptsinfo <- as.data.frame(apply(ptsinfo,2,str_to_upper))
    rownames(ptsinfo) <- NAMES
    
    #colorz <- c("yellow2", "grey50","black","white") 
    norm.data <-  data$normalized.data [4:ncol(data$normalized.data)]
    raw.data <-  data$raw.data [4:ncol(data$raw.data)]
    norm.data<- norm.data[,rownames(CartridgeID)]
    raw.data<- raw.data[,rownames(CartridgeID)]
    #togli gli 0
    matrix.for.plot = norm.data[ rowSums(norm.data) != 0 , ]
    #adesso fai la heatmap unsupervised per vedere se c'è del batch effect
    
    annot_col = data.frame( "CartridgeID"=factor(ptsinfo[,"CartridgeID"]))
    colnames(annot_col) <- "CartridgeID"
    rownames(annot_col) = colnames(matrix.for.plot)
    
    ######USE THIS CHUNK TO MAKE NA INTELLEGIBLE BY PHEATMAP
    NAMES <- rownames(annot_col)
    annot_col <- apply(annot_col, 2, as.character)
    annot_col[is.na(annot_col)] <- "NA"#change to character for label recognition
    annot_col <- as.data.frame(apply(annot_col, 2, as.factor))
    rownames(annot_col) <- NAMES
    
    pippo<- as.factor(as.vector(annot_col[,"CartridgeID"]))
    
    
    livelli<- levels(pippo)
    colori <- colorz[1:length(livelli)]
    names(colori) <- livelli
    
    annot_colors <- list( "CartridgeID"= colori)
    names(annot_colors) <- "CartridgeID"
    
    annot_col[] <- lapply(annot_col, as.factor)
    
    pheatmap(matrix.for.plot,
             clustering_distance_cols = "canberra",
             clustering_method = "ward.D2", 
             annotation_col = annot_col,
             annotation_colors = annot_colors,
             fontsize = 6,
             fontsize_row = 1/nrow(matrix.for.plot)*100,
             fontsize_col = 1/ncol(matrix.for.plot)*200,
             scale="row",
             cluster_rows = T,
             cluster_cols = T,
             main = "CartridgeID unsupervised clustering", #quindi prende dei valori filtrati 
             legend = F,
             breaks =  seq(-2.5, 2.5, length.out= 500),
             color = colorRampPalette(rev(brewer.pal(6,name="RdYlBu")))(500)
    )
    
    NumberOfPatients = ncol(data$normalized.data[4:ncol(data$normalized.data)])
    nopts = NumberOfPatients
    
    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,controlli,
         file = "saved_objects.RData")
    
    observe({
      print("entro observe")
      
      if (!all(outlier)) {
        print("vedi")
        updateTabsetPanel(session, "main_tabs", selected = "Overview")
        showTab("main_tabs", "Filtered Samples")
        ###
       
      } else {
        print("nascondi") #controllare perchè sono invertiti cioe che se non ha outlier mi parte questo
        hideTab("main_tabs", "Filtered Samples")
      }
    })
  })
  output$inputAnalpost1<- renderPlot({
    load("saved_objects.RData")
    ######plot filtrato 
    #plot con linee colorate 1 non filtrato e non aggiornato a questo codice 
    NumberOfPatients = ncol(data$normalized.data[4:ncol(data$normalized.data)])
    #png("./plots_StandardNorm/Controls.png")
    nin<- df$header[15,]
    nin <- names(nin)
    idx <- order(nin)
    #idx = order(df$header[15,])
    columnz = 3 + 1:NumberOfPatients
    CodeClass <- as.factor(df$x$CodeClass)
    ctrls<- which(CodeClass == "Negative" | CodeClass == "Positive")
    
    #palette of 16 colors in contrast
    Mypalette16 <- c("dodgerblue2", "#E31A1C", "yellow2", "#6A3D9A", "#FF7F00", "black", "gold1", "skyblue2", "palegreen2", "#FDBF6F", "gray70", "maroon", "orchid1", "darkturquoise", "darkorange4", "brown")
    
    colorz <- Mypalette16[1:length(ctrls)]
    
    for ( i in ctrls ) {
      plot (log2(as.numeric(data$raw.data [i, columnz ]
                            [idx]
      )),
      ylim = c(#log2(min(data$raw.data [ctrls, columnz])+0.1),
        ifelse(length(ctrls)>21,-3,
               ifelse(length(ctrls)>14,-2,-1)), #this line automatically gives you the correct amount of space for placing the legend
        log2(max(data$raw.data [ctrls, columnz]))),
      col = colorz[i-(ctrls[1]-1)],
      lwd=2,
      xlab = "",
      ylab = "Reads (log2 scaled)",
      xaxt = "n",
      type = "l",
      las = 2,
      main="Filtered Sample - Positive and Negative Controls")
      par (new=TRUE)
      
    }
    
    axis(side=1, at=1:NumberOfPatients, labels = colnames(data$raw.data)[columnz][idx], las = 2, cex = 1, cex.axis = 0.7)
    legend("bottom", legend = data$raw.data$Name [ctrls], bty = "n", cex= 0.6, lwd = 3, col = colorz,
           ncol= 6)#modify here if you want to add space to the legend

    save(df, data, Normdata, Rawdata, RawdataDGE, NormdataDGE, NumberOfPatients, 
         nin, idx, columnz, CodeClass, ctrls, outliers, mediana, treSUP, treMIN, out,
         Mypalette16, count, colorz, matrix.raw.HK,nopts,outlier,controlli,
         file = "saved_objects.RData")
  })
  
  output$inputAnalpost2<- renderPlot({
    
    load("saved_objects.RData")
    # secondo grafico  filtrato boxplot non normalizzate non filtrato
    treSUP<-median(as.matrix(Normdata))+ mad(as.matrix(Normdata))
    treMIN<-median(as.matrix(Normdata))- mad(as.matrix(Normdata))
    par (mfrow = c(1,1), mar =c(5, 4, 2, 2) + 0.1, las=2)
    
    boxplot(log2(data$raw.data[,3+idx]), ylim= c(0,round(max(log2(data$raw.data[,3+idx])), digits = -1)
                                                 
    ),
    col= ifelse(mediana>treMIN & mediana<treSUP,"green3","red"),
    las =2, cex.axis = 0.7)
    title ("Filtered samples - Signal distributions - Raw data")
    
  })
  
  output$inputAnalpost3<- renderPlot({
    #terzo grafico  filtrato
    load("saved_objects.RData")
    
    boxplot(data$normalized.data[,3+idx], ylim= c(0,round(max(log2(data$raw.data[,3+idx])), digits = -1)
                                                  
    ),
    col= ifelse(mediana>treMIN & mediana<treSUP,"green3","red"),
    las =2, cex.axis = 0.7)
    title ("Filtered samples - Signal distributions - Cartridge-corrected HK normalized data")
    abline(h = treSUP)
    abline(h = treMIN)
    text(x = -1, y=treSUP+0.4, label= round(treSUP, 3), col= "darkblue", font= 2, xpd = TRUE)
    text(x = -1, y=treMIN+0.4, label= round(treMIN, 3), col= "darkblue", font= 2, xpd = TRUE)
    
  })
  
  output$inputAnalpost4<- renderPlot({
    #quarto grafico  filtrato housekeeping 
    load("saved_objects.RData")
    
    nopts = NumberOfPatients
    
    count = grep ("Housekeeping", data$normalized.data$Code.Class )  ### nr housekeepings
    colorz<- Mypalette16[1:length(count)]
    ### plot bloxplot of housekeeping dsitributions
    matrix.raw.HK = data$raw.data [ count, ]
    
    boxplot(matrix.raw.HK [,columnz][idx],
            ylim = c(0,max(matrix.raw.HK[,columnz][idx])),
            names = colnames(data$raw.data)[columnz][idx] ,
            las=2,
            ylab = "Raw counts",
            cex = 0.8,
            cex.axis = 0.7,
            main= "Filtered samples - Housekeeping genes distributions",
            col= ifelse(mediana>treMIN & mediana<treSUP,"green3","red")
    )
    
    write.table(apply(matrix.raw.HK [,columnz],2,quantile), file="HK.matrix.txt", sep="\t", quote = F, )
  })
  
 
 ####?????######## 
  output$inputAnalpost5<- renderPlot({
    load("saved_objects.RData") 
    
    
    
    matplot (t(apply(matrix.raw.HK [,columnz][idx],2,quantile)),
             type = "l",
             pch =".",
             lwd=2,
             col = Mypalette16[1:5], #we are plotting the quantiles
             las = 2,
             xaxt = "n" ,
             ylab = "Raw counts",
             ylim= c(0,max(matrix.raw.HK[,columnz][idx])),
             cex.axis = 0.7,
             main= "Filtered samples - Housekeeping genes distributions quantiles")
    
    axis(side=1, at=1:nopts, labels = colnames(data$raw.data)[columnz][idx], las = 2, cex = 0.8, cex.axis = 0.7)
    legend("top", legend = colnames(t(apply(matrix.raw.HK [,columnz][idx],2,quantile))), bty = "n", cex= 0.6, lwd = 3, col = Mypalette16[1:5], ncol= 5)
  })
  
  output$inputAnalpost6<- renderPlot({
    load("saved_objects.RData")  
    
    columnz = 3 + 1:nopts ## nr experiments
    for ( i in count) {
      plot (as.numeric(data$raw.data [i, columnz ][idx]),
            ylim = c(min(data$raw.data [count, columnz]),max(data$raw.data [count, columnz])+5),
            col = colorz[which(count %in% i)],
            lwd=2,
            xlab = "",
            ylab = "Raw counts",
            xaxt = "n",
            type = "l",
            las = 2,
            cex.axis = 0.6,
            main= "Filtered samples - Housekeeping genes raw expression values")
      par (new=TRUE)
    }
    axis(side=1, at=1:nopts, labels = colnames(data$raw.data)[columnz][idx], las = 2, cex = 0.8,cex.axis = 0.7)
    legend("top", legend = data$raw.data$Name [count], bty = "n", cex= 0.6, lwd = 3, col = colorz, ncol= 5)
  })   
  
  output$inputAnalpost7<- renderPlot({
    load("saved_objects.RData")  
    
    columnz = 3 + 1:nopts ## nr experiments
    for ( i in count) {
      plot (as.numeric(data$normalized.data [i, columnz ][idx]),
            ylim = c(min(data$normalized.data [count, columnz]),max(data$normalized.data [count, columnz])+5),
            col = colorz[which(count %in% i)],
            lwd=2,
            xlab = "",
            ylab = "Normalized counts (log2 scaled)",
            xaxt = "n",
            type = "l",
            las = 2,
            cex.axis = 0.7,
            main= "Filtered samples - Housekeeping genes normalized expression values"
      )
      par (new=TRUE)
    }
    axis(side=1, at=1:nopts, labels = colnames(data$raw.data)[columnz][idx], las = 2, cex = 0.8, cex.axis = 0.7)
    legend("top", legend = data$raw.data$Name [count], bty = "n", cex= 0.6, lwd = 3, col = colorz, ncol= 5)
    
    NormDataFiltered = data$normalized.data
  })
  
  
  output$Heat<- renderPlot({ ###
    req(input$Btn_GetFolder)  # Assicura che l'utente abbia selezionato una cartella
    print("Heat")
    load("saved_objects.RData")
    
    
    CartridgeID<- as.data.frame(t(df$header["CartridgeID",]))
    
    Cartridges<- levels(as.factor(as.numeric(str_extract(CartridgeID$CartridgeID,"\\d+"))))
    
    for (n in Cartridges) {
      CartridgeID[grep(x = CartridgeID$CartridgeID, pattern = n),] <- paste("Batch", n)
    }
    CartridgeID$CartridgeID <- as.factor(CartridgeID$CartridgeID)
    
    ##ora mergi
    ptsinfo<- #cbind(# comment this line if no ptsinfo (xlsx) is available
      CartridgeID
    #,pheno) # comment this line if no ptsinfo (xlsx) is available
    NAMES <- rownames(ptsinfo)
    ptsinfo <- as.data.frame(apply(ptsinfo,2,str_to_upper))
    rownames(ptsinfo) <- NAMES
    
    #colorz <- c("yellow2", "grey50","black","white") 
    norm.data <-  data$normalized.data [4:ncol(data$normalized.data)]
    raw.data <-  data$raw.data [4:ncol(data$raw.data)]
    norm.data<- norm.data[,rownames(CartridgeID)]
    raw.data<- raw.data[,rownames(CartridgeID)]
    #togli gli 0
    matrix.for.plot = norm.data[ rowSums(norm.data) != 0 , ]
    #adesso fai la heatmap unsupervised per vedere se c'è del batch effect
    
    annot_col = data.frame( "CartridgeID"=factor(ptsinfo[,"CartridgeID"]))
    colnames(annot_col) <- "CartridgeID"
    rownames(annot_col) = colnames(matrix.for.plot)
    
    ######USE THIS CHUNK TO MAKE NA INTELLEGIBLE BY PHEATMAP
    NAMES <- rownames(annot_col)
    annot_col <- apply(annot_col, 2, as.character)
    annot_col[is.na(annot_col)] <- "NA"#change to character for label recognition
    annot_col <- as.data.frame(apply(annot_col, 2, as.factor))
    rownames(annot_col) <- NAMES
    
    pippo<- as.factor(as.vector(annot_col[,"CartridgeID"]))
    
    
    livelli<- levels(pippo)
    colori <- colorz[1:length(livelli)]
    names(colori) <- livelli
    
    annot_colors <- list( "CartridgeID"= colori)
    names(annot_colors) <- "CartridgeID"
    
    annot_col[] <- lapply(annot_col, as.factor)
    
    pheatmap(matrix.for.plot,
             clustering_distance_cols = "canberra",
             clustering_method = "ward.D2", 
             annotation_col = annot_col,
             annotation_colors = annot_colors,
             fontsize = 6,
             fontsize_row = 1/nrow(matrix.for.plot)*100,
             fontsize_col = 1/ncol(matrix.for.plot)*200,
             scale="row",
             cluster_rows = T,
             cluster_cols = T,
             main = "CartridgeID unsupervised clustering", #quindi prende dei valori filtrati 
             legend = F,
             breaks =  seq(-2.5, 2.5, length.out= 500),
             color = colorRampPalette(rev(brewer.pal(6,name="RdYlBu")))(500)
    )
    
    
    
  })
  #fare secondo panello poi si  oensa a if.
  print("prova per vedere se parte a prescindere il codice")
  load("saved_objects.RData")  
  
  
  
  
  
  
  
}

shinyApp(ui = ui, server = server)
