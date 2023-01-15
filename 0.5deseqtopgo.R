
####PREPARACI�N DE PAQUETERIAS####
library(topGO)
library(tximport)#Paqueter�a especifica para importar archivo de cuantificaci�n al formato de DESeq2
library(readr)#Paqueter�a para lectura de CSV
library(DESeq2)#Paqueter�a para el c�lculo de los DGEs
library(BiocParallel)#Paqueter�a para procesos en paralelizaci�n
library(IHW)  #Paqueter�a para filtar los DGEs por pruebas de hip�tesis
library(vsn) #Paqueter�a para los SD plots
library(pheatmap)  #Paqueteria para heatmaps
library(RColorBrewer)  #Paqueter�a para clusters
library( "genefilter" )
library( "gplots")
library("ggplot2")
library("vegan")
library("FactoMineR")
library("factoextra")
library("dunn.test")
library(readr)
register(SnowParam(6)) #Paralelizaci�n de nucleos
####iMPORTAR TSV's DE KALLISTO####
dir   <- "~/Documents/sctld_workspace/kal/of/k_of_host" #Define el directorio con los archivos CSVlist.files(dir)#Lista de archivos dentro de la carpeta
samples <- read.table(file.path("~/Documents/sctld_workspace/kal/of/of.txt"), header = TRUE)    #Objeto con los datos (factores) de las librer�as
files <- file.path(dir,samples$sample) #Objeto con una lista con los archivosobject with filepaths/samplename.tsv
names(files) <- samples$sample  #Designa el nombre de las librerias del archivo "samples" a los directorios "files"
all(file.exists(files)) #Doble check de concordancia entre los nombres las librer�as del objeto "files", TRUE indica concordancia
file.exists(files)#Final check, si alg�n nombre no concuerda se obtiene un "FALSE"
#tx2gene <- read.delim("~/Documents/DGE/symkal/sym_gene_to_trans.txt", header=FALSE)
#tx2gene[,3] = tx2gene$V1
#tx2gene = tx2gene[,2:3]
txi <- tximport(files, type = "kallisto", txOut = TRUE)#Sigue el nombre de los directorios y crea un objeto con los datos de expresi�n de todos los CSV
####DELIMITACI�N DE FACTORES####
#cond = factor(c("sctld","ctrl","sctld","ctrl","sctld","sctld","ctrl","sctld","ctrl","sctld","ctrl","ctrl","sctld","ctrl","ctrl","sctld","ctrl","sctld","ctrl","sctld","ctrl","sctld","ctrl","ctrl","ctrl", "sctld"))
#cond = factor(c("ctrl","sctld","ctrl","sctld","sctld","ctrl","sctld","ctrl","sctld","ctrl","sctld","sctld","sctld","ctrl","sctld","sctld"))
sampleTable<- data.frame(cond) #conjunta los objetos "spcode", "year", "site", "temp" y "cond" en un objeto
rownames(sampleTable) <- colnames(txi$counts) #Designa el nombre de los rows el nombre de los transcritos
as.integer(txi$counts) #Define como "integer" los datos de expresi�n
####PRE-DESEQ####
ddsmf <- DESeqDataSetFromTximport(txi, sampleTable, design = ~ cond)  #Crea un objeto tipo DESeq con los datos de expresi�n "txi", factores "sampletable" y un dise�o de acuerdo a los factores
#design(ddsmf) <- formula(~ 1)
design(ddsmf)
nrow(ddsmf)
keep <- rowSums(counts(ddsmf)) >1
ddsmf <- ddsmf[keep,]
nrow(ddsmf)
keep <- rowSums(counts(ddsmf) >=5) >= 3#Prefiltrado en el que se mantienen los rows (transcritos) con m�nimo 20 lecturas, en 3 librer�as como m�nimo
ddsmf<-ddsmf[keep,]   #Sustituye los rows del objeto dds con los rows que pasaron el pre-filtrado
nrow(ddsmf)
#N�mero de transcritos presentes
colData(ddsmf)   #Factores de la copia del nuevo objeto dds "ddsmf"
levels(ddsmf$cond) #Niveles de cada factor
design(ddsmf)   #Dise�o experimental del objeto "ddsmf"
#ddsmf$cond = relevel(ddsmf$cond, "ctrl")
####DESEQ####
ddsmf <- DESeq(ddsmf)
res = results(ddsmf, alpha = 0.05)
resultsNames(ddsmf)
#resmf<-results(ddsmf, alpha = 0.05, lfcThreshold = 2)#Crea un objeto con los resultados de la expresi�n diferencial, un dataFrame con 6 colums: "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"; Se puede modificar el argumento "alpha" (e.g. alpha=0.05)
#summary(resmf)
rld <-rlog(ddsmf)                                                                  #Transformaci�n de los DGEs a Log2
#vsd <-vst(ddsmf)                                                                   #Transformaci�n de los DGEs a Log2
#ntd <-normTransform(ddsmf)    
resmf_allcond_low_high <- results(ddsmf, filterFun=ihw, alpha = 0.05, contrast = c("cond", "sctld", "ctrl"))
resmf_allcond_low_high <- results(ddsmf, contrast = c("cond", "sctld", "ctrl"))
####Summary####
summary(resmf_allcond_low_high)
resOrdered <- resmf_allcond_low_high[order(resmf_allcond_low_high$padj),]
select_genes_low_high<-rownames(subset(resOrdered, padj < 0.05 ))
dbset  <- assay(rld)[ select_genes_low_high, ]
####MDS####
#dbdist=dist(t(dbset))
#mds=cmdscale(dbdist, eig = TRUE)
#ggordiplots::gg_ordiplot(mds, groups = cond, spiders = TRUE )
####PCA####
dbset=dbset[!duplicated(dbset), ]
X=t(dbset)
res.x<-PCA(X, scale.unit = TRUE, ncp = 5, graph = FALSE, col.w = TRUE)
eig.val <- get_eigenvalue(res.x)
var <- get_pca_var(res.x)
set.seed(123)
res.km <- kmeans(var$coord, centers = 3, nstart = 8)
grp <- as.factor(res.km$cluster)
fviz_pca_biplot(res.x, col.ind = cond,
                palette = c("#0073C2FF", "#EFC000FF", "#868686FF"),
                legend.title = "Cluster", addEllipses = TRUE, ellipse.level=0.95)


p <- fviz_pca_ind(res.x, label="none", habillage=cond,
                  addEllipses=TRUE, ellipse.level=0.95)
print(p)

#Atribución de varianza por PC
fviz_eig(res.x, addlabels = TRUE, ylim = c(0, 50))
# Contributions of variables to PC1
var <- get_pca_var(res.x)
fviz_contrib(res.x, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
fviz_contrib(res.x, choice = "var", axes = 2, top = 10)
# Contributions of variables to PC3
fviz_contrib(res.x, choice = "var", axes = 3, top = 10)
df=as.data.frame(res.x$ind$coord)[,c(1,2)]
df
df$tempcond=c("vc","vc","ic","vv","iv","iv","ic","vv")
#### PC1 Density Plot####
xdensity <- ggplot(df, aes(Dim.1, fill=cond)) + 
  geom_density(alpha=.95) + 
  xlim(-40,40) +
  theme(legend.position = "none")
xdensity
#### PC2 Density Plot####
ydensity <- ggplot(df, aes(Dim.2, fill=cond)) + 
  geom_density(alpha=.95) + 
  theme(legend.position = "none")
ydensity
####Heatmap####
heatmap.2( as.matrix(dbset),labRow=FALSE,labCol = cond, breaks=seq(from=0, to=10,length=25), 
           trace="none", dendrogram = "both", 
           col = colorRampPalette(c(
             "black",
             "blue",
             "yellow",
             "yellow2",
             "yellow3")) ,
           ColSideColors = c( Control="gray", DDD="blue" ,DPN="darkgreen", OHT="orange" )[
             cond ] )
####up down separation####
resmf_vent_temp_up                <- resmf_vent_temp              [resmf_vent_temp$log2FoldChange >0,]
resmf_vent_temp_down              <- resmf_vent_temp              [resmf_vent_temp$log2FoldChange <0,]
resmf_allcond_low_high_up         <- resmf_allcond_low_high       [resmf_allcond_low_high$log2FoldChange >0,]
resmf_allcond_low_high_down       <- resmf_allcond_low_high       [resmf_allcond_low_high$log2FoldChange <0,]
####p value####
x_vt_up             <- resmf_vent_temp_up$padj     # P-ajustado verano e invierno dentro de ventila
x_vt_down           <- resmf_vent_temp_down$padj     # P-ajustado verano e invierno dentro de ventila
x_ac_low_high_up    <- resmf_allcond_low_high_up$padj       # P-ajustado Condiciones (verano e invierno) en ventrila
x_ac_low_high_down  <- resmf_allcond_low_high_down$padj       # P-ajustado Condiciones (verano e invierno) en ventrila
#Crea un objeto con los nombres de los transcritos provenientes de DESeq2
names(x_vt_up)             <- row.names(resmf_vent_temp_up)       # Nombres verano e invierno dentro de ventila
names(x_vt_down)           <- row.names(resmf_vent_temp_down)       # Nombres verano e invierno dentro de ventila
names(x_ac_low_high_up)    <- row.names(resmf_allcond_low_high_up)         # Nombres Condiciones (verano e invierno) en ventrila
names(x_ac_low_high_down)  <- row.names(resmf_allcond_low_high_down)         # Nombres Condiciones (verano e invierno) en ventrila
####cleanNA####
x_vt_up    <- na.omit(x_vt_up)      
x_vt_down  <- na.omit(x_vt_down)  
x_ac_low_high_up    <- na.omit(x_ac_low_high_up)        
x_ac_low_high_down  <- na.omit(x_ac_low_high_down)        
####goterms####
geneID2GO <- readMappings("/home/erick/Documents/workspace/trinotate/go_mch.txt") #carga las categorias GO de trinotate
str(head(geneID2GO)) 
####Funci�n para separar los genes en TOPGO####
topDiffGenes <- function(allScore) { return(allScore < 0.05)}
####TopGO_Object####
#MF_Molecular function
GOdata_mf_ac_low_high_up    <- new("topGOdata", ontology = "MF", allGenes = x_ac_low_high_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_mf_ac_low_high_down  <- new("topGOdata", ontology = "MF", allGenes = x_ac_low_high_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
#BP_BIOLOGICAL PROCESS
GOdata_bp_ac_low_high_up    <- new("topGOdata", ontology = "BP", allGenes = x_ac_low_high_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_bp_ac_low_high_down  <- new("topGOdata", ontology = "BP", allGenes = x_ac_low_high_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
#CC_CELULAR COMPONENT
GOdata_cc_ac_low_high_up    <- new("topGOdata", ontology = "CC", allGenes = x_ac_low_high_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_cc_ac_low_high_down  <- new("topGOdata", ontology = "CC", allGenes = x_ac_low_high_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
####Fisher_parentchild####
resultFisher_mf_ac_low_high_up    <- runTest(GOdata_mf_ac_low_high_up,    algorithm = "parentchild", statistic = "fisher")  #MF
resultFisher_mf_ac_low_high_down  <- runTest(GOdata_mf_ac_low_high_down,  algorithm = "parentchild", statistic = "fisher")  #MF
resultFisher_bp_ac_low_high_up    <- runTest(GOdata_bp_ac_low_high_up,    algorithm = "parentchild", statistic = "fisher")  #BP
resultFisher_bp_ac_low_high_down  <- runTest(GOdata_bp_ac_low_high_down,  algorithm = "parentchild", statistic = "fisher")  #BP
resultFisher_cc_ac_low_high_up    <- runTest(GOdata_cc_ac_low_high_up,    algorithm = "parentchild", statistic = "fisher")  #CC
resultFisher_cc_ac_low_high_down  <- runTest(GOdata_cc_ac_low_high_down,  algorithm = "parentchild", statistic = "fisher")  #CC
####ks_weight01####
resultKS_mf_ac_low_high_up    <- runTest(GOdata_mf_ac_low_high_up,    algorithm = "weight01", statistic = "ks")  #MF
resultKS_mf_ac_low_high_down  <- runTest(GOdata_mf_ac_low_high_down,  algorithm = "weight01", statistic = "ks")  #MF
resultKS_bp_ac_low_high_up    <- runTest(GOdata_bp_ac_low_high_up,    algorithm = "weight01", statistic = "ks")  #BP
resultKS_bp_ac_low_high_down  <- runTest(GOdata_bp_ac_low_high_down,  algorithm = "weight01", statistic = "ks")  #BP
resultKS_cc_ac_low_high_up    <- runTest(GOdata_cc_ac_low_high_up,    algorithm = "weight01", statistic = "ks")  #CC
resultKS_cc_ac_low_high_down  <- runTest(GOdata_cc_ac_low_high_down,  algorithm = "weight01", statistic = "ks")  #CC
####t_weight01####
resultT_mf_ac_low_high_up    <- runTest(GOdata_mf_ac_low_high_up,    algorithm = "weight01", statistic = "t")  #MF
resultT_mf_ac_low_high_down  <- runTest(GOdata_mf_ac_low_high_down,  algorithm = "weight01", statistic = "t")  #MF
resultT_bp_ac_low_high_up    <- runTest(GOdata_bp_ac_low_high_up,    algorithm = "weight01", statistic = "t")  #BP
resultT_bp_ac_low_high_down  <- runTest(GOdata_bp_ac_low_high_down,  algorithm = "weight01", statistic = "t")  #BP
resultT_cc_ac_low_high_up    <- runTest(GOdata_cc_ac_low_high_up,    algorithm = "weight01", statistic = "t")  #CC
resultT_cc_ac_low_high_down  <- runTest(GOdata_cc_ac_low_high_down,  algorithm = "weight01", statistic = "t")  #CC
esultT_mf_wvc_mid_high_up   <- runTest(GOdata_mf_wvc_mid_high_up,   algorithm = "weight01", statistic = "t") #MF
####RESULTS####
allRes_mf_ac_low_high_up    <- GenTable(GOdata_mf_ac_low_high_up,    ParentChild_Fisher = resultFisher_mf_ac_low_high_up,    Weight01_T = resultT_mf_ac_low_high_up,    Weight01_KS = resultKS_mf_ac_low_high_up,   topNodes = 30)
allRes_mf_ac_low_high_down  <- GenTable(GOdata_mf_ac_low_high_down,  ParentChild_Fisher = resultFisher_mf_ac_low_high_down,  Weight01_T = resultT_mf_ac_low_high_down,  Weight01_KS = resultKS_mf_ac_low_high_down, topNodes = 30)
allRes_bp_ac_low_high_up    <- GenTable(GOdata_bp_ac_low_high_up,    ParentChild_Fisher = resultFisher_bp_ac_low_high_up,    Weight01_T = resultT_bp_ac_low_high_up,    Weight01_KS = resultKS_bp_ac_low_high_up,   topNodes = 30)
allRes_bp_ac_low_high_down  <- GenTable(GOdata_bp_ac_low_high_down,  ParentChild_Fisher = resultFisher_bp_ac_low_high_down,  Weight01_T = resultT_bp_ac_low_high_down,  Weight01_KS = resultKS_bp_ac_low_high_down, topNodes = 30)
allRes_cc_ac_low_high_up    <- GenTable(GOdata_cc_ac_low_high_up,    ParentChild_Fisher = resultFisher_cc_ac_low_high_up,    Weight01_KS = resultKS_cc_ac_low_high_up,   topNodes = 30)
allRes_cc_ac_low_high_down  <- GenTable(GOdata_cc_ac_low_high_down,  ParentChild_Fisher = resultFisher_cc_ac_low_high_down,  Weight01_KS = resultKS_cc_ac_low_high_down, topNodes = 30)

allRes_bp_ac_low_high_up$cat  <- "BP"; 
allRes_bp_ac_low_high_down$cat<- "BP"; 
allRes_cc_ac_low_high_up$cat  <- "CC"; 
allRes_cc_ac_low_high_down$cat<- "CC"; 
allRes_mf_ac_low_high_up$cat  <- "MF"; 
allRes_mf_ac_low_high_down$cat<- "MF"; 
#allRes_cc_ac_low_high_up$Weight01_T=NaN
#allRes_cc_ac_low_high_up$`Rank in Weight01_T`=NaN
#allRes_cc_ac_low_high_up$`Rank in Weight01_KS`=NULL
allRes_ac_low_high_up    <- data.frame(rbind(allRes_bp_ac_low_high_up,    allRes_cc_ac_low_high_up,    allRes_mf_ac_low_high_up))
#allRes_cc_ac_low_high_down$Weight01_T=NaN
#allRes_cc_ac_low_high_down$`Rank in Weight01_T`=NaN
#allRes_cc_ac_low_high_down$`Rank in Weight01_KS`=NULL
allRes_ac_low_high_down  <- data.frame(rbind(allRes_bp_ac_low_high_down,  allRes_cc_ac_low_high_down,  allRes_mf_ac_low_high_down))
####Lista de resultados
myResultsList_ac_low_high_up  <- list(  "Res" = allRes_ac_low_high_up,  "BP" = GOdata_bp_ac_low_high_up,    "BP_T" = resultT_bp_ac_low_high_up,    "BP_KS" = resultKS_bp_ac_low_high_up,    "BP_Fisher" = resultFisher_bp_ac_low_high_up,"CC" = GOdata_cc_ac_low_high_up,    "CC_T" = NULL,    "CC_KS" = resultKS_cc_ac_low_high_up,    "CC_Fisher" = resultFisher_cc_ac_low_high_up, "MF" = GOdata_mf_ac_low_high_up,    "MF_T" = resultT_mf_ac_low_high_up,    "MF_KS" = resultKS_mf_ac_low_high_up,    "MF_Fisher" = resultFisher_mf_ac_low_high_up)
myResultsList_ac_low_high_down  <- list("Res" = allRes_ac_low_high_down,  "BP" = GOdata_bp_ac_low_high_down,  "BP_T" = resultT_bp_ac_low_high_down,  "BP_KS" = resultKS_bp_ac_low_high_down,  "BP_Fisher" = resultFisher_bp_ac_low_high_down,  "CC" = GOdata_cc_ac_low_high_down,  "CC_T" = NULL,  "CC_KS" = resultKS_cc_ac_low_high_down,  "CC_Fisher" = resultFisher_cc_ac_low_high_down,  "MF" = GOdata_mf_ac_low_high_down,  "MF_T" = resultT_mf_ac_low_high_down,  "MF_KS" = resultKS_mf_ac_low_high_down,  "MF_Fisher" = resultFisher_mf_ac_low_high_down)
write.table(myResultsList_ac_low_high_up$Res, file = "sym_myResultsList_high_low_up",sep='\t')
write.table(myResultsList_ac_low_high_down$Res, file = "sym_myResultsList_high_low_down",sep='\t')
#####Graph structure####
showSigOfNodes(GOdata_bp_ac_low_high_up, score(resultFisher_bp_ac_low_high_up), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_bp_ac_low_high_down, score(resultFisher_bp_ac_low_high_down), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_cc_ac_low_high_up, score(resultFisher_cc_ac_low_high_up), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_cc_ac_low_high_down, score(resultFisher_cc_ac_low_high_down), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_mf_ac_low_high_up, score(resultFisher_mf_ac_low_high_up), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_mf_ac_low_high_down, score(resultFisher_mf_ac_low_high_down), firstSigNodes = 10, useInfo = 'all')

