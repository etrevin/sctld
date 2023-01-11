####invierno ventila 13 vs invierno requeson 7 - sitios_1
####vent ver(jun13) 13 vs vent (inene18) 13- temporada_2
####vent ver (low, mid y high) -ventila ver-condt 3
####vent in (low, mid y high) - ventil inv-cond 4
####vent (low, mid y high) -cond 5 _ Tesis
####req (nino) vs req (ver) - cond 6_ otro
#retirar 161, 163, v9, 560, 173 Y V12

####total_39 libs####
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
cond = factor(c(
  "sctld","ctrl","sctld","ctrl","sctld",
  "sctld","ctrl","sctld","ctrl","sctld",
  "ctrl","ctrl","sctld","ctrl","ctrl",
  "sctld","ctrl","sctld","ctrl","sctld",
  "ctrl","sctld","ctrl","ctrl","ctrl",
  "sctld"))





cond = factor(c(
  "ctrl","sctld","ctrl","sctld","sctld","ctrl","sctld","ctrl","sctld","ctrl","sctld","sctld","sctld","ctrl","sctld","sctld"))
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
#ddsmf$temp = relevel(ddsmf$temp, "win")
#ddsmf$cond = relevel(ddsmf$cond, "low")
####DESEQ####
#ddsmf$group <- factor(paste0(ddsmf$temp, ddsmf$cond))
#design(ddsmf) <- ~ group
ddsmf <- DESeq(ddsmf)
#res = results(ddsmf, alpha = 0.05)
resultsNames(ddsmf)
#resmf<-results(ddsmf, alpha = 0.05, lfcThreshold = 2)#Crea un objeto con los resultados de la expresi�n diferencial, un dataFrame con 6 colums: "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"; Se puede modificar el argumento "alpha" (e.g. alpha=0.05)
#summary(resmf)

rld <-rlog(ddsmf)                                                                  #Transformaci�n de los DGEs a Log2
#vsd <-vst(ddsmf)                                                                   #Transformaci�n de los DGEs a Log2
#ntd <-normTransform(ddsmf)    
resmf_allcond_low_high <- results(ddsmf, filterFun=ihw, alpha = 0.05, contrast = c("cond", "sctld", "ctrl"))
resmf_allcond_low_high <- results(ddsmf, contrast = c("cond", "sctld", "ctrl"))
####Summery####
summary(resmf_allcond_low_high)
#summary(resmf_vent_temp)
#Trinotate data for expre table
#sym_trinotate_annotation_report <- read.delim("~/Documents/DGE/symkal/sym_trinotate_annotation_report.xls")
#sym_trino2=distinct(sym_trinotate_annotation_report,X.gene_id,.keep_all = TRUE)         
#row.names(sym_trino2)<-sym_trino2$X.gene_id

#sym_all_low_high = as.data.frame(resmf_allcond_low_high)
#write.table(sym_all_low_high, file = "sym_all_low_high",sep='\t')

resOrdered <- resmf_allcond_low_high[order(resmf_allcond_low_high$padj),]

#resSig_low_high <- subset(resOrdered, padj < 0.05)
#write.table(resSig_low_high, file = "symSig_low_high",sep='\t')

select_genes_low_high<-rownames(subset(resOrdered, padj < 0.05 ))
mat_low_high  <- assay(rld)[ select_genes_low_high, ]
#select_genes_low_high_up<-rownames(subset(resOrdered, padj < 0.05 & log2FoldChange > 0))
#select_genes_low_high_down<-rownames(subset(resOrdered, padj < 0.05 & log2FoldChange < 0))
#mat_low_high_up  <- assay(rld)[ select_genes_low_high_up, ]
#mat_low_high_down  <- assay(rld)[ select_genes_low_high_down, ]


#trino_low_high  <- sym_trino2[ select_genes_low_high, ]
#write.table(trino_low_high, file = "sym_trino_low_high",sep='\t')

#trino_low_high_up  <- por_trino2[ select_genes_low_high_up, ]
#trino_low_high_down  <- por_trino2[ select_genes_low_high_down, ]

#write.table(select_genes_low_high_up, file = "por_select_genes_low_high_up",sep='\t')
#write.table(select_genes_low_high_down, file = "por_select_genes_low_high_down",sep='\t')
#write.table(mat_low_high_up, file = "por_mat_low_high_up",sep='\t')
#write.table(mat_low_high_down, file = "por_mat_low_high_down",sep='\t')
#write.table(trino_low_high_up, file = "por_trino_low_high_up",sep='\t')
#write.table(trino_low_high_down, file = "por_trino_low_high_down",sep='\t')




#sym_all_win_sum = as.data.frame(resmf_vent_temp)
#write.table(sym_all_win_sum, file = "sym_all_win_sum",sep='\t')

#resOrdered <- resmf_vent_temp[order(resmf_vent_temp$padj),]

#resSig_win_sum <- subset(resOrdered, padj < 0.05)
#write.table(resSig_win_sum, file = "sym_Sig_win_sum",sep='\t')

#select_genes_win_sum<-rownames(subset(resOrdered, padj < 0.05 ))
#mat_win_sum  <- assay(rld)[ select_genes_win_sum, ]

#select_genes_win_sum_up<-rownames(subset(resOrdered, padj < 0.05 & log2FoldChange > 0))
#select_genes_win_sum_down<-rownames(subset(resOrdered, padj < 0.05 & log2FoldChange < 0))
#mat_win_sum_up  <- assay(rld)[ select_genes_win_sum_up, ]
#mat_win_sum_down  <- assay(rld)[ select_genes_win_sum_down, ]

#trino_win_sum  <- sym_trino2[ select_genes_win_sum, ]
#write.table(trino_win_sum, file = "sym_trino_win_sum",sep='\t')

#trino_win_sum_up<- por_trino2[ select_genes_win_sum_up, ]
#trino_win_sum_down  <- por_trino2[ select_genes_win_sum_down, ]


#dbset=rbind(mat_low_high)
#dbset=dbset[!duplicated(dbset), ]
#tempcond<- interaction(temp,  cond)


####PCoA####
#dbset=rbind(mat_low_high)

#tempcond<- interaction(temp,  cond)

#dbdist=dist(t(mat_low_high))
#mds=cmdscale(dbdist, eig = TRUE)
#ggordiplots::gg_ordiplot(mds, groups = cond, spiders = TRUE )


####PCA####
#View(dbset)
#new_df <- dbset[ order(row.names(dbset)), ]
dbset[duplicated(dbset), ]
dbset[!duplicated(dbset), ]
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
#fviz_eig(res.x, addlabels = TRUE, ylim = c(0, 50))
# Contributions of variables to PC1
#var <- get_pca_var(res.x)
#fviz_contrib(res.x, choice = "var", axes = 1, top = 10)
# Contributions of variables to PC2
#fviz_contrib(res.x, choice = "var", axes = 2, top = 10)
# Contributions of variables to PC3
#fviz_contrib(res.x, choice = "var", axes = 3, top = 10)


#df=as.data.frame(res.x$ind$coord)[,c(1,2)]
#df
#df$tempcond=c("vc","vc","ic","vv","iv","iv","ic","vv")
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

####Results####
resultsNames(ddsmf)




#write.table(test, file='test.tsv', quote=FALSE, sep='\t', col.names = NA)
write.table(select_genes_win_sum_up, file = "por_select_genes_win_sum_up",sep='\t')
write.table(select_genes_win_sum_down, file = "por_select_genes_win_sum_down",sep='\t')
write.table(mat_win_sum_up, file = "por_mat_win_sum_up",sep='\t')
write.table(mat_win_sum_down, file = "por_mat_win_sum_down",sep='\t')
write.table(trino_win_sum_up, file = "por_trino_win_sum_up",sep='\t')
write.table(trino_win_sum_down, file = "por_trino_win_sum_down",sep='\t')





summary(resmf_allcond_low_high)
resLFC <- lfcShrink(ddsmf, coef="cond_sctld_vs_ctrl", type="apeglm")
plotMA(resLFC, ylim = c(-5, 5))



EnhancedVolcano(resLFC,
                lab = rownames(resLFC),
                x = 'log2FoldChange',
                y = 'pvalue')




resmf_allcond_low_mid <- results(ddsmf, filterFun=ihw, alpha = 0.05, contrast = c("cond", "low", "mid"))#Ponderaci�n de hip�tesis independiente; Se puede modificar el argumento "alpha" (e.g. alpha=0.05)
resmf_allcond_low_high <- results(ddsmf, filterFun=ihw, lfcThreshold = 2, alpha = 0.05, contrast = c("cond", "low", "high"))    #Ponderaci�n de hip�tesis independiente; Se puede modificar el argumento "alpha" (e.g. alpha=0.05)
resmf_allcond_mid_high <- results(ddsmf, filterFun=ihw, alpha = 0.05, contrast = c("cond", "mid", "high"))      
summary(resmf_allcond_low_mid)
summary(resmf_allcond_low_high)
summary(resmf_allcond_mid_high)

resmf_group_lowwin_lowsum<-results(ddsmf, alpha = 0.05, contrast=c("group","winlow","sumlow"))                                                                             #Crea un objeto con los resultados de la expresi�n diferencial, un dataFrame con 6 colums: "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"; Se puede modificar el argumento "alpha" (e.g. alpha=0.05)
resmf_group_midwin_midsum<-results(ddsmf, alpha = 0.05, contrast=c("group","winmid","summid"))                                                                             #Crea un objeto con los resultados de la expresi�n diferencial, un dataFrame con 6 colums: "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj"; Se puede modificar el argumento "alpha" (e.g. alpha=0.05)
resmf_group_highwin_highsum<-results(ddsmf, alpha = 0.05, contrast=c("group","winhigh","sumhigh"))  
summary(resmf_group_lowwin_lowsum)
summary(resmf_group_midwin_midsum)
summary(resmf_group_highwin_highsum)

resmf_wintvent_cond_low_mid<-results(ddsmf, alpha = 0.05, contrast=c("group","winlow","winmid"))  
resmf_wintvent_cond_low_high<-results(ddsmf, alpha = 0.05, contrast=c("group","winlow","winhigh"))  
resmf_wintvent_cond_mid_high<-results(ddsmf, alpha = 0.05, contrast=c("group","winmid","winhigh"))  
summary(resmf_wintvent_cond_low_mid)
summary(resmf_wintvent_cond_low_high)
summary(resmf_wintvent_cond_mid_high)

resmf_sumvent_cond_low_mid<-results(ddsmf, alpha = 0.05, contrast=c("group","sumlow","summid"))  
resmf_sumvent_cond_low_high<-results(ddsmf, alpha = 0.05, contrast=c("group","sumlow","sumhigh"))  
resmf_sumvent_cond_mid_high<-results(ddsmf, alpha = 0.05, contrast=c("group","summid","sumhigh")) 
summary(resmf_sumvent_cond_low_mid)
summary(resmf_sumvent_cond_low_high)
summary(resmf_sumvent_cond_mid_high)


####up down separation####
resmf_vent_temp_up                <- resmf_vent_temp              [resmf_vent_temp$log2FoldChange >0,]
resmf_vent_temp_down              <- resmf_vent_temp              [resmf_vent_temp$log2FoldChange <0,]

resmf_allcond_low_high_up         <- resmf_allcond_low_high       [resmf_allcond_low_high$log2FoldChange >0,]
resmf_allcond_low_high_down       <- resmf_allcond_low_high       [resmf_allcond_low_high$log2FoldChange <0,]
#resmf_allcond_low_mid_up          <- resmf_allcond_low_mid        [resmf_allcond_low_mid$log2FoldChange >0,]
#resmf_allcond_low_mid_down        <- resmf_allcond_low_mid        [resmf_allcond_low_mid$log2FoldChange <0,]
#resmf_allcond_mid_high_up         <- resmf_allcond_mid_high       [resmf_allcond_mid_high$log2FoldChange >0,]
#resmf_allcond_mid_high_down       <- resmf_allcond_mid_high       [resmf_allcond_mid_high$log2FoldChange <0,]

#resmf_condseas_lowwin_lowsum_up   <- resmf_group_lowwin_lowsum   [resmf_group_lowwin_lowsum$log2FoldChange >0,]
#resmf_condseas_lowwin_lowsum_down <- resmf_group_lowwin_lowsum   [resmf_group_lowwin_lowsum$log2FoldChange <0,]
#resmf_condseas_midwin_midsum_up   <- resmf_group_midwin_midsum   [resmf_group_midwin_midsum$log2FoldChange >0,]
#resmf_condseas_midwin_midsum_down <- resmf_group_midwin_midsum   [resmf_group_midwin_midsum$log2FoldChange <0,]
#resmf_condseas_highwin_highsum_up   <- resmf_group_highwin_highsum[resmf_group_highwin_highsum$log2FoldChange >0,]
#resmf_condseas_highwin_highsum_down <- resmf_group_highwin_highsum[resmf_group_highwin_highsum$log2FoldChange <0,]

#resmf_wintvent_cond_low_high_up   <- resmf_wintvent_cond_low_high [resmf_wintvent_cond_low_high$log2FoldChange >0,]
#resmf_wintvent_cond_low_high_down <- resmf_wintvent_cond_low_high [resmf_wintvent_cond_low_high$log2FoldChange <0,] 
#resmf_wintvent_cond_low_mid_up    <- resmf_wintvent_cond_low_mid  [resmf_wintvent_cond_low_mid$log2FoldChange  >0,]
#resmf_wintvent_cond_low_mid_down  <- resmf_wintvent_cond_low_mid  [resmf_wintvent_cond_low_mid$log2FoldChange  <0,]
#resmf_wintvent_cond_mid_high_up   <- resmf_wintvent_cond_mid_high [resmf_wintvent_cond_mid_high$log2FoldChange >0,]
#resmf_wintvent_cond_mid_high_down <- resmf_wintvent_cond_mid_high [resmf_wintvent_cond_mid_high$log2FoldChange <0,]

#resmf_sumvent_cond_low_high_up    <- resmf_sumvent_cond_low_high  [resmf_sumvent_cond_low_high$log2FoldChange >0,]
#resmf_sumvent_cond_low_high_down  <- resmf_sumvent_cond_low_high  [resmf_sumvent_cond_low_high$log2FoldChange <0,]
#resmf_sumvent_cond_low_mid_up     <- resmf_sumvent_cond_low_mid   [resmf_sumvent_cond_low_mid$log2FoldChange  >0,]
#resmf_sumvent_cond_low_mid_down   <- resmf_sumvent_cond_low_mid   [resmf_sumvent_cond_low_mid$log2FoldChange  <0,]
#resmf_sumvent_cond_mid_high_up    <- resmf_sumvent_cond_mid_high  [resmf_sumvent_cond_mid_high$log2FoldChange >0,]
#resmf_sumvent_cond_mid_high_down  <- resmf_sumvent_cond_mid_high  [resmf_sumvent_cond_mid_high$log2FoldChange <0,]

####p value####
x_vt_up             <- resmf_vent_temp_up$padj     # P-ajustado verano e invierno dentro de ventila
x_vt_down           <- resmf_vent_temp_down$padj     # P-ajustado verano e invierno dentro de ventila

x_ac_low_high_up    <- resmf_allcond_low_high_up$padj       # P-ajustado Condiciones (verano e invierno) en ventrila
x_ac_low_high_down  <- resmf_allcond_low_high_down$padj       # P-ajustado Condiciones (verano e invierno) en ventrila
#x_ac_low_mid_up     <- resmf_allcond_low_mid_up$padj       # P-ajustado Condiciones (verano e invierno) en ventrila
#x_ac_low_mid_down   <- resmf_allcond_low_mid_down$padj       # P-ajustado Condiciones (verano e invierno) en ventrila
#x_ac_mid_high_up    <- resmf_allcond_mid_high_up$padj       # P-ajustado Condiciones (verano e invierno) en ventrila
#x_ac_mid_high_down  <- resmf_allcond_mid_high_down$padj       # P-ajustado Condiciones (verano e invierno) en ventrila

#x_cs_lowwin_lowsum_up      <-resmf_condseas_lowwin_lowsum_up$padj
#x_cs_lowwin_lowsum_down    <-resmf_condseas_lowwin_lowsum_down$padj
#x_cs_midwin_midsum_up      <-resmf_condseas_midwin_midsum_up$padj
#x_cs_midwin_midsum_down    <-resmf_condseas_midwin_midsum_down$padj
#x_cs_highwin_highsum_up      <-resmf_condseas_highwin_highsum_up$padj
#x_cs_highwin_highsum_down    <-resmf_condseas_highwin_highsum_down$padj

#x_svc_low_high_up   <- resmf_sumvent_cond_low_high_up$padj  # P-ajustado Condiciones (verano) ventrila
#x_svc_low_high_down <- resmf_sumvent_cond_low_high_down$padj  # P-ajustado Condiciones (verano) ventrila
#x_svc_low_mid_up    <- resmf_sumvent_cond_low_mid_up$padj  # P-ajustado Condiciones (verano) ventrila
#x_svc_low_mid_down  <- resmf_sumvent_cond_low_mid_down$padj  # P-ajustado Condiciones (verano) ventrila
#x_svc_mid_high_up   <- resmf_sumvent_cond_mid_high_up$padj  # P-ajustado Condiciones (verano) ventrila
#x_svc_mid_high_down <- resmf_sumvent_cond_mid_high_down$padj  # P-ajustado Condiciones (verano) ventrila

#x_wvc_low_high_up   <- resmf_wintvent_cond_low_high_up$padj # P-ajustado Condiciones (invierno) ventrila
#x_wvc_low_high_down <- resmf_wintvent_cond_low_high_down$padj # P-ajustado Condiciones (invierno) ventrila
#x_wvc_low_mid_up    <- resmf_wintvent_cond_low_mid_up$padj # P-ajustado Condiciones (invierno) ventrila
#x_wvc_low_mid_down  <- resmf_wintvent_cond_low_mid_down$padj # P-ajustado Condiciones (invierno) ventrila
#x_wvc_mid_high_up   <- resmf_wintvent_cond_mid_high_up$padj # P-ajustado Condiciones (invierno) ventrila
#x_wvc_mid_high_down <- resmf_wintvent_cond_mid_high_down$padj # P-ajustado Condiciones (invierno) ventrila

#Crea un objeto con los nombres de los transcritos provenientes de DESeq2

names(x_vt_up)             <- row.names(resmf_vent_temp_up)       # Nombres verano e invierno dentro de ventila
names(x_vt_down)           <- row.names(resmf_vent_temp_down)       # Nombres verano e invierno dentro de ventila

names(x_ac_low_high_up)    <- row.names(resmf_allcond_low_high_up)         # Nombres Condiciones (verano e invierno) en ventrila
names(x_ac_low_high_down)  <- row.names(resmf_allcond_low_high_down)         # Nombres Condiciones (verano e invierno) en ventrila
#names(x_ac_low_mid_up)     <- row.names(resmf_allcond_low_mid_up)         # Nombres Condiciones (verano e invierno) en ventrila
#names(x_ac_low_mid_down)   <- row.names(resmf_allcond_low_mid_down)         # Nombres Condiciones (verano e invierno) en ventrila
#names(x_ac_mid_high_up)    <- row.names(resmf_allcond_mid_high_up)         # Nombres Condiciones (verano e invierno) en ventrila
#names(x_ac_mid_high_down)  <- row.names(resmf_allcond_mid_high_down)         # Nombres Condiciones (verano e invierno) en ventrila

#names(x_cs_lowwin_lowsum_up)   <-row.names(resmf_condseas_lowwin_lowsum_up)
#names(x_cs_lowwin_lowsum_down) <-row.names(resmf_condseas_lowwin_lowsum_down)
#names(x_cs_midwin_midsum_up)   <-row.names(resmf_condseas_midwin_midsum_up)
#names(x_cs_midwin_midsum_down) <-row.names(resmf_condseas_midwin_midsum_down)
#names(x_cs_highwin_highsum_up)   <-row.names(resmf_condseas_highwin_highsum_up)
#names(x_cs_highwin_highsum_down) <-row.names(resmf_condseas_highwin_highsum_down)

#names(x_svc_low_high_up)   <- row.names(resmf_sumvent_cond_low_high_up)    # Nombres Condiciones (verano) ventrila
#names(x_svc_low_high_down) <- row.names(resmf_sumvent_cond_low_high_down)    # Nombres Condiciones (verano) ventrila
#names(x_svc_low_mid_up)    <- row.names(resmf_sumvent_cond_low_mid_up)    # Nombres Condiciones (verano) ventrila
#names(x_svc_low_mid_down)  <- row.names(resmf_sumvent_cond_low_mid_down)    # Nombres Condiciones (verano) ventrila
#names(x_svc_mid_high_up)   <- row.names(resmf_sumvent_cond_mid_high_up)    # Nombres Condiciones (verano) ventrila
#names(x_svc_mid_high_down) <- row.names(resmf_sumvent_cond_mid_high_down)    # Nombres Condiciones (verano) ventrila

#names(x_wvc_low_high_up)   <- row.names(resmf_wintvent_cond_low_high_up)   # Nombres Condiciones (invierno) ventrila
#names(x_wvc_low_high_down) <- row.names(resmf_wintvent_cond_low_high_down)   # Nombres Condiciones (invierno) ventrila
#names(x_wvc_low_mid_up)    <- row.names(resmf_wintvent_cond_low_mid_up)   # Nombres Condiciones (invierno) ventrila
#names(x_wvc_low_mid_down)  <- row.names(resmf_wintvent_cond_low_mid_down)   # Nombres Condiciones (invierno) ventrila
#names(x_wvc_mid_high_up)   <- row.names(resmf_wintvent_cond_mid_high_up)   # Nombres Condiciones (invierno) ventrila
#names(x_wvc_mid_high_down) <- row.names(resmf_wintvent_cond_mid_high_down)   # Nombres Condiciones (invierno) ventrila

####cleanNA####

x_vt_up    <- na.omit(x_vt_up)      
x_vt_down  <- na.omit(x_vt_down)  
x_ac_low_high_up    <- na.omit(x_ac_low_high_up)        
x_ac_low_high_down  <- na.omit(x_ac_low_high_down)        
#x_ac_low_mid_up     <- na.omit(x_ac_low_mid_up)        
#x_ac_low_mid_down   <- na.omit(x_ac_low_mid_down)        
#x_ac_mid_high_up    <- na.omit(x_ac_mid_high_up)       
#x_ac_mid_high_down  <- na.omit(x_ac_mid_high_down)       

#x_cs_lowwin_lowsum_up     <- na.omit(x_cs_lowwin_lowsum_up)
#x_cs_lowwin_lowsum_down   <- na.omit(x_cs_lowwin_lowsum_down)
#x_cs_midwin_midsum_up     <- na.omit(x_cs_midwin_midsum_up)
#x_cs_midwin_midsum_down   <- na.omit(x_cs_midwin_midsum_down)
#x_cs_highwin_highsum_up   <- na.omit(x_cs_highwin_highsum_up)
#x_cs_highwin_highsum_down <- na.omit(x_cs_highwin_highsum_down)

#x_svc_low_high_up   <- na.omit(x_svc_low_high_up)  
#x_svc_low_high_down <- na.omit(x_svc_low_high_down)  
#x_svc_low_mid_up    <- na.omit(x_svc_low_mid_up)  
#x_svc_low_mid_down  <- na.omit(x_svc_low_mid_down)  
#x_svc_mid_high_up   <- na.omit(x_svc_mid_high_up)  
#x_svc_mid_high_down <- na.omit(x_svc_mid_high_down)  

#x_wvc_low_high_up   <- na.omit(x_wvc_low_high_up) 
#x_wvc_low_high_down <- na.omit(x_wvc_low_high_down) 
#x_wvc_low_mid_up    <- na.omit(x_wvc_low_mid_up) 
#x_wvc_low_mid_down  <- na.omit(x_wvc_low_mid_down) 
#x_wvc_mid_high_up   <- na.omit(x_wvc_mid_high_up)
#x_wvc_mid_high_down <- na.omit(x_wvc_mid_high_down)

####goterms####
geneID2GO <- readMappings("/home/erick/Documents/workspace/trinotate/go_mch.txt") #carga las categorias GO de trinotate
str(head(geneID2GO)) 

####Funci�n para separar los genes en TOPGO####
topDiffGenes <- function(allScore) { return(allScore < 0.05)}

####TopGO_Object####
#MF_Molecular function
GOdata_mf_vt_up             <- new("topGOdata", ontology = "MF", allGenes = x_vt_up,             geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_mf_vt_down           <- new("topGOdata", ontology = "MF", allGenes = x_vt_down,           geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)

GOdata_mf_ac_low_high_up    <- new("topGOdata", ontology = "MF", allGenes = x_ac_low_high_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_mf_ac_low_high_down  <- new("topGOdata", ontology = "MF", allGenes = x_ac_low_high_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
#GOdata_mf_ac_low_mid_up     <- new("topGOdata", ontology = "MF", allGenes = x_ac_low_mid_up,     geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_ac_low_mid_down   <- new("topGOdata", ontology = "MF", allGenes = x_ac_low_mid_down,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_ac_mid_high_up    <- new("topGOdata", ontology = "MF", allGenes = x_ac_mid_high_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_ac_mid_high_down  <- new("topGOdata", ontology = "MF", allGenes = x_ac_mid_high_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#GOdata_mf_cs_lowwin_lowsum_up     <- new("topGOdata", ontology = "MF", allGenes = x_cs_lowwin_lowsum_up,      geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_cs_lowwin_lowsum_down   <- new("topGOdata", ontology = "MF", allGenes = x_cs_lowwin_lowsum_down,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_cs_midwin_midsum_up     <- new("topGOdata", ontology = "MF", allGenes = x_cs_midwin_midsum_up,      geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_cs_midwin_midsum_down   <- new("topGOdata", ontology = "MF", allGenes = x_cs_midwin_midsum_down,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_cs_highwin_highsum_up   <- new("topGOdata", ontology = "MF", allGenes = x_cs_highwin_highsum_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_cs_highwin_highsum_down <- new("topGOdata", ontology = "MF", allGenes = x_cs_highwin_highsum_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#GOdata_mf_svc_low_high_up   <- new("topGOdata", ontology = "MF", allGenes = x_svc_low_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_svc_low_high_down <- new("topGOdata", ontology = "MF", allGenes = x_svc_low_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_svc_low_mid_up    <- new("topGOdata", ontology = "MF", allGenes = x_svc_low_mid_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_svc_low_mid_down  <- new("topGOdata", ontology = "MF", allGenes = x_svc_low_mid_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_svc_mid_high_up   <- new("topGOdata", ontology = "MF", allGenes = x_svc_mid_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_svc_mid_high_down <- new("topGOdata", ontology = "MF", allGenes = x_svc_mid_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#GOdata_mf_wvc_low_high_up   <- new("topGOdata", ontology = "MF", allGenes = x_wvc_low_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_wvc_low_high_down <- new("topGOdata", ontology = "MF", allGenes = x_wvc_low_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_wvc_low_mid_up    <- new("topGOdata", ontology = "MF", allGenes = x_wvc_low_mid_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_wvc_low_mid_down  <- new("topGOdata", ontology = "MF", allGenes = x_wvc_low_mid_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_wvc_mid_high_up   <- new("topGOdata", ontology = "MF", allGenes = x_wvc_mid_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_mf_wvc_mid_high_down <- new("topGOdata", ontology = "MF", allGenes = x_wvc_mid_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#BP_BIOLOGICAL PROCESS
GOdata_bp_vt_up             <- new("topGOdata", ontology = "BP", allGenes = x_vt_up,             geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_bp_vt_down           <- new("topGOdata", ontology = "BP", allGenes = x_vt_down,           geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)

GOdata_bp_ac_low_high_up    <- new("topGOdata", ontology = "BP", allGenes = x_ac_low_high_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_bp_ac_low_high_down  <- new("topGOdata", ontology = "BP", allGenes = x_ac_low_high_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
#GOdata_bp_ac_low_mid_up     <- new("topGOdata", ontology = "BP", allGenes = x_ac_low_mid_up,     geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_ac_low_mid_down   <- new("topGOdata", ontology = "BP", allGenes = x_ac_low_mid_down,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_ac_mid_high_up    <- new("topGOdata", ontology = "BP", allGenes = x_ac_mid_high_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_ac_mid_high_down  <- new("topGOdata", ontology = "BP", allGenes = x_ac_mid_high_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#GOdata_bp_cs_lowwin_lowsum_up     <- new("topGOdata", ontology = "BP",   allGenes = x_cs_lowwin_lowsum_up,      geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_cs_lowwin_lowsum_down   <- new("topGOdata", ontology = "BP",   allGenes = x_cs_lowwin_lowsum_down,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_cs_midwin_midsum_up     <- new("topGOdata", ontology = "BP",   allGenes = x_cs_midwin_midsum_up,      geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_cs_midwin_midsum_down   <- new("topGOdata", ontology = "BP",   allGenes = x_cs_midwin_midsum_down,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_cs_highwin_highsum_up   <- new("topGOdata", ontology = "BP",   allGenes = x_cs_highwin_highsum_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_cs_highwin_highsum_down <- new("topGOdata", ontology = "BP",   allGenes = x_cs_highwin_highsum_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#GOdata_bp_svc_low_high_up   <- new("topGOdata", ontology = "BP", allGenes = x_svc_low_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_svc_low_high_down <- new("topGOdata", ontology = "BP", allGenes = x_svc_low_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_svc_low_mid_up    <- new("topGOdata", ontology = "BP", allGenes = x_svc_low_mid_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_svc_low_mid_down  <- new("topGOdata", ontology = "BP", allGenes = x_svc_low_mid_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_svc_mid_high_up   <- new("topGOdata", ontology = "BP", allGenes = x_svc_mid_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_svc_mid_high_down <- new("topGOdata", ontology = "BP", allGenes = x_svc_mid_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#GOdata_bp_wvc_low_high_up   <- new("topGOdata", ontology = "BP", allGenes = x_wvc_low_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_wvc_low_high_down <- new("topGOdata", ontology = "BP", allGenes = x_wvc_low_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_wvc_low_mid_up    <- new("topGOdata", ontology = "BP", allGenes = x_wvc_low_mid_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_wvc_low_mid_down  <- new("topGOdata", ontology = "BP", allGenes = x_wvc_low_mid_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_wvc_mid_high_up   <- new("topGOdata", ontology = "BP", allGenes = x_wvc_mid_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_bp_wvc_mid_high_down <- new("topGOdata", ontology = "BP", allGenes = x_wvc_mid_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#CC_CELULAR COMPONENT
GOdata_cc_vt_up             <- new("topGOdata", ontology = "CC", allGenes = x_vt_up,             geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_cc_vt_down           <- new("topGOdata", ontology = "CC", allGenes = x_vt_down,           geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)

GOdata_cc_ac_low_high_up    <- new("topGOdata", ontology = "CC", allGenes = x_ac_low_high_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
GOdata_cc_ac_low_high_down  <- new("topGOdata", ontology = "CC", allGenes = x_ac_low_high_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO, nodeSize = 10)
#GOdata_cc_ac_low_mid_up     <- new("topGOdata", ontology = "CC", allGenes = x_ac_low_mid_up,     geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_ac_low_mid_down   <- new("topGOdata", ontology = "CC", allGenes = x_ac_low_mid_down,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_ac_mid_high_up    <- new("topGOdata", ontology = "CC", allGenes = x_ac_mid_high_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_ac_mid_high_down  <- new("topGOdata", ontology = "CC", allGenes = x_ac_mid_high_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#GOdata_cc_cs_lowwin_lowsum_up     <- new("topGOdata", ontology = "CC", allGenes = x_cs_lowwin_lowsum_up,      geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_cs_lowwin_lowsum_down   <- new("topGOdata", ontology = "CC", allGenes = x_cs_lowwin_lowsum_down,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_cs_midwin_midsum_up     <- new("topGOdata", ontology = "CC", allGenes = x_cs_midwin_midsum_up,      geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_cs_midwin_midsum_down   <- new("topGOdata", ontology = "CC", allGenes = x_cs_midwin_midsum_down,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_cs_highwin_highsum_up   <- new("topGOdata", ontology = "CC", allGenes = x_cs_highwin_highsum_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_cs_highwin_highsum_down <- new("topGOdata", ontology = "CC", allGenes = x_cs_highwin_highsum_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#GOdata_cc_svc_low_high_up   <- new("topGOdata", ontology = "CC", allGenes = x_svc_low_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_svc_low_high_down <- new("topGOdata", ontology = "CC", allGenes = x_svc_low_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_svc_low_mid_up    <- new("topGOdata", ontology = "CC", allGenes = x_svc_low_mid_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_svc_low_mid_down  <- new("topGOdata", ontology = "CC", allGenes = x_svc_low_mid_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_svc_mid_high_up   <- new("topGOdata", ontology = "CC", allGenes = x_svc_mid_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_svc_mid_high_down <- new("topGOdata", ontology = "CC", allGenes = x_svc_mid_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

#GOdata_cc_wvc_low_high_up   <- new("topGOdata", ontology = "CC", allGenes = x_wvc_low_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_wvc_low_high_down <- new("topGOdata", ontology = "CC", allGenes = x_wvc_low_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_wvc_low_mid_up    <- new("topGOdata", ontology = "CC", allGenes = x_wvc_low_mid_up,    geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_wvc_low_mid_down  <- new("topGOdata", ontology = "CC", allGenes = x_wvc_low_mid_down,  geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_wvc_mid_high_up   <- new("topGOdata", ontology = "CC", allGenes = x_wvc_mid_high_up,   geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)
#GOdata_cc_wvc_mid_high_down <- new("topGOdata", ontology = "CC", allGenes = x_wvc_mid_high_down, geneSel = topDiffGenes, annot =annFUN.gene2GO, gene2GO = geneID2GO)

####Fisher_parentchild####

#resultFisher_mf_vt_up             <- runTest(GOdata_mf_vt_up,             algorithm = "parentchild", statistic = "fisher")  #MF
#resultFisher_mf_vt_down           <- runTest(GOdata_mf_vt_down,           algorithm = "parentchild", statistic = "fisher")  #MF
#resultFisher_bp_vt_up             <- runTest(GOdata_bp_vt_up,             algorithm = "parentchild", statistic = "fisher")  #BP
#resultFisher_bp_vt_down           <- runTest(GOdata_bp_vt_down,           algorithm = "parentchild", statistic = "fisher")  #BP
#resultFisher_cc_vt_up             <- runTest(GOdata_cc_vt_up,             algorithm = "parentchild", statistic = "fisher")  #CC
#resultFisher_cc_vt_down           <- runTest(GOdata_cc_vt_down,           algorithm = "parentchild", statistic = "fisher")  #CC

resultFisher_mf_ac_low_high_up    <- runTest(GOdata_mf_ac_low_high_up,    algorithm = "parentchild", statistic = "fisher")  #MF
resultFisher_mf_ac_low_high_down  <- runTest(GOdata_mf_ac_low_high_down,  algorithm = "parentchild", statistic = "fisher")  #MF
resultFisher_bp_ac_low_high_up    <- runTest(GOdata_bp_ac_low_high_up,    algorithm = "parentchild", statistic = "fisher")  #BP
resultFisher_bp_ac_low_high_down  <- runTest(GOdata_bp_ac_low_high_down,  algorithm = "parentchild", statistic = "fisher")  #BP
resultFisher_cc_ac_low_high_up    <- runTest(GOdata_cc_ac_low_high_up,    algorithm = "parentchild", statistic = "fisher")  #CC
resultFisher_cc_ac_low_high_down  <- runTest(GOdata_cc_ac_low_high_down,  algorithm = "parentchild", statistic = "fisher")  #CC

#resultFisher_mf_ac_low_mid_up     <- runTest(GOdata_mf_ac_low_mid_up,     algorithm = "parentchild", statistic = "fisher")  #MF
#resultFisher_mf_ac_low_mid_down   <- runTest(GOdata_mf_ac_low_mid_down,   algorithm = "parentchild", statistic = "fisher")  #MF
#resultFisher_mf_ac_mid_high_up    <- runTest(GOdata_mf_ac_mid_high_up,    algorithm = "parentchild", statistic = "fisher")  #MF
#resultFisher_mf_ac_mid_high_down  <- runTest(GOdata_mf_ac_mid_high_down,  algorithm = "parentchild", statistic = "fisher")  #MF

#resultFisher_bp_ac_low_mid_up     <- runTest(GOdata_bp_ac_low_mid_up,     algorithm = "parentchild", statistic = "fisher")  #BP
#resultFisher_bp_ac_low_mid_down   <- runTest(GOdata_bp_ac_low_mid_down,   algorithm = "parentchild", statistic = "fisher")  #BP
#resultFisher_bp_ac_mid_high_up    <- runTest(GOdata_bp_ac_mid_high_up,    algorithm = "parentchild", statistic = "fisher")  #BP
#resultFisher_bp_ac_mid_high_down  <- runTest(GOdata_bp_ac_mid_high_down,  algorithm = "parentchild", statistic = "fisher")  #BP

#resultFisher_cc_ac_low_mid_up     <- runTest(GOdata_cc_ac_low_mid_up,     algorithm = "parentchild", statistic = "fisher")  #CC
#resultFisher_cc_ac_low_mid_down   <- runTest(GOdata_cc_ac_low_mid_down,   algorithm = "parentchild", statistic = "fisher")  #CC
#resultFisher_cc_ac_mid_high_up    <- runTest(GOdata_cc_ac_mid_high_up,    algorithm = "parentchild", statistic = "fisher")  #CC
#resultFisher_cc_ac_mid_high_down  <- runTest(GOdata_cc_ac_mid_high_down,  algorithm = "parentchild", statistic = "fisher")  #CC

#resultFisher_mf_cs_lowwin_lowsum_up       <- runTest(GOdata_mf_cs_lowwin_lowsum_up,     algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_cs_lowwin_lowsum_down     <- runTest(GOdata_mf_cs_lowwin_lowsum_down,   algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_cs_midwin_midsum_up       <- runTest(GOdata_mf_cs_midwin_midsum_up,     algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_cs_midwin_midsum_down     <- runTest(GOdata_mf_cs_midwin_midsum_down,   algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_cs_highwin_highsum_up     <- runTest(GOdata_mf_cs_highwin_highsum_up,   algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_cs_highwin_highsum_down   <- runTest(GOdata_mf_cs_highwin_highsum_down, algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_cc_cs_lowwin_lowsum_up       <- runTest(GOdata_cc_cs_lowwin_lowsum_up,     algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_cs_lowwin_lowsum_down     <- runTest(GOdata_cc_cs_lowwin_lowsum_down,   algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_cs_midwin_midsum_up       <- runTest(GOdata_cc_cs_midwin_midsum_up,     algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_cs_midwin_midsum_down     <- runTest(GOdata_cc_cs_midwin_midsum_down,   algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_cs_highwin_highsum_up     <- runTest(GOdata_cc_cs_highwin_highsum_up,   algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_cs_highwin_highsum_down   <- runTest(GOdata_cc_cs_highwin_highsum_down, algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_bp_cs_lowwin_lowsum_up       <- runTest(GOdata_bp_cs_lowwin_lowsum_up,     algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_cs_lowwin_lowsum_down     <- runTest(GOdata_bp_cs_lowwin_lowsum_down,   algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_cs_midwin_midsum_up       <- runTest(GOdata_bp_cs_midwin_midsum_up,     algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_cs_midwin_midsum_down     <- runTest(GOdata_bp_cs_midwin_midsum_down,   algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_cs_highwin_highsum_up     <- runTest(GOdata_bp_cs_highwin_highsum_up,   algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_cs_highwin_highsum_down   <- runTest(GOdata_bp_cs_highwin_highsum_down, algorithm = "parentchild", statistic = "fisher") #BP

#resultFisher_mf_svc_low_high_up   <- runTest(GOdata_mf_svc_low_high_up,   algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_svc_low_high_down <- runTest(GOdata_mf_svc_low_high_down, algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_svc_low_mid_up    <- runTest(GOdata_mf_svc_low_mid_up,    algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_svc_low_mid_down  <- runTest(GOdata_mf_svc_low_mid_down,  algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_svc_mid_high_up   <- runTest(GOdata_mf_svc_mid_high_up,   algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_svc_mid_high_down <- runTest(GOdata_mf_svc_mid_high_down, algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_cc_svc_low_high_up   <- runTest(GOdata_cc_svc_low_high_up,   algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_svc_low_high_down <- runTest(GOdata_cc_svc_low_high_down, algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_svc_low_mid_up    <- runTest(GOdata_cc_svc_low_mid_up,    algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_svc_low_mid_down  <- runTest(GOdata_cc_svc_low_mid_down,  algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_svc_mid_high_up   <- runTest(GOdata_cc_svc_mid_high_up,   algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_svc_mid_high_down <- runTest(GOdata_cc_svc_mid_high_down, algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_bp_svc_low_high_up   <- runTest(GOdata_bp_svc_low_high_up,   algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_svc_low_high_down <- runTest(GOdata_bp_svc_low_high_down, algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_svc_low_mid_up    <- runTest(GOdata_bp_svc_low_mid_up,    algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_svc_low_mid_down  <- runTest(GOdata_bp_svc_low_mid_down,  algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_svc_mid_high_up   <- runTest(GOdata_bp_svc_mid_high_up,   algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_svc_mid_high_down <- runTest(GOdata_bp_svc_mid_high_down, algorithm = "parentchild", statistic = "fisher") #BP

#resultFisher_mf_wvc_low_high_up   <- runTest(GOdata_mf_wvc_low_high_up,   algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_wvc_low_high_down <- runTest(GOdata_mf_wvc_low_high_down, algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_wvc_low_mid_up    <- runTest(GOdata_mf_wvc_low_mid_up,    algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_wvc_low_mid_down  <- runTest(GOdata_mf_wvc_low_mid_down,  algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_wvc_mid_high_up   <- runTest(GOdata_mf_wvc_mid_high_up,   algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_mf_wvc_mid_high_down <- runTest(GOdata_mf_wvc_mid_high_down, algorithm = "parentchild", statistic = "fisher") #MF
#resultFisher_bp_wvc_low_high_up   <- runTest(GOdata_bp_wvc_low_high_up,   algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_wvc_low_high_down <- runTest(GOdata_bp_wvc_low_high_down, algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_wvc_low_mid_up    <- runTest(GOdata_bp_wvc_low_mid_up,    algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_wvc_low_mid_down  <- runTest(GOdata_bp_wvc_low_mid_down,  algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_wvc_mid_high_up   <- runTest(GOdata_bp_wvc_mid_high_up,   algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_bp_wvc_mid_high_down <- runTest(GOdata_bp_wvc_mid_high_down, algorithm = "parentchild", statistic = "fisher") #BP
#resultFisher_cc_wvc_low_high_up   <- runTest(GOdata_cc_wvc_low_high_up,   algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_wvc_low_high_down <- runTest(GOdata_cc_wvc_low_high_down, algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_wvc_low_mid_up    <- runTest(GOdata_cc_wvc_low_mid_up,    algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_wvc_low_mid_down  <- runTest(GOdata_cc_wvc_low_mid_down,  algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_wvc_mid_high_up   <- runTest(GOdata_cc_wvc_mid_high_up,   algorithm = "parentchild", statistic = "fisher") #CC
#resultFisher_cc_wvc_mid_high_down <- runTest(GOdata_cc_wvc_mid_high_down, algorithm = "parentchild", statistic = "fisher") #CC


####ks_weight01####
resultKS_mf_vt_up             <- runTest(GOdata_mf_vt_up,             algorithm = "weight01", statistic = "ks")  #MF
resultKS_mf_vt_down           <- runTest(GOdata_mf_vt_down,           algorithm = "weight01", statistic = "ks")  #MF
resultKS_bp_vt_up             <- runTest(GOdata_bp_vt_up,             algorithm = "weight01", statistic = "ks")  #BP
resultKS_bp_vt_down           <- runTest(GOdata_bp_vt_down,           algorithm = "weight01", statistic = "ks")  #BP
resultKS_cc_vt_up             <- runTest(GOdata_cc_vt_up,             algorithm = "weight01", statistic = "ks")  #CC
resultKS_cc_vt_down           <- runTest(GOdata_cc_vt_down,           algorithm = "weight01", statistic = "ks")  #CC

resultKS_mf_ac_low_high_up    <- runTest(GOdata_mf_ac_low_high_up,    algorithm = "weight01", statistic = "ks")  #MF
resultKS_mf_ac_low_high_down  <- runTest(GOdata_mf_ac_low_high_down,  algorithm = "weight01", statistic = "ks")  #MF
resultKS_bp_ac_low_high_up    <- runTest(GOdata_bp_ac_low_high_up,    algorithm = "weight01", statistic = "ks")  #BP
resultKS_bp_ac_low_high_down  <- runTest(GOdata_bp_ac_low_high_down,  algorithm = "weight01", statistic = "ks")  #BP
resultKS_cc_ac_low_high_up    <- runTest(GOdata_cc_ac_low_high_up,    algorithm = "weight01", statistic = "ks")  #CC
resultKS_cc_ac_low_high_down  <- runTest(GOdata_cc_ac_low_high_down,  algorithm = "weight01", statistic = "ks")  #CC
#resultKS_mf_ac_low_mid_up     <- runTest(GOdata_mf_ac_low_mid_up,     algorithm = "weight01", statistic = "ks")  #MF
#resultKS_mf_ac_low_mid_down   <- runTest(GOdata_mf_ac_low_mid_down,   algorithm = "weight01", statistic = "ks")  #MF
#resultKS_mf_ac_mid_high_up    <- runTest(GOdata_mf_ac_mid_high_up,    algorithm = "weight01", statistic = "ks")  #MF
#resultKS_mf_ac_mid_high_down  <- runTest(GOdata_mf_ac_mid_high_down,  algorithm = "weight01", statistic = "ks")  #MF

#resultKS_bp_ac_low_mid_up     <- runTest(GOdata_bp_ac_low_mid_up,     algorithm = "weight01", statistic = "ks")  #BP
#resultKS_bp_ac_low_mid_down   <- runTest(GOdata_bp_ac_low_mid_down,   algorithm = "weight01", statistic = "ks")  #BP
#resultKS_bp_ac_mid_high_up    <- runTest(GOdata_bp_ac_mid_high_up,    algorithm = "weight01", statistic = "ks")  #BP
#resultKS_bp_ac_mid_high_down  <- runTest(GOdata_bp_ac_mid_high_down,  algorithm = "weight01", statistic = "ks")  #BP

#resultKS_cc_ac_low_mid_up     <- runTest(GOdata_cc_ac_low_mid_up,     algorithm = "weight01", statistic = "ks")  #CC
#resultKS_cc_ac_low_mid_down   <- runTest(GOdata_cc_ac_low_mid_down,   algorithm = "weight01", statistic = "ks")  #CC
#resultKS_cc_ac_mid_high_up    <- runTest(GOdata_cc_ac_mid_high_up,    algorithm = "weight01", statistic = "ks")  #CC
#resultKS_cc_ac_mid_high_down  <- runTest(GOdata_cc_ac_mid_high_down,  algorithm = "weight01", statistic = "ks")  #CC

#resultKS_mf_cs_lowwin_lowsum_up     <- runTest(GOdata_mf_cs_lowwin_lowsum_up,       algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_cs_lowwin_lowsum_down   <- runTest(GOdata_mf_cs_lowwin_lowsum_down,     algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_cs_midwin_midsum_up     <- runTest(GOdata_mf_cs_midwin_midsum_up,       algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_cs_midwin_midsum_down   <- runTest(GOdata_mf_cs_midwin_midsum_down,     algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_cs_highwin_highsum_up   <- runTest(GOdata_mf_cs_highwin_highsum_up,     algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_cs_highwin_highsum_down <- runTest(GOdata_mf_cs_highwin_highsum_down,   algorithm = "weight01", statistic = "ks") #MF
#resultKS_cc_cs_lowwin_lowsum_up     <- runTest(GOdata_cc_cs_lowwin_lowsum_up,       algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_cs_lowwin_lowsum_down   <- runTest(GOdata_cc_cs_lowwin_lowsum_down,     algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_cs_midwin_midsum_up     <- runTest(GOdata_cc_cs_midwin_midsum_up,       algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_cs_midwin_midsum_down   <- runTest(GOdata_cc_cs_midwin_midsum_down,     algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_cs_highwin_highsum_up   <- runTest(GOdata_cc_cs_highwin_highsum_up,     algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_cs_highwin_highsum_down <- runTest(GOdata_cc_cs_highwin_highsum_down,   algorithm = "weight01", statistic = "ks") #CC
#resultKS_bp_cs_lowwin_lowsum_up     <- runTest(GOdata_bp_cs_lowwin_lowsum_up,       algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_cs_lowwin_lowsum_down   <- runTest(GOdata_bp_cs_lowwin_lowsum_down,     algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_cs_midwin_midsum_up     <- runTest(GOdata_bp_cs_midwin_midsum_up,       algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_cs_midwin_midsum_down   <- runTest(GOdata_bp_cs_midwin_midsum_down,     algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_cs_highwin_highsum_up   <- runTest(GOdata_bp_cs_highwin_highsum_up,     algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_cs_highwin_highsum_down <- runTest(GOdata_bp_cs_highwin_highsum_down,   algorithm = "weight01", statistic = "ks") #BP

#resultKS_mf_svc_low_high_up   <- runTest(GOdata_mf_svc_low_high_up,   algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_svc_low_high_down <- runTest(GOdata_mf_svc_low_high_down, algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_svc_low_mid_up    <- runTest(GOdata_mf_svc_low_mid_up,    algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_svc_low_mid_down  <- runTest(GOdata_mf_svc_low_mid_down,  algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_svc_mid_high_up   <- runTest(GOdata_mf_svc_mid_high_up,   algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_svc_mid_high_down <- runTest(GOdata_mf_svc_mid_high_down, algorithm = "weight01", statistic = "ks") #MF
#resultKS_cc_svc_low_high_up   <- runTest(GOdata_cc_svc_low_high_up,   algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_svc_low_high_down <- runTest(GOdata_cc_svc_low_high_down, algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_svc_low_mid_up    <- runTest(GOdata_cc_svc_low_mid_up,    algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_svc_low_mid_down  <- runTest(GOdata_cc_svc_low_mid_down,  algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_svc_mid_high_up   <- runTest(GOdata_cc_svc_mid_high_up,   algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_svc_mid_high_down <- runTest(GOdata_cc_svc_mid_high_down, algorithm = "weight01", statistic = "ks") #CC
#resultKS_bp_svc_low_high_up   <- runTest(GOdata_bp_svc_low_high_up,   algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_svc_low_high_down <- runTest(GOdata_bp_svc_low_high_down, algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_svc_low_mid_up    <- runTest(GOdata_bp_svc_low_mid_up,    algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_svc_low_mid_down  <- runTest(GOdata_bp_svc_low_mid_down,  algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_svc_mid_high_up   <- runTest(GOdata_bp_svc_mid_high_up,   algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_svc_mid_high_down <- runTest(GOdata_bp_svc_mid_high_down, algorithm = "weight01", statistic = "ks") #BP

#resultKS_mf_wvc_low_high_up   <- runTest(GOdata_mf_wvc_low_high_up,   algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_wvc_low_high_down <- runTest(GOdata_mf_wvc_low_high_down, algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_wvc_low_mid_up    <- runTest(GOdata_mf_wvc_low_mid_up,    algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_wvc_low_mid_down  <- runTest(GOdata_mf_wvc_low_mid_down,  algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_wvc_mid_high_up   <- runTest(GOdata_mf_wvc_mid_high_up,   algorithm = "weight01", statistic = "ks") #MF
#resultKS_mf_wvc_mid_high_down <- runTest(GOdata_mf_wvc_mid_high_down, algorithm = "weight01", statistic = "ks") #MF
#resultKS_bp_wvc_low_high_up   <- runTest(GOdata_bp_wvc_low_high_up,   algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_wvc_low_high_down <- runTest(GOdata_bp_wvc_low_high_down, algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_wvc_low_mid_up    <- runTest(GOdata_bp_wvc_low_mid_up,    algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_wvc_low_mid_down  <- runTest(GOdata_bp_wvc_low_mid_down,  algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_wvc_mid_high_up   <- runTest(GOdata_bp_wvc_mid_high_up,   algorithm = "weight01", statistic = "ks") #BP
#resultKS_bp_wvc_mid_high_down <- runTest(GOdata_bp_wvc_mid_high_down, algorithm = "weight01", statistic = "ks") #BP
#resultKS_cc_wvc_low_high_up   <- runTest(GOdata_cc_wvc_low_high_up,   algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_wvc_low_high_down <- runTest(GOdata_cc_wvc_low_high_down, algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_wvc_low_mid_up    <- runTest(GOdata_cc_wvc_low_mid_up,    algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_wvc_low_mid_down  <- runTest(GOdata_cc_wvc_low_mid_down,  algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_wvc_mid_high_up   <- runTest(GOdata_cc_wvc_mid_high_up,   algorithm = "weight01", statistic = "ks") #CC
#resultKS_cc_wvc_mid_high_down <- runTest(GOdata_cc_wvc_mid_high_down, algorithm = "weight01", statistic = "ks") #CC

####t_weight01####

resultT_mf_vt_up             <- runTest(GOdata_mf_vt_up,             algorithm = "weight01", statistic = "t")  #MF
resultT_mf_vt_down           <- runTest(GOdata_mf_vt_down,           algorithm = "weight01", statistic = "t")  #MF
resultT_bp_vt_up             <- runTest(GOdata_bp_vt_up,             algorithm = "weight01", statistic = "t")  #BP
resultT_bp_vt_down           <- runTest(GOdata_bp_vt_down,           algorithm = "weight01", statistic = "t")  #BP
resultT_cc_vt_up             <- runTest(GOdata_cc_vt_up,             algorithm = "weight01", statistic = "t")  #CC
resultT_cc_vt_down           <- runTest(GOdata_cc_vt_down,           algorithm = "weight01", statistic = "t")  #CC

resultT_mf_ac_low_high_up    <- runTest(GOdata_mf_ac_low_high_up,    algorithm = "weight01", statistic = "t")  #MF
resultT_mf_ac_low_high_down  <- runTest(GOdata_mf_ac_low_high_down,  algorithm = "weight01", statistic = "t")  #MF
resultT_bp_ac_low_high_up    <- runTest(GOdata_bp_ac_low_high_up,    algorithm = "weight01", statistic = "t")  #BP
resultT_bp_ac_low_high_down  <- runTest(GOdata_bp_ac_low_high_down,  algorithm = "weight01", statistic = "t")  #BP
resultT_cc_ac_low_high_up    <- runTest(GOdata_cc_ac_low_high_up,    algorithm = "weight01", statistic = "t")  #CC
resultT_cc_ac_low_high_down  <- runTest(GOdata_cc_ac_low_high_down,  algorithm = "weight01", statistic = "t")  #CC

#resultT_mf_ac_low_mid_up     <- runTest(GOdata_mf_ac_low_mid_up,     algorithm = "weight01", statistic = "t")  #MF
#resultT_mf_ac_low_mid_down   <- runTest(GOdata_mf_ac_low_mid_down,   algorithm = "weight01", statistic = "t")  #MF
#resultT_mf_ac_mid_high_up    <- runTest(GOdata_mf_ac_mid_high_up,    algorithm = "weight01", statistic = "t")  #MF
#resultT_mf_ac_mid_high_down  <- runTest(GOdata_mf_ac_mid_high_down,  algorithm = "weight01", statistic = "t")  #MF

#resultT_bp_ac_low_mid_up     <- runTest(GOdata_bp_ac_low_mid_up,     algorithm = "weight01", statistic = "t")  #BP
#resultT_bp_ac_low_mid_down   <- runTest(GOdata_bp_ac_low_mid_down,   algorithm = "weight01", statistic = "t")  #BP
#resultT_bp_ac_mid_high_up    <- runTest(GOdata_bp_ac_mid_high_up,    algorithm = "weight01", statistic = "t")  #BP
#resultT_bp_ac_mid_high_down  <- runTest(GOdata_bp_ac_mid_high_down,  algorithm = "weight01", statistic = "t")  #BP

#resultT_cc_ac_low_mid_up     <- runTest(GOdata_cc_ac_low_mid_up,     algorithm = "weight01", statistic = "t")  #CC
#resultT_cc_ac_low_mid_down   <- runTest(GOdata_cc_ac_low_mid_down,   algorithm = "weight01", statistic = "t")  #CC
#resultT_cc_ac_mid_high_up    <- runTest(GOdata_cc_ac_mid_high_up,    algorithm = "weight01", statistic = "t")  #CC
#resultT_cc_ac_mid_high_down  <- runTest(GOdata_cc_ac_mid_high_down,  algorithm = "weight01", statistic = "t")  #CC

#resultT_mf_cs_lowwin_lowsum_up     <- runTest(GOdata_mf_cs_lowwin_lowsum_up,     algorithm = "weight01", statistic = "t") #MF
#resultT_mf_cs_lowwin_lowsum_down   <- runTest(GOdata_mf_cs_lowwin_lowsum_down,   algorithm = "weight01", statistic = "t") #MF
#resultT_mf_cs_midwin_midsum_up     <- runTest(GOdata_mf_cs_midwin_midsum_up,     algorithm = "weight01", statistic = "t") #MF
#resultT_mf_cs_midwin_midsum_down   <- runTest(GOdata_mf_cs_midwin_midsum_down,   algorithm = "weight01", statistic = "t") #MF
#resultT_mf_cs_highwin_highsum_up   <- runTest(GOdata_mf_cs_highwin_highsum_up,   algorithm = "weight01", statistic = "t") #MF
#resultT_mf_cs_highwin_highsum_down <- runTest(GOdata_mf_cs_highwin_highsum_down, algorithm = "weight01", statistic = "t") #MF
#resultT_cc_cs_lowwin_lowsum_up     <- runTest(GOdata_cc_cs_lowwin_lowsum_up,     algorithm = "weight01", statistic = "t") #CC
#resultT_cc_cs_lowwin_lowsum_down   <- runTest(GOdata_cc_cs_lowwin_lowsum_down,   algorithm = "weight01", statistic = "t") #CC
#resultT_cc_cs_midwin_midsum_up     <- runTest(GOdata_cc_cs_midwin_midsum_up,     algorithm = "weight01", statistic = "t") #CC
#resultT_cc_cs_midwin_midsum_down   <- runTest(GOdata_cc_cs_midwin_midsum_down,   algorithm = "weight01", statistic = "t") #CC
#resultT_cc_cs_highwin_highsum_up   <- runTest(GOdata_cc_cs_highwin_highsum_up,   algorithm = "weight01", statistic = "t") #CC
#resultT_cc_cs_highwin_highsum_down <- runTest(GOdata_cc_cs_highwin_highsum_down, algorithm = "weight01", statistic = "t") #CC
#resultT_bp_cs_lowwin_lowsum_up     <- runTest(GOdata_bp_cs_lowwin_lowsum_up,     algorithm = "weight01", statistic = "t") #BP
#resultT_bp_cs_lowwin_lowsum_down   <- runTest(GOdata_bp_cs_lowwin_lowsum_down,   algorithm = "weight01", statistic = "t") #BP
#resultT_bp_cs_midwin_midsum_up     <- runTest(GOdata_bp_cs_midwin_midsum_up,     algorithm = "weight01", statistic = "t") #BP
#resultT_bp_cs_midwin_midsum_down   <- runTest(GOdata_bp_cs_midwin_midsum_down,   algorithm = "weight01", statistic = "t") #BP
#resultT_bp_cs_highwin_highsum_up   <- runTest(GOdata_bp_cs_highwin_highsum_up,   algorithm = "weight01", statistic = "t") #BP
#resultT_bp_cs_highwin_highsum_down <- runTest(GOdata_bp_cs_highwin_highsum_down, algorithm = "weight01", statistic = "t") #BP

#resultT_mf_svc_low_high_up   <- runTest(GOdata_mf_svc_low_high_up,   algorithm = "weight01", statistic = "t") #MF
#resultT_mf_svc_low_high_down <- runTest(GOdata_mf_svc_low_high_down, algorithm = "weight01", statistic = "t") #MF
#resultT_mf_svc_low_mid_up    <- runTest(GOdata_mf_svc_low_mid_up,    algorithm = "weight01", statistic = "t") #MF
#resultT_mf_svc_low_mid_down  <- runTest(GOdata_mf_svc_low_mid_down,  algorithm = "weight01", statistic = "t") #MF
#resultT_mf_svc_mid_high_up   <- runTest(GOdata_mf_svc_mid_high_up,   algorithm = "weight01", statistic = "t") #MF
#resultT_mf_svc_mid_high_down <- runTest(GOdata_mf_svc_mid_high_down, algorithm = "weight01", statistic = "t") #MF
#resultT_cc_svc_low_high_up   <- runTest(GOdata_cc_svc_low_high_up,   algorithm = "weight01", statistic = "t") #CC
#resultT_cc_svc_low_high_down <- runTest(GOdata_cc_svc_low_high_down, algorithm = "weight01", statistic = "t") #CC
#resultT_cc_svc_low_mid_up    <- runTest(GOdata_cc_svc_low_mid_up,    algorithm = "weight01", statistic = "t") #CC
#resultT_cc_svc_low_mid_down  <- runTest(GOdata_cc_svc_low_mid_down,  algorithm = "weight01", statistic = "t") #CC
#resultT_cc_svc_mid_high_up   <- runTest(GOdata_cc_svc_mid_high_up,   algorithm = "weight01", statistic = "t") #CC
#resultT_cc_svc_mid_high_down <- runTest(GOdata_cc_svc_mid_high_down, algorithm = "weight01", statistic = "t") #CC
#resultT_bp_svc_low_high_up   <- runTest(GOdata_bp_svc_low_high_up,   algorithm = "weight01", statistic = "t") #BP
#resultT_bp_svc_low_high_down <- runTest(GOdata_bp_svc_low_high_down, algorithm = "weight01", statistic = "t") #BP
#resultT_bp_svc_low_mid_up    <- runTest(GOdata_bp_svc_low_mid_up,    algorithm = "weight01", statistic = "t") #BP
#resultT_bp_svc_low_mid_down  <- runTest(GOdata_bp_svc_low_mid_down,  algorithm = "weight01", statistic = "t") #BP
#resultT_bp_svc_mid_high_up   <- runTest(GOdata_bp_svc_mid_high_up,   algorithm = "weight01", statistic = "t") #BP
#resultT_bp_svc_mid_high_down <- runTest(GOdata_bp_svc_mid_high_down, algorithm = "weight01", statistic = "t") #BP

#resultT_mf_wvc_low_high_up   <- runTest(GOdata_mf_wvc_low_high_up,   algorithm = "weight01", statistic = "t") #MF
#resultT_mf_wvc_low_high_down <- runTest(GOdata_mf_wvc_low_high_down, algorithm = "weight01", statistic = "t") #MF
#resultT_mf_wvc_low_mid_up    <- runTest(GOdata_mf_wvc_low_mid_up,    algorithm = "weight01", statistic = "t") #MF
#resultT_mf_wvc_low_mid_down  <- runTest(GOdata_mf_wvc_low_mid_down,  algorithm = "weight01", statistic = "t") #MF
#resultT_mf_wvc_mid_high_up   <- runTest(GOdata_mf_wvc_mid_high_up,   algorithm = "weight01", statistic = "t") #MF
#resultT_mf_wvc_mid_high_down <- runTest(GOdata_mf_wvc_mid_high_down, algorithm = "weight01", statistic = "t") #MF
#resultT_bp_wvc_low_high_up   <- runTest(GOdata_bp_wvc_low_high_up,   algorithm = "weight01", statistic = "t") #BP
#resultT_bp_wvc_low_high_down <- runTest(GOdata_bp_wvc_low_high_down, algorithm = "weight01", statistic = "t") #BP
#resultT_bp_wvc_low_mid_up    <- runTest(GOdata_bp_wvc_low_mid_up,    algorithm = "weight01", statistic = "t") #BP
#resultT_bp_wvc_low_mid_down  <- runTest(GOdata_bp_wvc_low_mid_down,  algorithm = "weight01", statistic = "t") #BP
#resultT_bp_wvc_mid_high_up   <- runTest(GOdata_bp_wvc_mid_high_up,   algorithm = "weight01", statistic = "t") #BP
#resultT_bp_wvc_mid_high_down <- runTest(GOdata_bp_wvc_mid_high_down, algorithm = "weight01", statistic = "t") #BP
#resultT_cc_wvc_low_high_up   <- runTest(GOdata_cc_wvc_low_high_up,   algorithm = "weight01", statistic = "t") #CC
#resultT_cc_wvc_low_high_down <- runTest(GOdata_cc_wvc_low_high_down, algorithm = "weight01", statistic = "t") #CC
#resultT_cc_wvc_low_mid_up    <- runTest(GOdata_cc_wvc_low_mid_up,    algorithm = "weight01", statistic = "t") #CC
#resultT_cc_wvc_low_mid_down  <- runTest(GOdata_cc_wvc_low_mid_down,  algorithm = "weight01", statistic = "t") #CC
#resultT_cc_wvc_mid_high_up   <- runTest(GOdata_cc_wvc_mid_high_up,   algorithm = "weight01", statistic = "t") #CC
#resultT_cc_wvc_mid_high_down <- runTest(GOdata_cc_wvc_mid_high_down, algorithm = "weight01", statistic = "t") #CC

####RESULTS####

allRes_mf_vt_up    <- GenTable(GOdata_mf_vt_up,    ParentChild_Fisher = resultFisher_mf_vt_up,    Weight01_T = resultT_mf_vt_up,       Weight01_KS = resultKS_mf_vt_up,   topNodes = 30)
allRes_mf_vt_down  <- GenTable(GOdata_mf_vt_down,  ParentChild_Fisher = resultFisher_mf_vt_down,  Weight01_T = resultT_mf_vt_down,     Weight01_KS = resultKS_mf_vt_down, topNodes = 30)
allRes_bp_vt_up    <- GenTable(GOdata_bp_vt_up,    ParentChild_Fisher = resultFisher_bp_vt_up,    Weight01_T = resultT_bp_vt_up,       Weight01_KS = resultKS_bp_vt_up,   topNodes = 30)
allRes_bp_vt_down  <- GenTable(GOdata_bp_vt_down,  ParentChild_Fisher = resultFisher_bp_vt_down,  Weight01_T = resultT_bp_vt_down,     Weight01_KS = resultKS_bp_vt_down, topNodes = 30)
allRes_cc_vt_up    <- GenTable(GOdata_cc_vt_up,    ParentChild_Fisher = resultFisher_cc_vt_up,    Weight01_T = resultT_cc_vt_up,       Weight01_KS = resultKS_cc_vt_up,   topNodes = 30)
allRes_cc_vt_down  <- GenTable(GOdata_cc_vt_down,  ParentChild_Fisher = resultFisher_cc_vt_down,  Weight01_T = resultT_cc_vt_down,     Weight01_KS = resultKS_cc_vt_down, topNodes = 30)

allRes_mf_ac_low_high_up    <- GenTable(GOdata_mf_ac_low_high_up,    ParentChild_Fisher = resultFisher_mf_ac_low_high_up,    Weight01_T = resultT_mf_ac_low_high_up,    Weight01_KS = resultKS_mf_ac_low_high_up,   topNodes = 30)
allRes_mf_ac_low_high_down  <- GenTable(GOdata_mf_ac_low_high_down,  ParentChild_Fisher = resultFisher_mf_ac_low_high_down,  Weight01_T = resultT_mf_ac_low_high_down,  Weight01_KS = resultKS_mf_ac_low_high_down, topNodes = 30)
allRes_bp_ac_low_high_up    <- GenTable(GOdata_bp_ac_low_high_up,    ParentChild_Fisher = resultFisher_bp_ac_low_high_up,    Weight01_T = resultT_bp_ac_low_high_up,    Weight01_KS = resultKS_bp_ac_low_high_up,   topNodes = 30)
allRes_bp_ac_low_high_down  <- GenTable(GOdata_bp_ac_low_high_down,  ParentChild_Fisher = resultFisher_bp_ac_low_high_down,  Weight01_T = resultT_bp_ac_low_high_down,  Weight01_KS = resultKS_bp_ac_low_high_down, topNodes = 30)
allRes_cc_ac_low_high_up    <- GenTable(GOdata_cc_ac_low_high_up,    ParentChild_Fisher = resultFisher_cc_ac_low_high_up,    Weight01_KS = resultKS_cc_ac_low_high_up,   topNodes = 30)
allRes_cc_ac_low_high_down  <- GenTable(GOdata_cc_ac_low_high_down,  ParentChild_Fisher = resultFisher_cc_ac_low_high_down,  Weight01_KS = resultKS_cc_ac_low_high_down, topNodes = 30)

#allRes_mf_ac_low_mid_up     <- GenTable(GOdata_mf_ac_low_mid_up,     ParentChild_Fisher = resultFisher_mf_ac_low_mid_up,     Weight01_T = resultT_mf_ac_low_mid_up,     Weight01_KS = resultKS_mf_ac_low_mid_up,    topNodes = 30)
#allRes_mf_ac_low_mid_down   <- GenTable(GOdata_mf_ac_low_mid_down,   ParentChild_Fisher = resultFisher_mf_ac_low_mid_down,   Weight01_T = resultT_mf_ac_low_mid_down,   Weight01_KS = resultKS_mf_ac_low_mid_down,  topNodes = 30)
#allRes_mf_ac_mid_high_up    <- GenTable(GOdata_mf_ac_mid_high_up,    ParentChild_Fisher = resultFisher_mf_ac_mid_high_up,    Weight01_T = resultT_mf_ac_mid_high_up,    Weight01_KS = resultKS_mf_ac_mid_high_up,   topNodes = 30)
#allRes_mf_ac_mid_high_down  <- GenTable(GOdata_mf_ac_mid_high_down,  ParentChild_Fisher = resultFisher_mf_ac_mid_high_down,  Weight01_T = resultT_mf_ac_mid_high_down,  Weight01_KS = resultKS_mf_ac_mid_high_down, topNodes = 30)

#allRes_bp_ac_low_mid_up     <- GenTable(GOdata_bp_ac_low_mid_up,     ParentChild_Fisher = resultFisher_bp_ac_low_mid_up,     Weight01_T = resultT_bp_ac_low_mid_up,     Weight01_KS = resultKS_bp_ac_low_mid_up,    topNodes = 30)
#allRes_bp_ac_low_mid_down   <- GenTable(GOdata_bp_ac_low_mid_down,   ParentChild_Fisher = resultFisher_bp_ac_low_mid_down,   Weight01_T = resultT_bp_ac_low_mid_down,   Weight01_KS = resultKS_bp_ac_low_mid_down,  topNodes = 30)
#allRes_bp_ac_mid_high_up    <- GenTable(GOdata_bp_ac_mid_high_up,    ParentChild_Fisher = resultFisher_bp_ac_mid_high_up,    Weight01_T = resultT_bp_ac_mid_high_up,    Weight01_KS = resultKS_bp_ac_mid_high_up,   topNodes = 30)
#allRes_bp_ac_mid_high_down  <- GenTable(GOdata_bp_ac_mid_high_down,  ParentChild_Fisher = resultFisher_bp_ac_mid_high_down,  Weight01_T = resultT_bp_ac_mid_high_down,  Weight01_KS = resultKS_bp_ac_mid_high_down, topNodes = 30)

#allRes_cc_ac_low_mid_up     <- GenTable(GOdata_cc_ac_low_mid_up,     ParentChild_Fisher = resultFisher_cc_ac_low_mid_up,     Weight01_T = resultT_cc_ac_low_mid_up,     Weight01_KS = resultKS_cc_ac_low_mid_up,    topNodes = 30)
#allRes_cc_ac_low_mid_down   <- GenTable(GOdata_cc_ac_low_mid_down,   ParentChild_Fisher = resultFisher_cc_ac_low_mid_down,   Weight01_T = resultT_cc_ac_low_mid_down,   Weight01_KS = resultKS_cc_ac_low_mid_down,  topNodes = 30)
#allRes_cc_ac_mid_high_up    <- GenTable(GOdata_cc_ac_mid_high_up,    ParentChild_Fisher = resultFisher_cc_ac_mid_high_up,    Weight01_T = resultT_cc_ac_mid_high_up,    Weight01_KS = resultKS_cc_ac_mid_high_up,   topNodes = 30)
#allRes_cc_ac_mid_high_down  <- GenTable(GOdata_cc_ac_mid_high_down,  ParentChild_Fisher = resultFisher_cc_ac_mid_high_down,  Weight01_T = resultT_cc_ac_mid_high_down,  Weight01_KS = resultKS_cc_ac_mid_high_down, topNodes = 30)

#allRes_mf_cs_lowwin_lowsum_up       <-GenTable(GOdata_mf_cs_lowwin_lowsum_up,      ParentChild_Fisher = resultFisher_mf_cs_lowwin_lowsum_up,      Weight01_T = resultT_mf_cs_lowwin_lowsum_up,      Weight01_KS = resultKS_mf_cs_lowwin_lowsum_up,     topNodes = 30)
#allRes_mf_cs_lowwin_lowsum_down     <-GenTable(GOdata_mf_cs_lowwin_lowsum_down,    ParentChild_Fisher = resultFisher_mf_cs_lowwin_lowsum_down,    Weight01_T = resultT_mf_cs_lowwin_lowsum_down,    Weight01_KS = resultKS_mf_cs_lowwin_lowsum_down,   topNodes = 30)
#allRes_mf_cs_midwin_midsum_up       <-GenTable(GOdata_mf_cs_midwin_midsum_up,      ParentChild_Fisher = resultFisher_mf_cs_midwin_midsum_up,      Weight01_T = resultT_mf_cs_midwin_midsum_up,      Weight01_KS = resultKS_mf_cs_midwin_midsum_up,     topNodes = 30)
#allRes_mf_cs_midwin_midsum_down     <-GenTable(GOdata_mf_cs_midwin_midsum_down,    ParentChild_Fisher = resultFisher_mf_cs_midwin_midsum_down,    Weight01_T = resultT_mf_cs_midwin_midsum_down,    Weight01_KS = resultKS_mf_cs_midwin_midsum_down,   topNodes = 30)
#allRes_mf_cs_highwin_highsum_up     <-GenTable(GOdata_mf_cs_highwin_highsum_up,    ParentChild_Fisher = resultFisher_mf_cs_highwin_highsum_up,    Weight01_T = resultT_mf_cs_highwin_highsum_up,    Weight01_KS = resultKS_mf_cs_highwin_highsum_up,   topNodes = 30)
#allRes_mf_cs_highwin_highsum_down   <-GenTable(GOdata_mf_cs_highwin_highsum_down,  ParentChild_Fisher = resultFisher_mf_cs_highwin_highsum_down,  Weight01_T = resultT_mf_cs_highwin_highsum_down,  Weight01_KS = resultKS_mf_cs_highwin_highsum_down, topNodes = 30)
#allRes_bp_cs_lowwin_lowsum_up       <-GenTable(GOdata_bp_cs_lowwin_lowsum_up,      ParentChild_Fisher = resultFisher_bp_cs_lowwin_lowsum_up,      Weight01_T = resultT_bp_cs_lowwin_lowsum_up,      Weight01_KS = resultKS_bp_cs_lowwin_lowsum_up,     topNodes = 30)
#allRes_bp_cs_lowwin_lowsum_down     <-GenTable(GOdata_bp_cs_lowwin_lowsum_down,    ParentChild_Fisher = resultFisher_bp_cs_lowwin_lowsum_down,    Weight01_T = resultT_bp_cs_lowwin_lowsum_down,    Weight01_KS = resultKS_bp_cs_lowwin_lowsum_down,   topNodes = 30)
#allRes_bp_cs_midwin_midsum_up       <-GenTable(GOdata_bp_cs_midwin_midsum_up,      ParentChild_Fisher = resultFisher_bp_cs_midwin_midsum_up,      Weight01_T = resultT_bp_cs_midwin_midsum_up,      Weight01_KS = resultKS_bp_cs_midwin_midsum_up,     topNodes = 30)
#allRes_bp_cs_midwin_midsum_down     <-GenTable(GOdata_bp_cs_midwin_midsum_down,    ParentChild_Fisher = resultFisher_bp_cs_midwin_midsum_down,    Weight01_T = resultT_bp_cs_midwin_midsum_down,    Weight01_KS = resultKS_bp_cs_midwin_midsum_down,   topNodes = 30)
#allRes_bp_cs_highwin_highsum_up     <-GenTable(GOdata_bp_cs_highwin_highsum_up,    ParentChild_Fisher = resultFisher_bp_cs_highwin_highsum_up,    Weight01_T = resultT_bp_cs_highwin_highsum_up,    Weight01_KS = resultKS_bp_cs_highwin_highsum_up,   topNodes = 30)
#allRes_bp_cs_highwin_highsum_down   <-GenTable(GOdata_bp_cs_highwin_highsum_down,  ParentChild_Fisher = resultFisher_bp_cs_highwin_highsum_down,  Weight01_T = resultT_bp_cs_highwin_highsum_down,  Weight01_KS = resultKS_bp_cs_highwin_highsum_down, topNodes = 30)
#allRes_cc_cs_lowwin_lowsum_up       <-GenTable(GOdata_cc_cs_lowwin_lowsum_up,      ParentChild_Fisher = resultFisher_cc_cs_lowwin_lowsum_up,      Weight01_T = resultT_cc_cs_lowwin_lowsum_up,      Weight01_KS = resultKS_cc_cs_lowwin_lowsum_up,     topNodes = 30)
#allRes_cc_cs_lowwin_lowsum_down     <-GenTable(GOdata_cc_cs_lowwin_lowsum_down,    ParentChild_Fisher = resultFisher_cc_cs_lowwin_lowsum_down,    Weight01_T = resultT_cc_cs_lowwin_lowsum_down,    Weight01_KS = resultKS_cc_cs_lowwin_lowsum_down,   topNodes = 30)
#allRes_cc_cs_midwin_midsum_up       <-GenTable(GOdata_cc_cs_midwin_midsum_up,      ParentChild_Fisher = resultFisher_cc_cs_midwin_midsum_up,      Weight01_T = resultT_cc_cs_midwin_midsum_up,      Weight01_KS = resultKS_cc_cs_midwin_midsum_up,     topNodes = 30)
#allRes_cc_cs_midwin_midsum_down     <-GenTable(GOdata_cc_cs_midwin_midsum_down,    ParentChild_Fisher = resultFisher_cc_cs_midwin_midsum_down,    Weight01_T = resultT_cc_cs_midwin_midsum_down,    Weight01_KS = resultKS_cc_cs_midwin_midsum_down,   topNodes = 30)
#allRes_cc_cs_highwin_highsum_up     <-GenTable(GOdata_cc_cs_highwin_highsum_up,    ParentChild_Fisher = resultFisher_cc_cs_highwin_highsum_up,    Weight01_T = resultT_cc_cs_highwin_highsum_up,    Weight01_KS = resultKS_cc_cs_highwin_highsum_up,   topNodes = 30)
#allRes_cc_cs_highwin_highsum_down   <-GenTable(GOdata_cc_cs_highwin_highsum_down,  ParentChild_Fisher = resultFisher_cc_cs_highwin_highsum_down,  Weight01_T = resultT_cc_cs_highwin_highsum_down,  Weight01_KS = resultKS_cc_cs_highwin_highsum_down, topNodes = 30)

#allRes_mf_svc_low_high_up   <- GenTable(GOdata_mf_svc_low_high_up,   ParentChild_Fisher = resultFisher_mf_svc_low_high_up,   Weight01_T = resultT_mf_svc_low_high_up,   Weight01_KS = resultKS_mf_svc_low_high_up,  topNodes = 30)
#allRes_mf_svc_low_high_down <- GenTable(GOdata_mf_svc_low_high_down, ParentChild_Fisher = resultFisher_mf_svc_low_high_down, Weight01_T = resultT_mf_svc_low_high_down, Weight01_KS = resultKS_mf_svc_low_high_down,topNodes = 30)
#allRes_mf_svc_low_mid_up    <- GenTable(GOdata_mf_svc_low_mid_up,    ParentChild_Fisher = resultFisher_mf_svc_low_mid_up,    Weight01_T = resultT_mf_svc_low_mid_up,    Weight01_KS = resultKS_mf_svc_low_mid_up,   topNodes = 30)
#allRes_mf_svc_low_mid_down  <- GenTable(GOdata_mf_svc_low_mid_down,  ParentChild_Fisher = resultFisher_mf_svc_low_mid_down,  Weight01_T = resultT_mf_svc_low_mid_down,  Weight01_KS = resultKS_mf_svc_low_mid_down, topNodes = 30)
#allRes_mf_svc_mid_high_up   <- GenTable(GOdata_mf_svc_mid_high_up,   ParentChild_Fisher = resultFisher_mf_svc_mid_high_up,   Weight01_T = resultT_mf_svc_mid_high_up,   Weight01_KS = resultKS_mf_svc_mid_high_up,  topNodes = 30)
#allRes_mf_svc_mid_high_down <- GenTable(GOdata_mf_svc_mid_high_down, ParentChild_Fisher = resultFisher_mf_svc_mid_high_down, Weight01_T = resultT_mf_svc_mid_high_down, Weight01_KS = resultKS_mf_svc_mid_high_down,topNodes = 30)
#allRes_bp_svc_low_high_up   <- GenTable(GOdata_bp_svc_low_high_up,   ParentChild_Fisher = resultFisher_bp_svc_low_high_up,   Weight01_T = resultT_bp_svc_low_high_up,   Weight01_KS = resultKS_bp_svc_low_high_up,  topNodes = 30)
#allRes_bp_svc_low_high_down <- GenTable(GOdata_bp_svc_low_high_down, ParentChild_Fisher = resultFisher_bp_svc_low_high_down, Weight01_T = resultT_bp_svc_low_high_down, Weight01_KS = resultKS_bp_svc_low_high_down,topNodes = 30)
#allRes_bp_svc_low_mid_up    <- GenTable(GOdata_bp_svc_low_mid_up,    ParentChild_Fisher = resultFisher_bp_svc_low_mid_up,    Weight01_T = resultT_bp_svc_low_mid_up,    Weight01_KS = resultKS_bp_svc_low_mid_up,   topNodes = 30)
#allRes_bp_svc_low_mid_down  <- GenTable(GOdata_bp_svc_low_mid_down,  ParentChild_Fisher = resultFisher_bp_svc_low_mid_down,  Weight01_T = resultT_bp_svc_low_mid_down,  Weight01_KS = resultKS_bp_svc_low_mid_down, topNodes = 30)
#allRes_bp_svc_mid_high_up   <- GenTable(GOdata_bp_svc_mid_high_up,   ParentChild_Fisher = resultFisher_bp_svc_mid_high_up,   Weight01_T = resultT_bp_svc_mid_high_up,   Weight01_KS = resultKS_bp_svc_mid_high_up,  topNodes = 30)
#allRes_bp_svc_mid_high_down <- GenTable(GOdata_bp_svc_mid_high_down, ParentChild_Fisher = resultFisher_bp_svc_mid_high_down, Weight01_T = resultT_bp_svc_mid_high_down, Weight01_KS = resultKS_bp_svc_mid_high_down,topNodes = 30)
#allRes_cc_svc_low_high_up   <- GenTable(GOdata_cc_svc_low_high_up,   ParentChild_Fisher = resultFisher_cc_svc_low_high_up,   Weight01_T = resultT_cc_svc_low_high_up,   Weight01_KS = resultKS_cc_svc_low_high_up,  topNodes = 30)
#allRes_cc_svc_low_high_down <- GenTable(GOdata_cc_svc_low_high_down, ParentChild_Fisher = resultFisher_cc_svc_low_high_down, Weight01_T = resultT_cc_svc_low_high_down, Weight01_KS = resultKS_cc_svc_low_high_down,topNodes = 30)
#allRes_cc_svc_low_mid_up    <- GenTable(GOdata_cc_svc_low_mid_up,    ParentChild_Fisher = resultFisher_cc_svc_low_mid_up,    Weight01_T = resultT_cc_svc_low_mid_up,    Weight01_KS = resultKS_cc_svc_low_mid_up,   topNodes = 30)
#allRes_cc_svc_low_mid_down  <- GenTable(GOdata_cc_svc_low_mid_down,  ParentChild_Fisher = resultFisher_cc_svc_low_mid_down,  Weight01_T = resultT_cc_svc_low_mid_down,  Weight01_KS = resultKS_cc_svc_low_mid_down, topNodes = 30)
#allRes_cc_svc_mid_high_up   <- GenTable(GOdata_cc_svc_mid_high_up,   ParentChild_Fisher = resultFisher_cc_svc_mid_high_up,   Weight01_T = resultT_cc_svc_mid_high_up,   Weight01_KS = resultKS_cc_svc_mid_high_up,  topNodes = 30)
#allRes_cc_svc_mid_high_down <- GenTable(GOdata_cc_svc_mid_high_down, ParentChild_Fisher = resultFisher_cc_svc_mid_high_down, Weight01_T = resultT_cc_svc_mid_high_down, Weight01_KS = resultKS_cc_svc_mid_high_down,topNodes = 30)

#allRes_mf_wvc_low_high_up   <- GenTable(GOdata_mf_wvc_low_high_up,   ParentChild_Fisher = resultFisher_mf_wvc_low_high_up,   Weight01_T = resultT_mf_wvc_low_high_up,   Weight01_KS = resultKS_mf_wvc_low_high_up,  topNodes = 30)
#allRes_mf_wvc_low_high_down <- GenTable(GOdata_mf_wvc_low_high_down, ParentChild_Fisher = resultFisher_mf_wvc_low_high_down, Weight01_T = resultT_mf_wvc_low_high_down, Weight01_KS = resultKS_mf_wvc_low_high_down,topNodes = 30)
#allRes_mf_wvc_low_mid_up    <- GenTable(GOdata_mf_wvc_low_mid_up,    ParentChild_Fisher = resultFisher_mf_wvc_low_mid_up,    Weight01_T = resultT_mf_wvc_low_mid_up,    Weight01_KS = resultKS_mf_wvc_low_mid_up,   topNodes = 30)
#allRes_mf_wvc_low_mid_down  <- GenTable(GOdata_mf_wvc_low_mid_down,  ParentChild_Fisher = resultFisher_mf_wvc_low_mid_down,  Weight01_T = resultT_mf_wvc_low_mid_down,  Weight01_KS = resultKS_mf_wvc_low_mid_down, topNodes = 30)
#allRes_mf_wvc_mid_high_up   <- GenTable(GOdata_mf_wvc_mid_high_up,   ParentChild_Fisher = resultFisher_mf_wvc_mid_high_up,   Weight01_T = resultT_mf_wvc_mid_high_up,   Weight01_KS = resultKS_mf_wvc_mid_high_up,  topNodes = 30)
#allRes_mf_wvc_mid_high_down <- GenTable(GOdata_mf_wvc_mid_high_down, ParentChild_Fisher = resultFisher_mf_wvc_mid_high_down, Weight01_T = resultT_mf_wvc_mid_high_down, Weight01_KS = resultKS_mf_wvc_mid_high_down,topNodes = 30)
#allRes_bp_wvc_low_high_up   <- GenTable(GOdata_bp_wvc_low_high_up,   ParentChild_Fisher = resultFisher_bp_wvc_low_high_up,   Weight01_T = resultT_bp_wvc_low_high_up,   Weight01_KS = resultKS_bp_wvc_low_high_up,  topNodes = 30)
#allRes_bp_wvc_low_high_down <- GenTable(GOdata_bp_wvc_low_high_down, ParentChild_Fisher = resultFisher_bp_wvc_low_high_down, Weight01_T = resultT_bp_wvc_low_high_down, Weight01_KS = resultKS_bp_wvc_low_high_down,topNodes = 30)
#allRes_bp_wvc_low_mid_up    <- GenTable(GOdata_bp_wvc_low_mid_up,    ParentChild_Fisher = resultFisher_bp_wvc_low_mid_up,    Weight01_T = resultT_bp_wvc_low_mid_up,    Weight01_KS = resultKS_bp_wvc_low_mid_up,   topNodes = 30)
#allRes_bp_wvc_low_mid_down  <- GenTable(GOdata_bp_wvc_low_mid_down,  ParentChild_Fisher = resultFisher_bp_wvc_low_mid_down,  Weight01_T = resultT_bp_wvc_low_mid_down,  Weight01_KS = resultKS_bp_wvc_low_mid_down, topNodes = 30)
#allRes_bp_wvc_mid_high_up   <- GenTable(GOdata_bp_wvc_mid_high_up,   ParentChild_Fisher = resultFisher_bp_wvc_mid_high_up,   Weight01_T = resultT_bp_wvc_mid_high_up,   Weight01_KS = resultKS_bp_wvc_mid_high_up,  topNodes = 30)
#allRes_bp_wvc_mid_high_down <- GenTable(GOdata_bp_wvc_mid_high_down, ParentChild_Fisher = resultFisher_bp_wvc_mid_high_down, Weight01_T = resultT_bp_wvc_mid_high_down, Weight01_KS = resultKS_bp_wvc_mid_high_down,topNodes = 30)
#allRes_cc_wvc_low_high_up   <- GenTable(GOdata_cc_wvc_low_high_up,   ParentChild_Fisher = resultFisher_cc_wvc_low_high_up,   Weight01_T = resultT_cc_wvc_low_high_up,   Weight01_KS = resultKS_cc_wvc_low_high_up,  topNodes = 30)
#allRes_cc_wvc_low_high_down <- GenTable(GOdata_cc_wvc_low_high_down, ParentChild_Fisher = resultFisher_cc_wvc_low_high_down, Weight01_T = resultT_cc_wvc_low_high_down, Weight01_KS = resultKS_cc_wvc_low_high_down,topNodes = 30)
#allRes_cc_wvc_low_mid_up    <- GenTable(GOdata_cc_wvc_low_mid_up,    ParentChild_Fisher = resultFisher_cc_wvc_low_mid_up,    Weight01_T = resultT_cc_wvc_low_mid_up,    Weight01_KS = resultKS_cc_wvc_low_mid_up,   topNodes = 30)
#allRes_cc_wvc_low_mid_down  <- GenTable(GOdata_cc_wvc_low_mid_down,  ParentChild_Fisher = resultFisher_cc_wvc_low_mid_down,  Weight01_T = resultT_cc_wvc_low_mid_down,  Weight01_KS = resultKS_cc_wvc_low_mid_down, topNodes = 30)
#allRes_cc_wvc_mid_high_up   <- GenTable(GOdata_cc_wvc_mid_high_up,   ParentChild_Fisher = resultFisher_cc_wvc_mid_high_up,   Weight01_T = resultT_cc_wvc_mid_high_up,   Weight01_KS = resultKS_cc_wvc_mid_high_up,  topNodes = 30)
#allRes_cc_wvc_mid_high_down <- GenTable(GOdata_cc_wvc_mid_high_down, ParentChild_Fisher = resultFisher_cc_wvc_mid_high_down, Weight01_T = resultT_cc_wvc_mid_high_down, Weight01_KS = resultKS_cc_wvc_mid_high_down,topNodes = 30)

allRes_bp_vt_up$cat   <- "BP";
allRes_bp_vt_down$cat   <- "BP";
allRes_bp_ac_low_high_up$cat   <- "BP"; 
allRes_bp_ac_low_high_down$cat   <- "BP"; 

allRes_cc_vt_down$cat <- "CC"; 
allRes_cc_vt_up$cat   <- "CC"; 
allRes_cc_ac_low_high_up$cat   <- "CC"; 
allRes_cc_ac_low_high_down$cat   <- "CC"; 

allRes_mf_vt_up$cat   <- "MF"; 
allRes_mf_vt_down$cat <- "MF"; 
allRes_mf_ac_low_high_up$cat   <- "MF"; 
allRes_mf_ac_low_high_down$cat <- "MF"; 



#allRes_bp_ac_low_mid_up$cat   <- "BP";
#allRes_bp_ac_mid_high_up$cat   <- "BP"; 


#allRes_bp_svc_low_high_up$cat   <- "BP";
#allRes_bp_svc_low_mid_up$cat   <- "BP";
#allRes_bp_svc_mid_high_up$cat   <- "BP"; 
#allRes_bp_wvc_low_high_up$cat   <- "BP"; 
#allRes_bp_wvc_low_mid_up$cat   <- "BP"; 
#allRes_bp_wvc_mid_high_up$cat      <- "BP" ; 
#allRes_bp_cs_lowwin_lowsum_up$cat      <- "BP" ; 
#allRes_bp_cs_midwin_midsum_up$cat      <- "BP" ; 
#allRes_bp_cs_highwin_highsum_up$cat      <- "BP"
#allRes_bp_vt_down$cat <- "BP"; 
#allRes_bp_ac_low_high_down$cat <- "BP"; 
#allRes_bp_ac_low_mid_down$cat <- "BP";
#allRes_bp_ac_mid_high_down$cat <- "BP"; 
#allRes_bp_svc_low_high_down$cat <- "BP";
#allRes_bp_svc_low_mid_down$cat <- "BP";
#allRes_bp_svc_mid_high_down$cat <- "BP"; 
#allRes_bp_wvc_low_high_down$cat <- "BP"; 
#allRes_bp_wvc_low_mid_down$cat <- "BP"; 
#allRes_bp_wvc_mid_high_down$cat    <- "BP" ; 
#allRes_bp_cs_lowwin_lowsum_down$cat    <- "BP" ; 
#allRes_bp_cs_midwin_midsum_down$cat    <- "BP" ; 
#allRes_bp_cs_highwin_highsum_up$cat      <- "BP"
allRes_cc_vt_up$cat   <- "CC"; 
allRes_cc_ac_low_high_up$cat   <- "CC"; 
#allRes_cc_ac_low_mid_up$cat   <- "CC";
#allRes_cc_ac_mid_high_up$cat   <- "CC"; 
#allRes_cc_svc_low_high_up$cat   <- "CC";
#allRes_cc_svc_low_mid_up$cat   <- "CC";
#allRes_cc_svc_mid_high_up$cat   <- "CC"; 
#allRes_cc_wvc_low_high_up$cat   <- "CC"; 
#allRes_cc_wvc_low_mid_up$cat   <- "CC"; 
#allRes_cc_wvc_mid_high_up$cat      <- "CC" ; 
#allRes_cc_cs_lowwin_lowsum_up$cat      <- "CC" ; 
#allRes_cc_cs_midwin_midsum_up$cat      <- "CC" ; 
#allRes_cc_cs_highwin_highsum_up$cat      <- "CC"

allRes_cc_ac_low_high_down$cat <- "CC"; 
#allRes_cc_ac_low_mid_down$cat <- "CC";
#allRes_cc_ac_mid_high_down$cat <- "CC"; 
#allRes_cc_svc_low_high_down$cat <- "CC";
#allRes_cc_svc_low_mid_down$cat <- "CC";
#allRes_cc_svc_mid_high_down$cat <- "CC"; 
#allRes_cc_wvc_low_high_down$cat <- "CC"; 
#allRes_cc_wvc_low_mid_down$cat <- "CC"; 
#allRes_cc_wvc_mid_high_down$cat    <- "CC" ; 
#allRes_cc_cs_lowwin_lowsum_down$cat    <- "CC" ; 
#allRes_cc_cs_midwin_midsum_down$cat    <- "CC" ; 
#allRes_cc_cs_highwin_highsum_up$cat      <- "CC"

#allRes_mf_ac_low_mid_up$cat   <- "MF";


#allRes_mf_ac_mid_high_up$cat   <- "MF"; 
#allRes_mf_svc_low_high_up$cat   <- "MF";
#allRes_mf_svc_low_mid_up$cat   <- "MF";
#allRes_mf_svc_mid_high_up$cat   <- "MF"; 
#allRes_mf_wvc_low_high_up$cat   <- "MF"; 
#allRes_mf_wvc_low_mid_up$cat   <- "MF"; 
#allRes_mf_wvc_mid_high_up$cat      <- "MF" ; 
#allRes_mf_cs_lowwin_lowsum_up$cat      <- "MF" ; 
#allRes_mf_cs_midwin_midsum_up$cat      <- "MF" ; 
#allRes_mf_cs_highwin_highsum_up$cat      <- "MF"
#allRes_mf_ac_low_mid_down$cat <- "MF";
#allRes_mf_ac_mid_high_down$cat <- "MF"; 
#allRes_mf_svc_low_high_down$cat <- "MF";
#allRes_mf_svc_low_mid_down$cat <- "MF";
#allRes_mf_svc_mid_high_down$cat <- "MF"; 
#allRes_mf_wvc_low_high_down$cat <- "MF"; 
#allRes_mf_wvc_low_mid_down$cat <- "MF"; 
#allRes_mf_wvc_mid_high_down$cat    <- "MF" ; 
#allRes_mf_cs_lowwin_lowsum_down$cat    <- "MF" ; 
#allRes_mf_cs_midwin_midsum_down$cat    <- "MF" ; 
#allRes_mf_cs_highwin_highsum_up$cat      <- "MF"


allRes_vt_up    <- data.frame(rbind(allRes_bp_vt_up,    allRes_cc_vt_up,    allRes_mf_vt_up))
allRes_vt_down  <- data.frame(rbind(allRes_bp_vt_down,  allRes_cc_vt_down,  allRes_mf_vt_down))

allRes_cc_ac_low_high_up$Weight01_T=NaN
allRes_cc_ac_low_high_up$`Rank in Weight01_T`=NaN
allRes_cc_ac_low_high_up$`Rank in Weight01_KS`=NULL
allRes_ac_low_high_up    <- data.frame(rbind(allRes_bp_ac_low_high_up,    allRes_cc_ac_low_high_up,    allRes_mf_ac_low_high_up))

allRes_cc_ac_low_high_down$Weight01_T=NaN
allRes_cc_ac_low_high_down$`Rank in Weight01_T`=NaN
allRes_cc_ac_low_high_down$`Rank in Weight01_KS`=NULL
allRes_ac_low_high_down  <- data.frame(rbind(allRes_bp_ac_low_high_down,  allRes_cc_ac_low_high_down,  allRes_mf_ac_low_high_down))


#allRes_ac_low_mid_up     <- data.frame(rbind(allRes_bp_ac_low_mid_up,     allRes_cc_ac_low_mid_up,     allRes_mf_ac_low_mid_up))
#allRes_ac_low_mid_down   <- data.frame(rbind(allRes_bp_ac_low_mid_down,   allRes_cc_ac_low_mid_down,   allRes_mf_ac_low_mid_down))
#allRes_ac_mid_high_up    <- data.frame(rbind(allRes_bp_ac_mid_high_up,    allRes_cc_ac_mid_high_up,    allRes_mf_ac_mid_high_up))
#allRes_ac_mid_high_down  <- data.frame(rbind(allRes_bp_ac_mid_high_down,  allRes_cc_ac_mid_high_down,  allRes_mf_ac_mid_high_down))

#allRes_cs_lowwin_lowsum_up     <- data.frame(rbind(allRes_bp_cs_lowwin_lowsum_up,     allRes_cc_cs_lowwin_lowsum_up,     allRes_mf_cs_lowwin_lowsum_up))
#allRes_cs_lowwin_lowsum_down   <- data.frame(rbind(allRes_bp_cs_lowwin_lowsum_down,   allRes_cc_cs_lowwin_lowsum_down,   allRes_mf_cs_lowwin_lowsum_down))
#allRes_cs_midwin_midsum_up     <- data.frame(rbind(allRes_bp_cs_midwin_midsum_up,     allRes_cc_cs_midwin_midsum_up,     allRes_mf_cs_midwin_midsum_up))
#allRes_cs_midwin_midsum_down   <- data.frame(rbind(allRes_bp_cs_midwin_midsum_down,   allRes_cc_cs_midwin_midsum_down,   allRes_mf_cs_midwin_midsum_down))
#allRes_cs_highwin_highsum_up   <- data.frame(rbind(allRes_bp_cs_highwin_highsum_up,   allRes_cc_cs_highwin_highsum_up,   allRes_mf_cs_highwin_highsum_up))
#allRes_cs_highwin_highsum_down <- data.frame(rbind(allRes_bp_cs_highwin_highsum_down, allRes_cc_cs_highwin_highsum_down, allRes_mf_cs_highwin_highsum_down))

#allRes_svc_low_high_up   <- data.frame(rbind(allRes_bp_svc_low_high_up,   allRes_cc_svc_low_high_up,   allRes_mf_svc_low_high_up))
#allRes_svc_low_high_down <- data.frame(rbind(allRes_bp_svc_low_high_down, allRes_cc_svc_low_high_down, allRes_mf_svc_low_high_down))
#allRes_svc_low_mid_up    <- data.frame(rbind(allRes_bp_svc_low_mid_up,    allRes_cc_svc_low_mid_up,    allRes_mf_svc_low_mid_up))
#allRes_svc_low_mid_down  <- data.frame(rbind(allRes_bp_svc_low_mid_down,  allRes_cc_svc_low_mid_down,  allRes_mf_svc_low_mid_down))
#allRes_svc_mid_high_up   <- data.frame(rbind(allRes_bp_svc_mid_high_up,   allRes_cc_svc_mid_high_up,   allRes_mf_svc_mid_high_up))
#allRes_svc_mid_high_down <- data.frame(rbind(allRes_bp_svc_mid_high_down, allRes_cc_svc_mid_high_down, allRes_mf_svc_mid_high_down))

#allRes_wvc_low_high_up   <- data.frame(rbind(allRes_bp_wvc_low_high_up,   allRes_cc_wvc_low_high_up,   allRes_mf_wvc_low_high_up))
#allRes_wvc_low_high_down <- data.frame(rbind(allRes_bp_wvc_low_high_down, allRes_cc_wvc_low_high_down, allRes_mf_wvc_low_high_down))
#allRes_wvc_low_mid_up    <- data.frame(rbind(allRes_bp_wvc_low_mid_up,    allRes_cc_wvc_low_mid_up,    allRes_mf_wvc_low_mid_up))
#allRes_wvc_low_mid_down  <- data.frame(rbind(allRes_bp_wvc_low_mid_down,  allRes_cc_wvc_low_mid_down,  allRes_mf_wvc_low_mid_down))
#allRes_wvc_mid_high_up   <- data.frame(rbind(allRes_bp_wvc_mid_high_up,   allRes_cc_wvc_mid_high_up,   allRes_mf_wvc_mid_high_up))
#allRes_wvc_mid_high_down <- data.frame(rbind(allRes_bp_wvc_mid_high_down, allRes_cc_wvc_mid_high_down, allRes_mf_wvc_mid_high_down))
#######################################################AGREGAR CS

myResultsList_vt_up    <- list("Res" = allRes_vt_up,  "BP" = GOdata_bp_vt_up,      "BP_T" = resultT_bp_vt_up,    "BP_KS" = resultKS_bp_vt_up,    "BP_Fisher" = resultFisher_bp_vt_up,  "CC" = GOdata_cc_vt_up,      "CC_T" = resultT_cc_vt_up,    "CC_KS" = resultKS_cc_vt_up,    "CC_Fisher" = resultFisher_cc_vt_up,  "MF" = GOdata_mf_vt_up,      "MF_T" = resultT_mf_vt_up,    "MF_KS" = resultKS_mf_vt_up,    "MF_Fisher" = resultFisher_mf_vt_up)
myResultsList_vt_down  <- list("Res" = allRes_vt_down,  "BP" = GOdata_bp_vt_down,  "BP_T" = resultT_bp_vt_down,  "BP_KS" = resultKS_bp_vt_down,  "BP_Fisher" = resultFisher_bp_vt_down,  "CC" = GOdata_cc_vt_down,  "CC_T" = resultT_cc_vt_down,  "CC_KS" = resultKS_cc_vt_down,  "CC_Fisher" = resultFisher_cc_vt_down,  "MF" = GOdata_mf_vt_down,  "MF_T" = resultT_mf_vt_down,  "MF_KS" = resultKS_mf_vt_down,  "MF_Fisher" = resultFisher_mf_vt_down)


myResultsList_ac_low_high_up  <- list(  "Res" = allRes_ac_low_high_up,  "BP" = GOdata_bp_ac_low_high_up,    "BP_T" = resultT_bp_ac_low_high_up,    "BP_KS" = resultKS_bp_ac_low_high_up,    "BP_Fisher" = resultFisher_bp_ac_low_high_up,"CC" = GOdata_cc_ac_low_high_up,    "CC_T" = NULL,    "CC_KS" = resultKS_cc_ac_low_high_up,    "CC_Fisher" = resultFisher_cc_ac_low_high_up, "MF" = GOdata_mf_ac_low_high_up,    "MF_T" = resultT_mf_ac_low_high_up,    "MF_KS" = resultKS_mf_ac_low_high_up,    "MF_Fisher" = resultFisher_mf_ac_low_high_up)
myResultsList_ac_low_high_down  <- list("Res" = allRes_ac_low_high_down,  "BP" = GOdata_bp_ac_low_high_down,  "BP_T" = resultT_bp_ac_low_high_down,  "BP_KS" = resultKS_bp_ac_low_high_down,  "BP_Fisher" = resultFisher_bp_ac_low_high_down,  "CC" = GOdata_cc_ac_low_high_down,  "CC_T" = NULL,  "CC_KS" = resultKS_cc_ac_low_high_down,  "CC_Fisher" = resultFisher_cc_ac_low_high_down,  "MF" = GOdata_mf_ac_low_high_down,  "MF_T" = resultT_mf_ac_low_high_down,  "MF_KS" = resultKS_mf_ac_low_high_down,  "MF_Fisher" = resultFisher_mf_ac_low_high_down)
######

#myResultsList_ac_low_mid_up  <- list(  "Res" = allRes_ac_low_mid_up, #"BP" = GOdata_bp_ac_low_mid_up,  "BP_T" = resultT_bp_ac_low_mid_up,  "BP_KS" = #resultKS_bp_ac_low_mid_up,  "BP_Fisher" = resultFisher_bp_ac_low_mid_up, #"CC" = GOdata_cc_ac_low_mid_up,  "CC_T" = resultT_cc_ac_low_mid_up,  "CC_KS" = #resultKS_cc_ac_low_mid_up,  "CC_Fisher" = resultFisher_cc_ac_low_mid_up,   #"MF" = GOdata_bp_ac_low_mid_up,  "MF_T" = resultT_mf_ac_low_mid_up,  "MF_KS" = #resultKS_mf_ac_low_mid_up,  "MF_Fisher" = resultFisher_mf_ac_low_mid_up)
#myResultsList_ac_low_mid_down  <- list("Res" = allRes_ac_low_mid_down,  "BP" = GOdata_bp_ac_low_mid_down,  "BP_T" = resultT_bp_ac_low_mid_down,  "BP_KS" = resultKS_bp_ac_low_mid_down,  "BP_Fisher" = resultFisher_bp_ac_low_mid_down, "CC" = GOdata_cc_ac_low_mid_down,  "CC_T" = resultT_cc_ac_low_mid_down,  "CC_KS" = resultKS_cc_ac_low_mid_down,  "CC_Fisher" = resultFisher_cc_ac_low_mid_down, "MF" = GOdata_bp_ac_low_mid_down,  "MF_T" = resultT_mf_ac_low_mid_down,  "MF_KS" = resultKS_mf_ac_low_mid_down,  "MF_Fisher" = resultFisher_mf_ac_low_mid_down)


#myResultsList_ac_mid_high_up  <- list("Res" = allRes_ac_mid_high_up, "BP" = GOdata_bp_ac_mid_high_up,  "BP_T" = resultT_bp_ac_mid_high_up,  "BP_KS" = resultKS_bp_ac_mid_high_up,  "BP_Fisher" = resultFisher_bp_ac_mid_high_up, "CC" = GOdata_cc_ac_mid_high_up,  "CC_T" = resultT_cc_ac_mid_high_up,  "CC_KS" = resultKS_cc_ac_mid_high_up,  "CC_Fisher" = resultFisher_cc_ac_mid_high_up, "MF" = GOdata_bp_ac_mid_high_up,  "MF_T" = resultT_mf_ac_mid_high_up,  "MF_KS" = resultKS_mf_ac_mid_high_up,  "MF_Fisher" = resultFisher_mf_ac_mid_high_up)
#myResultsList_ac_mid_high_down  <- list("Res" = allRes_ac_mid_high_down,  "BP" = GOdata_bp_ac_mid_high_down,  "BP_T" = resultT_bp_ac_mid_high_down,  "BP_KS" = resultKS_bp_ac_mid_high_down,  "BP_Fisher" = resultFisher_bp_ac_mid_high_down, "CC" = GOdata_cc_ac_mid_high_down,  "CC_T" = resultT_cc_ac_mid_high_down,  "CC_KS" = resultKS_cc_ac_mid_high_down,  "CC_Fisher" = resultFisher_cc_ac_mid_high_down, "MF" = GOdata_bp_ac_mid_high_down,  "MF_T" = resultT_mf_ac_mid_high_down,  "MF_KS" = resultKS_mf_ac_mid_high_down,  "MF_Fisher" = resultFisher_mf_ac_mid_high_down)



#myResultsList_cs_lowwin_lowsum_up  <- list("Res" = allRes_cs_lowwin_lowsum_up,  "BP" = GOdata_bp_cs_lowwin_lowsum_up,  "BP_T" = resultT_bp_cs_lowwin_lowsum_up,  "BP_KS" = resultKS_bp_cs_lowwin_lowsum_up,  "BP_Fisher" = resultFisher_bp_cs_lowwin_lowsum_up,"CC" = GOdata_cc_cs_lowwin_lowsum_up,  "CC_T" = resultT_cc_cs_lowwin_lowsum_up,  "CC_KS" = resultKS_cc_cs_lowwin_lowsum_up,  "CC_Fisher" = resultFisher_cc_cs_lowwin_lowsum_up, "MF" = GOdata_bp_cs_lowwin_lowsum_up,  "MF_T" = resultT_mf_cs_lowwin_lowsum_up,  "MF_KS" = resultKS_mf_cs_lowwin_lowsum_up,  "MF_Fisher" = resultFisher_mf_cs_lowwin_lowsum_up)
#myResultsList_cs_lowwin_lowsum_down  <- list("Res" = allRes_cs_lowwin_lowsum_down,  "BP" = GOdata_bp_cs_lowwin_lowsum_down,  "BP_T" = resultT_bp_cs_lowwin_lowsum_down,  "BP_KS" = resultKS_bp_cs_lowwin_lowsum_down,  "BP_Fisher" = resultFisher_bp_cs_lowwin_lowsum_down,  "CC" = GOdata_cc_cs_lowwin_lowsum_down,  "CC_T" = resultT_cc_cs_lowwin_lowsum_down,  "CC_KS" = resultKS_cc_cs_lowwin_lowsum_down,  "CC_Fisher" = resultFisher_cc_cs_lowwin_lowsum_down,  "MF" = GOdata_bp_cs_lowwin_lowsum_down,  "MF_T" = resultT_mf_cs_lowwin_lowsum_down,  "MF_KS" = resultKS_mf_cs_lowwin_lowsum_down,  "MF_Fisher" = resultFisher_mf_cs_lowwin_lowsum_down)


#myResultsList_cs_midwin_midsum_up  <- list("Res" = allRes_cs_midwin_midsum_up,  "BP" = GOdata_bp_cs_midwin_midsum_up,  "BP_T" = resultT_bp_cs_midwin_midsum_up,  "BP_KS" = resultKS_bp_cs_midwin_midsum_up,  "BP_Fisher" = resultFisher_bp_cs_midwin_midsum_up,  "CC" = GOdata_cc_cs_midwin_midsum_up,  "CC_T" = resultT_cc_cs_midwin_midsum_up,  "CC_KS" = resultKS_cc_cs_midwin_midsum_up,  "CC_Fisher" = resultFisher_cc_cs_midwin_midsum_up,  "MF" = GOdata_bp_cs_midwin_midsum_up,  "MF_T" = resultT_mf_cs_midwin_midsum_up,  "MF_KS" = resultKS_mf_cs_midwin_midsum_up,  "MF_Fisher" = resultFisher_mf_cs_midwin_midsum_up)
#myResultsList_cs_midwin_midsum_down  <- list("Res" = allRes_cs_midwin_midsum_down,  "BP" = GOdata_bp_cs_midwin_midsum_down,  "BP_T" = resultT_bp_cs_midwin_midsum_down,  "BP_KS" = resultKS_bp_cs_midwin_midsum_down,  "BP_Fisher" = resultFisher_bp_cs_midwin_midsum_down,  "CC" = GOdata_cc_cs_midwin_midsum_down,  "CC_T" = resultT_cc_cs_midwin_midsum_down,  "CC_KS" = resultKS_cc_cs_midwin_midsum_down,  "CC_Fisher" = resultFisher_cc_cs_midwin_midsum_down, "MF" = GOdata_bp_cs_midwin_midsum_down,  "MF_T" = resultT_mf_cs_midwin_midsum_down,  "MF_KS" = resultKS_mf_cs_midwin_midsum_down,  "MF_Fisher" = resultFisher_mf_cs_midwin_midsum_down)


#myResultsList_cs_highwin_highsum_up  <- list("Res" = allRes_cs_highwin_highsum_up,  "BP" = GOdata_bp_cs_highwin_highsum_up,  "BP_T" = resultT_bp_cs_highwin_highsum_up,  "BP_KS" = resultKS_bp_cs_highwin_highsum_up,  "BP_Fisher" = resultFisher_bp_cs_highwin_highsum_up,"CC" = GOdata_cc_cs_highwin_highsum_up,  "CC_T" = resultT_cc_cs_highwin_highsum_up,  "CC_KS" = resultKS_cc_cs_highwin_highsum_up,  "CC_Fisher" = resultFisher_cc_cs_highwin_highsum_up, "MF" = GOdata_bp_cs_highwin_highsum_up,  "MF_T" = resultT_mf_cs_highwin_highsum_up,  "MF_KS" = resultKS_mf_cs_highwin_highsum_up,  "MF_Fisher" = resultFisher_mf_cs_highwin_highsum_up)
#myResultsList_cs_highwin_highsum_down  <- list("Res" = allRes_cs_highwin_highsum_down, "BP" = GOdata_bp_cs_highwin_highsum_down,  "BP_T" = resultT_bp_cs_highwin_highsum_down,  "BP_KS" = resultKS_bp_cs_highwin_highsum_down,  "BP_Fisher" = resultFisher_bp_cs_highwin_highsum_down,"CC" = GOdata_cc_cs_highwin_highsum_down,  "CC_T" = resultT_cc_cs_highwin_highsum_down,  "CC_KS" = resultKS_cc_cs_highwin_highsum_down,  "CC_Fisher" = resultFisher_cc_cs_highwin_highsum_down,"MF" = GOdata_bp_cs_highwin_highsum_down,  "MF_T" = resultT_mf_cs_highwin_highsum_down,  "MF_KS" = resultKS_mf_cs_highwin_highsum_down,  "MF_Fisher" = resultFisher_mf_cs_highwin_highsum_down)

myResultsListsvc_low_high_up <- list("Res" = allRes_svc_low_high_up, "BP" = GOdata_bp_svc_low_high_up, "BP_T" = resultT_bp_svc_low_high_up, "BP_KS" = resultKS_bp_svc_low_high_up, "BP_Fisher" = resultFisher_bp_svc_low_high_up,"CC" = GOdata_cc_svc_low_high_up, "CC_T" = resultT_cc_svc_low_high_up, "CC_KS" = resultKS_cc_svc_low_high_up, "CC_Fisher" = resultFisher_cc_svc_low_high_up,"MF" = GOdata_bp_svc_low_high_up, "MF_T" = resultT_mf_svc_low_high_up, "MF_KS" = resultKS_mf_svc_low_high_up, "MF_Fisher" = resultFisher_mf_svc_low_high_up)
myResultsListsvc_low_high_down <- list("Res" = allRes_svc_low_high_down, "BP" = GOdata_bp_svc_low_high_down, "BP_T" = resultT_bp_svc_low_high_down, "BP_KS" = resultKS_bp_svc_low_high_down, "BP_Fisher" = resultFisher_bp_svc_low_high_down,"CC" = GOdata_cc_svc_low_high_down, "CC_T" = resultT_cc_svc_low_high_down, "CC_KS" = resultKS_cc_svc_low_high_down, "CC_Fisher" = resultFisher_cc_svc_low_high_down, "MF" = GOdata_bp_svc_low_high_down, "MF_T" = resultT_mf_svc_low_high_down, "MF_KS" = resultKS_mf_svc_low_high_down, "MF_Fisher" = resultFisher_mf_svc_low_high_down)


myResultsList_svc_low_mid_up <- list("Res" = allRes_svc_low_mid_up, "BP" = GOdata_bp_svc_low_mid_up, "BP_T" = resultT_bp_svc_low_mid_up, "BP_KS" = resultKS_bp_svc_low_mid_up, "BP_Fisher" = resultFisher_bp_svc_low_mid_up, "CC" = GOdata_cc_svc_low_mid_up, "CC_T" = resultT_cc_svc_low_mid_up, "CC_KS" = resultKS_cc_svc_low_mid_up, "CC_Fisher" = resultFisher_cc_svc_low_mid_up, "MF" = GOdata_bp_svc_low_mid_up, "MF_T" = resultT_mf_svc_low_mid_up, "MF_KS" = resultKS_mf_svc_low_mid_up, "MF_Fisher" = resultFisher_mf_svc_low_mid_up)
myResultsList_svc_low_mid_down <- list("Res" = allRes_svc_low_mid_down, "BP" = GOdata_bp_svc_low_mid_down, "BP_T" = resultT_bp_svc_low_mid_down, "BP_KS" = resultKS_bp_svc_low_mid_down, "BP_Fisher" = resultFisher_bp_svc_low_mid_down, "CC" = GOdata_cc_svc_low_mid_down, "CC_T" = resultT_cc_svc_low_mid_down, "CC_KS" = resultKS_cc_svc_low_mid_down, "CC_Fisher" = resultFisher_cc_svc_low_mid_down,"MF" = GOdata_bp_svc_low_mid_down, "MF_T" = resultT_mf_svc_low_mid_down, "MF_KS" = resultKS_mf_svc_low_mid_down, "MF_Fisher" = resultFisher_mf_svc_low_mid_down)


myResultsList_svc_mid_high_up <- list("Res" = allRes_svc_mid_high_up, "BP" = GOdata_bp_svc_mid_high_up, "BP_T" = resultT_bp_svc_mid_high_up, "BP_KS" = resultKS_bp_svc_mid_high_up, "BP_Fisher" = resultFisher_bp_svc_mid_high_up, "CC" = GOdata_cc_svc_mid_high_up, "CC_T" = resultT_cc_svc_mid_high_up, "CC_KS" = resultKS_cc_svc_mid_high_up, "CC_Fisher" = resultFisher_cc_svc_mid_high_up,"MF" = GOdata_bp_svc_mid_high_up, "MF_T" = resultT_mf_svc_mid_high_up, "MF_KS" = resultKS_mf_svc_mid_high_up, "MF_Fisher" = resultFisher_mf_svc_mid_high_up)
myResultsList_svc_mid_high_down <- list("Res" = allRes_svc_mid_high_down, "BP" = GOdata_bp_svc_mid_high_down, "BP_T" = resultT_bp_svc_mid_high_down, "BP_KS" = resultKS_bp_svc_mid_high_down, "BP_Fisher" = resultFisher_bp_svc_mid_high_down,"CC" = GOdata_cc_svc_mid_high_down, "CC_T" = resultT_cc_svc_mid_high_down, "CC_KS" = resultKS_cc_svc_mid_high_down, "CC_Fisher" = resultFisher_cc_svc_mid_high_down, "MF" = GOdata_bp_svc_mid_high_down, "MF_T" = resultT_mf_svc_mid_high_down, "MF_KS" = resultKS_mf_svc_mid_high_down, "MF_Fisher" = resultFisher_mf_svc_mid_high_down)


myResultsList_wvc_low_high_up <- list("Res" = allRes_wvc_low_high_up, "BP" = GOdata_bp_wvc_low_high_up, "BP_T" = resultT_bp_wvc_low_high_up, "BP_KS" = resultKS_bp_wvc_low_high_up, "BP_Fisher" = resultFisher_bp_wvc_low_high_up,"CC" = GOdata_cc_wvc_low_high_up, "CC_T" = resultT_cc_wvc_low_high_up, "CC_KS" = resultKS_cc_wvc_low_high_up, "CC_Fisher" = resultFisher_cc_wvc_low_high_up,"MF" = GOdata_bp_wvc_low_high_up, "MF_T" = resultT_mf_wvc_low_high_up, "MF_KS" = resultKS_mf_wvc_low_high_up, "MF_Fisher" = resultFisher_mf_wvc_low_high_up)
myResultsList_wvc_low_high_down <- list("Res" = allRes_wvc_low_high_down, "BP" = GOdata_bp_wvc_low_high_down, "BP_T" = resultT_bp_wvc_low_high_down, "BP_KS" = resultKS_bp_wvc_low_high_down, "BP_Fisher" = resultFisher_bp_wvc_low_high_down, "CC" = GOdata_cc_wvc_low_high_down, "CC_T" = resultT_cc_wvc_low_high_down, "CC_KS" = resultKS_cc_wvc_low_high_down, "CC_Fisher" = resultFisher_cc_wvc_low_high_down, "MF" = GOdata_bp_wvc_low_high_down, "MF_T" = resultT_mf_wvc_low_high_down, "MF_KS" = resultKS_mf_wvc_low_high_down, "MF_Fisher" = resultFisher_mf_wvc_low_high_down)


myResultsList_wvc_low_mid_up <- list("Res" = allRes_wvc_low_mid_up, "BP" = GOdata_bp_wvc_low_mid_up, "BP_T" = resultT_bp_wvc_low_mid_up, "BP_KS" = resultKS_bp_wvc_low_mid_up, "BP_Fisher" = resultFisher_bp_wvc_low_mid_up,"CC" = GOdata_cc_wvc_low_mid_up, "CC_T" = resultT_cc_wvc_low_mid_up, "CC_KS" = resultKS_cc_wvc_low_mid_up, "CC_Fisher" = resultFisher_cc_wvc_low_mid_up, "MF" = GOdata_bp_wvc_low_mid_up, "MF_T" = resultT_mf_wvc_low_mid_up, "MF_KS" = resultKS_mf_wvc_low_mid_up, "MF_Fisher" = resultFisher_mf_wvc_low_mid_up)

myResultsList_wvc_low_mid_down <- list("Res" = allRes_wvc_low_mid_down, "BP" = GOdata_bp_wvc_low_mid_down, "BP_T" = resultT_bp_wvc_low_mid_down, "BP_KS" = resultKS_bp_wvc_low_mid_down, "BP_Fisher" = resultFisher_bp_wvc_low_mid_down, "CC" = GOdata_cc_wvc_low_mid_down, "CC_T" = resultT_cc_wvc_low_mid_down, "CC_KS" = resultKS_cc_wvc_low_mid_down, "CC_Fisher" = resultFisher_cc_wvc_low_mid_down, "MF" = GOdata_bp_wvc_low_mid_down, "MF_T" = resultT_mf_wvc_low_mid_down, "MF_KS" = resultKS_mf_wvc_low_mid_down, "MF_Fisher" = resultFisher_mf_wvc_low_mid_down)


myResultsList_wvc_mid_high_up <- list("Res" = allRes_wvc_mid_high_up, "BP" = GOdata_bp_wvc_mid_high_up, "BP_T" = resultT_bp_wvc_mid_high_up, "BP_KS" = resultKS_bp_wvc_mid_high_up, "BP_Fisher" = resultFisher_bp_wvc_mid_high_up,"CC" = GOdata_cc_wvc_mid_high_up, "CC_T" = resultT_cc_wvc_mid_high_up, "CC_KS" = resultKS_cc_wvc_mid_high_up, "CC_Fisher" = resultFisher_cc_wvc_mid_high_up, "MF" = GOdata_bp_wvc_mid_high_up, "MF_T" = resultT_mf_wvc_mid_high_up, "MF_KS" = resultKS_mf_wvc_mid_high_up, "MF_Fisher" = resultFisher_mf_wvc_mid_high_up)
myResultsList_wvc_mid_high_down <- list("Res" = allRes_wvc_mid_high_down, "BP" = GOdata_bp_wvc_mid_high_down, "BP_T" = resultT_bp_wvc_mid_high_down, "BP_KS" = resultKS_bp_wvc_mid_high_down, "BP_Fisher" = resultFisher_bp_wvc_mid_high_down,"CC" = GOdata_cc_wvc_mid_high_down, "CC_T" = resultT_cc_wvc_mid_high_down, "CC_KS" = resultKS_cc_wvc_mid_high_down, "CC_Fisher" = resultFisher_cc_wvc_mid_high_down, "MF" = GOdata_bp_wvc_mid_high_down, "MF_T" = resultT_mf_wvc_mid_high_down, "MF_KS" = resultKS_mf_wvc_mid_high_down, "MF_Fisher" = resultFisher_mf_wvc_mid_high_down)

write.table(myResultsList_vt_up$Res, file = "sym_myResultsList_vt_up",sep='\t')
write.table(myResultsList_vt_down$Res, file = "sym_myResultsList_vt_down",sep='\t')
write.table(myResultsList_ac_low_high_up$Res, file = "sym_myResultsList_high_low_up",sep='\t')
write.table(myResultsList_ac_low_high_down$Res, file = "sym_myResultsList_high_low_down",sep='\t')



#####Graph structure####
showSigOfNodes(GOdata_bp_vt_up, score(resultKS_bp_vt_up), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_bp_vt_down, score(resultFisher_bp_vt_down), firstSigNodes = 10, useInfo = 'all')

showSigOfNodes(GOdata_cc_vt_up, score(resultFisher_cc_vt_up), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_cc_vt_down, score(resultFisher_cc_vt_down), firstSigNodes = 10, useInfo = 'all')

showSigOfNodes(GOdata_mf_vt_up, score(resultFisher_mf_vt_up), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_mf_vt_down, score(resultFisher_mf_vt_up), firstSigNodes = 10, useInfo = 'all')

showSigOfNodes(GOdata_bp_ac_low_high_up, score(resultFisher_bp_ac_low_high_up), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_bp_ac_low_high_down, score(resultFisher_bp_ac_low_high_down), firstSigNodes = 10, useInfo = 'all')

showSigOfNodes(GOdata_cc_ac_low_high_up, score(resultFisher_cc_ac_low_high_up), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_cc_ac_low_high_down, score(resultFisher_cc_ac_low_high_down), firstSigNodes = 10, useInfo = 'all')

showSigOfNodes(GOdata_mf_ac_low_high_up, score(resultFisher_mf_ac_low_high_up), firstSigNodes = 10, useInfo = 'all')
showSigOfNodes(GOdata_mf_ac_low_high_down, score(resultFisher_mf_ac_low_high_down), firstSigNodes = 10, useInfo = 'all')

