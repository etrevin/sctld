####WGCNA####
####1.Dendro####
library(WGCNA);
library(dynamicTreeCut);
#setwd("~/Documents/DGE"); #Establecer directorio
#options(stringsAsFactors = FALSE);       #tratar las datos todos los datos como strings
#ppanvst = read.csv("~/Documents/DGE/symvst_filter_mad_0.5_8485.csv", row.names=1); #lectura de la base de datos de expresion y creación de un objeto en R
ppanvst=vsd2
ppanvst$mad=NULL
#ppanvst$k559.tsv=NULL
#ppanvst$k559.tsv=NULL

datExpr0 = as.data.frame(t(ppanvst));           #Crea una matriz con formato sin la primera columna de ppanvst
names(datExpr0) = rownames(ppanvst);                             #extrae la columna de los nombres y los agrega al vector data
rownames(datExpr0) = names(ppanvst);              #define como nombres de columna los nombres de la primera columna de ppanvst
gsg = goodSamplesGenes(datExpr0, verbose = 3);           #prueba para retirar transcripts con valores de "0"
gsg$allOK                                                #Compureba cada row para determinar si existen valores "0" como False
if (!gsg$allOK)
{
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));          #Optionally, print the gene names that were removed
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));   #Optionally, print the sample names that were removed
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]    #Remove the offending genes and samples from the data
}
sampleTree = hclust(dist(datExpr0), method = "ward.D2"); #Hierarchical Clustering de los datos datEprx0
x = as.matrix(dist(datExpr0, method = "euclidian"))
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
abline(h = 330, col = "red");                           #Linea de corte en el plot que se realizará posteriormente  #concluye con el proceso del PDF
clust = cutreeDynamic(sampleTree, method = "hybrid", distM = x, minClusterSize = 1)
#clust = cutreeDynamicTree(sampleTree, minModuleSize = 1)
#clust = cutreeStatic(sampleTree, cutHeight = 200, minSize = 1)                                             #Altura del cluster donde se realizará el corte para los modulos con tamaño mínimo de 1
table(clust)
keepSamples = (clust>0)                                 #mantener las muestras segregadas por el corte
datExpr = datExpr0[keepSamples, ]                        #Crea un objeto a partir de las muestras del corte con todos los transcritos de la primera base de datos
datExpr = datExpr0                       #Crea un objeto a partir de las muestras del corte con todos los transcritos de la primera base de datos
nGenes = ncol(datExpr)                                   #crea un objeto con el número de columnas de datExpr
nSamples = nrow(datExpr)                                 #crea un objeto con el número de transcritos de datExpr
ventraits <- read.csv("~/Documents/DGE/ventraits.csv")
allTraits=ventraits
dim(allTraits)                                           #Dimensiones (Col y Rows) de la base de datos
names(allTraits)                                         #Obtiene los nombres de los factores
pporamples = rownames(datExpr);                          #Extrae los nombres de las muestras del corte datExpr
traitRows = match(pporamples, allTraits$sample);         #Empareja los nombres de datExpr con la base de datos del ambiente
datTraits = allTraits[traitRows, -1];                    #extrae los valores ambientales de ciertas librerías unicamente
rownames(datTraits) = allTraits[traitRows, 1];                            
collectGarbage();                                        #Performs garbage collection until free memory idicators show no change.
sampleTree2 = hclust(dist(datExpr), method = "average")  #Segundo dendogram con base a la base de datos procesada 
traitColors = numbers2colors(datTraits, signed = FALSE); #Paleta de colores
plotDendroAndColors(sampleTree2, traitColors, groupLabels = names(datTraits), main = "Sample dendrogram and trait heatmap") #Dendograma con representación por colores de los factores ambientales con el dendograma:blanco->poco, rojo->mucho, gris->sin informacion
#save(datExpr, datTraits, file = "/LUSTRE/usuario/etrevino/r_stuff/wgcna/pordatainput.RData")   #Guarda el ambiente de trabajo de los datos de expreción y los traits

allowWGCNAThreads() 
#enableWGCNAThreads(4)                                     #permite el uso de multiples nucleos
#lnames1 = load(file = "/LUSTRE/usuario/etrevino/r_stuff/wgcna/pordatainput.RData");            #extrae el nombre de la bases de datos
#lnames1                                                   #muestra los nombres
powers = c(c(1:10), seq(from = 12, to=30, by=2))         #Define el umbral de la funcion "power" para definir los modulos
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)                         #Análisis de topología de redes tipo "weigthed" para definir "power" más óptimo
#save(sft, file = "sft.Rdata")
#load("/LUSTRE/usuario/etrevino/r_stuff/wgcna/sft.Rdata");
sizeGrWindow(9, 5);                                       # open a graphics window
par(mfrow = c(1,2));                                     #Características del plot
cex1 = 0.9;                                              #Características del plot
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n", main = paste("Scale independence")); #scale-free fit index en función del power
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers,cex=cex1,col="red");   
abline(h=0.90,col="red")                                 #Linea de corte con un R2 al 90%
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))  # mean connectivity en función del power
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
gc()
net = blockwiseModules(datExpr, power = 8,TOMType = "signed Nowick 2", corType ="pearson", networkType = "signed", deepsplit = 2 , minModuleSize = 75, numericLabels = FALSE, pamStage = FALSE,pamRespectsDendro = FALSE, saveTOMs = TRUE, saveTOMFileBase = "pporTOM", verbose = 5)  #construcción de la red de genes e identificación de módulos; Power se selecciona de acuerdo al R2=90% en la escala de independencia; minModuleSize tamaño mínimo del modulo en genes; 
table(net$colors)                                       #Tamaños de los modulos
sizeGrWindow(12, 9)                                     # open a graphics window
mergedColors = labels2colors(net$colors)                #Convierte las etiquetas a colores antes del plotting
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]], "Module colors", dendroLabels = FALSE, hang = 0.03, addGuide = TRUE, guideHang = 0.05) # Plot the dendrogram and the module colors underneath
moduleLabels = net$colors                               #objeto con las etiquetas de los módulos

write.table(moduleLabels, file = "Por_module_col",sep='\t')


moduleColors = labels2colors(net$colors)                #objeto con los colores de los modulos



MEs = net$MEs;                                          #Module eigengenes

write.table(MEs, file = "Por_MEs",sep='\t')


geneTree = net$dendrograms[[1]];                        #Objeto con el dendograma de los módulos
#save(MEs, moduleLabels, moduleColors, geneTree, file = "/LUSTRE/usuario/etrevino/r_stuff/wgcna/por-networkConstruction-auto.RData")              #Guarda los objetos con los módulos en el ambiente global 


####2.Module Epr-heatmaps####
colorh1=moduleColors
datME=moduleEigengenes(datExpr, colorh1)$eigengenes
signif(cor(datME, use="p"), 2)
dissimME=(1-t(cor(datME, method="p")))/2
hclustdatME=hclust(as.dist(dissimME), method="average" )
# Plot the eigengene dendrogram
par(mfrow=c(1,1))
plot(hclustdatME, main="Clustering tree based of the module eigengenes")
sizeGrWindow(8,9)
plotMEpairs(datME)

sizeGrWindow(8,9)
par(mfrow=c(3,1), mar=c(1, 2, 4, 1))
which.module="red";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
# for the second (blue) module we use
which.module="blue";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
which.module="brown";
plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=T,
        clabels=T,rcols=which.module,
        title=which.module )
sizeGrWindow(8,7);

which.module="black"
ME=datME[, paste("ME",which.module, sep="")]
#par(mfrow=c(2,1), mar=c(0, 5.5, 3, 2))
#plotMat(t(scale(datExpr[,colorh1==which.module ]) ),nrgcols=30,rlabels=F, clabels=T,rcols=which.module, main=which.module, cex.main=2)
#par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", cex.main=2,
        ylab="eigengene expression",xlab="array sample")

####3.MODULE-TRAIT relationship heatmap####
#lnames1 = load(file = "symdatainput.RData");            #Carga los datos de expreción y los traits
#lnames2 = load(file = "sym-networkConstruction-auto.RData");            #Carga los datos con los módulos
nGenes = ncol(datExpr);                                  #Objeto con el número de genes
nSamples = nrow(datExpr);                                #Objeto con el número muestras
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes                                                   #Calcula los Eigengenes
MEs = orderMEs(MEs0)                                     #Ordena los Eigengenes de acuerdo a la relevancia del componente
moduleTraitCor = cor(MEs, datTraits, use = "p");         #Correlaciones entre module-Eigengenes y trais externos
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);                                             #P valor de las correlaciones de Student
sizeGrWindow(10,6)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "");          #Muestra las correlaciones y los P valores
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(3, 8.5, 1, 1));                                       #Genera un archivo tipo PDF 
labeledHeatmap.multiPage(Matrix = moduleTraitCor, xLabels = names(datTraits), yLabels = names(MEs), ySymbols = names(MEs), 
      colorLabels = TRUE, colors = blueWhiteRed(100, gamma = 1, endSaturation = 1), textMatrix = textMatrix, setStdMargins = FALSE, cex.text = 0.8, zlim = c(-1,1),
      main = paste("Module-trait relationships")) # valores de la correlacion en un heatmap plot
 
dev.off()
####4.GS-MSmodule-trait Relationship####
# Define variable weight containing the weight column of datTrait
weight = as.data.frame(datTraits$salinity);
names(weight) = "weight"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, weight, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(weight), sep="");
names(GSPvalue) = paste("p.GS.", names(weight), sep="");

module = "blue"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for tempcol",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)
####5.TOMPLOT & expression patters####
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(datExpr, power = 6);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")
nSelect = 400
# For reproducibility, we set the random seed
set.seed(10);
select = sample(nGenes, size = nSelect);
selectTOM = dissTOM[select, select];
# There’s no simple way of restricting a clustering tree to a subset of genes, so we must re-cluster.
selectTree = hclust(as.dist(selectTOM), method = "average")
selectColors = moduleColors[select];
# Open a graphical window
sizeGrWindow(9,9)
# Taking the dissimilarity to a power, say 10, makes the plot more informative by effectively changing
# the color palette; setting the diagonal to NA also improves the clarity of the plot
plotDiss = selectTOM^7;
diag(plotDiss) = NA;
TOMplot(plotDiss, selectTree, selectColors, main = "Network heatmap plot, selected genes")
# Recalculate module eigengenes
MEs = moduleEigengenes(datExpr, moduleColors)$eigengenes
##### Isolate weight from the clinical traits#####
weight = as.data.frame(datTraits$weight_g);
names(weight) = "weight"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, weight))
# Plot the relationships among the eigengenes and the trait
sizeGrWindow(5,7.5);
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
####EXTRA_???####

names(datExpr)
names(datExpr)[moduleColors=="turquoise"]
annot = read.csv(file = "GeneAnnotation.csv"); dim(annot)
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$substanceBXH)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

annot <- read.delim("~/Documents/DGE/porkal/por_trinotate_annotation_report.xls", header=FALSE, comment.char="#")
names(annot)
probes = names(datExpr)
probes2annot = match(probes, annot$V1)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0.

# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       geneSymbol = NA,
                       LocusLinkID = NA,
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, weight, use = "p")));
# Add module membership information in the chosen order for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.weight)); geneInfo = geneInfo0[geneOrder, ]



