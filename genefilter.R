vsd2=(assay(vsd))
vsd2=as.data.frame(vsd2)

#keep <- rowSums(vsd) >10
#vsd <- vsd[keep,]
#keep <- rowSums(vsd2 >=10) >= 3
#vsd2 <- vsd2[keep,]
#x2 <- vsd2

vsd2$mad <-apply(vsd2,1,mad)
keep <- (vsd2$mad) >= .05
vsd2 <- vsd2[keep,]
View(as.data.frame(vsd2))

vsd2$variance <- NULL
vsd=log2(vsd+1)
write.csv(vsd2,"por_mad2_log2.csv")

vsd3 <-vst(vsd2) 

#SYM

# mad >= 0.5 8485
# mad >= 1.0 1498

#host

# mad >= 0.5 
# mad >= 1.0 

vsd$variance = apply(vsd, 1, var)
vsd2 = vsd[vsd$variance >= quantile(vsd$variance, c(.50)), ]

vsd2<-varFilter(assay(vsd), var.func=IQR, var.cutoff=0.5, filterByQuantile=TRUE)
vsd2<-as.data.frame(vsd2)
