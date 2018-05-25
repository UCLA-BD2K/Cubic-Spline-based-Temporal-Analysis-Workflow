#####################################################################################################
### Implementing Cross Validation to select value of λ and Implement Smoothing Splines
#####################################################################################################
MB3ratio4 <- MB3ratio3
for (i in 1:nrow(MB3ratio4)){
  s2 <- smooth.spline(t(rbind(c(1,3,5,7,10,14),MB3ratio4[i,])), cv = TRUE)
  MB3ratio4[i,] <- predict(s2, c(1,3,5,7,10,14))$y
}

#####################################################################################################
### Perform kmean clustering
#####################################################################################################
km1 <- kmeans(MB3ratio4,centers = 8)$cluster
MB3ratio4$cluster <- km1
MB3ratio4$cys <- MB3ratio2$Group.1
MB3ratio5 <- MB3ratio4[order(MB3ratio4$cluster),]
colors <- c(seq(-2,-0.5,length=50),seq(-0.49,0.49,length=50),seq(0.5,2,length=50))
my_palette <- colorRampPalette(c("black", "grey", "white"))(n = 149)

#####################################################################################################
### Display a result in Heatmap
#####################################################################################################
#install.packages("gplots")
library(gplots)
pdf("Result/Figure 3/4-2-18 Ratio and mean heatmap significant black white.pdf",width = 7, height = 7)
heatmap.2(as.matrix(MB3ratio5[,-c(7:8)]), col=my_palette, scale="none", breaks = colors, symbreaks=TRUE, 
          Rowv = "none", Colv = "none", dendrogram = "none",
          trace="none", margin=c(3,1), density.info= "none",
          #Colv = dendcompletemc, 
          keysize=1, 
          #( "bottom.margin", "left.margin", "top.margin", "left.margin" )
          key.par=list(mar=c(3.5,0,3,0)),
          # lmat -- added 2 lattice sections (5 and 6) for padding
          lmat=rbind(c(5, 4, 2), c(6, 1, 3)), lhei=c(2.5, 5), lwid=c(1, 10, 1))
dev.off()
write.csv(MB3ratio5,"Result/Figure 3/4-2-18 Ratio and mean spline smoothing.csv",row.names = F)
