
p<-dim(G.final)[1]
name<-colnames(y)
cutoff=0.03

adjacency<-diag(1,p,p)
for (i in 2:(p-1)) {
  for (j in 1:(i-1)) {
    ifelse(abs(G.final[i,j])<=cutoff & abs(G.final[j,i])<=cutoff &
           abs(H1.final[i,j])<=cutoff & abs(H1.final[j,i])<=cutoff &
           abs(H2.final[i,j])<=cutoff & abs(H2.final[j,i])<=cutoff &
           abs(K.final[i,j])<=cutoff & abs(K.final[j,i])<=cutoff, 
           adjacency[i,j]<-adjacency[j,i]<-NA,
           adjacency[i,j]<-adjacency[j,i]<-1)
  }
} 

count.na=sum(is.na(adjacency))
p.nzero=1-count.na/(p*(p-1))
p.nzero


library("WGCNA")

TOM = TOMsimilarity(adjacency)
tom_connectivity = apply(TOM,1,sum)
dissTOM = 1-TOM
colnames(dissTOM) = name

NEITree = hclust(as.dist(dissTOM), method = "average")

minModuleSize = 5
dynamicMods = cutreeDynamic(dendro = NEITree, distM = dissTOM, method = "hybrid", deepSplit = 4, pamRespectsDendro = TRUE, minClusterSize = minModuleSize)
table(dynamicMods)

plotDendroAndColors(NEITree,colors=dynamicMods)
plotTOM = dissTOM
diag(plotTOM) = NA
TOMplot(plotTOM, NEITree, Colors =dynamicMods,
        main = "Network Heatmap Plot")
connectivity_matrix<-intramodularConnectivity(adjacency, dynamicMods)
ID = seq(1,nrow(TOM))
label = name
weight = tom_connectivity
module=dynamicMods
connectivity<-connectivity_matrix$kTotal
intra_connectivity<-connectivity_matrix$kWithin
net<-cbind(ID,label,module,weight,connectivity,intra_connectivity)
net=as.data.frame(net)
write.csv(net,file = paste("vinfo.csv", sep=""),fileEncoding = "GBK")

edges = as.data.frame(t(combn(nrow(adjacency),2)))
edges$weight=NA

for (i in 1:nrow(edges)) {
  a=edges[i,1]
  b=edges[i,2]
  edges[i,3]=adjacency[a,b]
}

edges=edges[-c(which(is.na(edges$weight))),]
colnames(edges)[1:2]<-c("source","target")
write.csv(edges,file = paste("einfo.csv", sep=""),fileEncoding = "GBK")





