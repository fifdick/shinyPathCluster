#library(networkly)
library(plotly)




##############################################################

######################FUNCTIONS##################################
get_gene_info<-function(pathwayAnalyses,db)
{
	genes<-pathwayAnalyses[[db]]$genesets
	genes<-do.call(rbind,genes)
	genes<-subset(genes,!duplicated(gene_id) & !duplicated(gene_name))
#	print(sum(duplicated(genes)))
#	print(which(duplicated(genes)))
	return(genes)
}
#############################################################################
#############################################################################

get_GP_mat<-function(pathwayAnalyses,db,parsed)
{
library(PathCluster)
#make Occ matrix
#NEED TO TAKE THIS OUT OF HERE BECAUSE IT IS CALCULATED FOR EACH CLUSTER AGAIN...DOESNT MAKE SENS ITS A WASTE OF TIME
genesets<-parsed$geneset
gp <- gene_pathway_matrix(genesets)
print("idx of duplicates in global gp matrix")
print(which(duplicated(rownames(gp))))
#exclude all non significnt genes
info<-get_gene_info(pathwayAnalyses,db)
info<-subset(info,pvalue<0.05)
#print(info)
#print(head(rownames(gp)))
gp<-gp[info$gene_name,]
print("idx of duplicates in gloval gp after subsetting for significant genes")
print(which(duplicated(rownames(gp))))
return(gp)
}

#############################################################################
#############################################################################

modify_GP_mat <- function(gp,pathwayList) {

#subset by only pathways that are in the cluster
cat("dim of gene-pathway mat:",dim(gp),"\n")
#print(length(pathwayList))
cat("pathwayList: ",pathwayList,"\n")
print(pathwayList %in% colnames(gp))
gp_subset <- subset(gp,select=pathwayList)
#cat("dim of gp_mat subset by items in cluster: ",dim(gp_subset),"\n")

#discard all genes with less than 1 occurance in row
#by using apply on each row in matrix and taking its sum, returning indices whose some is greater than one and subseting rownames of the matrix by those indices
if(ncol(gp_subset)==1)
{#we normally only want to plot genes that are occuring in at least 2 pathways, but if a cluster only exist of on pathway and we want to plot that pathway we want of corse all genes in that pathway
genes<-rownames(gp_subset)[which(apply(gp_subset,1,sum)>0) ]
print("single pathway genes:")
print(head(genes))
print(sum(duplicated(genes)))
}
else
{
genes<-rownames(gp_subset)[which(apply(gp_subset,1,sum)>1) ]
print("numer of overlapping genes between pathways:")
print(length(genes))
print(head(genes))
}
print("dim of gp_subset")
print(dim(gp_subset))
print("head of gp_subset")
print(head(gp_subset))
gp_mat <- gp_subset[genes,]
if(length(genes)==1)
{

gp_mat<-t(gp_mat)
gp_mat<-as.data.frame(gp_mat)
colnames(gp_mat)<-pathwayList
rownames(gp_mat)<-genes
print("gp_mat after subsetting gp_subset by genes")
cat("genes:",genes,"\n")
#cat("gp_mat: \n", gp_mat, "\n")
#cat("dim gp_mat: ",dim(gp_mat),"\n")
cat("type: ",typeof(gp_mat),"\n")
}
else
{
gp_mat<-as.data.frame(gp_mat)
print("same thing next two lines")
print(dim(gp_mat))
print(gp_mat)
print("###")
rownames(gp_mat)<-genes
print(head(gp_mat))
colnames(gp_mat)<-pathwayList
cat("Number of genes in net:",length(genes),"\n")
}
#print(dim(gp_mat))
return(gp_mat)
}

#MAKE
#### edge list
#############################################################################
#############################################################################

#TODO check if colname is there when iterating
get_edges <- function(occmat)
{
abschnitte <- lapply(seq(1,ncol(occmat)),function (idx)
       {       #col = pathway
	       col<-occmat[,idx]  
	       #number of genes in this pathway equals the numbe rof occurances of this pathway as a node in the source list of the edge list
	       no_genes <- sum(col)
		cat("in func get_edges(), print no_genes ",no_genes,"\n")
	       source <- rep(as.character(colnames(occmat)[idx]),no_genes)
		cat("edgelist head source: ",head(source),"\n")	
	       #this pathway has edges to all its genes
		print("rownames from this modified occmat")
		print(rownames(occmat))	       
               target <- rownames(occmat)[col==1]
	       cat("edgelist, head target: ", head(target),"\n")
	       abschnitt <- data.frame("source"=source,"target"=target)
	      # print(head(abschnitt))
	       return(abschnitt)
	})
edge_list <- do.call(rbind,abschnitte)
#cat("after rbind: ",print(head(edge_list)),"\n")
colnames(edge_list)<-c("from","to")
edge_list$color<-rep("grey",nrow(edge_list))
#edge_list$color.highlight<-rep("pink",nrow(edge_list))
return(edge_list)
}

#############################################################################
#############################################################################

get_nodes <- function(edge_list,pvalues,genesets,ensembl,func,external=FALSE)
{

library(visNetwork)
#collect all pathways as nodes 
sources<-unique(as.character(edge_list[,1]))
print("#########DEUBUG############")
print(length(sources))
# collect all genes as nodes
print("checking if name was last before get nodes function:")
print("printing whole edge_list to")
print(edge_list[,2])
print(edge_list$to)
targets <- unique(as.character(edge_list[,2]))
print(length(targets))
#modified to only have target labels and pathway names as title in tooltip hover
#labels <- c(sources,targets)
labels<-c(sources,targets)
# use raw geneset obj to find out fold change value for each gene node
genesets <- do.call(rbind,genesets)
#there are of course multiple duplicates comng from different pathways
#it doesnt matter from which row to select the values because they are the same everywhere
genesets <- subset(genesets, !duplicated(gene_name)) 
print("func getnodes(), printing head of genesets after removal of dupliactes")
print(head(genesets))
print("func getnodes(),printing targets")
print(head(targets))
if(external ==TRUE)
{
	
	total_len <- length(sources)+length(targets)
	gene_pvals <- sapply(targets,function(x){
		return(genesets[genesets$gene_name==x,"pvalue"])
	})
	pvalue<-c(pvalues[sources],gene_pvals)
	groups <- c(rep("pathway",length(sources)),rep("gene",length(targets)))
	foldChange <- sapply(targets,function(x){
		return(genesets[genesets$gene_name==x,"log2FoldChange"])
	})
	foldChange <- c(rep(NA,length(sources)),foldChange)
	gene_ids <-  sapply(targets,function(x){
		return(genesets[genesets$gene_name==x,"gene_id"])
	})
	nodes<-data.frame("name"=labels,"ids"=c(rep(NA,length(sources)),gene_ids),"log2FoldChange"=foldChange,"type"=groups,"pvalue"=pvalue)

}
#if fold change is positive, return green
else
{
	gene_col <- sapply(targets,function(x){
#print(x)
		ifelse(genesets[genesets$gene_name==x,"log2FoldChange"]>0,"green","red")
	})
	print("head gene_col")
	print(head(gene_col))
	#print("table of gene_col")
	#print(table(gene_col))
	#pathways all have the same colour
	colours <- c(rep("darkblue",length(sources)),gene_col)
	print("head colours")
	print(head(colours))
	total_len <- length(sources)+length(targets)
	#pathways should be entities and have a squared shape, genes circles
	shapes <- c(rep("diamond",length(sources)),rep("dot",length(targets)))
	#values set the size of the shape, all genes having same size pathways sizes correspond to their pvalue
	print(pvalues)
	print(sources)
	pvals<-pvalues[sources]
	#check what is the smalles value and make all zero pvalues even smaller
	#then negate all pvalues to flip the value such that bigger shapes indicate more significant pvals
	#cat("in func get_nodes(), print pvals per cluster: ",pvals,"\n")
	if(length(pvals[pvals!=0])==0)
	{
	min <- 0.0001
	}
	else
	{
	min<- min(na.omit(pvals[pvals!=0]))
	}
	#adding a pseudocount for zero values
	#taking min pvalue and adding a 0 behind the comma to make the 0 pvalues not zero but still smaller than the smallest
	pseudocount <- min/10
	pvals <- sapply(pvals,function(x){ifelse(x==0,pseudocount,x)})
	pvals<- -(pvals)
	pvals<- pvals+2
	connectivity_table<-table(edge_list$to)
	print("connectivity")
	print(head(connectivity_table))
	gene_shape_size <-sapply(targets,function(t){
		size<-connectivity_table[as.character(t)]
		})
	values <- c(pvals,gene_shape_size)#rep(gene_shape_size,length(targets )))
	groups <- c(rep("pathway",length(sources)),rep("gene",length(targets)))
	nlevels <- sapply(seq(1,total_len),function(i) {
		if(groups[i]=="gene")
		{ 
			if(colours[i]=="red")
			{
				if((i %% 2)==0)
					level<-1
				else
					level<-2
			}
			else if(colours[i]=="green")
			{
				if((i %% 2)==0)
					level<-5
				else
					level<-6

			}
			else
			{print("colour of failing gene which causes level bug:")
			print(colours[i])}
		}
		else
		{level<-4}
		return(level)
		})
	print("Head of levels for cluster")
	print(head(nlevels))

	# first get ensemblidi
	ensembl<-as.data.table(ensembl)
	cat("dim targets",length(targets),"\n")
	targets_as_ensembl<-sapply(targets,function(t){
	id<-ensembl[ gene_name == t]$gene_id
	return(id)
	})
	cat("dim targetsasensembl",length(targets_as_ensembl),"\n")
	# get funciton
	func<-as.data.table(func)
	target_titles<-lapply(targets_as_ensembl,function(eid){
	des<-func[ ensembl_gene_id == eid]$description
	print(des)
	#type<-func[ensembl_gene_id == eid]$gene_biotype
	link<-paste0("http://www.ensembl.org/id/",eid)
	print(link)
	return(list("d"=des,"l"=link))
	})
	cat("print dim target_titles:",length(target_titles),"\n")

	titles<-c(sources,target_titles)
	titles<-sapply(titles,function(t){
	if(is.list(t))
	{
	#<p><b>1</b><br>Node !</p> 
	str<-paste0("<p><b>", t$d,"</b><br>","<a href=\"",t$l,"\" target=\"_blank\">",t$l,"</a></p>")
	print(str)
	}
	else
	{	
	str<-paste0("<p><b>", t,"</b></p>")
	}
	return(str)
	})
	print("################DEEEEBUUGGGGG#################")
	print(length(sources))
	print(length(targets))
	ids<-seq(1,total_len)
	cat("print dim titles:",length(titles),"\n")
	nodes<-data.frame(id=ids,label=labels,title=titles,level=nlevels,shape=shapes,value=values,color=colours,group=groups,font.color=rep("black",total_len))

}

return(nodes)


}
#############################################################################
#############################################################################

rekey_labels <- function(edge_list,nodes)
{
#that should work since nodes are unique in node list
rownames(nodes) <- nodes$label
#source label is already idi
edge_list_raw<-edge_list$from
edge_list$to<- sapply(edge_list$to, function(t) {
	#for each source look up its id in nodes list
		id<-nodes[as.character(t),]$id
		return(id)})
edge_list$from<- sapply(edge_list$from, function(s) {
	#for each source look up its id in nodes list
		id<-nodes[as.character(s),]$id
		return(id)})

len_sources<-length(unique(edge_list$from))
head(nodes$label)
nodes$label<-sapply(seq(1,length(nodes$label)),function(i){
if(as.character(nodes$label[i]) %in% as.character(unique(edge_list_raw)))
{
label<-i
print(label)
}
else
{
label<-nodes$label[i]
print("ELSE")
print(label)
}
return(as.character(label))
})
return(list("edge_list"=edge_list,"nodes"=nodes))
}
#############################################################################
#############################################################################
get_func_annot<-function(ensembl_ids)
{
library(biomaRt)
library(EnsDb.Hsapiens.v75)
print(sum(duplicated(ensembl_ids)))
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
function_map = getBM(attributes = c('ensembl_gene_id', "description"), 
	                    filters = 'ensembl_gene_id', 
			                  values = ensembl_ids, 
			                  mart = ensembl)
rownames(function_map)<-as.character(function_map[,1])
function_map$description<-sapply(function_map$description,function(d){
subbed_d<- gsub("\\[.*","",d)
return(subbed_d)
})
return(function_map)
}
#############################################################################
#############################################################################
get_nets<-function(pathwayAnalyses,dbclusters,db,external=FALSE)
{
#get all pathways in this db
require(visNetwork, quietly = TRUE)
options(browser="chromium-browser")
library(PathCluster)
#db<-"NBB.c2.cp.kegg"
parsed_obj <- parse_input(pathwayAnalyses[[db]],length(pathwayAnalyses[[db]]$genesets))
#generate key map for interpro function annotation of each gene
ensemblGenes<-get_gene_info(pathwayAnalyses,db)[,c("gene_id","gene_name")]
gene_func<-get_func_annot(ensemblGenes$gene_id)
#retrieve clusters df associated with this db
clusters<-dbclusters
#split df by cluster
library(dplyr)
#cat("clusters title \n",unique(clusters$title), "\n")
##TO BE DEBUGGED FROM HERE##
library(data.table)
clist<-split(as.data.table(clusters),f=as.data.table(clusters)$title)
clist<-clist[unique(clusters$title)]
print(length(clist))
#names(clist)<-unique(clusters$title
#clist<-list(clist[unique(clusters$title)])
#IS THIS CORRECT? CHECK
#print(head(clist))
print(head(clusters))
#names(clist)<-unique(clusters$title)
gp<-get_GP_mat(pathwayAnalyses,db,parsed_obj)

nets<- lapply(clist,function(cluster_df) {
print("#####################################################")
cat("generating net for clutser: ",cluster_df$title,"\n")
cat("We are in db: ",db,"\n")
print("##########################################################")
test_cluster<-cluster_df$items
title<-cluster_df$title[1]
print("calc occmat")
occmat<-modify_GP_mat(gp,test_cluster)
print("getting edges")
if(nrow(occmat)==0)
{
	print("found no overlapping genes for the pathways in this cluster, all genes that were common among the pathways in this cluster (causing the high similarity) have been filtered out due to non-significance in the data")
edge_list<-NULL
edge_list_key<-NULL
nodes<-NULL
node_list_key<-NULL
net<-NULL
return(list("nodes"=node_list_key,"edges"=edge_list_key,"network"=net))

}
else
{
edge_list <- get_edges(occmat)
print(head(edge_list))
genesets<-pathwayAnalyses[[db]]$genesets[test_cluster]
pvalues <-pathwayAnalyses[[db]]$p.values[test_cluster]
cat("in func lapply clist, pvalues are: ",pvalues,"\n")
if(external==FALSE)
{
nodes <- get_nodes(edge_list,pvalues,genesets,ensemblGenes,gene_func)
data_with_keys<-rekey_labels(edge_list,nodes)
node_list_key <-data_with_keys$nodes
edge_list_key<-data_with_keys$edge_list

net<-visNetwork(node_list_key, edge_list_key,main=title,width="200%") %>% visHierarchicalLayout(levelSeparation=200,blockShifting=F,sortMethod="directed",direction="DU",nodeSpacing=500)  %>% visEdges(width=0.5,selectionWidth=2,color=list(color="black",highlight="yellow")) %>% visInteraction(selectable = TRUE, multiselect = TRUE, hover=TRUE) 
return(list("nodes"=node_list_key,"edges"=edge_list_key,"network"=net))

}
else
{
	nodes <- get_nodes(edge_list,pvalues,genesets,ensemblGenes,gene_func,external)
	return(list("nodes"=nodes,"edges"=edge_list))
}
}
#TODO also return network parameters for plotting

})
print("###############################################")
names(nets)<-names(clist)
return(nets)
}


########################################################################
#############################TESTS######################################

#######################DATA#####################################
#use selfmade package to generate pathway gene binary occurance matrix (firstupdate packages)
library(readr)
load("./RData/pathwayAnalyses.Rda")

#list RData files in clusters folder
files <- list.files("/data/content/RNAseq/Pathway-sim/RData/clusters/",full.names=F)
cluster_names <- sapply(files,function(f){gsub(".RData","",f)})
files <- list.files("/data/content/RNAseq/Pathway-sim/RData/clusters/",full.names=T)
names(files) <- cluster_names

#iterate through names of pathwayAnalyses object to keep order, load cluster object to environment, make nets object and return to list
#iterativley rm cluster (overwritten anyways)
all_nets<-lapply(names(pathwayAnalyses),function(db){
file <- files[db]
cat("Now loading: ",file, "\n")
load(paste0("/data/content/RNAseq/Pathway-sim/RData/clusters/",file))
this_nets<-get_nets(pathwayAnalyses,cluster,db)
rm(cluster)
return(this_nets)
})
names(all_nets)<-names(pathwayAnalyses)
save(all_nets,file="/data/content/RNAseq/Pathway-sim/RData/nets/networks_all.RData")


################################################################
#######TODO######################################################
#move above functions to PathCluster package
#move import and call and net generation to app.R

#######################Save node and edge list as csv for external use (e.g. cytoscape)#############################
library(dplyr)
all_nets<-lapply(names(pathwayAnalyses),function(db){
file <- files[db]
cat("Now loading: ",file, "\n")
load(paste0("/data/content/RNAseq/Pathway-sim/RData/clusters/",file))
external=TRUE
this_nets<-get_nets(pathwayAnalyses,cluster,db,external)
unlisted_nodes <- lapply(this_nets,function(x)
       {
       		return(x$nodes)	
       })
names(unlisted_nodes) <- names(this_nets)
unlisted_nodes <- bind_rows(unlisted_nodes,.id="cluster_title")
unlisted_edges <- lapply(this_nets,function(x)
       {
       		return(x$edges)	
       })
names(unlisted_edges) <- names(this_nets)
unlisted_edges <- bind_rows(unlisted_edges,.id="cluster_title")

write.table(unlisted_nodes,file=paste0("/data/content/RNAseq/Pathway-sim/RData/txtFiles/",db,"_nodes.txt"),col.names=TRUE,quote=FALSE,sep="\t",row.names=FALSE)
write.table(unlisted_edges[,-c(4)],file=paste0("/data/content/RNAseq/Pathway-sim/RData/txtFiles/",db,"_edges.txt"),col.names=TRUE,quote=FALSE,sep="\t",row.names=FALSE)
rm(cluster)
return(this_nets)
})

