############################################
######Function definitions##################
############################################

###########################################################

#input: obj:two slots:list of pathways (genesets) and list of corrected p-values of pathways
#first slot is named list of data.frames. The name corresponds to the pathway ID.
#The data.frame contains the genes of that pathway (rows) and p.values, ids, and log fold-chances (columns).
#ouput:  list[vector[char]]. Each entry in list is a named vector, where name is pathway name and values of vector are gene_names

parse_input_go<-function(df,subsetsize)
{
  library(sets)
  library(stringr)
  input<-lapply(seq(1,subsetsize),function(y){
    names<-df$genesets[[y]]$gene_name
    vec<-ifelse(df$genesets[[y]]$pvalue>0.05,0,1)
    names(vec)<-names
    return(vec)
  })
  names(input)<-names(df$genesets)[1:subsetsize]
  names(input)<-sapply(seq(1:length(names(input))),function(x){
    gsub("^GO:[0-9]*","",names(input)[x])
  })
  pvalues<-df$p.values[1:subsetsize]
  names(pvalues)<-names(input)
  return(list("geneset"=input,"pvals"=pvalues))
}

parse_input_kegg<-function(df)
{
  library(sets)
  library(stringr)
  input<-lapply(seq(1,length(df$genesets)),function(y){
    names<-df$genesets[[y]]$gene_name
    vec<-ifelse(df$genesets[[y]]$pvalue>0.05,0,1)
    names(vec)<-names
    return(vec)
  })
  names(input)<-names(df$genesets)
  names(input)<-sapply(seq(1:length(names(input))),function(x){
    gsub("^hsa[0-9]*","",names(input)[x])
  })
  pvalues<-df$p.values
  names(pvalues)<-names(input)
  return(list("geneset"=input,"pvals"=pvalues))
}

parse_input<-function(df,subsetsize)
{
  library(sets)
  library(stringr)
  input<-lapply(seq(1,subsetsize),function(y){
    names<-df$genesets[[y]]$gene_name
    vec<-ifelse(df$genesets[[y]]$pvalue>0.05,0,1)
    names(vec)<-names
    return(vec)
  })
  names(input)<-names(df$genesets)[1:subsetsize]
  pvalues<-df$p.values[1:subsetsize]
  names(pvalues)<-names(input)[1:subsetsize]
  return(list("geneset"=input,"pvals"=pvalues))
}


###########################################################

#input: list[vector[char]]. Each entry in list is a named vector, where name is pathway name and values of vector are gene_names
#output: data.frame[numeric]]. Binary occurance matrix. Cols are pathways, rows are genes. 1 indicates occurance 0 indicates absence

gene_pathway_matrix<-function(alist){
  
 gene_names<-unique(unlist(lapply(alist,names)))
 pathway_names<-names(alist)
 
 df<-sapply(seq(1:length(pathway_names)),function(y){
   vec<-alist[[y]]
   df_col<-sapply(seq(1:length(gene_names)),function(x){
     if(gene_names[x] %in% names(vec) )
        1
     else
       0
   })
   return(df_col)
 })
 
 colnames(df)<-pathway_names
 rownames(df)<-gene_names
 
 return(df)
}

###########################################################

#input: data.frame[numeric]. Binary occurance matrix. Cols are pathways, rows are genes. 1 indicates occurance 0 indicates absence
#output: data.frame[double]. Similarity matrix. Cols and rows are pathways. All versus all similarities. Values are Kappa scores.

sim_matrix<-function (alist){
  amatrix<-gene_pathway_matrix(alist)
  require(irr)
#  require(parallel)
#  cl <- makeCluster(8)
  df<-sapply(seq(1:ncol(amatrix)),function(y){
    one<-as.data.frame(amatrix[,y])
    col_kappas<-sapply(seq(1:ncol(amatrix)),function(x) {
      kappa2(data.frame(one,amatrix[,x]))$value
      })
    col_kappas  
    })
  rownames(df)<-colnames(amatrix)
  colnames(df)<-colnames(amatrix)
  return(df)
}

###########################################################

#input: list[vector[char]]. Each entry in list is a named vector, where name is pathway name and values of vector are gene_names
#output: "void" plots distribution of pathway sizes (number of genes in pathway) and gives summary stats

plot_size_dist<-function(alist)
{
  require(ggplot2)
  bin_mat<-gene_pathway_matrix(alist)
  gene_count<-sapply(seq(1,ncol(bin_mat)),function(x) {
   sum(bin_mat[,x])
  })
  filter<-sapply(gene_count,function(x){if(x>50 && x<1000) return("passed") else return("failed")})
  gene_count<-data.frame(count=gene_count,filter=filter)
  cols <- c("passed" = "palegreen2", "failed" = "grey68")
  ggplot(data=gene_count,aes(count,fill=filter)) +
    geom_histogram(aes(y = ..density..)) +
    labs( title="Size distribution of pathways", x="Number of genes per pathway",y="density",fill="hard filter") +
    geom_vline(xintercept=50,col="red") +
    geom_vline(xintercept=1000,col="red") +
    scale_fill_manual(labels = c("failed", "50<x<1000"),values=cols) +
    stat_function(fun = dnorm, 
                  args = list(mean = mean(gene_count$count), sd = sd(gene_count$count)), 
                  lwd = 1, 
                  col = 'red')
  summary(gene_count$count)
}

###########################################################

#construction of netwrok from enriched pathways
#evaluation of significance of all non-empoty intersections between two pathways

intersect_signif_matrix<-function(alist){
  df<-sapply(seq(1:length(alist)),function(y){
    one<-alist[[y]]
    sig_values<-sapply(seq(1:length(alist)),function(x) {
      int<-intersect(names(one),names(alist[[x]]))
      N_i<-length(int)
      S_i<-sum(one[int])
      S_u<-sum(one)+sum(alist[[x]])
      N_u<-length(one)+length(alist[[x]])
      pvalue<-phyper(S_i,S_u,N_u-S_u,N_i)
      #pvalue<-fisher.test(matrix(c(S_i,S_u-S_i,N_i-S_i,(N_u-S_i)-S_u),2,2),alternative="less")$value
      return(pvalue)
    })
    sig_values
  })
  rownames(df)<-names(alist)
  colnames(df)<-names(alist)
  return(df)
}

adj_matrix<-function(alist){
  sig_mat<-intersect_signif_matrix(alist)
  adj_matrix<-apply(sig_mat,c(1,2),function(x){ifelse(x<0.05,1,0)})
  rownames(adj_matrix)<-rownames(sig_mat)
  colnames(adj_matrix)<-colnames(sig_mat)
  return(adj_matrix)
}

###########DEBUGGING########################

#test_input1<-list("pathway1"=c("gene1"=1,"gene2"=1,"gene4"=0),"pathway2"=c("gene1"=1,"gene2"=1,"gene5"=1),"pathway3"=c("gene1"=1,"gene2"=1,"gene3"=1))
#gene_count<-sample(seq(1:10000),1000,replace=TRUE)
#gene_count<-rnorm(1000,mean=500,sd=100)

#############################################

# Acutual algorithm loop #

#############################################




redundant_filter<-function (threshold,alist,pvalues) {
print("calc adj matrix")
adjmat<-adj_matrix(alist)
print("calc similarity matrix")
  simmat<-sim_matrix(alist)
  terms<-rownames(simmat)
  len<-length(terms)
  cluster<-matrix(rep(0,len*len),nrow=len,ncol=len)
  cluster<-as.data.frame(cluster)
  diag(cluster)<-1
#  print(length(terms))
#  print(sum(is.na(terms)))
  colnames(cluster)<-terms
  rownames(cluster)<-terms
  diag(simmat)<-NA
  
  update_simmat<-function(failed,passed,simmat,alist){
    diag(simmat)<-1
    gene_pathway_mat<-gene_pathway_matrix(alist)
    genes_in_cluster<-rep(0,length(rownames(gene_pathway_mat)))
    names(genes_in_cluster)<-rownames(gene_pathway_mat)
    merged<-unique(c(names(alist[[failed]]),names(alist[[passed]])))
    genes_in_cluster<-sapply(names(genes_in_cluster),function(z){
      if( z %in% merged)
      { return(1)
      }
      else
        return(0)
    })
    names(genes_in_cluster)<-rownames(gene_pathway_mat)
    new_kappa_col<-sapply(seq(1,nrow(simmat)),function(x){
      new_sim<-kappa2(data.frame(gene_pathway_mat[,x],genes_in_cluster))$value
      return(new_sim)
    })
    simmat[,passed]<-new_kappa_col
    simmat<-simmat[-which(rownames(simmat)==failed),-which(rownames(simmat)==failed)]
    diag(simmat)<-NA
    return(simmat)
  }
 # print("starting algorithm while loop")
 # print(simmat)
  while ( max(simmat[upper.tri(simmat, diag = FALSE)])>=threshold)
  {
    
    index<-which(simmat == max(simmat[upper.tri(simmat, diag = FALSE)]), arr.ind = TRUE)
    t_i=index[1,1] #row
    t_j=index[1,2] #col
    t<-c(t_i,t_j)
 #   print(t)
    #if both pathways are either too general or too specific throw out the too general
    if( (length(alist[[rownames(simmat)[t_i]]])<=50 || length(alist[[rownames(simmat)[t_i]]])>=1000) && (length(alist[[rownames(simmat)[t_j]]])<=50 || length(alist[[rownames(simmat)[t_i]]])>=1000) )
    {
      #t_i1, t_i=2

      pair_size<-c(length(alist[[rownames(simmat)[t_i]]]),length(alist[[rownames(simmat)[t_j]]]))
      failed<-rownames(simmat)[t[which.max(pair_size)]]
      passed<-rownames(simmat)[t[-which.max(pair_size)]]
      simmat<-update_simmat(failed,passed,simmat,alist)
      cluster[failed,passed]<- 1
      cluster[failed,failed]<- -1
      next
    }
  #if one of the pathways is either too general or too specific throw it
  else if (  (length(alist[[rownames(simmat)[t_i]]])<=50 || length(alist[[rownames(simmat)[t_i]]])>=1000) )
  {
    failed<-rownames(simmat)[t_i]
    passed<-rownames(simmat)[t_j]
    simmat<-update_simmat(failed,passed,simmat,alist)
    cluster[failed,passed]<- 1
    cluster[failed,failed]<- -1
    next
    
  }
  #if one of the pathways is either too general or too specific throw it
  else if((length(alist[[rownames(simmat)[t_j]]])<=50 || length(alist[[rownames(simmat)[t_j]]])>=1000))
  {
    failed<-rownames(simmat)[t_j]
    passed<-rownames(simmat)[t_i]
    simmat<-update_simmat(failed,passed,simmat,alist)
    cluster[failed,passed]<- 1
    cluster[failed,failed]<- -1
    next
  }
  #if both pathways fulfill size requirements 
  else
  {# print(pvalues[rownames(simmat)[t_i]])
  	print(pvalues)
  #print(rownames(simmat)[t_i])
  #print(rownames(simmat))
    if ( pvalues[rownames(simmat)[t_i]]!=pvalues[rownames(simmat)[t_j]] )
    {
      
      less_sig<-which.max(c(pvalues[rownames(simmat)[t_i]], pvalues[rownames(simmat)[t_j]]))
      failed<-rownames(simmat)[t[less_sig]]
      passed<-rownames(simmat)[t[-less_sig]]
      simmat<-update_simmat(failed,passed,simmat,alist)
      cluster[failed,passed]<- 1
      cluster[failed,failed]<- -1
      next
    }
    else if (adjmat[rownames(simmat)[t_i],rownames(simmat)[t_j]]!=0)
    {
      #reject childterm (more specific)
     
      child<-which.min(c(sum(alist[[rownames(simmat)[t_i]]]),sum(alist[[rownames(simmat)[t_j]]])))
      failed<-rownames(simmat)[t[child]]
      passed<-rownames(simmat)[t[-child]]
      simmat<-update_simmat(failed,passed,simmat,alist)

      next
    }
    #non of the above filter characteristics applied-> random choice
    else
    {
      failed<-rownames(simmat)[t_i]
      passed<-rownames(simmat)[t_j]
      simmat<-update_simmat(failed,passed,simmat,alist)
      cluster[failed,passed]<- 1
      cluster[failed,failed]<- -1
      next
    }
  }
    
    
    
  }#while
  return(list(cluster,simmat))
}#function


process_result<-function(amatrix,alist,pvalues){
  library(qdapTools)
print("processing result")
  cluster_titles<-names(which(diag(as.matrix(amatrix))!=-1, arr.ind=TRUE))
  clusters<-lapply(seq(1,length(cluster_titles)),function (y) {
    title<-cluster_titles[y]
    belongings<-rownames(amatrix)[which(amatrix[title]==1)]
    return(belongings)
  })
  names(clusters)<-cluster_titles
  clusters<-list2df(clusters,col1="items",col2="title")
  clusters$pvalue<-pvalues[clusters$items]
  clusters$abslog10<-abs(log10(clusters$pvalue))
 # clusters$items_long<-
  #cluster$title_long<-
  return(clusters)
  }
  

print_map<-function(clusters,title){
  library(treemap)
  tmPlot(
    clusters,
    index = c("title","items"),
    vSize = "abslog10",
    type = "categorical",
    vColor = "title",
    title = title,
    inflate.labels =FALSE,      # set this to TRUE for space-filling group labels - good for posters
    lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
    bg.labels = "#CCCCCCAA",     # define background color of group labels
    # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
    position.legend = "none",
    palette="PiYG",
    fontcolor.labels=c("black","white"),
    border.col=c("black","white"),             # Color of borders of groups, of subgroups, of subsubgroups ....
    border.lwds=c(5,2)               
  )
  
}

rekey_go<-function(clusters,names){
  clusters$item_names<-sapply(seq(1,length(clusters$items)),function(x){
   
    pattern<-clusters$items[x]
    names[grep(pattern,names)]
  })
  clusters$title_names<-sapply(seq(1,length(clusters$title)),function(x){
    pattern<-clusters$title[x]
    names[grep(pattern,names)]
  })
  clusters$item_names<-sapply(seq(1,length(clusters$items)),function(x){
    gsub(clusters$items[x],"",clusters$item_names[x])
    })
  clusters$title_names<-sapply(seq(1,length(clusters$title)),function(x){
    gsub(clusters$title[x],"",clusters$title_names[x])
  })
  return(clusters)
}

###########DEBUGGING########################



#test3<-parse_input_go(go_cc_NBB)
#test_result3<-redundant_filter(0.4,test3$geneset,test3$pvals)
#clusters<-process_result(test_result3[[1]],test3$geneset,test3$pvals)
#clusters<-rekey(clusters,names(go_cc_NBB$genesets))
#print_map(clusters)


###################TESTING###########################################

allDup <- function (value)
{
  duplicated(value) | duplicated(value, fromLast = TRUE)
}

run_pathway_vis<-function (input,title,subsetsize=length(input$genesets)) {
if ( subsetsize > 1000)
{
	subsetsize <- 1000
}

test4<-parse_input(input,subsetsize)
#print(head(test4$geneset))
#print(test4$pvals)
start_time <- Sys.time()
filtered<-redundant_filter(0.4,test4$geneset,test4$pvals)
end_time <- Sys.time()
end_time - start_time
clusters<-process_result(filtered[[1]],test4$geneset,test4$pvals)
#morethan2only<-clusters[allDup(clusters$title),]
#pdf(paste0("./",title,".pdf"))
#print(print_map(clusters,title))
#dev.off()
return(clusters)
}
###############################################################################################
###########################################DO#####################################################
#library(PathCluster)
library(readr)
load("/data/content/RNAseq/Pathway-sim/RData/pathwayAnalyses.Rda")
do.call(rbind,lapply(pathwayAnalyses,function(x){return(c(length(x$genesets),length(x$p.values)))}))

#intersect to exclude NAs
pathwayAnalyses<-lapply(pathwayAnalyses, function(obj) {
obj$genesets<-obj$genesets[names(obj$p.values)]
return(obj)
  })


library(parallel) 
#library(mc2d)

# First the small ones:
lengths <- sapply(pathwayAnalyses,function(x){length(x$genesets)}) 
names(lengths) <- names(pathwayAnalyses)
small <- names(which(lengths<1000))
big<- names(which(lengths>=1000)) 
big
small
save.image("./image22-05.RData")  

clusters <- mclapply(missing, function(x) {
library(irr)
print("Now processing:")
input <- pathwayAnalyses[[x]]
cat("Title: ",x,"  Length: ",lengths[x],"\n")		
print("#############START#####################")
 title <- x
cluster <- run_pathway_vis(input,title)
print("saving")
save(cluster,file=paste0("./RData/clusters/",title,".RData"))
print("saved")
print("##############DONE###################")
return(cluster)
  },mc.cores=12)
Sys.time()
save(clusters,file="./RData/clusters/clusters_missing.RData")

parent_drivingGene <- lapply(seq(1,length(clusters)),function (x) {
  
  parent_ids <-unique(clusters[[x]]$title)

  gene_lists <- sapply(seq(1,length(parent_ids)),function (y) {
			 name <- parent_ids[y]
			       genes <-pathwayAnalyses[[x]]$genesets[[name]]
			sig_genes <- subset(genes,pvalue < 0.05)
return(sig_genes)
   
    })
  per_pathway_info <- list("parent_ids"=parent_ids,"driving_genes"=gene_lists)
  
  
  
  })


clusters[[10]]<-run_pathway_vis(pathwayAnalyses[[10]],"NBB.c5.mf")
save(clusters,file="./RData/clusters/clusters_all<1000.RData") 


#TODO !!! TO MAYBE ADJUST FILTER
plot_size_dist(test4$geneset)
library("qusage")
#TODO map names of msig to description
#TODO check kappa paper
#TODO MAYBE
#for each group add up pvalues (with apropritae statistic) and then plot top50: metap sumlog, take extra df with only title and combined pvalue, sort by pvalue and subset result clutsers by ttop 50 titles
#TABLE(clusters$title) also as list for HAris

hallmark_nbb_cl<-run_pathway_vis(hallmark_NBB,"Hallmark (msigDB) NBB")
hallmark_parkvest_cl<-run_pathway_vis(hallmark_parkvest,"Hallmark (msigDB) parkvest")
kegg_disease_NBB_cl<-run_pathway_vis(kegg_disease_NBB,"KEGG disease NBB")
kegg_disease_parkvest_cl<-run_pathway_vis(kegg_disease_parkvest,"KEGG disease Parkvest")
kegg_sigmet_NBB_cl<-run_pathway_vis(kegg_sigmet_NBB,"KEGG sigmet NBB")
kegg_sigmet_parkvest_cl<-run_pathway_vis(kegg_sigmet_parkvest,"KEGG sigmet Parkvest")
go_bp_parkvest_cl<-run_pathway_vis(go_bp_parkvest,"GO biological processes ParkVest",100)
go_mf_parkvest_cl<-run_pathway_vis(go_mf_parkvest,"GO molecular functions ParkVest",100)
go_cc_parkvest_cl<-run_pathway_vis(go_cc_parkvest,"GO cellular components ParkVest",100)


save(kegg_disease_NBB_cl,kegg_disease_parkvest_cl,kegg_sigmet_NBB_cl,kegg_sigmet_parkvest_cl,file="/data/content/RNAseq/kegg_clusters.RData")
save(go_bp_parkvest_cl,go_mf_parkvest_cl,file="/data/content/RNAseq/go_clusters.RData")


install.packages("PathCluster")
library(PathCluster)
start_time <- Sys.time()
msig_parkvest_cl<-run_pathway_vis(msig_parkvest,"MsigDB parkvest",100)
save(msig_parkvest_cl,file="/data/content/RNAseq/msig_parkvest_clusters.RData")
end_time <- Sys.time()
end_time - start_time


install.packages("wordcloud") # word-cloud generator 
library(wordcloud)
install.packages("RColorBrewer") # color palettes
