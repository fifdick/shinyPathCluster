

print("loading/installing dependencies")
list.of.packages <- c("shiny", "visNetwork","DT","shinyalert","metap")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)



require(shiny)
library(visNetwork)
library(DT)
library(shinyalert)
library(metap)
load("./data/networks_all.RData")
load("./data/pathwayAnalyses.Rda")


print("preprocessing data")
nets<-all_nets
nets<-lapply(nets,function(net){
idx<-lapply(net,function(x){!(is.null(x$edges))})
net<-net[unlist(idx)]
return(net)})

allgenes<-lapply(pathwayAnalyses[names(nets)],function(db){
genes<-db$genesets
genes<-do.call(rbind,genes) 
genes<-subset(genes,!duplicated(gene_id) & !duplicated(gene_name))
return(genes)
})
#allgenes <- as.data.table(allgenes)
allgenes<-lapply(allgenes,function(g){
#as.data.table(allgenes)
g<-g[,-c(1)]
colnames(g)<-c("ENSEMBL_ID","gene_name","pvalue","log2FoldChange")
return(g)
})
names(allgenes) <- names(pathwayAnalyses[names(nets)])


#################################################################
server<-function(input, output, session) {
#cat(file=stderr(),input$db,input$type(),"\n")

	     

observeEvent(c(input$db,input$type),{
#cat(file=stderr(),"choices",choices())
noPwInNet<-lapply(nets[[input$db]],function(n){
	table <- table(n$nodes$group)
	return(table["pathway"])
	})

single_nets<-names(nets[[input$db]])[noPwInNet==1]
single_nets_order<-lapply(single_nets,function(title){
	pathways_in_net<-rownames(subset(nets[[input$db]][[title]]$nodes,group=="pathway"))
	pval <- pathwayAnalyses[[input$db]]$p.values[pathways_in_net]
	return(list("names"=title ,"value"=pval))
	})
single_nets_order <- as.data.frame(do.call(rbind,single_nets_order))
single_nets_order<-single_nets_order[order(as.numeric(single_nets_order$value)),1]








multi_nets <- names(nets[[input$db]])[noPwInNet>1]
multi_nets_order<-lapply(multi_nets,function(title){
	pathways_in_net<-rownames(subset(nets[[input$db]][[title]]$nodes,group=="pathway"))
	pvals <- pathwayAnalyses[[input$db]]$p.values[pathways_in_net]
	if((length(pvals)-sum(pvals==0))<2  || length(pvals)<2)
	{	
		combined_pvalue <- 0
	}
	else
	{
		combined_pvalue <- sumlog(pvals)$p
	}
	return(list("names"=title ,"value"=combined_pvalue))
	})
multi_nets_order <- as.data.frame(do.call(rbind,multi_nets_order))
print(multi_nets_order)
multi_nets_order<-multi_nets_order[order(as.numeric(multi_nets_order$value)),1]
print(multi_nets_order)
if(length(multi_nets_order)==0)	
{
	multi_nets_order <- c("NA")
}
choices<-list("Both"=names(nets[[input$db]]),"Single pathway networks"=single_nets_order,"Multi pathway networks"=multi_nets_order)

updateSelectizeInput(session=session,inputId='title',choices=choices[[input$type]])#choices=names(choices()[[input$type]]))
})    
	     


output$network <- renderVisNetwork({
validate(
	 need(input$title!="NA","Pathway similarities (geneset intersections) did not reach the required threshold. No clusters found. All pathways are listed as single network pathways")
)
	

	net<-nets[[input$db]]

	nodes <-net[[input$title]]$nodes

		 nodes$title <-sapply(nodes$title,function(l){
		 label<-gsub("KEGG_|GO_","",l)
		 return(label)}) 

		 edges<-net[[input$title]]$edges
	  
		legend_df <- subset(nodes,group=="pathway")
		legend_df$title<-sapply(legend_df$title,function(t){
			stripped <- gsub("<p><b>", "", t)
			stripped <- gsub("</b></p>", "", stripped)
			return(stripped)
			})   
		legend_df$label  <- sapply(seq(1,nrow(legend_df)),function(i){
			title<-legend_df$title[i]
			label<-legend_df$label[i]
			return(paste0(label,": ",title))
			})
		legend_df <- subset(legend_df, select = -c(6))

	visNetwork(nodes, edges,main=input$title,width="100%",height="300%") %>% visHierarchicalLayout(levelSeparation=100,blockShifting=F,sortMethod="directed",direction="DU",nodeSpacing=500)  %>% visEdges(width=0.5,selectionWidth=4,color=list(highlight="yellow"),dashes=F) %>% visInteraction(selectable = TRUE, multiselect = TRUE, hover=TRUE) %>% visLegend(addNodes=legend_df,useGroups=FALSE,width=0.2)
	  
})

output$mytable =renderDataTable({
validate(
	 need(input$title!="NA","Select \"Single pathway networks\" in 2nd dropdown menu above")
)
	
	net<-nets[[input$db]]
	nodes <-net[[input$title]]$nodes
	nodes<-subset(nodes,group=="gene")
	df<-subset(allgenes[[input$db]], gene_name %in% nodes$label)
	rownames(df)<-NULL
	df$absLogFold<-abs(as.numeric(df[,4]))
	DT::datatable(df,caption=paste0("Pathway database:     ",input$db,"       Pathway cluster: ",input$title),extensions='Buttons', options = list(dom = 'Bfrtip', buttons = c('copy', 'excel', 'pdf', 'print', 'colvis'),paging = FALSE,columnDefs = list(list(searchable = FALSE, targets = c(1,3)))),filter='top')


})#renderDataTable

observeEvent(input$info,{   shinyalert(
    title = "Help",
    text = paste0("1) Choose desired database (predefined lists of pathways with specific genelists per pathway) <br> 2) Choose to explore networks that consist of more than one pathway to see what genes they share, or choose to explore single pathway networks, which are pathways that have no similarity with any other pathway in the database and are thus not clustered. <br> 3) Choose a cluster <br> Pathways are shown in blue diamonds. <br> Size of diamonds represents significance of pathway in dataset.<br> Genes are shown in red or green dots, indicating under or overexpression in dataset. <br> Size of genes represent connectivity of the node. Bigger (smaller) dots represent higher (lower) degree (number of edges, number of pathways they belong to). <br> Edges indicate which genes belong to which pathway. <br> Pathway titles are extended on the left-hand side legend or by hovering above the pathway node. <br> Gene description is visible from hovering above gene node. <br> More detailed ensembl information is linked. <br> Nodes and edges can be selected to highlight. <br> Selecting nodes, highlights all belonging edges. <br> Multiselect by pressing \"Ctrl\" . <br> Data table on the bottom lists all genes that are represented in the plot. <br> List can be downloaded to different data formats by pressing one of the buttons."), 
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = TRUE,
    type = "info",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "Thanks!",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ) })

  observeEvent(input$about,{   shinyalert(
    title = "Info",
    text = "This App displays gene-pathway networks for all pathways in all available databases. <br> Those pathways are predefined genesets, according to their sources (e.g. KEGG, GO, etc.) <br> All pathways were significantly enriched in our dataset. <br> For each database, significant pathways were clustered based on similarity (Kappa score, based on number of shared genes). <br> Plots show only genes, which are shared among the pathways, and which came up significant in the DGE analysis. <br> Cluster titles were selected according to multiple criteria including pvalue, size of pathway, etc.. <br> Single pathway clusters are listed separately and are sorted in the dropdown by significance. <br> Multi pathway clusters are sorted in the dropdown by significance, where significance per cluster was calculated by combining all pathway pvalues from this cluster using Fisher's method.",
    closeOnEsc = TRUE,
    closeOnClickOutside = TRUE,
    html = TRUE,
    type = "info",
    showConfirmButton = TRUE,
    showCancelButton = FALSE,
    confirmButtonText = "Thanks!",
    confirmButtonCol = "#AEDEF4",
    timer = 0,
    imageUrl = "",
    animation = TRUE
  ) })

  
  }#server


dbs<-names(nets)
ui <- fluidPage(
titlePanel("RNAseq gene-pathway networks"),
 useShinyalert(),
actionButton("info","Help"),
actionButton("about","Info"),
hr(),
fluidRow(
column(4,
selectInput('db','Database',names(nets),selected=names(nets)[1])),
	 column(4,
selectInput('type',"Network type",choices=c("Both","Single pathway networks","Multi pathway networks"),selected="Multi pathway networks")),
	 column(4,
selectInput('title', 'Cluster',names(nets[[1]])) )
),
hr(),
visNetworkOutput("network",width="100%",height="600px"),
dataTableOutput("mytable")
)



#########################CALLING###################################
print("setting chrome as browser")
print("Have fun exploring!")
options(browser="google-chrome")
shinyApp(ui=ui,server=server)
#Sys.sleep(10)
#stopApp(returnValue = invisible())

#library(rsconnect)
#rsconnect::deployApp('./app.R')
