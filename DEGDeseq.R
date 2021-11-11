DEG.DESeq2 <- function (rawCounts,
                        maxP=.05, 
                        minFC=2, 
                        selectedComparisons=NULL, 
                        sampleInfo = NULL,
                        modelFactors=NULL, 
                        blockFactor = NULL, 
                        referenceLevels=NULL,
                        groups = NULL){

g = unique(groups)# order is reversed
  
  
  # check for replicates, removes samples without replicates
  reps = as.matrix(table(groups)) # number of replicates per biological sample
  if ( sum( reps[,1] >= 2) <2 ) # if less than 2 samples with replicates
    return( list(results= NULL, comparisons = NULL, Exp.type="Failed to parse sample names to define groups. 
		Cannot perform DEGs and pathway analysis. Please double check column names! Use WT_Rep1, WT_Rep2 etc. ", topGenes=NULL)) 
  # remove samples without replicates
  g <- rownames(reps)[which(reps[,1] >1)]
  ix <- which( groups %in% g)  
  groups <- groups[ix]   
  rawCounts <- rawCounts[,ix] 
  
  
  Exp.type = paste(length(g)," sample groups detected.")
  comparisons = ""
  for( i in 1:(length(g)-1) )
    for (j in (i+1):length(g)) 
      comparisons = c(comparisons,paste(g[j],"-",g[i],sep="" ) )
  comparisons <- comparisons[-1]
  
  colData = cbind(colnames(rawCounts), groups )
  
  # no sample file, but user selected comparisons using column names
  if( is.null(modelFactors) & length( selectedComparisons) >0  ) 	
    comparisons = selectedComparisons
  
  comparisons2 = comparisons	 # this is for showing comparison names, which might be different from internally	
  # Set up the DESeqDataSet Object and run the DESeq pipeline
  dds = DESeqDataSetFromMatrix(countData=rawCounts,
                               colData=colData,
                               design=~groups)								
  
  if( is.null(modelFactors)  ) 
    dds = DESeq(dds)  	else  
    {    # using selected factors and comparisons
      # build model
      modelFactors = c(modelFactors,blockFactor) # block factor is just added in. 
      
      factors = modelFactors   # selected factors and interactions: c( "strain", "treatment",  "strain:treatment")
      factors = factors[ !grepl(":",factors )]   # non-interaction terms
      # interaction terms like strain:treatment
      Interactions = modelFactors[ grepl(":",modelFactors )]
      
      colData = sampleInfo  
      factorsCoded = toupper(letters )[1: dim(colData)[2] ]   # Factors are encoded as "A", "B", "C"; this avoid illigal letters
      names(factorsCoded) =  colnames(colData)  # this is for look up; each column of sampleInfo  
      colnames(colData) = factorsCoded # all columns named A B C D 
      
      colData = as.data.frame(colData)
      
      # set reference levels for factors
      if(! is.null( referenceLevels) ) {   # c("genotype:wt", "treatment:control" )
        # first factor
        for ( refs in referenceLevels)
          if(! is.null( refs) ) {
            ix = match(gsub(":.*","",refs), colnames(sampleInfo) ) # corresponding column id for factor
            colData[,ix] = as.factor( colData[,ix] )
            colData[,ix] = relevel(colData[,ix],gsub(".*:","",refs)  )
          }
      }
      
      # base model
      DESeq2.Object= paste("dds = DESeqDataSetFromMatrix(countData=rawCounts, colData=colData, design=~ ", 
                           paste( factorsCoded[factors],collapse="+")) # only use selected factors		
      Exp.type = paste("Model: ~", paste(modelFactors,collapse=" + ") )
      
      
      # create model
      if( length(Interactions)>0 ) { # if there is interaction
        for( interactionTerms in Interactions) {
          interactingFactors = unlist(strsplit(interactionTerms,":" ) )  # split strain:treatment as "strain" and "mutant"
          tem = paste(factorsCoded [ interactingFactors ],collapse=":")   # convert "strain:mutant" to "A:B"
          DESeq2.Object = paste(DESeq2.Object, " + ",tem)
        }			
      }
      DESeq2.Object= paste( DESeq2.Object, ")") # ends the model
      
      eval(parse(text = DESeq2.Object) )
      dds = DESeq(dds)  # main function		
      
      # comparisons 
      # "group: control vs. mutant"
      comparisons = gsub(".*: ","",selectedComparisons)
      comparisons = gsub(" vs\\. ","-",comparisons)
      factorsVector= gsub(":.*","",selectedComparisons) # corresponding factors for each comparison
      
      # comparison2 holds names for display with real factor names
      # comparison  is used in calculation it is A, B, C for factors
      comparisons2 = comparisons    
      #comparisons2 = gsub(" vs\\. ","-",selectedComparisons) 
      #comparisons2 = gsub(":","_",comparisons2)		
      # Note that with interaction terms, not all meaningful comparisons is listed for selection. 
      # this is complex. Only under reference level.
      
      # comparisons due to interaction terms
      if( length(Interactions)>0 ) { # if there is interaction
        interactionComparisons = resultsNames(dds)
        interactionComparisons = interactionComparisons[ grepl("\\.",interactionComparisons )   ]
        
        comparisons = c(comparisons,interactionComparisons )
        
        # translate comparisons generated in interaction terms back to real factor names
        interactionComparisons2 = interactionComparisons
        for ( i in 1:length(interactionComparisons2 ) ) {
          tem = unlist(strsplit(interactionComparisons2[i],"\\." ) )
          tem_factors = substr(tem,1,1) 
          
          tem_factors[1] = names(factorsCoded)[factorsCoded == tem_factors[1]]  # get the first letter and translate into real factor names
          tem_factors[2] = names(factorsCoded)[factorsCoded == tem_factors[2]]  # get the 2nd letters and translate into real factor names
          
          interactionComparisons2[i] <- paste0( "I:",tem_factors[1], "_",substr(tem[1],2,nchar(tem[1]) ),".",
                                                tem_factors[2], "_",substr(tem[2],2,nchar(tem[2]) ) 
          )				
        }
        comparisons2 = c(comparisons2,interactionComparisons2 )
      }
    } # if selected factors	
  
  # extract contrasts according to comprisons defined above
  result1 = NULL; allCalls = NULL;
  topGenes = list(); pk = 1 # counter
  pp=0 # first results?
  for( kk in 1:length(comparisons) ) {
    tem = unlist( strsplit(comparisons[kk],"-") )
    
    if(is.null(modelFactors)) # if just group comparison using sample names
      selected = results(dds, contrast=c("groups", tem[1], tem[2]) )   else {
        if(!grepl("\\.", comparisons[kk] ) )    # if not interaction term: they contain .  interaction term
          selected = results(dds, contrast=c( factorsCoded[ factorsVector[kk] ],tem[1], tem[2]) ) else # either A, B, C ...
            selected = results(dds, name=comparisons[kk] ) # interaction term
      }
    
    selected$calls =0   
    selected$calls [which( selected$log2FoldChange > log2(minFC) & selected$padj < maxP ) ]  <-  1
    selected$calls [ which( selected$log2FoldChange <  -log2(minFC) & selected$padj < maxP ) ] <-  -1
    colnames(selected)= paste( as.character(comparisons2[kk]), "___",colnames(selected),sep="" )
    selected = as.data.frame(selected)
    if (pp==0){  # if first one with significant genes, collect gene list and Pval+ fold
      result1 = selected; pp = 1; 
      # selected[,2] <- -1 * selected[,2] # reverse fold change direction
      topGenes[[1]] = selected; 
      names(topGenes)[1] = comparisons2[kk]; } else 
      { result1 = merge(result1,selected,by="row.names"); 
      rownames(result1) = result1[,1]; 
      result1 <- result1[,-1]
      pk= pk+1; 
      # selected[,2] <- -1 * selected[,2] # reverse fold change direction
      topGenes[[pk]] = selected; 
      names(topGenes)[pk] = comparisons2[kk];  # assign name to comprison
      }
  }
  
  Interactions = c()
  if( !is.null(modelFactors) )
    Interactions = modelFactors[ grepl(":",modelFactors )]
  
  #---  add comprisons for non-reference levels. It adds to the results1 object.	
  if( length(Interactions)>0 ) { # if there is interaction
    factorLookup=c() # a factor whose values are factors and names are factor and level combination conditionTreated, genotypeWT
    levelLookup = c()
    
    for( i in 1:dim(sampleInfo)[2]) {
      sampleInfo2 = unique(sampleInfo)
      tem = rep(toupper(letters)[i],dim(sampleInfo2)[1]  )
      names(tem) = paste0(toupper(letters)[i],sampleInfo2[,i])
      factorLookup = c(factorLookup,tem)  
      
      tem = as.character( sampleInfo2[,i] )
      names(tem) = paste0(toupper(letters)[i],sampleInfo2[,i])
      levelLookup = c(levelLookup, tem)
    }
    
    # split  genotypeI.conditionTrt --> c("genotype","I","conditoin","Trt")
    splitInteractionTerms <- function (term) {
      if(!grepl("\\.",term) ) return(NULL)
      terms2 = unlist(strsplit(term,"\\.") )
      # factor1, level1, factor2, level2
      return(c(factorLookup[terms2[1]], levelLookup[terms2[1]],factorLookup[terms2[2]], levelLookup[terms2[2]]   ) )
    }
    # none interaction terms 
    NoneInterTerms = resultsNames(dds)[ !grepl( "\\.", resultsNames(dds)) ]
    NoneInterTerms=NoneInterTerms[-1]
    allInteractionTerms = resultsNames(dds)[ grepl( "\\.", resultsNames(dds)) ]
    
    
    for( kk in 1:length(NoneInterTerms) ) { # for each none interaction term
      if(!is.null(modelFactors) ) {# if not just group comparison using sample names
        #current factor
        cFactor = gsub("_.*","",NoneInterTerms[kk] )
        
        for(interactionTerm in allInteractionTerms ) {
          
          splited = splitInteractionTerms (interactionTerm)  # 4 components
          if (cFactor != splited[1] & cFactor != splited[3]  ) 
            next;						
          
          selected = results(dds, list(c( NoneInterTerms[kk],interactionTerm ) ) ) 
          comparisonName = paste0( NoneInterTerms[kk],"__", gsub("\\.","",interactionTerm) )
          
          if( cFactor == splited[1] )
            otherLevel = splited[4] else otherLevel = splited[2]
          
          comparisonName = paste0(#names(factorsCoded)[which(factorsCoded==cFactor)], # real factor name
            gsub("_vs_","-", substr(NoneInterTerms[kk], 3, nchar(NoneInterTerms[kk]  )  )), # the comparison
            "_for_",otherLevel)
          comparisons2 = c(comparisons2, comparisonName)
          selected$calls =0   
          selected$calls [which( selected$log2FoldChange > log2(minFC) & selected$padj < maxP ) ]  <-  1
          selected$calls [ which( selected$log2FoldChange <  -log2(minFC) & selected$padj < maxP ) ] <-  -1
    
          colnames(selected)= paste( comparisonName, "___",colnames(selected),sep="" )
          selected = as.data.frame(selected)
          if (pp==0){  # if first one with significant genes, collect gene list and Pval+ fold
            result1 = selected; pp = 1; 
            # selected[,2] <- -1 * selected[,2] # reverse fold change direction
            topGenes[[1]] = selected; 
            names(topGenes)[1] = comparisonName; } else 
            { result1 = merge(result1,selected,by="row.names"); 
            rownames(result1) = result1[,1]; 
            result1 <- result1[,-1]
            pk= pk+1; 
            # selected[,2] <- -1 * selected[,2] # reverse fold change direction
            topGenes[[pk]] = selected; 
            names(topGenes)[pk] = comparisonName;  # assign name to comprison
            }
        } #for	
        
      } #if
    } #for
    
    
  } #if
  
  
  
  
  #---
  #if( length(comparisons) == 1) topGenes <- topGenes[[1]] # if only one comparison, topGenes is not a list, just a data frame itself.
  if(! is.null(result1)) { 
    # note that when you only select 1 column from a data frame it automatically converts to a vector. drop =FALSE prevents that.
    allCalls = as.matrix( result1[,grep("calls",colnames(result1)), drop = FALSE  ] )
    colnames(allCalls)= gsub("___.*","", colnames(allCalls))
    colnames(allCalls)= gsub("\\.","-", colnames(allCalls)) # note that samples names should have no "."
    colnames(allCalls)= gsub("^I-","I:", colnames(allCalls))
  }
  return( list(results= allCalls, comparisons = comparisons2, Exp.type=Exp.type, topGenes=topGenes)) 
}

writeDEGxlsx <- function(){
  
  wb <- createWorkbook()
  wb2 <- createWorkbook()
  wb3 <- createWorkbook()
  
  
  for(i in dds$comparisons)
  {
    annotation=annot %>% rownames_to_column("ID")
    factor_labeling=as.data.frame(dds$topGenes[c(paste(i))]) %>% rownames_to_column("ID") %>% left_join(annotation,"ID")
    colnames(factor_labeling)=gsub("\\S+___","",colnames(factor_labeling))
      
    sname<-paste(i,sep="")
    ifelse(i==i,app<-"FALSE",app<-"TRUE")
    addWorksheet(wb, sname)
    writeData(wb, sname, factor_labeling)
    
  }
  
  for(i in dds$comparisons)
  {
    annotation=annot %>% rownames_to_column("ID")
    factor_labeling=as.data.frame(dds$topGenes[c(paste(i))]) %>% rownames_to_column("ID")
    colnames(factor_labeling)=gsub("\\S+___","",colnames(factor_labeling))
    factor_labeling=filter(factor_labeling,padj <= as.numeric(0.01))
    factor_labeling=filter(factor_labeling,calls != 0) %>%
      left_join(annotation,"ID")
    sname<-paste(i,sep="")
    ifelse(i==i,app<-"FALSE",app<-"TRUE")
    addWorksheet(wb2, sname)
    writeData(wb2, sname, factor_labeling)
    
  }
  
  for(i in dds$comparisons)
  {
    annotation=annot %>% rownames_to_column("ID")
    factor_labeling=as.data.frame(dds$topGenes[c(paste(i))]) %>% rownames_to_column("ID")
    colnames(factor_labeling)=gsub("\\S+___","",colnames(factor_labeling))
    factor_labeling=filter(factor_labeling,padj <= as.numeric(0.05))
    factor_labeling=filter(factor_labeling,calls != 0) %>%
      left_join(annotation,"ID")
    sname<-paste(i,sep="")
    ifelse(i==i,app<-"FALSE",app<-"TRUE")
    addWorksheet(wb3, sname)
    writeData(wb3, sname, factor_labeling)
    
  }
  
  saveWorkbook(wb, file = "DEG.xlsx", overwrite = TRUE)
  saveWorkbook(wb2, file = "DEG_0.01.xlsx", overwrite = TRUE)
  saveWorkbook(wb3, file = "DEG_0.05.xlsx", overwrite = TRUE)
}


detectGroups <- function (x){  # x are col names
  tem <- gsub("[0-9]*$","",x) # Remove all numbers from end
  #tem = gsub("_Rep|_rep|_REP","",tem)
  tem <- gsub("_$","",tem); # remove "_" from end
  tem <- gsub("_Rep$","",tem); # remove "_Rep" from end
  tem <- gsub("_rep$","",tem); # remove "_rep" from end
  tem <- gsub("_REP$","",tem)  # remove "_REP" from end
  return( tem )
}

GO <- function(contrast,gene_lengths,background,Go_info,pcuttof = 0.01){
gowb<- createWorkbook()
  for (i in contrast) {
    
  factor_labeling=as.data.frame(dds$topGenes[c(paste(i))]) 
  colnames(factor_labeling)=gsub("\\S+___","",colnames(factor_labeling))
  factor_labeling=filter(factor_labeling,padj <= as.numeric(pcuttof))
  factor_labeling=filter(factor_labeling,calls != 0)
  
  # capture list of genes for functional enrichment testing
  factor_labeling = factor_labeling
  factor_labeling[,1] = rep('custom_list', dim(factor_labeling)[1])
  factor_labeling = factor_labeling[,1,drop=F]
  colnames(factor_labeling) = c('type')
  factor_list = unique(factor_labeling[,1])
  DE_genes = rownames(factor_labeling)
  
  
  # get gene lengths
  gene_lengths = gene_lengths
  gene_lengths = as.matrix(gene_lengths[,1,drop=F])
  
  
  # get background gene list
  background = background
  background.gene_ids = rownames(background)
  background.gene_ids = unique(c(background.gene_ids, DE_genes))
  sample_set_gene_ids = background.gene_ids
  
  
  # parse GO assignments
  GO_info = GO_info
  GO_info_listed = apply(GO_info, 1, function(x) unlist(strsplit(x,',')))
  names(GO_info_listed) = rownames(GO_info)
  get_GO_term_descr =  function(x) {
    d = 'none';
    go_info = GOTERM[[x]];
    if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
    return(d);
  }
  
  
  #organize go_id -> list of genes
  GO_to_gene_list = list()
  for (gene_id in intersect(names(GO_info_listed), sample_set_gene_ids)) {
    go_list = GO_info_listed[[gene_id]]
    for (go_id in go_list) {
      GO_to_gene_list[[go_id]] = c(GO_to_gene_list[[go_id]], gene_id)
    }
  }
  
  
  # GO-Seq protocol: build pwf based on ALL DE features
  missing_gene_lengths = sample_set_gene_ids[! sample_set_gene_ids %in% rownames(gene_lengths)]
  if (length(missing_gene_lengths) > 0) {
    stop("Error, missing gene lengths for features: ", paste(missing_gene_lengths, collapse=', '))
  }
  sample_set_gene_lengths = gene_lengths[sample_set_gene_ids,]
  GO_info_listed = GO_info_listed[ names(GO_info_listed) %in% sample_set_gene_ids ]
  cat_genes_vec = as.integer(sample_set_gene_ids %in% rownames(factor_labeling))
  pwf=nullp(cat_genes_vec, bias.data=sample_set_gene_lengths,plot.fit = FALSE)
  rownames(pwf) = sample_set_gene_ids
  
  
  # perform functional enrichment testing for each category.
  for (feature_cat in factor_list) {
    message('Processing category: ', feature_cat)
    gene_ids_in_feature_cat = rownames(factor_labeling)[factor_labeling$type == feature_cat]
    cat_genes_vec = as.integer(sample_set_gene_ids %in% gene_ids_in_feature_cat)
    pwf$DEgenes = cat_genes_vec
    res = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=FALSE,method ="Wallenius")
    ## over-represented categories:
    pvals = res$over_represented_pvalue
    pvals[pvals > 1 - 1e-10] = 1 - 1e-10
    q = qvalue(pvals)
    res$over_represented_FDR = q$qvalues
    go_enrich_filename = paste(i,sep="")
    result_table = res[res$over_represented_pvalue<=0.05,]
    descr = unlist(lapply(result_table$category, get_GO_term_descr))
    result_table$go_term = descr;
    result_table$list = paste(i)
    result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) { 
      gene_list = GO_to_gene_list[[x]]
      gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
      paste(gene_list, collapse=', ');
      
      
      
    }) )
    
    result_table=filter(result_table ,over_represented_pvalue <= 0.05)
    top=order(result_table$over_represented_pvalue,decreasing = FALSE)[1:10]
    result_table=result_table[top,]
    result=result_table[order(result_table$over_represented_pvalue),]
    sname=go_enrich_filename
  }
  addWorksheet(gowb, sname)
  writeData(gowb, sname, result)
  }

saveWorkbook(gowb, file = "Goseq.xlsx", overwrite = TRUE) 
}

Goplot <- function(){
  sheets <- openxlsx::getSheetNames("Goseq.xlsx")
  sigGO=do.call(rbind, lapply(sheets,openxlsx::read.xlsx,xlsxFile="Goseq.xlsx"))

  sigGO$`-logp`=-log10(sigGO$over_represented_pvalue)
ggdotchart(sigGO, 
             x = "term", 
             y = "-logp",
             ylab = "-log10(pvalue)",
             xlab = "Ontology Terms",
             title=paste("GO Analysis:","\nP < 0.01"),
             color = "list",
             fill = "list",
             palette = "jco", 
             #sorting = "descending",                       
             rotate = TRUE,
             yscale = "none" ,    
             group = "ontology"   ,                     
             ggtheme = theme_pubr(legend = "right"),
             dot.size = "numDEInCat"
  )+theme_cleveland() +facet_wrap(ontology~., drop=TRUE,strip.position = "top",scales = "free_x")+theme(text = element_text(family="serif", face="bold")) }

Ma <- function(contrast=contrast,fc=2,fdr=0.01,top=5){
  res <- read.xlsx(xlsxFile = "DEG.xlsx",sheet = paste(contrast))

ggmaplot(res, main = paste(contrast),
         fdr = fdr, fc = fc, size = 0.4,
         palette = c("#B31B21", "#1465AC", "darkgray"),
         genenames = as.vector(res$SYMPOL),
         legend = "top", 
         top = top,
         font.label = c("bold", 11), label.rectangle = TRUE,
         font.legend = "bold",repel =TRUE,
         font.main = "bold",ggtheme = ggplot2::theme_minimal())
}

VOL <- function(contrast, selectLab = "",fdr=0.05,fc=2 ){
  
  res <- read.xlsx(xlsxFile = "DEG.xlsx",sheet = paste(contrast))
  EnhancedVolcano(res,
                  FCcutoff = fc,
                  pCutoff = fdr,
                  drawConnectors = TRUE,
                  typeConnectors ="open",
                  endsConnectors = "last",
                  widthConnectors = 0.75,
                  lab = res$SYMPOL,
                  pointSize = 1,
                  x = 'log2FoldChange',
                  boxedLabels = TRUE,
                  title = paste(contrast),
                  subtitle = "",
                  y = 'padj',
  )
}


NETscat <- function(Trait){
for (i in Trait) {
  #i="Rice"
module=row.names(corr[order(corr[,paste(i)],decreasing = TRUE),])[1]
Traitdat = as.data.frame(datTraits[paste(i)]);
names(Traitdat) = paste(i)
modNames = names(MEs)
geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datExpr, Traitdat, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(Traitdat), sep="");
names(GSPvalue) = paste("p.GS.", names(Traitdat), sep="");

column = match(module, modNames);
moduleGenes = moduleLabels==module;
gmm=data.frame(MM=c(abs(geneModuleMembership[moduleGenes, column])),GT=c(abs(geneTraitSignificance[moduleGenes, 1])))

gg=ggscatter(gmm,x = "MM",y="GT",xlab=paste("Module Membership in", "module", module),ylab =paste( "Gene significance for",i) ,conf.int = TRUE,add = "reg.line",add.params = list(color = "blue", fill = "lightgray"))+theme_get()+theme(text = element_text(family="serif", face="bold"))
cat("### ",paste("Module Membership of Module", module, "with Trait" ,i),"\n")
par(margin(t=4))
print(gg)
plot.new()
#print(paste("Figure",i,":",GOcap,sep = ""))
dev.off()
cat('\n\n')
}}

dectoutlier <- function(countData,rep){
  
  pca <- PcaHubert(t(vstcounts[,c(grep(paste(rep),colnames(vstcounts)))]),alpha = 0.7) 
  which(pca@flag=='FALSE')
  
}

outlier <- function(x,y){
  outliers <- foreach(m=unique(y)) %do%
  dectoutlier(x,m)
  outliers<-as.data.frame(unlist(outliers))
  row.names(outliers)
}

KeepDrop = function(data=df,cols="var",newdata=df2,drop=1) {
  
  # Double Quote Output Dataset Name
  t = deparse(substitute(newdata))
  
  
  # Drop Columns
  if(drop == 1){
    newdata = data [ , !(names(data) %in% scan(textConnection(cols), what="", sep=" ",quiet = TRUE))]}
  
  # Keep Columns
  else {
    newdata = data [ , names(data) %in% scan(textConnection(cols), what="", sep=" ",quiet = TRUE)]}
  assign(t, newdata, .GlobalEnv)
  
}

NETGO <- function(module,gene_lengths,background,Go_info,pcuttof = 0.05){
  
  
  gowb<- createWorkbook()
  for (i in module) {
  
  
    # Select modules
   
    # Select module probes
    probes=names(as.data.frame(datExpr))
    inModule=is.finite(match(moduleLabels,i));
    modProbes=probes[inModule];
    modGenes=modProbes[inModule]
    # Select the corresponding Topological Overlap
    modTOM =TOM[inModule, inModule];
    nTop =30;
    IMConn =softConnectivity(datExpr[,modProbes],power = power,verbose = 0,corFnc = "bicor",type = networktype);
    top = (rank(-IMConn) <= nTop)
    ntop=modTOM[top, top]
    dimnames(modTOM) = list(modProbes, modProbes)
    dir.create("cytoscape")
    # Export the network into edge and node list files Cytoscape can read
    cyt = exportNetworkToCytoscape(ntop,
                                   edgeFile = paste("cytoscape/CytoscapeInput-edges", paste(i, collapse="-"), ".txt", sep=""),
                                   nodeFile = paste("cytoscape/CytoscapeInput-nodes", paste(i, collapse="-"), ".txt", sep=""),
                                   weighted = TRUE,includeColNames = TRUE,
                                   threshold = 0,
                                   nodeNames = modProbes,
                                   altNodeNames = annot[modProbes,],
                                   nodeAttr = labels2colors(annot[modProbes,]))
    
    
    
  # capture list of genes for functional enrichment testing
  factor_labeling = data.frame(modProbes,modProbes,row.names = TRUE)
  factor_labeling[,1] = rep('custom_list', dim(factor_labeling)[1])
  factor_labeling = factor_labeling[,1,drop=F]
  colnames(factor_labeling) = c('type')
  factor_list = unique(factor_labeling[,1])
  DE_genes = rownames(factor_labeling)
  
  
  # get gene lengths
  gene_lengths = gene_lengths
  gene_lengths = as.matrix(gene_lengths[,1,drop=F])
  
  
  # get background gene list
  background = background
  background.gene_ids = rownames(background)
  background.gene_ids = unique(c(background.gene_ids, DE_genes))
  sample_set_gene_ids = background.gene_ids
  
  
  # parse GO assignments
  GO_info = GO_info
  GO_info_listed = apply(GO_info, 1, function(x) unlist(strsplit(x,',')))
  names(GO_info_listed) = rownames(GO_info)
  get_GO_term_descr =  function(x) {
    d = 'none';
    go_info = GOTERM[[x]];
    if (length(go_info) >0) { d = paste(Ontology(go_info), Term(go_info), sep=' ');}
    return(d);
  }
  
  
  #organize go_id -> list of genes
  GO_to_gene_list = list()
  for (gene_id in intersect(names(GO_info_listed), sample_set_gene_ids)) {
    go_list = GO_info_listed[[gene_id]]
    for (go_id in go_list) {
      GO_to_gene_list[[go_id]] = c(GO_to_gene_list[[go_id]], gene_id)
    }
  }
  
  
  # GO-Seq protocol: build pwf based on ALL DE features
  missing_gene_lengths = sample_set_gene_ids[! sample_set_gene_ids %in% rownames(gene_lengths)]
  if (length(missing_gene_lengths) > 0) {
    stop("Error, missing gene lengths for features: ", paste(missing_gene_lengths, collapse=', '))
  }
  sample_set_gene_lengths = gene_lengths[sample_set_gene_ids,]
  GO_info_listed = GO_info_listed[ names(GO_info_listed) %in% sample_set_gene_ids ]
  cat_genes_vec = as.integer(sample_set_gene_ids %in% rownames(factor_labeling))
  pwf=nullp(cat_genes_vec, bias.data=sample_set_gene_lengths,plot.fit = FALSE)
  rownames(pwf) = sample_set_gene_ids
  
  
  # perform functional enrichment testing for each category.
  for (feature_cat in factor_list) {
    message('Processing category: ', feature_cat)
    gene_ids_in_feature_cat = rownames(factor_labeling)[factor_labeling$type == feature_cat]
    cat_genes_vec = as.integer(sample_set_gene_ids %in% gene_ids_in_feature_cat)
    pwf$DEgenes = cat_genes_vec
    res = goseq(pwf,gene2cat=GO_info_listed, use_genes_without_cat=FALSE,method ="Wallenius",)
    ## over-represented categories:
    pvals = res$over_represented_pvalue
    pvals[pvals > 1 - 1e-10] = 1 - 1e-10
    q = qvalue(pvals)
    res$over_represented_FDR = q$qvalues
    go_enrich_filename = paste("Module",i, sep='')
    result_table = res[res$over_represented_pvalue<=0.05,]
    descr = unlist(lapply(result_table$category, get_GO_term_descr))
    result_table$go_term = descr;
    result_table$list = paste(i)
    result_table$gene_ids = do.call(rbind, lapply(result_table$category, function(x) { 
      gene_list = GO_to_gene_list[[x]]
      gene_list = gene_list[gene_list %in% gene_ids_in_feature_cat]
      paste(gene_list, collapse=', ');
      
      
      
    }) )
    result_table=filter(result_table ,over_represented_pvalue <= 0.05)
    top=order(result_table$over_represented_pvalue,decreasing = FALSE)[1:10]
    result_table=result_table[top,]
    result=result_table[order(result_table$over_represented_pvalue),]
    sname=go_enrich_filename
  }
  addWorksheet(gowb, sname)
  writeData(gowb, sname, result)
  }
  
  
  saveWorkbook(gowb, file = "ModulesGoseq.xlsx", overwrite = TRUE) 
}

Goplot2 <- function(){
  
  sheets <- openxlsx::getSheetNames("ModulesGoseq.xlsx")
  sigGO=do.call(rbind, lapply(sheets,openxlsx::read.xlsx,xlsxFile="ModulesGoseq.xlsx"))
  
  sigGO$`-logp`=-log10(sigGO$over_represented_pvalue)
  ggdotchart(sigGO, 
             x = "term", 
             y = "-logp",
             ylab = "-log10(pvalue)",
             xlab = "Ontology Terms",
             title=paste("GO Analysis:","\nP < 0.01"),
             color = "list",
             fill = "list",
             palette = "jco", 
             #sorting = "descending",                       
             rotate = TRUE,
             yscale = "none" ,    
             group = "ontology"   ,                     
             ggtheme = theme_pubr(legend = "right"),
             dot.size = "numDEInCat"
  )+theme_cleveland() +facet_wrap(ontology~., drop=TRUE,strip.position = "top",scales = "free_x") }
   
  

cureannot <- function(){
  nokeep <- data.frame(x=annot[grep(paste("hypothetical"),annot$Name,ignore.case = TRUE),])
  nokeep2 <- data.frame(x=annot[grep(paste("uncharacteri"),annot$Name,ignore.case = TRUE),])
nokeep <- rbind(nokeep,nokeep2)  
  
  
  return(nokeep$x)
  }
  
Plotnetwork <- function(module){
   
   for (i in module ) {
      
     result_table=read.xlsx(xlsxFile = "ModulesGoseq.xlsx",sheet = paste("Module",i,sep = ""))
     sigGO=filter(result_table ,over_represented_pvalue <= 0.05)
     top=order(sigGO$over_represented_pvalue,decreasing = FALSE)[1:30]
     sigGO=sigGO[top,]
     sigGO$`-logp`=-log10(sigGO$over_represented_pvalue) 
     sigGO=sigGO[c(1,2,6,7,11)] %>% 
       separate_rows(gene_ids, sep = ",") %>%  
       group_by(gene_ids) 
     sigGO=sigGO[order(sigGO$category,decreasing = TRUE), ]
     sigGO=sigGO[complete.cases(sigGO),] %>% rename(ID=gene_ids)
     sigGO$ID=gsub(" ","",sigGO$ID)
     sigGO=sigGO[!duplicated(sigGO$ID),] 
    anno=annot %>% rownames_to_column("ID")
    edge=read.delim(paste("cytoscape/CytoscapeInput-edges",paste(i,collapse="-"),".txt",sep = ""),)
    edge2=order(edge$fromNode,decreasing = TRUE)["CSG00000008064"]
    edge=as.matrix(edge[edge2,])
    edge=as.data.frame(edge[complete.cases(edge), ]) 
    node=data.frame(ID=unique(unlist(edge[c(1,2)],use.names =FALSE))) %>% left_join(sigGO,by="ID") %>% left_join(anno,by="ID")
    node[is.na(node)]="no record"
    Net2=tbl_graph(edges = edge, nodes = node ,directed = FALSE)
    Net2=as_tbl_graph(Net2)
    
    
      gg=Net2 %>% activate(nodes) %>% 
      ggraph(layout = "graphopt") + geom_edge_link(width = 1, colour = "lightgray") +
      geom_node_point(size = 8 ,aes(colour = term,shape = ontology,na.rm = TRUE))+
      
      geom_node_text(aes(label = Name), repel = TRUE)+
      ggtitle(paste("Module",i))+theme_graph()
      
    cat("### Module",i,"\n") 
    print(gg)
    plot.new()
    #print(paste("Figure",i,":",GOcap,sep = ""))
    dev.off()
    cat('\n\n')
    }
    
   
    
}  

mostVar <- function(data, n, i_want_most_var = TRUE) {
  data.var <- apply(data, 1, stats::var)
  data[order(data.var, decreasing = i_want_most_var)[1:n],] 
}











  


