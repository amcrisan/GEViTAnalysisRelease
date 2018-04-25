#A little helper function for this part of data cleaning
maphow<-function(x=NULL,howList = NULL){
  return(str_extract(howList[as.numeric(x)+1,]$name,pattern="_[A-Za-z]+") %>% gsub("_","",.))
}


# NOTE TO SELF - NO ERROR CORRECTING YET, EXPECTS PERFECT OR WILL BREAK

cleanCodes<-function(tmp=NULL,img=NULL){
  # seperate high level whats
  what <- filter(tmp,grepl("what",name))
  how1 <- filter(tmp,grepl("how",name))
  
  
  #parsing and redefining the what structure
  what<-mutate(what,name=gsub("what_","",name))
  
  tmp <- data.table::setDT(what)[, list(Values = unlist(strsplit(name, "\\/"))), by = id] # splitting the values by `,` for each `ID`
  tmp[, Var := paste0("v", seq_len(.N)), by = id] # Adding the `Var` variable
  tmp<-dcast(tmp, id~ Var, fill = NA_character_, value.var = "Values") # decasting
  
  colnames(tmp)<-c("id",paste("WhatLvl",(1:(ncol(tmp)-1)),sep=""))
  what<-merge(x=what,y=tmp,by="id")
  
  #structure into multiple columns, adding the attributes as lower level how
  how1<-how1 %>%
    mutate(howIndx = strsplit(as.character(parts), ",")) %>% 
    tidyr::unnest(howIndx) %>%
    mutate(name = gsub("how_","",name))
  
  
  what$HowLvl1<- plyr::mapvalues(from=how1$howIndx,to=how1$name,what$id)
  
  topo<-what %>%
    mutate(basename=rep(fName,nrow(what))) %>%
    mutate(HowLvl2 = strsplit(as.character(attributes), ",")) %>% 
    tidyr::unnest(HowLvl2) %>%
    select(contains("What"),contains("How"),contains("basename"),contains("submitterID"))
  
  return(topo)
  
}


#retrieve and format data
formatData<-function(ids = NULL){
  #limitation : extract resuts as sets of 100. There's a limit to how many
  #ids can be queried. 1000 is too many, so I used 100 because it makes for nice numbers
  
  if(length(ids) < 200){
    #process all at once
    numVals<-c(1,length(ids))
  }else{
    #break it up into chunks
    numVals<-seq(from=0,to=length(ids),by=200)
    
    if(tail(numVals, n=1)< length(ids)){
      numVals<-c(numVals,tail(numVals, n=1) + (length(ids) %% 200))
    }
  }
  allData<-data.frame(PMID=NULL,
                      YearPub=NULL,
                      Journal=NULL,
                      Authors=NULL,
                      Title=NULL,
                      Abstract=NULL,
                      meshTerms=NULL)
  
  for(i in 1:(length(numVals)-1)){
    #to make this faster, form new query on ID run in parallel
    start = numVals[i]+1
    end = numVals[i + 1]
    pubResults<-EUtilsGet(paste0(ids[start:end],collapse = ","),type="efetch",db="pubmed")
    
    #collaposing the author list for the data frame
    authors<-sapply(pubResults@Author,function(x){
      apply(x,1,function(y){sprintf("%s, %s",y[1],y[2])}) %>% paste0(.,collapse=";")
    })
    
    #putting all the data together into a data frame
    temp<-data.frame(PMID=pubResults@PMID,
                     YearPub=pubResults@YearPubmed,
                     Journal=pubResults@Title,
                     Authors=authors,
                     Title=pubResults@ArticleTitle,
                     Abstract=pubResults@AbstractText)
    
    
    temp$meshTerms<-pubResults@Mesh
    
    allData<-rbind(allData,temp)
  }
  
  return(allData)
}


# Re-computing tsne clusters - WILL WORK FOR SHINY APP
tsneClusters<-function(tsneCord = NULL, genEpitext_df = NULL, minPtsVal = 100){
  
  #get the clusters for the tsne perplexity 100 plot using HDBSCAN
  cl <- hdbscan(tsneCord[,c("comp1","comp2")], minPts = minPtsVal) #minimum cluster size of 50 documents
  tsneCord$tsneCluster<-cl$cluster
  
  print(cl)
  
  #name the clusters according to the most commonly used term in each cluster
  clusterMap <- c()
  
  for(cluster in unique(cl$cluster)){
    
    if(cluster == 0){
      clusterMap<-rbind(clusterMap,c(cluster,"Noise"))
      next()
    }
    clusterMembers<-tsneCord %>%
      filter(tsneCluster == cluster) %>%
      select(PMID)
    
    topWord<-genEpitext_df %>%
      filter(PMID %in% clusterMembers$PMID) %>%
      ungroup() %>%
      group_by(wordStemmed) %>%
      tally() %>%
      arrange(-nn) %>%
      top_n(3)
    
    clusterMap<-rbind(clusterMap,c(cluster,paste0(topWord$wordStemmed,collapse="-")))
  }
  
  #create new variable of re-assigned names
  tsneCord$tsneClusterNames <- plyr::mapvalues(from=clusterMap[,1],to=clusterMap[,2],tsneCord$tsneCluster)
  return(tsneCord)
}


nameClustersNice<-function(tsneCord = NULL, clustVals = NULL,genEpitext_df=NULL,top=3,returnList= TRUE){
  #name the clusters according to the most commonly used term in each cluster
  clusterMap <- c()
  
  for(cluster in unique(clustVals)){
    
    addWord<-NA
    #for numerical or text clusters
    if(cluster == 0 | cluster == "0"){
      #clusterMap<-rbind(clusterMap,c(cluster,"Noise"))
      #next()
      addWord = "Noise"
    }
    clusterMembers<-tsneCord %>%
      filter(tsneCluster == cluster) %>%
      select(PMID)
    
    topWord<-genEpitext_df %>%
      filter(PMID %in% clusterMembers$PMID) %>%
      ungroup() %>%
      group_by(wordStemmed) %>%
      tally() %>%
      arrange(-nn) %>%
      top_n(top)
    
    if(returnList){
        clusterMap<-rbind(clusterMap,c(cluster,list(bigrams=topWord$wordStemmed)))
    }else{
      if(!is.na(addWord)){
        clusterMap<-rbind(clusterMap,c(cluster,paste0(c(addWord,topWord$wordStemmed),collapse="-")))
      }else{
        clusterMap<-rbind(clusterMap,c(cluster,paste0(topWord$wordStemmed,collapse="-")))
      }
    }
  }
  
  #create new variable of re-assigned names
  tsneClusterNames <- plyr::mapvalues(from=clusterMap[,1],to=clusterMap[,2],tsneCord$tsneCluster)
  return(tsneClusterNames)
}

# Took this code from :https://github.com/cstubben/pmcXML/blob/master/R/pmcText.R
pmcText <-function(doc, sentence=TRUE ){
  
  z <- vector("list")
  
  ## Also in metadata..
  z[["Main title"]] <- xpathSApply(doc, "//front//article-title", xmlValue)
  
  # ABSTRACT 
  
  x <- paste( xpathSApply(doc, "//abstract[not(contains(@abstract-type,'summary'))]//p", xmlValue), collapse=" ")
  z[["Abstract"]] <- fixText(x) 
  
  ## author or executive summary...
  x<- paste( xpathSApply(doc, "//abstract[contains(@abstract-type,'summary')]//p", xmlValue), collapse=" ")  
  if(x != "") z[["Summary"]] <- fixText(x) 
  
  ## check for //body/p  instead of //body/sec/p
  
  intro <- xpathSApply(doc, "//body/p", xmlValue)
  intro <- gsub("\n", " ", intro)
  
  ## figures - run before removing any nested fig nodes
  f1 <- suppressMessages( pmcFigure(doc) )
  
  ## BODY - split into sections
  x <- getNodeSet(doc, "//body//sec")
  
  ## IF no  sections? - EID journal
  if(length(x) == 0){
    message("NOTE: No sections found, using Main")
    z[["Main"]] <-   fixText(intro)
  }else{
    if(length(intro) > 0){
      message("NOTE: adding untitled Introduction")
      intro <- gsub("\n", " ", intro)
      z[["Introduction"]] <-   fixText(intro)
    }   
    
    ## check for tables and figs within p tags  
    
    x1 <- xpathSApply(doc, "//sec/p/table-wrap", removeNodes)
    
    if(length(x1) > 0) message(paste("WARNING: removed", length(x1), "nested table tag"))
    x1 <- xpathSApply(doc, "//sec/p/fig", removeNodes)
    if(length(x1) > 0) message(paste("WARNING: removed", length(x1), "nested fig tag"))
    
    sec <- xpathSApply(doc, "//body//sec/title", xmlValue)
    n <- xpathSApply(doc, "//body//sec/title", function(y) length(xmlAncestors(y) ))
    path <- path.string(sec, n)
    
    y <- lapply(x, function(y) xpathSApply(y, "./p", xmlValue))
    
    ##LOOP through subsections
    for(i in 1: length(y) ){
      subT <- path[i]
      subT <- gsub("\\.$", "", subT)
      subT <- gsub("[; ]{3,}", "; ", subT)  # in case of "; ; ; "
      if(length(y[[i]]) > 0)  z[[ subT ]] <- fixText( y[[i]] )
    }
  } 
  
  
  ## ACKNOWLEDGEMENTS
  ack <- xpathSApply(doc, "//back//sec/title[starts-with(text(), 'Acknow')]/../p", xmlValue)
  
  if (length(ack) == 0)  ack <- xpathSApply(doc, "//back/ack/p", xmlValue)  # scientificWorldJournal
  if (length(ack) > 0)  z[["Acknowledgements"]]<- fixText(ack)
  
  # Funding
  funds <-  xpathSApply(doc, "//funding-statement", xmlValue)
  if (length(funds) > 0)  z[["Funding"]]<- fixText(funds )
  
  sec <- names(z)
  sec<- sec[!sec %in% c("Main title", "Abstract", "Summary", "Acknowledgements", "Funding")] 
  
  # Figures
  if(!is.null(f1)){
    z<- c(z, f1)
    z[["Figure caption"]] <- names(f1)
  }
  
  ##   # SPLIT sections (not title)
  if(sentence) z[-1] <- lapply(z[-1], splitP)
  
  z[["Section title"]] <- sec
  
  # add attributes
  attr(z, "id") <- attr(doc, "id")
  z
}

fixText <-function(x){
  
  x <- gsub("“|”|″", "\"", x)
  x <- gsub("’|′", "'", x)
  x <- gsub("−|–", "-", x)
  x <- gsub("∼", "~", x)
  x <- gsub("×", "x", x)
  x <- gsub("‥", "..", x)
  x <- gsub("ö", "o", x)
  x <- gsub("NA; ", "", x)  # section titles only
  x <- gsub("\n", " ", x) 
  x <- gsub("  *", " ", x) 
  x <- gsub("^ *", "", x) 
  x <- gsub(" *$", "", x) 
  x
}

pmcFigure <- function(doc,  attr = FALSE ){
  
  x <-  getNodeSet(doc, "//fig" )
  if (length(x) > 0) {
    
    ## should have label and caption (and caption with title and p)
    f1 <- sapply(x, xpathSApply, "./label", xmlValue)
    f1 <- gsub("[ .]+$", "", f1)
    
    # get caption title and paragraphs together since some caption titles are missing, in bold tags or have long caption title names that should be split 
    f2 <- sapply(x, xpathSApply, "./caption", xmlValue)
    
    ## get caption title to first :; or . (check if some caption/title missing a period?)
    
    z <-  lapply(f2, splitP, "[;:.]")
    cap <-  sapply(z, "[", 1)
    cap <- gsub("\\.$", "", cap) # drop period ??
    
    z <- lapply(z, function(x) paste(x[-1], collapse=" "))
    
    names(z) <- paste(f1, cap, sep=". ")
    
    if(attr) 
    {
      ## URL  
      ids <- sapply(x, xpathSApply, ".", xmlGetAttr, "id")
      
      ids <- paste("http://www.ncbi.nlm.nih.gov/pmc/articles", attr(doc, "id"), "figure", ids, sep="/")     
      txt <- pmcText2(doc)
      
      for(i in 1:length(x) ){
        message(paste(" ",  f1[i], ". ", cap[i], sep="" ))
        attr(z[[i]], "label") <-  f1[i]
        attr(z[[i]], "caption") <-  cap[i]
        attr(z[[i]], "file") <-  ids[i]
        # cite
        fig <- f1[i]  
        # how to match "Fig. 1, 3"
        
        fig <- gsub("Fig[^ ]*", "Fig[.ure]*", fig)
        cs <- searchP(txt, fig)
        if(!is.null(cs)){
          attr(z[[i]], "cite") <- fixText( as.vector( apply(cs, 1, paste, collapse=". ") )) 
        }else{
          message("No sentences citing " , fig)   
        }
      }
    }else{
      for(i in 1:length(z)) message(" ", names(z)[i])
      
    }
    z
  }else{
    NULL
  }
}


# split paragraphs into sentences... 
# could use sentDetect in old openNLP package (now Maxent_Sent_Token_Annotator)
## but that does not split on sentences ending with numbers like table 2 or  roman numbers I 

splitP <- function( z,  split= "[.?]"){
  # check if empty list (returned by xpathSApply)
  if(length(z) == 0){
    NULL
  }else{
    z <- gsub("\n", " ", z)
    # non-breaking spaces in figure 1 or table 1
    z <- gsub("\u00A0", " ", z)
    # some abbrevations may be after parentheses (so drop space if no other words end in that abbrev.)
    z <- gsub(".^",       ". ^",        z, fixed = TRUE)   # reference superscripts after period 
    z <- gsub("Suppl. ",   "SupplX.X ",     z, fixed = TRUE)   # Supplement
    z <- gsub("Supp. ",    "SuppX.X ",     z, fixed = TRUE)   # Supplement 
    z <- gsub("et al. ",  "et alX.X ",   z, fixed = TRUE)   # et al. 
    z <- gsub("et. al. ",  "et alX.X ",   z, fixed = TRUE)   # et. al. 
    z <- gsub(" no. ",    " noX.X ",     z, fixed = TRUE)   # acc no. 
    z <- gsub(" nos. ",   " nosX.X ",    z, fixed = TRUE)   # acc nos.
    z <- gsub(" No. ",    " NoX.X ",     z, fixed = TRUE) 
    
    z <- gsub("e.g. ",    "e.gX.X ",     z, fixed = TRUE)   # e.g.
    z <- gsub("(eg. ",    "(egX.X ",     z, fixed = TRUE)   # eg.  .. add paren to avoid words ending in eg.
    z <- gsub("i.e. ",    "i.eX.X ",     z, fixed = TRUE)   # i.e. 
    z <- gsub(" spp. ",   " sppX.X ",    z, fixed = TRUE)   # species
    z <- gsub(" sp. ",    " spX.X ",     z, fixed = TRUE) 
    z <- gsub("subsp. ", "subspX.X ",  z, fixed = TRUE)   # subspecies
    z <- gsub(" var. ",   " varX.X ",    z, fixed = TRUE)   # varieties
    z <- gsub(" bv. ",   " bvX.X ",    z, fixed = TRUE)   # biovars
    z <- gsub(" sv. ",   " svX.X ",    z, fixed = TRUE)   # serovars
    z <- gsub(" hr. ",   " hrX.X ",    z, fixed = TRUE)   # hours
    z <- gsub(" hrs. ",  " hrsX.X ",    z, fixed = TRUE)   # hours
    z <- gsub("i.n. ",   "i.nX.X ",    z, fixed = TRUE)  #intranasal i.n.
    z <- gsub("i.p. ",   "i.pX.X ",    z, fixed = TRUE) 
    
    z <- gsub(" ca. ",    " caX.X ",     z, fixed = TRUE)  # approx 
    z <- gsub("approx. ", "approxX.X ",  z, fixed = TRUE)
    z <- gsub(" vs. ",    " vsX.X ",     z, fixed = TRUE)  # vs.
    z <- gsub(" Dr. ",    " DrX.X ",     z, fixed = TRUE)  # Dr.
    z <- gsub(" Drs. ",    " DrsX.X ",     z, fixed = TRUE)  # Dr.
    z <- gsub(" Mr. ",    " MrX.X ",     z, fixed = TRUE)  # Mr.
    z <- gsub(" Mrs. ",    " MrsX.X ",     z, fixed = TRUE)  # Mrs.
    z <- gsub("cfu. ",    "cfuX.X ",     z, fixed = TRUE)  # cfu
    z <- gsub("c.f.u. ",    "c.f.uX.X ",     z, fixed = TRUE)  #
    z <-  gsub(" Ref. ", " RefX.X ", z, fixed = TRUE)
    
    ## Don't split if sentence starts with table or fig labels
    z <-  gsub("^(Table S?[0-9I]+)\\. ", "\\1X.X ", z)
    z <-  gsub("^(Figure S?[0-9I]+)\\. ", "\\1X.X ", z)
    z <-  gsub("^(Fig.? S?[0-9I]+)\\. ", "\\1X.X ", z)
    z <-  gsub("Fig. ", "FigX.X ", z, fixed = TRUE)
    
    ## single letter and period.  
    # dont split on single characters like B. subtilis
    z <- gsub("\\b([A-Za-z])\\. ([a-z])", "\\1X.X \\2", z)
    
    # or letters at start of sentence,  list of things.  A. stuff. B. more stuff
    z <- gsub("\\. ([A-Za-z])\\. ([A-Z])", ". \\1X.X \\2", z)
    
    # split sentences, 
    # z2 <- unlist( strsplit(z, split) )
    
    #  change ? or . to ?~ and .~ and split on "~"
    z2 <- strsplit( gsub( paste0("(", split, ") "),  "\\1~", z), "~")  
    z2 <- unlist(z2)
    
    
    ## remove placeholders
    z2 <- gsub("X.X",    ".", z2, fixed = TRUE)
    z2
  }  
}

path.string<-function(x,n){
  n2 <-length(n)
  z <- vector("list", n2)
  if(min(n) > 1) n <- n - min(n) + 1
  path<-""
  for(i in 1: n2){
    path[n[i] ] <- x[i]
    path <- path[1:n[i]]
    z[[i]] <- paste(path, collapse="; ")
  }
  unlist(z)
}
