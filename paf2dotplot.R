#!/usr/bin/env Rscript

## Make Dot Plot with Percent Divergence on color scale
rm(list=ls())
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(plotly))
library(gridExtra)
library(viridis)

plotdotplot<-function(alignments, opt, xlab, ylab){
  # Fixes for PAF
  # Some measure of similarity - need to check on this
  alignments$percentID = alignments$numResidueMatches / alignments$lenAln
  queryStartTemp = alignments$queryStart
  # Flip starts, ends for negative strand alignments
  alignments$queryStart[which(alignments$strand == "-")] = alignments$queryEnd[which(alignments$strand == "-")]
  alignments$queryEnd[which(alignments$strand == "-")] = queryStartTemp[which(alignments$strand == "-")]
  
  #Only keep alignments in which queryLen and refLen are at least least minseqlength 
  alignments=alignments[which(alignments$queryLen > opt$minseqlength & alignments$refLen > opt$minseqlength),]
  rm(queryStartTemp)
  cat(paste0("paffile:",opt$paf_filename,"\n"))
  cat(paste0("\nNumber of alignments: ", nrow(alignments),"\n"))
  cat(paste0("Number of query sequences: ", length(unique(alignments$queryID)),"\n"))
  cat(paste0("Number of Reference sequences: ", length(unique(alignments$refID)),"\n"))
  
  # sort by ref chromosome sizes, keep top X chromosomes OR keep specified IDs
  if(is.null(opt$refIDs)){
    chromMax = tapply(alignments$refEnd, alignments$refID, max)
    if(is.null(opt$keep_ref)){
      opt$keep_ref = length(chromMax)
    }
    refIDsToKeepOrdered = names(sort(chromMax, decreasing = T)[1:opt$keep_ref])
    alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered),]
    
  } else {
    refIDsToKeepOrdered = unlist(strsplit(opt$refIDs, ","))
    alignments = alignments[which(alignments$refID %in% refIDsToKeepOrdered),]
  }
  
  # filter queries by alignment length, for now include overlapping intervals
  queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
  alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]
  
  # filter alignment by length
  alignments = alignments[which(alignments$lenAln > opt$min_align),]
  
  # re-filter queries by alignment length, for now include overlapping intervals
  queryLenAgg = tapply(alignments$lenAln, alignments$queryID, sum)
  alignments = alignments[which(alignments$queryID %in% names(queryLenAgg)[which(queryLenAgg > opt$min_query_aln)]),]
  
  cat(paste0("\nAfter filtering... Number of alignments: ", nrow(alignments),"\n"))
  cat(paste0("After filtering... Number of query sequences: ", length(unique(alignments$queryID)),"\n\n"))
  cat(paste0("Number of Reference sequences: ", length(unique(alignments$refID)),"\n"))
  summary(alignments$queryLen)

  # sort df on ref
  alignments$refID = factor(alignments$refID, levels = refIDsToKeepOrdered) # set order of refID
  alignments = alignments[with(alignments,order(refID,refStart)),]
  chromMax = tapply(alignments$refEnd, alignments$refID, max)

  # make new ref alignments for dot plot
  if(length(levels(alignments$refID)) > 1){
    alignments$refStart2 = alignments$refStart + sapply(as.character(alignments$refID), function(x) ifelse(x == names((chromMax))[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
    alignments$refEnd2 = alignments$refEnd +     sapply(as.character(alignments$refID), function(x) ifelse(x == names((chromMax))[1], 0, cumsum(as.numeric(chromMax))[match(x, names(chromMax)) - 1]) )
  } else {
    alignments$refStart2 = alignments$refStart
    alignments$refEnd2 = alignments$refEnd
  }

  ## queryID sorting step 1/2
  # sort levels of factor 'queryID' based on longest alignment
  alignments$queryID = factor(alignments$queryID, levels=unique(as.character(alignments$queryID))) 
  queryMaxAlnIndex = tapply(alignments$lenAln,
                          alignments$queryID,
                          which.max,
                          simplify = F)
  alignments$queryID = factor(alignments$queryID, levels = unique(as.character(alignments$queryID))[order(mapply(
    function(x, i)
      alignments$refStart2[which(i == alignments$queryID)][x],
    queryMaxAlnIndex,
    names(queryMaxAlnIndex)
  ))])

  ## queryID sorting step 2/2
  ## sort levels of factor 'queryID' based on longest aggregrate alignmentst to refID's
  # per query ID, get aggregrate alignment length to each refID 
  queryLenAggPerRef = sapply((levels(alignments$queryID)), function(x) tapply(alignments$lenAln[which(alignments$queryID == x)], alignments$refID[which(alignments$queryID == x)], sum) )
  if(length(levels(alignments$refID)) > 1){
    queryID_Ref = apply(queryLenAggPerRef, 2, function(x) rownames(queryLenAggPerRef)[which.max(x)])
  } else {queryID_Ref = sapply(queryLenAggPerRef, function(x) names(queryLenAggPerRef)[which.max(x)])}
  # set order for queryID
  alignments$queryID = factor(alignments$queryID, levels = (levels(alignments$queryID))[order(match(queryID_Ref, levels(alignments$refID)))])

  #  flip query starts stops to forward if most align are in reverse complement
  queryRevComp = tapply(alignments$queryEnd - alignments$queryStart, alignments$queryID, function(x) sum(x)) < 0
  queryRevComp = names(queryRevComp)[which(queryRevComp)]
  queryMax = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), max)
  names(queryMax) = levels(alignments$queryID)
  alignments$queryStart[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryStart[which(alignments$queryID %in% queryRevComp)] + 1
  alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] = queryMax[match(as.character(alignments$queryID[which(alignments$queryID %in% queryRevComp)]), names(queryMax))] - alignments$queryEnd[which(alignments$queryID %in% queryRevComp)] + 1

  ## make new query alignments for dot plot
  # subtract queryStart and Ends by the minimum alignment coordinate + 1
  queryMin = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), min)
  names(queryMin) = levels(alignments$queryID)
  alignments$queryStart = as.numeric(alignments$queryStart - queryMin[match(as.character(alignments$queryID),names(queryMin))] + 1)
  alignments$queryEnd = as.numeric(alignments$queryEnd - queryMin[match(as.character(alignments$queryID),names(queryMin))] + 1)

  queryMax = tapply(c(alignments$queryEnd, alignments$queryStart), c(alignments$queryID,alignments$queryID), max)
  names(queryMax) = levels(alignments$queryID)
  alignments$queryStart2 = alignments$queryStart + sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )
  alignments$queryEnd2 = alignments$queryEnd +     sapply(as.character(alignments$queryID), function(x) ifelse(x == names(queryMax)[1], 0, cumsum(queryMax)[match(x, names(queryMax)) - 1]) )

  # get mean percent ID per contig
  #   calc percent ID based on on-target alignments only
  if(opt$on_target & length(levels(alignments$refID)) > 1){
    alignments$queryTarget = queryID_Ref[match(as.character(alignments$queryID), names(queryID_Ref))]
    alignmentsOnTarget = alignments[which(as.character(alignments$refID) == alignments$queryTarget),]
    scaffoldIDmean = tapply(alignmentsOnTarget$percentID, alignmentsOnTarget$queryID, mean)
    alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
    alignments$percentIDmean[which(as.character(alignments$refID) != alignments$queryTarget)] = NA
  } else{
    scaffoldIDmean = tapply(alignments$percentID, alignments$queryID, mean)
    alignments$percentIDmean = as.numeric(scaffoldIDmean[match(as.character(alignments$queryID), names(scaffoldIDmean))])
  }
  
  print(dim(alignments))
  # plot
  yTickMarks = tapply(alignments$queryEnd2, alignments$queryID, max)
  options(warn = -1) # turn off warnings
  if (opt$similarity) {
    gp = ggplot(alignments) +
      geom_point(
        mapping = aes(x = refStart2, y = queryStart2, color = percentIDmean),
        size = 1
      ) +
      geom_point(
        mapping = aes(x = refEnd2, y = queryEnd2, color = percentIDmean),
        size = 1
      ) +
      geom_segment(
        aes(
          x = refStart2,
          xend = refEnd2,
          y = queryStart2,
          yend = queryEnd2,
          color = percentIDmean,
          text = sprintf(
            'Query ID: %s<br>Query Start Pos: %s<br>Query End Pos: %s<br>Target ID: %s<br>Target Start Pos: %s<br>Target End Pos: %s<br>Length: %s kb',
            queryID,
            queryStart,
            queryEnd,
            refID,
            refStart,
            refEnd,
            round(lenAln / 1000, 1)
          )
        ),
        size = 2
      ) +
      scale_x_continuous(breaks = cumsum(as.numeric(chromMax)),
                        labels = levels(alignments$refID)) +
      theme_bw() +
      theme(text = element_text(size = 8)) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = 10, angle = 15),
        axis.text.x = element_text(size = 10, angle = 45, hjust=1)
      ) +
      scale_y_continuous(breaks = yTickMarks, labels = substr(levels(alignments$queryID), start = 1, stop = 35)) +
      #scale_y_continuous(breaks = yTickMarks, labels = NULL) +
      { if(opt$h_lines){ geom_hline(yintercept = yTickMarks,
                                    color = "grey60",
                                    size = .1) }} +
      scale_color_viridis_c(option = "D", limits = c(0, 0.3)) +
      #scale_color_distiller(palette = "Spectral", limits = c(0, 0.5)) +
      #scale_color_gradient(palette = "Spectral", low = "blue", high = "red", limits = c(0, 1)) +
      labs(color = "Mean Percent Identity (per query)", 
          title = paste0(   paste0("Post-filtering number of alignments: ", nrow(alignments),"\t\t\t\t"),
                            paste0("minimum alignment length (-m): ", opt$min_align,"\n"),
                            paste0("Post-filtering number of queries: ", length(unique(alignments$queryID)),"\t\t\t\t\t\t\t\t"),
                            paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln)
          )) +
      xlab(xlab) +
      ylab(ylab)
  } else {
    gp = ggplot(alignments) +
      geom_point(mapping = aes(x = refStart2, y = queryStart2, color = refChr),
                size = 1) +
      geom_point(mapping = aes(x = refEnd2, y = queryEnd2, color = refChr),
                size = 1) +
      geom_segment(aes(
        x = refStart2,
        xend = refEnd2,
        y = queryStart2,
        yend = queryEnd2,
        text = sprintf(
          'Query ID: %s<br>Query Start Pos: %s<br>Query End Pos: %s<br>Target ID: %s<br>Target Start Pos: %s<br>Target End Pos: %s<br>Length: %s kb',
          queryID,
          queryStart,
          queryEnd,
          refID,
          refStart,
          refEnd,
          round(lenAln / 1000, 1)
        )
      ),
      size = 2) +
      scale_x_continuous(breaks = cumsum(as.numeric(chromMax)),
                        labels = levels(alignments$refID)) +
      theme_bw() +
      theme(text = element_text(size = 8)) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.text.y = element_text(size = 10, angle = 15),
        axis.text.x = element_text(size = 10, angle = 45, hjust=1)
      ) +
      #scale_y_continuous(breaks = yTickMarks, labels = NULL) +
      scale_y_continuous(breaks = yTickMarks, labels = substr(levels(alignments$queryID), start = 1, stop = 35)) + 
      { if(opt$h_lines){ geom_hline(yintercept = yTickMarks,
                                    color = "grey60",
                                    size = .1) }} +
      labs(color = "Target Chromosome", 
         title = paste0(   paste0("Post-filtering number of alignments: ", nrow(alignments),"\t\t\t\t"),
                           paste0("minimum alignment length (-m): ", opt$min_align,"\n"),
                           paste0("Post-filtering number of queries: ", length(unique(alignments$queryID)),"\t\t\t\t\t\t\t\t"),
                           paste0("minimum query aggregate alignment length (-q): ", opt$min_query_aln)
         )) +
      xlab(xlab) +
      ylab(ylab)
  }
  gp
}


opt=list()
opt$keep_ref=1000
opt$min_query_aln=900000
opt$min_align=100000
opt$minseqlength=900000
opt$on_target=FALSE
opt$similarity=FALSE
opt$h_lines=TRUE
opt$groupsfile=NULL
#opt$paf_filename="paffiles/SOFF_LAPpurple__vs__SSPO_AP85441/SOFF_LAPpurple__vs__SSPO_AP85441.paf"
opt$paf_filename="paffiles/SOFF_LAPpurple__vs__SSPO_NPX/SOFF_LAPpurple__vs__SSPO_NPX.paf"
opt$chr_names="chrnames.txt"
opt$workdir="/data/diriano/sugarcaneGenome/"

if (!is.null(opt$groupsfile)){
 chrgroups=read.table(opt$groupsfile,header = FALSE)
 colnames(chrgroups)<-c('chr','sp','seqid')
 for (chr in unique(chrgroups$chr)){
  for (sp in unique(chrgroups[which(chrgroups$chr == chr), 'sp'])){
    nseqs=nrow(chrgroups[which(chrgroups$chr==chr & chrgroups$sp == sp),])
    print(paste(chr,sp,nseqs,sep=" "))
  } 
 }
} else{
  aligns = read.table(paste(opt$workdir,opt$paf_filename,sep=''), stringsAsFactors = F, fill = T)
  colnames(aligns)[1:12] = c("queryID","queryLen","queryStart","queryEnd","strand","refID","refLen","refStart","refEnd","numResidueMatches","lenAln","mapQ")
  aligns$queryID <- toupper(aligns$queryID)
  aligns$refID <- toupper(aligns$refID)
  aligns2=aligns[which(aligns$queryLen >= opt$minseqlength),]
  dim(aligns2)
  if (!is.null(opt$chr_names)){
    chrnames=read.delim(paste(opt$workdir,opt$chr_names,sep=""), stringsAsFactors = FALSE,header=FALSE)
    colnames(chrnames)<-c('contig','chrname1','chrname2')
    chrnames$contig<-toupper(chrnames$contig)
    merged_df <- merge(aligns2, chrnames[, c("contig", "chrname2")], by.x = "queryID", by.y='contig', all.x = TRUE)
    dim(merged_df)
    aligns2 = merged_df
    colnames(aligns2)[19]=c("queryChr")
    rm(merged_df)
    merged_df <- merge(aligns2, chrnames[, c("contig", "chrname2")], by.x = "refID", by.y='contig', all.x = TRUE)
    dim(merged_df)
    aligns2 = merged_df
    colnames(aligns2)[20]=c("refChr")
    rm(merged_df)
    aligns2$queryID2=paste(aligns2$queryID,aligns2$queryChr,sep="__")
    aligns2$refID2=paste(aligns2$refID,aligns2$refChr,sep="__")
    aligns2$queryID=aligns2$queryID2
    aligns2$queryID2=NULL
    aligns2$refID=aligns2$refID2
    aligns2$refID2=NULL
    head(aligns2)
  }
  
  if(nrow(aligns2) > 0 ){
    
    # Remove the ".paf" extension
    compseqs <- gsub("\\.paf$", "", basename(opt$paf_filename))
    # Split the string at "__vs__"
    compseqs <- strsplit(compseqs, "__vs__")[[1]]
    compseqs[2]
    gp=plotdotplot(aligns2, opt,xlab=paste('Target:',compseqs[1]), ylab=paste('Query:',compseqs[2]))
    dotplotfile_id=paste(opt$workdir,opt$paf_filename,".pdf",sep='')
    ggsave(filename = dotplotfile_id, plot= gp, width = 30, height = 20, units = "in", dpi = 600)
  }
}

