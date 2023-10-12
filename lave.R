library(KernSmooth)
library(LSD)	
library(lattice)
library(raster)
#library(rasterVis)	#commented out 2023-02-07; rasterVis is no longer supported it appears.

## TODO: use $stop instad of $end to better interface w/ other bio suites like Sushi

### LAVE <--> EVAL
### Lave = wash. 
### Some tools to help EVALuate and wash assemblies (or explore pairwise interactions).


### EVALUATION MATRICES
### - Tools to interpret the evaluation tables from the SLURM-GEAR tools
### - TO ADD

### PAIRWISE ASSEMBLY ALIGNMENT VISUALIZATION
###   - USE: typcally used to visualize assembly alignments, discordant PE reads, or split alignments from long reads
###   - assembalign
###   - gapaln (really bed aligner)
###   - revcomppaf
###   - revcompbed
###   - revcomptignames

### PAIRED LOCUS PLOTS
###   - USE: typically used in a Hi-C like fashion, but w/ normal PE or split read data (using BED like inputs here)
###       - Can visualize hot spots loci for SVs such as those created during UR
###       - Can use to spot misassemblies and potential joins (new scaffold events)
###       - can also easily filter contact matrix given some cutoff to then visualize and make the most confident decision
###       - the contact matrix technically specifies an assembly graph of sorts, so filtering can really help viz high conf scaffolds
###       - for assembly improvement recommend viewing different SV-types separately DEL, DUP, INV, TRA
###       - one needs to be more careful when interpreting DEL/DUP/INV signals as misassembled areas than TRA (which given high MAPQ is relatively easy to interpret -- but still be careful)
###       - high mapq will be an important guide for misassembly generally speaking -- however looking at lower mapq and even 0 can help understand the asm structure as well
###   - pairedLociLevelPlot
###   - pairedLociLevelPlot2
###   - splitreadlevelplot (antiquated and now just wrapper for pairedLociLevelPlot)
###   - pairedLociStackedLevelPlot
###   - insertplot (TOFINISH)
###   - pairedlociscatters (Antiquated mostly)
###   - BELOW: PART OF NAMESPACE BUT NOT NEC INTENDED TO BE USED BY END-USER
###   - get_target_chr
###   - get_target_region
###   - get_target_insert_sizes
###   - prep_for_full_matrix
###   - get_xyz_vars
###   - get_color_code
###   - get_levelplot
###   - prep_paired_loci_values
###   - get_stacked_levelplot
###   - addpoly




### PAIRWISE ASSEMBLY ALIGNMENT VISUALIZATION
### CAN ALSO BE USED WITH DISCORDANT READS (RECOMMEND USING HIGH TRANSPARENCY TO LET HOTSPOTS BUILD UP AS SOLIDS)
### NOTE TO JMU: SEE EXAMPLES IN /Users/johnurban/Gerbi_Lab/Projects/Sciara-Genome/NoteBook/2018/08-Aug/genome/minimap2/falcScafs-v-canuScafs-all/miniasm/canu-falc-scaf-paf-explore....R on MBP

## DEFINE ASSEMBLY ALIGNMENT FUNCTION
assembalign <- function(falctigs=NA, canutigs=NA, paf=NA, canugaps=NA, falcgaps=NA, rc.falclist=NA, rc.canulist=NA, cadd=c(1,2,3), step=50e6, goncol=NA, bgoncol="black", tigcol=NA, btigcol="black"){
  ## paf need not be PAF object, but should be dataframe w/ following indexes: query, qstart, qend, target, tstart, tend, strand, and mapq --- mapq not currently used, but will have option soon
  ##  -- if do not have genome files (falctigs, canutigs) -- then they are learned from PAF in which case you need the DataFrame to be set up as a PAF: q,qlen,qend,strand,t,tlen,tstart,tend,...
  ## falctigs and canutigs should be renamed -- they specify data frames w/ 2 columns: seqnames, seqlengths
  ## canugaps and falcgaps would be better described as BED-like dataframes for the query and target -- i.e. where the BED ticks are put in the plt (top or bottom)
  ##  In some cases when the assemblies are the same and the links are from a dataset (e.g. discordant PE or splitreads), 
  ##    the same DataFrame can be given so the top and bottom mirror each other
  ##    Or one can highlight different features on the top v bottom.
  ##    TODO: In the future, perhaps diff features could be identified by colors
  ##      -- actually I can already do that modularly with the "gapaln" function which would be better called bedaln or featurealn
  ## the rc.lists are lists of sequences to orient as the revcomp (which specifies operations to do on paired alignments)
  ##    This is useful to untangle knots, twists, and turns in assembly comparisons
  ## cadd is just an integer offset to how to color polygons -- it specifies the offset for the top, bottom, and connecting polygons
  ## step is a tick mark spacer basically 
  ## goncol and bgoncol -- can tell it what to color the polygons and borders (overriding the cadd approach)
  ## tigcol and btigcol -- same but for contig/seq polygons
  
  ## If a paf is not given, use this as default
  if(sum(is.na(paf))>0){paf <- read.table("aln.filtered-pysorted.rnd2.scaf.forR", as.is = TRUE, col.names = c("query", "qlen", "qstart", "qend", "strand", "target", "tlen", "tstart", "tend", "nmatch","alnlen","mapq", "group"))}
  
  ## If falctigs not given, use all from longest to shortest
  if(sum(is.na(falctigs))>0){falctigs <- unique(paf[order(paf$qlen, decreasing = TRUE),1:2])}
  colnames(falctigs) <- c("chr", "len")
  falctigs$cumsum <- cumsum(falctigs$len)
  if(length(falctigs$cumsum) == 1){falctigs$starts <- 0}
  else {falctigs$starts <- cumsum(c(0,falctigs$len[1:(length(falctigs$len)-1)]))}
  
  ## If falctigs not given, use all from longest to shortest
  if(sum(is.na(canutigs))>0){canutigs <- unique(paf[order(paf$tlen, decreasing = TRUE),6:7])}
  colnames(canutigs) <- c("chr", "len")
  canutigs$cumsum <- cumsum(canutigs$len)
  if(length(canutigs$cumsum)==1){canutigs$starts <- 0}
  else{canutigs$starts <- cumsum(c(0,canutigs$len[1:(length(canutigs$len)-1)]))}
  
  ## If revcomp lists are given, revcomp the PAF entries -- will return input PAF if both rc.lists are NA
  paf <- revcomppaf(paf, falctigs, canutigs, rc.falclist, rc.canulist)
  
  ## If gaps and revcomp lists are given, revcomp the BED entries -- will return input BED is rclist is NA
  falcgaps <- revcompbed(falcgaps, falctigs, rc.falclist)
  canugaps <- revcompbed(canugaps, canutigs, rc.canulist)
  
  ## If revcomp lists given, change names to revcomp_name -- THIS STEP NEEDS TO BE AFTER ABOVE 2
  falctigs <- revcomptignames(falctigs, rc.falclist)
  canutigs <- revcomptignames(canutigs, rc.canulist)
  
  ## Define plotting parameters
  xlim <- c(0,max(falctigs$cumsum, canutigs$cumsum))
  
  ## Begin plotting
  plot(seq(0,xlim[2], length.out = 16), rep(0, 16), xlim=xlim, ylim=c(0,1), type="n", yaxt="n", xaxt="n", xlab="Pos (Mb)", ylab="")
  axis(side=1, seq(0,ceiling(xlim[2]/1e8)*1e8, step), seq(0,ceiling(xlim[2]/1e8)*1e8, step)/1e6)
  axis(side=2, c(0,1), c("canu","falcon"), las=1)
  
  ## Add Falcon contig/scaffold boxes
  start <- 1
  for (i in 1:dim(falctigs)[1]){
    end <- falctigs$cumsum[i]
    if(sum(is.na(tigcol))>0){ptigcol <- i+cadd[1]}else{ptigcol<-tigcol}
    polygon(x = c(start, start, end, end, start), y = c(0.95,1.05,1.05,0.95,0.95), col = ptigcol)
    axis(side=3, tick=FALSE, at = mean(c(falctigs$cumsum[i],falctigs$starts[i])), labels = falctigs$chr[i], las=2, cex.axis=0.6)
    start <- end+1
  }
  ## Add Canu contig/scaffold boxes
  start <- 1
  for (i in 1:dim(canutigs)[1]){
    end <- canutigs$cumsum[i]
    if(sum(is.na(tigcol))>0){ptigcol <- i+cadd[2]}else{ptigcol<-tigcol}
    polygon(x = c(start, start, end, end, start), y = c(-0.05,0.05,0.05,-0.05,-0.05), col = ptigcol)
    axis(side=1, tick=FALSE, at = mean(c(canutigs$cumsum[i],canutigs$starts[i])), labels = canutigs$chr[i], las=2, cex.axis=0.6)
    start <- end+1
  }
  ## Draw polygon connections between falcon and canu assemblies
  for (i in 1:dim(paf)[1]){
    faltig <- paf$query[i]
    cantig <- paf$target[i]
    fstart <- falctigs$starts[falctigs$chr == faltig]
    cstart <- canutigs$starts[canutigs$chr == cantig]
    if(sum(is.na(goncol))>0){pgoncol <- i+cadd[3]}else{pgoncol<-goncol}
    polygon(x = c(paf$qstart[i]+fstart, paf$tstart[i]+cstart, paf$tend[i]+cstart, paf$qend[i]+fstart, paf$qstart[i]+fstart), y = c(1,0,0,1,1), col = pgoncol, border=bgoncol)
  }  
  ## Add Scaffold Gap info if provided - THERE CAN BE NO NA VALUES IN GAPS OBJECT
  ## Falcon gaps 
  if(sum(is.na(falcgaps)) == 0){gapaln(tigs=falctigs, gaps=falcgaps, top=1.01, bottom=0.99)}
  if(sum(is.na(canugaps)) == 0){gapaln(tigs=canutigs, gaps=canugaps, top=0.01, bottom=-0.01)}
}

## DEFINE FUNCTION FOR ADDING GAP INFO TO ASSEMBLALIGN PLOT
gapaln <- function(tigs=NA, gaps=NA, top=NA, bottom=NA, col="white", border="white"){
  print("OLD GAPALN")
  ## tigs (e.g. faltigs object), gaps (e.g. falcgaps object), top (e.g.1), bottom (e.g. 0.9)
  ## Require all arguments be provided explicitly to this fxn
  if(sum(is.na(c(tigs, gaps, top, bottom))) > 0){return("Need to provide all arguments.")}
  ## Ensure tigs has 4 colums
  if(!(dim(tigs)[2] == 4)){return("Tigs should be a 4-column dataframe")}
  ## Ensure tigs named correctly -- ASSUMES COLS IN CORRECT ORDER
  colnames(tigs) <- c("chr", "len", "cumsum", "starts")
  
  ## ADD GAPS
  for (i in 1:dim(gaps)[1]){
    tig <- gaps$chr[i]
    if (tig %in% tigs$chr){
      start <- tigs$starts[tigs$chr == tig]
      left <- gaps$start[i]+start
      right <- gaps$end[i]+start
      polygon(x = c(left, left, right, right, left), y = c(bottom, top, top, bottom, bottom), col = col, border = border)
      
    }
  }  
}

## DEFINE REVCOMP FUNCTION FOR PAF (AND/OR GAPS)
revcomppaf <- function(paf=NA, falctigs=NA, canutigs=NA, rc.falclist=NA, rc.canulist=NA){
  ## Require all arguments be provided explicitly to this fxn
  if(sum(is.na(c(paf))) > 0){return("Need to provide PAF.")}
  newpaf <- paf
  if(sum(is.na(rc.falclist))==0){
    for (tig in rc.falclist){
      len <- falctigs$len[falctigs$chr == tig]
      newpaf$qstart[paf$query == tig] <- len - paf$qend[paf$query == tig] + 1 ## adjust pos and make the end = start
      newpaf$qend[paf$query == tig] <- len - paf$qstart[paf$query == tig] + 1 ## adjust pos and make the start = end
      newpaf$query[paf$query == tig] <- paste0('rc_',tig)
    }
  }
  if(sum(is.na(rc.canulist))==0){
    for (tig in rc.canulist){
      len <- canutigs$len[canutigs$chr == tig]
      newpaf$tstart[paf$target == tig] <- len - paf$tend[paf$target == tig] + 1 ## adjust pos and make the end = start
      newpaf$tend[paf$target == tig] <- len - paf$tstart[paf$target == tig] + 1 ## adjust pos and make the start = end
      newpaf$target[paf$target == tig] <- paste0('rc_',tig)
    }
  }
  #RETURN
  newpaf
}

revcompbed <- function(bed=NA, tigs=NA, rc.list=NA){
  ## Require all arguments be provided explicitly to this fxn
  if(sum(is.na(c(bed))) > 0){return("Need to provide BED.")}
  newbed <- bed
  if(sum(is.na(rc.list))==0){
    for (tig in rc.list){
      len <- tigs$len[tigs$chr == tig]
      newbed$start[bed$chr == tig] <- len - bed$end[bed$chr == tig] ## adjust pos and make the end = start ## no +1 b/c BED 0-based
      newbed$end[bed$chr == tig] <- len - bed$start[bed$chr == tig] ## adjust pos and make the start = end ## no +1 b/c BED 0-based
      newbed$chr[bed$chr == tig] <- paste0('rc_',tig)
    }
  }
  #RETURN
  newbed
}


revcomptignames <- function(tigs, rc.list){
  if(sum(is.na(rc.list))==0){
    new <- tigs$chr[tigs$chr %in% rc.list]
    for (i in 1:length(new)){new[i] <- paste0('rc_',new[i])}
    tigs$chr[tigs$chr %in% rc.list] <- new 
  }
  tigs
}



####################################################################################
read_bedpe_sushi <- function(bedpe){
  ##naming as required for Sushi.R
  read.table(bedpe, as.is=TRUE, col.names = c("chromosome1", "start1", "stop1", "chromosome2", "start2", "stop2","name", "score", "strand1", "strand2", "sv"))
}

read_bedpe <- function(fname, as.is = TRUE, col.names=c("chr1", "start1","end1","chr2","start2","end2", "read", "mapq","X1","X2","sv"), colClasses=c(rep(c("character","numeric","numeric"), 2),"character","numeric", rep("character",3)), ...){
  read.table(fname, as.is = as.is, col.names = col.names, colClasses = colClasses, ...)
}

read_bed <- function(bed, chrvarname="chromosome", endvarname='stop'){
  ## naming as required for Sushi.R
  colnames <- c(chrvarname, "start", endvarname,"name", "score", "strand")
  tmp <- read.table(bed, as.is = TRUE, nrows = 2)
  nc <- ncol(tmp)
  read.table(bed, as.is=TRUE, col.names = colnames[1:nc])
}

read_bed4 <- function(bed, colnames=c("chr","start","end","label"), colClasses=c("character","numeric","numeric","character")){
  ## naming as required for Sushi.R
  read.table(bed, as.is=TRUE, col.names = colnames, colClasses = colClasses)
}


read_bed_pe_as_paf <- function(fname, as.is = TRUE, col.names=c("query", "qstart","qend","target","tstart","tend", "read", "mapq","X1","X2","sv"), colClasses=c(rep(c("character","numeric","numeric"), 2),"character","numeric", rep("character",3)), ...){
  read.table(fname, as.is = as.is, col.names = col.names, colClasses = colClasses, ...)
}

read_bed_pe_as_paf_alt <- function(fname, as.is = TRUE, col.names=c("query", "qstart","tstart","target","qend","tend", "read", "mapq","X1","X2","sv"), colClasses=c(rep(c("character","numeric","numeric"), 2),"character","numeric", rep("character",3)), ...){
  read.table(fname, as.is = as.is, col.names = col.names, colClasses = colClasses, ...)
}

read_genome_file <- function(fname, as.is = TRUE, col.names=c("chr","len"), colClasses=c("character","numeric"), ...){
  read.table(fname, as.is = as.is, col.names = col.names, colClasses = colClasses, ...)
}












####OTHER INDEVO OR ABANDONED#################################################
###############################################################################
## Try to find solution that makes everything positive
adjust <- function(solpaf, falc=NA, canu=NA){
  if(!(is.na(falc))){
    solpaf$strand[solpaf$query == falc] = -1*solpaf$strand[solpaf$query == falc]
  }
  if(!(is.na(canu))){
    solpaf$strand[solpaf$target == canu] = -1*solpaf$strand[solpaf$target == canu]
  }
  print(100*sum(solpaf$strand == 1)/length(solpaf$strand))
  solpaf
}
showgrp <- function(solpaf, group='group_1'){
  solpaf[solpaf$group == group,]
}
solvedgrp <- function(solpaf, group='group_1'){
  100*sum(solpaf[solpaf$group == group,]$strand == 1)/length(solpaf[solpaf$group == group,]$strand)
}
showrows <- function(solpaf, falc=NA, canu=NA){
  if(!(is.na(falc))){
    solpaf[solpaf$query == falc, ]
  }else if(!(is.na(canu))){
    solpaf[solpaf$target == canu, ]
  }
}
getsolpaf <- function(paf){
  solpaf <- paf[, c(1,6,5,13)]
  solpaf$strand[solpaf$strand == '-'] = -1
  solpaf$strand[solpaf$strand == '+'] = 1
  class(solpaf$strand) <-"numeric"
  print(100*sum(solpaf$strand == 1)/length(solpaf$strand))
  solpaf
}




#########################################
### pairWiseLociPlotting

get_target_chr_len <- function(genome, target=NA){
  if(!(is.na(target))){d <-genome[genome$chr == target,]}
  d
}

get_target_chr <- function(d, target=NA){
  if(!(is.na(target))){
    d.out <- d[d$chr == target,]
  } else {
    d.out <- d
  }
  d.out
}

get_target_region <- function(d, target=NA, start=NA, end=NA){
  if(!(is.na(target))){
    d <- get_target_chr(d, target)
  } 
  if(!(is.na(start))){
    d <- d[d$start >= start,]
  }
  if(!(is.na(end))){
    d <- d[d$end <= end,]
  }
  d
}

get_midpoints <- function(d){
  rowMeans(data.frame(start=d$start, end=d$end))
}
get_target_insert_sizes <-function(d,shortest, longest){
  d[abs(d$end-d$start) >= shortest & abs(d$end-d$start) <= longest,]
}

prep_for_full_matrix <- function(filtpairs){
  revpairs <- filtpairs
  revpairs$start <- filtpairs$end
  revpairs$end <- filtpairs$start
  xy <- rbind(filtpairs, revpairs)
  filtpairs <- xy[order(xy$start, xy$end),]
  filtpairs
}


get_xyz_vars <- function(filtpairs, bx, by, gridsize, maxpair, maxz, xlim, ylim, znorm=1){
  print(maxpair)
  ## GET SMOOTHED MATRIX
  ans <- bkde2D(x = filtpairs[,2:3], bandwidth = c(bx,by), gridsize = c(gridsize, gridsize), range.x = list(x1=c(0,maxpair),x2=c(0,maxpair)))
  print(sum(ans$fhat))
  z<-znorm*ans$fhat
  #### z<- z/sum(z) ## I wonder if just normalizing to self makes more sense than total mappable reads....
  x <- ans$x1
  y <- ans$x2
  
  ## PRE-PROCESS MAX VALUE PRIOR TO SUBSETTING MATRIX -- THIS VALUE CAN BE PROVIDED BY USER
  #calculate the following before possible subsetting so that coloring stays same closer up
  if(is.na(maxz)){
    print("HELLO")
    maxz <- max(z)
  } 
  
  ## TELL USER WHAT MAX MATRIX VALUE WAS SO THEY CAN USE IT IN OTHER PLOTS WITH MAXZ VARIABLE
  print(c("maxz",maxz))
  
  
  ## MATRIX SUBSETTING - X AND Z VECTORS
  if(!(is.na(xlim[1]))){
    if(length(xlim) != 2){print("Xlim must have 2 floats: min and max.")}
    rows <- x >= xlim[1] & x <= xlim[2]
    x <- x[rows]
    z <- z[rows,]
    ## Even though you'd think the levelplot visualization x-axis corresponds to the columns - I think the matrix is tilted 90* backwards in vis -- so need to swap logic to get right vis result
  }
  
  ## MATRIX SUBSETTING - Y AND Z VECTORS
  if(!(is.na(ylim[1]))){
    if(length(ylim) != 2){print("Ylim must have 2 floats: min and max.")}
    columns <- y >= ylim[1] & y <= ylim[2]
    y <- y[columns]
    z <- z[,columns]
  }
  
  return(list(x=x, y=y, z=z, maxz=maxz, bx=bx, xmax=maxpair, gridsize=gridsize))
}

get_color_code <- function(matcols, matcolreps){
  if(sum(is.na(matcolreps))>0){matcolreps <- rep(1,length(matcols))}
  colorcode <- c()
  for (i in 1:length(matcols)){
    colorcode <- c(colorcode, rep(matcols[i], matcolreps[i]))
  }
  colorcode
}

get_lave_colorcode <- function(x=1){
  if(x==1){colorcode <- c("dark blue", "white", "red")}
  else if (x==2){colorcode <- c("dark blue", "cyan", "white", "red", "yellow")}
  else if (x==3){colorcode <- c("white","black","red","orange","yellow")}
  else if (x==4){colorcode <- c("dark blue", "white","red","yellow")}
  else if (x==5){colorcode <- c("black","dark blue","white","red","yellow")}
  else if (x==6){colorcode <- c("black", "white")}
  else if (x==7){colorcode <- c("black", "white", "red")}
  else if (x==8){colorcode <- c("white","black")}
  else if (x==9){colorcode <- c("white","black","red")}
  else if (x==10){colorcode <- c("white","black","red","yellow")}
  else if (x==11){colorcode <- c("white","black","red","orange", "yellow")}
  else if (x==12){colorcode <- c("violet","blue","dark blue","dark green","green","orange","red")}
  else if (x==13){colorcode <- c("black","dark blue","cyan","white","red","orange", "yellow")}
  else if (x==14){colorcode <- c("white", "cyan","dark blue","black","red","orange", "yellow")}
  else if (x==15){colorcode <- c("white","dark blue","black","red","orange", "yellow")}
  else if (x==16){colorcode <- c("cyan","blue","dark blue","black","red","orange", "yellow")}
  else if (x==17){colorcode <- c("cyan","blue","dark blue","grey","red","orange", "yellow")}
  else if (x==18){colorcode <- c("white","dark blue","dark grey","red","orange", "yellow")}
  else if (x==19){colorcode <- c("cyan","blue","dark blue","white","red","orange", "yellow")}
  else if (x==20){colorcode <- c("white", "cyan","blue","black","red","orange", "yellow")}
  else if (x==21){colorcode <- c("white", "cyan","blue", "dark blue","black", "dark red","red","orange", "yellow")}
  else if (x==22){colorcode <- c("blue","cyan","white","black","yellow","orange","red")}
  else if(x==23){colorcode <- c("white","red","black")}
  else if(x==24){colorcode <- c("white","blue","black")}
  else if(x==25){colorcode <- c("white","green","black")}
  else if(x==26){colorcode <- c("white","cyan","black")}
  else if(x==27){colorcode <- c("white","orange","black")}
  return(colorcode)
}

get_levelplot <- function(z, x, y, colorkey, coltrio, seqmaxz, ust, uen, seglwd, segcol){
  lattice.options(axis.padding=list(factor=0.05))
  levelplot(z, row.values = x, column.values = y, 
            colorkey=colorkey, 
            xlab="", ylab="", scales=list(x=list(rot=90), 
                                          pretty=TRUE, ylab="", xlab="", tck = c(1,0)),
            col.regions=colorRampPalette(coltrio), at=seqmaxz,
            # panel = function(...){panel.levelplot(...); panel.segments(x0=c(ust, ust, uen, uen), x1=c(ust, uen, uen, ust), y0=c(ust, ust, uen, uen), y1=c(uen, ust, ust, uen), lwd=seglwd, col=segcol)})
            ## Mar 15, 2018 - Trying suggestion from here:http://r.789695.n4.nabble.com/Removing-cell-borders-from-svg-or-eps-in-levelplot-td4684807.html
            panel = function(...){panel.levelplot.raster(...); panel.segments(x0=c(ust, ust, uen, uen), x1=c(ust, uen, uen, ust), y0=c(ust, ust, uen, uen), y1=c(uen, ust, ust, uen), lwd=seglwd, col=segcol)})
            ## IN FUTURE: perhaps can add multiple BED files by (i) concatenating all to ust and uen appropriately, and (ii) providing a segcol list matching the number of feature in each file
  
}


prep_paired_loci_values <- function(data, target, genome, shortest, longest, makefullmatrix, bx, by, gridsize, maxz, xlim, ylim, znorm){
  ## PREPARE DATA SET
  pairs <- get_target_chr(data, target)
  filtpairs <- get_target_insert_sizes(pairs, shortest,longest)
  if(makefullmatrix){filtpairs <- prep_for_full_matrix(filtpairs)}
  
  ## GET END COORD
  if(sum(is.na(genome))==0){
    if(!is.na(target)){maxpair<-genome$size[genome$chr == target]}
    else{maxpair <- genome$size[genome$chr == data$chr[1]]}
  }else{
    maxpair <- max(filtpairs[,2:3], na.rm = TRUE) ## 2018Aug29-note: this seems to determine chrom end via right-most coord of read mappings -- should allow a .genome file to be passed for actual chrom lens
  }
  
  ## GET AND RETURN SMOOTHED XYZ MATRIX VARIABLES
  # xyz <- get_xyz_vars(filtpairs, bx, by, gridsize, maxpair, maxz, xlim, ylim, znorm)
  get_xyz_vars(filtpairs, bx, by, gridsize, maxpair, maxz, xlim, ylim, znorm)
}


pairedLociLevelPlot <- function(data, UR, target=NA, genome=NA, shortest=0, longest=1e10, bx=1e4, by=1e4, gridsize=200, makefullmatrix=FALSE, seglwd=2, segcol="black", matcols=c("dark blue","white","red"), matcolreps=NA, xlim=c(NA), ylim=c(NA), maxz=NA, znorm=1, cnnormdata=c(NA), finalnormbymaxz=FALSE, returnlevelplot=TRUE, plotcn=TRUE){
  
  ## BED start and end
  # print("A")
  target.urs <- get_target_chr(UR, target)
  ust <- target.urs$start
  uen <- target.urs$end
  
  ## PREPARE PAIRED LOCI VALUES
  # print("B")
  xyz <- prep_paired_loci_values(data, target, genome, shortest, longest, makefullmatrix, bx, by, gridsize, maxz, xlim, ylim, znorm)
  
  ## POTENTIALLY NORMALIZE BY COPY NUMBER OR OTHER
  ## For now it also uses znorm -- but I have considered not -- all in all it is a scalar so will not hurt
  # print("C")
  if(sum(is.na(cnnormdata)==0)){
    ## PREPARE PAIRED LOCI VALUES FOR COVERAGE DATASET
    # print("D")
    cn.xyz <- prep_paired_loci_values(cnnormdata, target=target, genome=genome, shortest=0, longest=1e12, makefullmatrix=makefullmatrix, bx=bx, by=by, gridsize=gridsize, maxz=NA, xlim=xlim, ylim=ylim, znorm=znorm)
    

    ## GET COVERAGE VALUES FOR EACH BIN BY TAKING THE SUM OF SIGNAL FOR THAT BIN WITH ALL OTHER BINS ()
    cncov <- colSums(cn.xyz$z, na.rm = TRUE)
    n.cncov <- length(cncov)
    
    ## NORMALIZE COVERAGE BY THE MEDIAN COVERAGE VALUE (ASSUMED TO BE CLOSE TO CN=1) [[NOTE THIS NEGATES ANY GIVEN ZNORM CONSTANT]]
    cncov.med <- median(cncov)

    cncov.mednorm <- cncov/cncov.med
    

    ## SHOW USER A PLOT OF THE COVERAGE ALONG WITH PROVIDED BED LOCATIONS
    if(is.na(target)){covxlab<-"Position (bp)"}else{covxlab<-paste0(target," Position (bp)")}
    if(plotcn){
      print("plotcn")
      plot(cn.xyz$x, cncov.mednorm, type='l', xlab=covxlab, ylab="RCN", main="Coverage analysis given random sample\n and binning/smoothing parameters.", las=1)
      covy <- rep(1, dim(target.urs)[1])
      segments(x0 = target.urs$start, x1=target.urs$end, y0=covy, y1=covy, lwd=2, col="blue")
    }
    # print("E")
    ## NORMALIZE EACH PAIRWISE LOCI INTERACTION TO THE MINIMUM COVERAGE VALUE BETWEEN THE TWO BINS (AS THE MIN SETS THE LIMIT TO HOW MANY TIMES THAT INTERACTION COULD BE DETECTED)
    ## xyz$z <- xyz$z/cncov.mednorm  ## For now does not take min val for 2 loci -- It is just normalizing matrix in on orientation
    ## # FIRST CREATE CN MATRIX
    cnmat <- matrix(nrow = n.cncov, ncol = n.cncov)
    for(i in 1:n.cncov){
      for(j in i:n.cncov){
        mv <- min(cncov.mednorm[i], cncov.mednorm[j])
        cnmat[i,j] <- mv
        cnmat[j,i] <- mv
      }
    }

    ## # DIVIDE CURRENT PAIRWISE CONTACT MATRIX BY CN MATRIX ELEMENT-WISE
    print(c("Number CN=0", sum(cnmat == 0)))
    xyz$z <- xyz$z/cnmat

    #print(sum(xyz$z == 0))
    #print(sum(is.na(xyz$z)))
    
    ## REPLACE MAXZ ONLY IF MAXZ WAS NOT SPECIFIED BY USER
    if(is.na(maxz)){
      xyz$maxz <- max(xyz$z, na.rm = TRUE)
      print(c("maxz.update", xyz$maxz))
    }
  }
  
  ## OPTIONALLY DIVIDE MATRIX BY MAXZ TO PUT BETWEEN 0 AND 1 (NOT NEC BTWN 0-1 WHEN USER PROVIDES MAXZ)
  if(finalnormbymaxz){
    print(c("maxz.used-for-final-norm", xyz$maxz))
    xyz$z <- xyz$z/xyz$maxz
    xyz$maxz <- 1
  }
  
  ## DESIGNING THE COLOR VECTOR TO FEED LEVELPLOT
  coltrio <- get_color_code(matcols, matcolreps)
  
  ## PREPARE SOME INPUTS FOR LEVEL PLOT
  ## aug2018-tck controls tick mark length for colorkey -- , tick.number=3 seemed to do nothing -- replaced with labels list.
  seqmaxz <- seq(0,xyz$maxz,xyz$maxz/1000)
  seqmaxz2 <- seq(0,xyz$maxz,xyz$maxz/3)
  colorkey <- list(space="right", col=colorRampPalette(coltrio), at=seqmaxz, raster=TRUE, tck=1, labels=list(labels=seqmaxz2, at=seqmaxz2))
  
  ## LEVEL PLOT
  # return(list(z=z, x=x, y=y, colorkey=colorkey))
  #print(xyz$x[2:length(xyz$x)]-xyz$x[1:(length(xyz$x)-1)])
  if(returnlevelplot){
    get_levelplot(xyz$z, xyz$x, xyz$y, colorkey, coltrio, seqmaxz, ust, uen, seglwd, segcol) 
  } else {
    return(list(z=xyz$z, x=xyz$x, y=xyz$y, colorkey=colorkey, coltrio=coltrio, maxz=xyz$maxz, seqmaxz=seqmaxz, ust=ust, uen=uen, seglwd=seglwd, segcol=segcol, bx=xyz$bx, xmax=xyz$xmax, gridsize=xyz$gridsize))
  }
}

pairedLociLevelPlot2 <- function(data1, data2, UR, target=NA, genome=NA, shortest=0, longest=1e10, bx=1e4, by=1e4, gridsize=200, makefullmatrix=FALSE, seglwd=2, segcol="black", matcols=c("dark blue","white","red"), matcolreps=NA, xlim=c(NA), ylim=c(NA), maxz=NA, znorm1=1, znorm2=1, cnnormdata1=c(NA), cnnormdata2=c(NA), finalnormbymaxz=FALSE, setfinalto1=TRUE, returnlevelplot=TRUE, plotcn=TRUE){
  ## No need to normalize by some constant maxz normed by read and copy number -- norming by constant_maxz for both matrices is just a constant that is weeded out when matrices compared
  ## Normalizing each independently to their own maxz is not necessarily appropriate either
  # print(1)
  info1 <- pairedLociLevelPlot(data1, UR, target, genome, shortest, longest, bx, by, gridsize, makefullmatrix, seglwd, segcol, matcols, matcolreps, xlim, ylim, maxz, znorm1, cnnormdata1, finalnormbymaxz, returnlevelplot=FALSE, plotcn=plotcn)
  # print(2)
  info2 <- pairedLociLevelPlot(data2, UR, target, genome, shortest, longest, bx, by, gridsize, makefullmatrix, seglwd, segcol, matcols, matcolreps, xlim, ylim, maxz, znorm2, cnnormdata2, finalnormbymaxz, returnlevelplot=FALSE, plotcn=plotcn)
  # print(3)
  info.fe <- info1
  # print(4)
  l1 <- sum(info1$x==info2$x)/length(info1$x)
  l2 <- sum(info1$y==info2$y)/length(info1$y)
  print(c(l1,l2))
  ##info.fe$z <- info1$z/info2$z ## PRODUCES NANS THAT ARE AUTOMATICALLY WHITE -- WOULD NEED PSEUDOCOUNTS SOMEWHERE ALONG THE WAY
  ## info.fe$seqmaxz <- seq(0,info.fe$maxz,info.fe$maxz/1000)
  info.fe$z <- info1$z-info2$z
  
  info.fe$maxz <- max(info.fe$z)
  info.fe$minz <- min(info.fe$z)
  mostextreme <- max(abs(c(info.fe$minz, info.fe$maxz)))
  if(setfinalto1){
    info.fe$z <- info.fe$z/mostextreme
    info.fe$maxz <- max(info.fe$z)
    info.fe$minz <- min(info.fe$z) 
    mostextreme <- max(abs(c(info.fe$minz, info.fe$maxz)))
  }
  seqrange <- c(-1*mostextreme, mostextreme)
  # info.fe$seqmaxz <- seq(info.fe$minz,info.fe$maxz,(info.fe$maxz-info.fe$minz)/1000)
  # seqmaxz2 <- seq(info.fe$minz,info.fe$maxz,(info.fe$maxz-info.fe$minz)/3)
  info.fe$seqmaxz <- seq(seqrange[1],seqrange[2], (seqrange[2]-seqrange[1])/1000)
  seqmaxz2 <- seq(seqrange[1],seqrange[2], (seqrange[2]-seqrange[1])/3)
  info.fe$colorkey <- list(space="right", col=colorRampPalette(info.fe$coltrio), at=info.fe$seqmaxz, raster=TRUE, tck=1, labels=list(labels=seqmaxz2, at=seqmaxz2))
  ## LEVEL PLOT
  # return(list(z=z, x=x, y=y, colorkey=colorkey))
  if(returnlevelplot){
    get_levelplot(info.fe$z, info.fe$x, info.fe$y, info.fe$colorkey, info.fe$coltrio, info.fe$seqmaxz, info.fe$ust, info.fe$uen, info.fe$seglwd, info.fe$segcol) 
  } else {
    return(info.fe)
  }
}


splitreadlevelplot <- function(data, UR, target=NA, genome=NA, shortest=0, longest=1e10, bx=1e4, by=1e4, gridsize=200, makefullmatrix=FALSE, seglwd=2, segcol="black", matcols=c("dark blue","white","red"), matcolreps=c(1,3,1), xlim=c(NA), ylim=c(NA), maxz=NA, znorm=1, cnnormdata=c(NA), finalnormbymaxz=FALSE, returnlevelplot=TRUE){
  ## splitreadlevelplot now just wraps over pairedLociLevelPlot
  ## pairedLociLevelPlot is a more suitable name since it processes both split-read and PE data
  pairedLociLevelPlot(data, UR, target, genome, shortest, longest, bx, by, gridsize, makefullmatrix, seglwd, segcol, matcols, matcolreps, xlim, ylim, maxz, znorm, cnnormdata, finalnormbymaxz, returnlevelplot)
}

get_stacked_levelplot <- function(s, x, y, colorkey, coltrio, seqmaxz, ust, uen, seglwd, segcol){
  lattice.options(axis.padding=list(factor=0.05))
  levelplot(s, layout=c(1,2), 
            row.values = x, column.values = y, 
            colorkey=colorkey, 
            xlab="", ylab="", scales=list(x=list(rot=90), 
                                          pretty=TRUE, ylab="", xlab="", tck = c(1,0)),
            col.regions=colorRampPalette(coltrio), at=seqmaxz,
            # panel = function(...){panel.levelplot(...); panel.segments(x0=c(ust, ust, uen, uen), x1=c(ust, uen, uen, ust), y0=c(ust, ust, uen, uen), y1=c(uen, ust, ust, uen), lwd=seglwd, col=segcol)})
            ## Mar 15, 2018 - Trying suggestion from here:http://r.789695.n4.nabble.com/Removing-cell-borders-from-svg-or-eps-in-levelplot-td4684807.html
            panel = function(...){panel.levelplot.raster(...); panel.segments(x0=c(ust, ust, uen, uen), x1=c(ust, uen, uen, ust), y0=c(ust, ust, uen, uen), y1=c(uen, ust, ust, uen), lwd=seglwd, col=segcol)})
  
}

pairedLociStackedLevelPlot <- function(data1, data2, UR, target=NA, shortest=0, longest=1e10, bx=1e4, by=1e4, gridsize=200, makefullmatrix=FALSE, seglwd=2, segcol="black", matcols=c("dark blue","white","red"), matcolreps=c(1,3,1), xlim=c(NA), ylim=c(NA), maxz=NA, znorm=1){
  ## copy/poasted pairedLociLevelPlot on Aug29,2018 to begin dev
  ## Want to be able to pass 2 (or more) datasets and have them plotted together
  ## PREPARE PAIRED LOCI VALUES
  xyz1 <- prep_paired_loci_values(data1, target, shortest, longest, makefullmatrix, bx, by, gridsize, maxz, xlim, ylim, znorm)
  xyz2 <- prep_paired_loci_values(data2, target, shortest, longest, makefullmatrix, bx, by, gridsize, maxz, xlim, ylim, znorm)
  
  ## DESIGNING THE COLOR VECTOR TO FEED LEVELPLOT
  coltrio <- get_color_code(matcols, matcolreps)
  
  ## BED start and end
  target.urs <- get_target_chr(UR, target)
  ust <- target.urs$start
  uen <- target.urs$end
  
  ## PREPARE SOME INPUTS FOR LEVEL PLOT
  ## aug2018-tck controls tick mark length for colorkey -- , tick.number=3 seemed to do nothing -- replaced with labels list.
  maxz <- max(xyz1$maxz, xyz2$maxz)
  seqmaxz <- seq(0,maxz,maxz/1000)
  seqmaxz2 <- seq(0,maxz,maxz/3)
  colorkey <- list(space="right", col=colorRampPalette(coltrio), at=seqmaxz, raster=TRUE, tck=1, labels=list(labels=seqmaxz2, at=seqmaxz2))
  
  ## LEVEL PLOT
  # return(list(z=z, x=x, y=y, colorkey=colorkey))
  s <- stack(raster(xyz1$z), raster(xyz2$z))
  print(sum( xyz1$x == xyz2$x)/length(xyz1$x))
  print(sum( xyz1$y == xyz2$y)/length(xyz1$y))
  get_stacked_levelplot(s,  xyz1$x, xyz1$y, colorkey, coltrio, seqmaxz, ust, uen, seglwd, segcol)
}

get_row_plot <- function(data, UR, target=NA, genome=NA, shortest=0, longest=1e10, bx=1e4, by=1e4, gridsize=200, makefullmatrix=FALSE, seglwd=2, segcol="black", matcols=c("dark blue","white","red"), matcolreps=NA, xylim=c(NA), maxz=NA, znorm=1, cnnormdata=c(NA), finalnormbymaxz=FALSE, returnlevelplot=TRUE, plotcn=TRUE, rowslice=c(0,1e100), fxn=identity, presentation="sum", marginalize=FALSE, xlab="Position", ylab="Relative Interactions", xlim=c(NA)){
  ## Take a given min and max locus index (given by xylim) and plot as regular line plot
  ## When it is 1 row -- as is
  ## When it is more than 1 row, either as Sum or each is its own line
  ## presentation can be "sum" or "all"
  
  
  # if(sum(is.na(xylim))>0){xylim <- c(0,1e100)}
  xyz <- pairedLociLevelPlot(data, UR, target, genome, shortest, longest, bx, by, gridsize, makefullmatrix, seglwd, segcol, matcols, matcolreps, xylim, xylim, maxz, znorm, cnnormdata, finalnormbymaxz, returnlevelplot=FALSE, plotcn=plotcn)

  xbool <- xyz$x >= rowslice[1] & xyz$x <= rowslice[2]
  
  if(sum(is.na(xlim))> 0){xlim=c(0,xyz$xmax)}
  if(presentation == "return"){
    y <- colSums(xyz$z[xbool,])
    if(marginalize){y <- y/sum(y)}
    return(list(x=xyz$x,y=fxn(y),type="l", las=1, xlab=xlab, ylab=ylab))
  } else if(presentation == "sum"){
    y <- colSums(xyz$z[xbool,])
    if(marginalize){y <- y/sum(y)}
    plot(xyz$x,fxn(y),type="l", las=1, xlab=xlab, ylab=ylab, xlim=xlim)
  } else if (presentation == "all"){
    ylim <- range(xyz$z[xbool,])
    plot(xyz$x, xyz$x, ylim=ylim, type="n", las=1, xlab=xlab, ylab=ylab, xlim=xlim)
    for (i in (1:length(xyz$x))[xbool]){
      y <- xyz$z[i,]
      if(marginalize){y <- y/sum(y)}
      lines(xyz$x, fxn(y), col=i)
    }
  }
}  

get_row_plot_2_samples <- function(data1, data2, UR, target=NA, genome=NA, shortest=0, longest=1e10, bx=1e4, by=1e4, gridsize=200, makefullmatrix=FALSE, seglwd=2, segcol="black", matcols=c("dark blue","white","red"), matcolreps=NA, xylim=c(NA), maxz=NA, znorm1=1, znorm2=1, cnnormdata1=c(NA), cnnormdata2=c(NA), finalnormbymaxz=FALSE, returnlevelplot=TRUE, plotcn=TRUE, rowslice=c(0,1e100), fxn=identity, presentation="subtract", marginalize=FALSE, xlab="Position", ylab="Relative Interactions", col1="black", col2="grey", xlim=c(NA)){
  ## Take a given min and max locus index (given by xylim) and plot as regular line plot
  ## When it is 1 row -- as is
  ## When it is more than 1 row, either as Sum or each is its own line
  ## presentation can be "sum" or "all"
  
  
  # if(sum(is.na(xylim))>0){xylim <- c(0,1e100)}
  xyz1 <- pairedLociLevelPlot(data1, UR, target, genome, shortest, longest, bx, by, gridsize, makefullmatrix, seglwd, segcol, matcols, matcolreps, xylim, xylim, maxz, znorm1, cnnormdata1, finalnormbymaxz, returnlevelplot=FALSE, plotcn=plotcn)
  xyz2 <- pairedLociLevelPlot(data2, UR, target, genome, shortest, longest, bx, by, gridsize, makefullmatrix, seglwd, segcol, matcols, matcolreps, xylim, xylim, maxz, znorm2, cnnormdata2, finalnormbymaxz, returnlevelplot=FALSE, plotcn=plotcn)
  
  xbool1 <- xyz1$x >= rowslice[1] & xyz1$x <= rowslice[2]
  xbool2 <- xyz2$x >= rowslice[1] & xyz2$x <= rowslice[2]
  print(c("Subsetting agrees:",sum(xbool1==xbool2)/length(xbool1)))
  
  y1 <- colSums(xyz1$z[xbool1,])
  y2 <- colSums(xyz2$z[xbool2,])
  if(marginalize){y1 <- y1/sum(y1); y2 <- y2/sum(y2)}
  print(xlim)
  if(sum(is.na(xlim))> 0){xlim=c(0,xyz1$xmax)}
  print(xlim)
  if(presentation == "return"){
    return(list(x1=xyz1$x,x2=xyz2$x, y1=fxn(y1), y2=fxn(y2), type="l", las=1, xlab=xlab, ylab=ylab))
  } else if(presentation == "both"){
    ylim <- c(0,max(y1,y2))
    plot(xyz1$x,fxn(y1),type="n", las=1, xlab=xlab, ylab=ylab, ylim=ylim, xlim=xlim)
    lines(xyz1$x,fxn(y1), col=col1)
    lines(xyz2$x,fxn(y2), col=col2)
  } else if (presentation == "subtract"){
    y <- y1-y2
    plot(xyz1$x, y, type="l", las=1, xlab=xlab, ylab=ylab, xlim=xlim)
  } else {return(list(x1=xyz1$x,x2=xyz2$x, y1=fxn(y1), y2=fxn(y2), type="l", las=1, xlab=xlab, ylab=ylab))}
}  

addpoly <- function(d,t=1000,b=-100, bcol=rgb(0,0,1,0.2), col=rgb(0,0,1,0.2)){
  for(i in 1:length(d$start)){
    l <- d$start[i]
    r <- d$end[i]
    polygon(x = c(l,l,r,r,l), y=c(b,t,t,b,b), border = bcol, col = col)
  }
}

##TODO: add to pairedloci triangle -- also have an addsquares to reg matrix
# add_bed_diamonds <- function(d,t=1000,b=-100, bcol=rgb(0,0,1,0.2), col=rgb(0,0,1,0.2)){
#   for(i in 1:length(d$start)){
#     l <- d$start[i]
#     r <- d$end[i]
#     polygon(x = c(l,l,r,r,l), y=c(b,t,t,b,b), border = bcol, col = col)
#   }
# }

insertplot <- function(insertstats, UR, target){
  target.chr <- insertstats[insertstats$chr == target, ]
  target.urs <- UR[UR$chr == target,]
  
  ## number of PE reads with insert sizes between 10kb-500kb that overlap 5kb bin
  plot(target.chr$start, target.chr$count, cex=0.5)
  segments(x0 = target.urs$start, x1 = target.urs$end, y0 = rep(0, length(target.urs$start)), y1 = rep(0, length(target.urs$start)), lwd=5, col="blue")
  addpoly(target.urs)
}

pairedlociscatters <- function(data, UR, target=NA, shortest=0, longest=1e10, bx=1e4, by=1e4, gridsize=200, heatgrid=1e3, makefullmatrix=FALSE, segcol="blue", segcol2="black", segcol3="blue"){
  ## Formerly known as splitreadplot
  ## This function is basically abandoned for the pairedLociLevel and pairedLociImage series...
  ## TODO - adjust counts by dividing by copy number.... (may also want to adjust to a control sample -- else perhaps GC bias and/or mappability can cause artifacts)
  if(!(is.na(target))){
    pairs <- data[data$chr == target,]
    target.urs <- UR[UR$chr == target,]
  } else {
    pairs <- data
    target.urs <- UR
  }
  
  filtpairs <- pairs[abs(pairs$end-pairs$start) >= shortest & abs(pairs$end-pairs$start) <= longest,]
  
  if(makefullmatrix){
    revpairs <- filtpairs
    revpairs$start <- filtpairs$end
    revpairs$end <- filtpairs$start
    xy <- rbind(filtpairs, revpairs)
    filtpairs <- xy[order(xy$start, xy$end),]
  }
  
  
  par(mfrow=c(3,2))
  plot(filtpairs$start, filtpairs$end)
  segments(x0=target.urs$start, x1=target.urs$end, y0=target.urs$start, y1=target.urs$start, col=segcol, lwd=5)
  segments(x0=target.urs$end, x1=target.urs$end, y0=target.urs$end, y1=target.urs$start, col=segcol, lwd=5)
  
  # heatscatter(filtpairs$start, filtpairs$end, grid=heatgrid)
  heatscatter(filtpairs$start, filtpairs$end)
  for(i in 1:length(target.urs$start)){
    s <- target.urs$start[i]
    e <- target.urs$end[i]
    x <- c(s,s,e,e,s)
    y <- c(0, s, e, 0, 0)
    polygon(x = x, y = y, col = "grey")
  }
  # abline(v=(target.urs$start+target.urs$end)/2, lty=3)
  # abline(h=(target.urs$start+target.urs$end)/2, lty=3)
  
  smoothScatter(filtpairs$start, filtpairs$end)
  segments(x0=target.urs$start, x1=target.urs$end, y0=target.urs$start, y1=target.urs$start, col=segcol2, lwd=5)
  segments(x0=target.urs$end, x1=target.urs$end, y0=target.urs$end, y1=target.urs$start, col=segcol2, lwd=5)
  
  ## LATTICE DOESNT PLAY WELL WITH THESE PLOTS
  # maxpair <- max(filtpairs[,2:3])
  # ans <- bkde2D(x = filtpairs[,2:3], bandwidth = c(bx,by), gridsize = c(gridsize, gridsize), range.x = list(x1=c(0,maxpair),x2=c(0,maxpair)))
  # lattice.options(axis.padding=list(factor=0.5))
  # z<-ans$fhat
  # x <- ans$x1
  # y <- ans$x2
  # maxz <- max(z)
  # seqmaxz <- seq(0,maxz,maxz/1000)
  # coltrio <- c("dark blue", "white", "red")
  # ust <- target.urs$start
  # uen <- target.urs$end
  # levelplot(z, row.values = x, column.values = y, 
  #           colorkey=list(space="right", col=colorRampPalette(coltrio), at=seqmaxz, raster=TRUE, tck=1, tick.number=3), 
  #           xlab="", ylab="", scales=list(x=list(rot=90), 
  #                                         pretty=TRUE, ylab="", xlab="", tck = c(1,0)),
  #           col.regions=colorRampPalette(coltrio), at=seqmaxz,
  #           panel = function(...){panel.levelplot(...); panel.segments(x0=c(ust, ust, uen, uen), x1=c(ust, uen, uen, ust), y0=c(ust, ust, uen, uen), y1=c(uen, ust, ust, uen))})
  
  
  
  ans <- bkde2D(x = filtpairs[,2:3], bandwidth = c(bx,by))
  contour(ans$x1, ans$x2, ans$fhat)
  segments(x0=target.urs$start, x1=target.urs$end, y0=target.urs$start, y1=target.urs$start, col=segcol3, lwd=5)
  segments(x0=target.urs$end, x1=target.urs$end, y0=target.urs$end, y1=target.urs$start, col=segcol3, lwd=5)
  
  # filled.contour(ans$x1, ans$x2, ans$fhat, color.palette = heat.colors) ## cant use inside panel b/c starts new layout
  persp(ans$fhat, theta=235)
  
  par(mfrow=c(1,1))
}





pairedLociTriangle <- function (hicdata, chrom, chromstart, chromend, max_y = 30, zrange = NULL, palette = SushiColors(7), flip = FALSE) {
  ## Derived from Sushi.R's plothic function
  rows = as.numeric(rownames(hicdata))
  cols = as.numeric(colnames(hicdata))
  hicregion = hicdata[which(rows >= chromstart & rows <= chromend), 
                      which(cols >= chromstart & cols <= chromend)]
  nbins = nrow(hicregion)
  stepsize = abs(chromstart - chromend)/(2 * nbins)
  hicm = as.matrix(hicregion)
  if (is.null(zrange) == TRUE) {
    max_z = max(hicm)
    min_z = min(hicm)
  }
  if (is.null(zrange) == FALSE) {
    min_z = zrange[1]
    max_z = zrange[2]
  }
  hicmcol = matrix(maptocolors(hicm, palette, num = 100, range = zrange), 
                   nrow = nrow(hicm))
  if (flip == FALSE) {
    plot(1, 1, xlim = c(chromstart, chromend), ylim = c(0, 
                                                        max_y), type = "n", xaxs = "i", yaxs = "i", bty = "n", 
         xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  }
  if (flip == TRUE) {
    plot(1, 1, xlim = c(chromstart, chromend), ylim = c(-max_y, 
                                                        0), type = "n", xaxs = "i", yaxs = "i", bty = "n", 
         xaxt = "n", yaxt = "n", xlab = "", ylab = "")
  }
  for (rownum in (1:nrow(hicm))) {
    y = -0.5
    if (flip == TRUE) {
      y = 0.5
    }
    x = chromstart + (rownum * 2 * stepsize) - (stepsize * 
                                                  2)
    for (colnum in (rownum:ncol(hicm))) {
      x = x + stepsize
      if (flip == FALSE) {
        y = y + 0.5
      }
      if (flip == TRUE) {
        y = y - 0.5
      }
      xs = c(x - stepsize, x, x + stepsize, x, x - stepsize)
      ys = c(y, y + 0.5, y, y - 0.5, y)
      polygon(xs, ys, border = NA, col = hicmcol[colnum, 
                                                 rownum])
    }
  }
  return(list(c(min_z, max_z), palette))
}

draw_triangle_link <- function(x0,x1, y0=0, y1 =NA, lwd=1, lty=1, col='dark grey'){
  hyp <- sqrt(x0^2 + x1^2)
  h <- hyp/2
  xmid <- mean(c(x0,x1))
  if(is.na(y1)){y1 <- h}
  segments(x0 = x0, x1 = xmid, y0 = y0, y1 = y1, lwd=lwd, lty=lty, col=col)
  segments(x0 = xmid, x1 = x1, y0 = y1, y1 = y0, lwd=lwd, lty=lty, col=col)
}


triangle_links_from_bed <- function(bed, y0=0, y1 =NA, lwd=1, lty=1, col='dark grey'){
    for(i in 1:nrow(bed)){
      draw_triangle_link(x0 = bed$start[i], x1 = bed$end[i], y0 = y0, y1 = y1, lwd=lwd, lty=lty, col=col)
    }
}
