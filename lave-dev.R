### THIS SET OF FUNCTIONS ARE FOR DEV PURPOSES ONLY
### ONCE WORKED OUT THEY CAN BE PORTED INTO lave.R 
###   - either as is or as different names
### This can be sourced after lave.R is sourced s.t. functions of the same name are over-written, using the dev versions


## I'm assuming paired plots are 0-based from now on
## PE and split read inputs should specify 0-based starts --- it should be BED/python like

get_genome_offsets <- function(genome){
  ## genome is read in using read_genome_file() fxn.
  ## assumes colnames = chr, len
  ## if len col is named size, just do g$len = g$size
  ## Assuming 0-based indexing (starts at 0)
  ## Thus the i-th chr/seq in a genome file will start at the cumulative sum of lengths up to and including the (i-1)-th chr/seq
  #genome$cumsum <- cumsum(genome$len)
  genome$cumsum <- cumsum(genome$size)
  N <- length(genome$chr)
  G <- sum(genome$size)
  genome$offset <- c(0, genome$cumsum[1:(N-1)])
  allchr <- data.frame(chr="all", size=G, cumsum=G, offset=0)
  genome <- rbind(genome,allchr)
  return(genome)
}

read_genome_file <- function(fname, as.is = TRUE, col.names=c("chr","size"), colClasses=c("character","numeric"), ...){
  ## 1. reads in file
  ## 2. computes cumsums and offsets
  ## 3. invents new chr called "all" for indexing if need be
  ## 4. returns 4 col df -- chr, size, cumsum, offset
  get_genome_offsets(read.table(fname, as.is = as.is, col.names = col.names, colClasses = colClasses, ...))
}


get_bed_offsets <- function(d, genome){
  # d is BED-like dataframe that has cols start, end
  # genome is a product of get_genome_offsets()
  targets <- as.character(unique(d$chr))
  for (target in targets){
    offset <- genome$offset[genome$chr == target]
    dsubset <- d$chr == target
    d$start[dsubset] <- d$start[dsubset] + offset
    d$end[dsubset] <- d$end[dsubset] + offset
  }
  return(d)
}

get_bedpe_offsets <- function(d, genome){
  # d is BEDPE-like dataframe that has cols: chr1, start1, end1, chr2, start2, end2
  # genome is a product of get_genome_offsets()
  targets <- as.character(unique(c(d$chr1,d$chr2)))
  for (target in targets){
    offset <- genome$offset[genome$chr == target]
    dsubset1 <- d$chr1 == target
    dsubset2 <- d$chr2 == target
    d$start1[dsubset1] <- d$start1[dsubset1] + offset
    d$end1[dsubset1] <- d$end1[dsubset1] + offset
    d$start2[dsubset2] <- d$start2[dsubset2] + offset
    d$end2[dsubset2] <- d$end2[dsubset2] + offset
  }
  return(d)
}

convert_bedpe_to_all_bed <- function(d, genome, return.mapq=FALSE){
  # Takes in bedpe-like dataframe (d) and genome dataframe output fromn get_genome_offsets()
  # adds appropriate offset to each bed1 and bed2 feature
  # takes the mid-point for each bed1 and bed2 feature (if start==end, then start==end==midpoint)
  #   assigns those midpoints to start and end respectively
  # returns new dataframe: chr, start, end
  # Wher chr is arbitrarily set to "all" for each entry
  # For input bedpe-like dataframe, 
  #   expects first 6 columns to be: chr1, start1, end1, chr2, start2, end2
  #   expects a col called sv as well (typicall filled with DEL, DUP, INV, TRA)
  #   can have any columns after that
  
  d <- get_bedpe_offsets(d, genome)
  if(return.mapq){
    return(data.frame(chr=rep("all", length(d$chr1)), start=(d$start1+d$end1)/2, end=(d$start2+d$end2)/2, sv=d$sv, mapq=d$mapq))
  } else {
    return(data.frame(chr=rep("all", length(d$chr1)), start=(d$start1+d$end1)/2, end=(d$start2+d$end2)/2, sv=d$sv))
    }
}



convert_bedpe_to_paflike <- function(d, genome, g.col1="chr", g.col2="size"){
  # Takes in bedpe-like dataframe (d) and genome dataframe that has chr/contig and size/len columns - can also be output from get_genome_offsets()
  # 
  nc <- ncol(d)
  nr <- nrow(d)
  # nr <- 1
  paf <- data.frame( matrix(rep(0, 12*nr), nrow=nr) )
  colnames(paf) <- c("query", "qlen","qstart","qend","strand","target","tlen","tstart","tend","tstrand","type","nreads")
  for (i in 1:nr){
    row <- d[i,]
    size1 <- genome[[g.col2]][genome[[g.col1]] == row$chr1]
    size2 <- genome[[g.col2]][genome[[g.col1]] == row$chr2]
    # print(row)
    # print( size1 ) #, row$start1, row$end1, row$X1, row$chr2, size2, row$start2, row$end2, row$X2, row$sv, row$mapq))
    # print("")
    paf[i,] <- c(row$chr1, size1, row$start1, row$end1, row$X1, row$chr2, size2, row$start2, row$end2, row$X2, row$sv, row$mapq)
  }
  # print("HI")
  class(paf$query) <- "character"
  # print("HI")
  class(paf$target) <- "character"
  # print("HI")
  class(paf$strand) <- "character"
  # print("HI")
  class(paf$tstrand) <- "character"
  # print("HI")
  return(paf)
}


get_target_chr <- function(d, target=NA, genome=NA, return_offset=FALSE){
  # d is BED-like dataframe that includes target value in $chr
  # target is a chr/seq name found in d$chr
  # genome is a dataframe object from genome file - 2 cols need be present: chr, len
  # return_offset - by default if a target is provided, it will be retured w/ coords from 0 to chr_length.
  #               - when a genome file is given, IF the user needs to see the offset values of the target set this to TRUE
  #               - this is not a typical need so it is not default
  ## Needs target="all"
  ## If target and genome both NA, return d
  ## If target specified but genome NA or genome present but return_offset FALSE, return d[d$chr == target,]
  ## If target NA but genome specified, return modified d s.t. that chromosome starts and ends are offset by cumulative length values up to their occurence in genome file
  ## If target specified and genome file specified, return only the target from modified d: mod_d[mod_d$chr == target,] -- this option is not expected to be used much. 
  
  if(is.na(target) & sum(is.na(genome))>0){ ## If target and genome both NA, return d
    print(10)
    d.out <- d
  } else if(!(is.na(target)) & target!="all" & (sum(is.na(genome))>0 | !(return_offset))){
    print(11)
    ## If target specified but genome NA or genome present but return_offset FALSE, return d[d$chr == target,]
    d.out <- d[d$chr == target,]
  } else if(sum(is.na(genome))==0){
    print(12)
    # ##get genome offsets and make mod_d
    # genome <- get_genome_offsets(genome)
    # genome should have already come with offsets
    d <- get_bed_offsets(d, genome)
    if (is.na(target) | target=="all"){
      print(121)
      d.out <- d
    } else if (is.na(target) & return_offset){## target is specified, user wants offset values
      print(122)
      d.out <- d[d$chr == target,]}
  }
  
  else { ##catch-all -- be careful! (same behavior as both NA)
    print(13)
    d.out <- d
  }
  
  return(d.out)
}


prep_paired_loci_values <- function(data, target, genome, shortest, longest, makefullmatrix, bx, by, gridsize, maxz, xlim, ylim, znorm, return_offset){
  ## PREPARE DATA SET
  pairs <- get_target_chr(data, target, genome, return_offset)
  filtpairs <- get_target_insert_sizes(pairs, shortest,longest)
  if(makefullmatrix){filtpairs <- prep_for_full_matrix(filtpairs)}
  
  ## GET END COORD
  if(sum(is.na(genome))==0){
    if(!is.na(target)){
      print(target)
      maxpair <- genome$size[genome$chr == target]
    } else{ 
      maxpair <- genome$size[genome$chr == data$chr[1]]
    }
  }else{
    maxpair <- max(filtpairs[,2:3], na.rm = TRUE) ## 2018Aug29-note: this seems to determine chrom end via right-most coord of read mappings -- should allow a .genome file to be passed for actual chrom lens
  }
  
  ## GET AND RETURN SMOOTHED XYZ MATRIX VARIABLES
  # xyz <- get_xyz_vars(filtpairs, bx, by, gridsize, maxpair, maxz, xlim, ylim, znorm)
  get_xyz_vars(filtpairs, bx, by, gridsize, maxpair, maxz, xlim, ylim, znorm)
}


pairedLociLevelPlot <- function(data, UR, target=NA, genome=NA, shortest=0, longest=1e10, bx=1e4, by=1e4, gridsize=200, makefullmatrix=FALSE, seglwd=2, segcol="black", matcols=c("dark blue","white","red"), matcolreps=NA, xlim=c(NA), ylim=c(NA), maxz=NA, znorm=1, cnnormdata=c(NA), finalnormbymaxz=FALSE, returnlevelplot=TRUE, plotcn=TRUE, return_offset=FALSE){
  
  ## BED start and end
  # print("A")
  target.urs <- get_target_chr(UR, target, genome, return_offset)
  print(target.urs)
  ust <- target.urs$start
  uen <- target.urs$end
  
  ## PREPARE PAIRED LOCI VALUES
  # print("B")
  xyz <- prep_paired_loci_values(data, target, genome, shortest, longest, makefullmatrix, bx, by, gridsize, maxz, xlim, ylim, znorm, return_offset)
  
  ## POTENTIALLY NORMALIZE BY COPY NUMBER OR OTHER
  ## For now it also uses znorm -- but I have considered not -- all in all it is a scalar so will not hurt
  # print("C")
  if(sum(is.na(cnnormdata)==0)){
    ## PREPARE PAIRED LOCI VALUES FOR COVERAGE DATASET
    # print("D")
    cn.xyz <- prep_paired_loci_values(cnnormdata, target=target, genome=genome, shortest=0, longest=1e12, makefullmatrix=makefullmatrix, bx=bx, by=by, gridsize=gridsize, maxz=NA, xlim=xlim, ylim=ylim, znorm=znorm, return_offset)
    
    
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



