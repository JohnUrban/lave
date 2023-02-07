#lave-assembalign-dev.R




##### V02: REDO ASSEMBALIGN WITH MORE GENERIC VARIABLE NAMES AND PLOT LABELING, RDNA UPDATES INCL ARROWS, COV, ETC
##### V03: MODULARIZATION, DOTPLOT, ETC
assembalign <- function(querylengths=NA, targetlengths=NA, paf=NA, targetgaps=NA, querygaps=NA, 
                        rc.querylist=NA, rc.targetlist=NA, cadd=c(1,2,3), step=NA, goncol=NA, bgoncol=NA, 
                        pos.goncol=NA, pos.bgoncol="black", neg.goncol=NA, neg.bgoncol="black", tigcol=NA, 
                        btigcol=NA, qtigcol=NA, qbtigcol="black", ttigcol=NA, tbtigcol="black", querylengthnorm=1, 
                        targetlengthnorm=1,querytiglabels=TRUE, targettiglabels=TRUE, ylabels=c("target","query"), 
                        gapcol="white", gapbcol="white", xticknorm.target=1e6, xticknorm.query=1e6, 
                        plotqueryticks=FALSE, xlab="Pos (Mb)", tstep=50e6, qstep=50e6, ttickline=0, qtickline=0, 
                        t.padj=0, q.padj=0, xlabline=1.5, font=1, ygonpars=NA, ylim=c(-0.01,1.01), qgonlimits="small",
                        tgonlimits="small", qstart=NA, qend=NA, tstart=NA, tend=NA, qoffset=0, toffset=0, 
                        longest2shortest=TRUE, qspecificticks=FALSE, tspecificticks=FALSE, xlim=NA, force.int=TRUE, 
                        noyaxisticks=FALSE, inner.qlabels=FALSE, inner.tlabels=FALSE, reverse.qlabels=FALSE, 
                        xtext.line=1.5, xtext.cex=0.75, orderOfAppearance=FALSE, targetOrder=NA, queryOrder=NA, 
                        changestrand=TRUE,  gappy=NA, gap.scale=0.1, gapborderwidth=0.1, y.qgons=1, y.tgons=0, qgon.halfwidth=NA, 
                        tgon.halfwidth=NA, add.new.level=NA, useQueryLengths.start=FALSE, useTargetLengths.start=FALSE, 
                        pairwiseAlnGons=TRUE, ...){
  ## VERSION 03: modularization.
  ##    Trying to modularize chunks of code into callable functions that can be shared between 
  ##    assembalign, assembalign.dotplot, and future functions that expand upon both.
  
  ## paf need not be PAF object, but should be dataframe w/ following indexes: query, qstart, qend, target, tstart, tend, strand, and mapq --- mapq not currently used, but will have option soon
  ##  -- if do not have genome files (querylengths, targetlengths) -- then they are learned from PAF in which case you need the DataFrame to be set up as a PAF: q,qlen,qend,strand,t,tlen,tstart,tend,...
  ## querylengths and targetlengths should be renamed -- they specify data frames w/ 2 columns: seqnames, seqlengths
  ## targetgaps and querygaps would be better described as BED-like dataframes for the query and target -- i.e. where the BED ticks are put in the plt (top or bottom)
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
  ##    This will color both query and target tigs the same color
  ## qtigcol=NA, qbtigcol="black", ttigcol=NA, tbtigcol="black", 
  ##    assign different colors to query and target
  ##    If these are left NA they default to what is set for tigcol and btigcol
  ##querylengthnorm=1, targetlengthnorm=1,
  ##    These will multiply the sequence lengths by given number
  ##    TODO:::: THIS DOES NOT YET AFFECT ADD-ONS SUCH AS BED LOCATIONS
  ##querytiglabels=TRUE,targettiglabels=TRUE,ylabels=c("query","target")
  ##  Whether or not to use tiglabels for query or target
  ## ylabels=c("query","target")
  ##    set to whatever you want (e.g. c("","") for no ylabels)
  ## gapcol="white", gapbcol="white" -- colors of gaps (BED features) and their borders
  ## xticknorm.target=1e6, xticknorm.query=1e6, plotqueryticks=FALSE
  ## orderOfAppearance -- the contigs will appear in the order they appear in the PAF, not shortToLong, LongToShort
  ##    Overrides longestToShortest TRUE/FALSE.
  ## changestrand -- when revcomping, this changes the strand of the output so colors match accordingly. If something mapped to the neg strand of fwd, then it maps to pos of rev.
  ##
  ## add.new.level: NA by default. Use with "query" or "target". Make sure qgon/tgon params set right to add the new tigGons and alnGons to the right area of plot. Make sure your first call to assembalin provided big enough space.
  ##
  ## 2019-05-04: useQueryLength.start=FALSE, useTargetLength.start=FALSE
  ##  These are band-aids for the building multiple level plots
  ##  I don't know how other changes to the code would affect previous applications, so this option is just added to override way of doing things for current application (multiple level assembly alignments)
  
  ## COLOR CODE: tigcol=NA, btigcol=NA, qtigcol=NA, qbtigcol="black", ttigcol=NA, tbtigcol="black"
  ## Phasing out tigcol and btigcol in favor of sep cols for target and query
  ## The next 2 lines will help older analyses continue to work in a legacy fashion
  if(!(is.na(tigcol))){qtigcol <- tigcol; ttigcol <- tigcol}
  if(!(is.na(btigcol))){qbtigcol <- btigcol; tbtigcol <- btigcol}
  
  ## Phasing out step in favor of qstep and tstep
  if(!(is.na(step))){tstep <- step; qstep <- step}
  if(!(is.na(goncol))){pos.goncol <- goncol; neg.goncol <- goncol}
  if(!(is.na(bgoncol))){pos.bgoncol <- bgoncol; neg.bgoncol <- bgoncol}
  
  ## GON LIMITS: SETTING UPPER AND LOWER BOUNDARIES ON CONTIG POLYGONS
  ## IF HALFWIDTH PROVIDED, IT WILL DEFAULT TO USING Y.GON +/- HALFWIDTH
  ## IF NOT, IT WILL EXPLORE [Q/T]GONLIMITS VARABLE BY DEFAULT
  ## IF GONLIMITS VAR IS LEN=1 AND IN [SMALL, MEDIUM, LARGE, XLARGE], IT WILL USE Y.GON +/- PREDIFINED UP/DOWN SIZES
  ## IF GONLIMITS VAR IS LEN=2, THEN IT WILL USE THOSE AS BOTTOM AND TOP VALUES IN VECTOR: c(b,t,t,b,b)
  ## IF GONLIMITS VAR IS LEN=5, IT IS ASSUMED THAT IT IS THE VECTOR: c(b,t,t,b,b)
  qgonlimits <- assembalign.get.gonlimits(gonlimits=qgonlimits, y.gons=y.qgons, gon.halfwidth=qgon.halfwidth, method="query", k=0.05)
  tgonlimits <- assembalign.get.gonlimits(gonlimits=tgonlimits, y.gons=y.tgons, gon.halfwidth=tgon.halfwidth, method="target", k=0.05)
  qgon.span <- abs(qgonlimits[2]-qgonlimits[1])
  tgon.span <- abs(tgonlimits[2]-tgonlimits[1])
  qgon.mid <- (qgonlimits[2]+qgonlimits[1])/2
  tgon.mid <- (tgonlimits[2]+tgonlimits[1])/2
  
  ## YGONPARS: PARAMETERS FOR TOP AND BOTTOM LIMITS OF ALNGONS
  if(is.na(ygonpars)){
    if(tgonlimits[2] < qgonlimits[1]){
      ygonpars <- c(tgonlimits[2],qgonlimits[1])
    } else {
      ygonpars <- c(tgonlimits[1],qgonlimits[2])
    }
  }
  
  
  ## Normalize PAF given instructions
  paf <- update.paf.given.norms(paf=paf, querylengthnorm = querylengthnorm, targetlengthnorm = targetlengthnorm)
  
  
  ## If querylengths not given, use all from longest to shortest
  ## CAUTION: querylengthnorm=1 on purpose, since that was already applied to PAF in above PAf-norm section
  if(sum(is.na(querylengths))>0){
    querylengths <- get_query_lengths_from_paf(paf, querylengthnorm=1, 
                                               longest2shortest = longest2shortest, 
                                               orderOfAppearance = orderOfAppearance,
                                               orderByGiven = queryOrder)
   } else {
    ## ensure this variable only has contigs from paf
    querylengths <- querylengths[querylengths$chr %in% unique(paf$query),]
  }
  
  
  ## If targetlengths not given, use all from longest to shortest
  ## CAUTION: targetlengthnorm=1 on purpose, since that was already applied to PAF in above PAF-norm section
  if(sum(is.na(targetlengths))>0){
    targetlengths <- get_query_lengths_from_paf(paf, querylengthnorm=1, 
                                               longest2shortest = longest2shortest, 
                                               orderOfAppearance = orderOfAppearance,
                                               orderByGiven = targetOrder)
  } else {
    ## ensure this variable only has contigs from paf
    targetlengths <- targetlengths[targetlengths$chr %in% unique(paf$target),]
  }
  
  
 
  
  
  ## REVCOMP_01: If revcomp lists are given, revcomp the PAF entries -- will return input PAF if both rc.lists are NA
  paf <- revcomppaf(paf, querylengths, targetlengths, rc.querylist, rc.targetlist, changestrand = changestrand)
  
  ## REVCOMP_02: If gaps and revcomp lists are given, revcomp the BED entries -- will return input BED is rclist is NA
  if(sum(is.na(querygaps)) == 0){querygaps <- revcompbed(querygaps, querylengths, rc.querylist)}
  if(sum(is.na(targetgaps)) == 0){targetgaps <- revcompbed(targetgaps, targetlengths, rc.targetlist)}
  
  ## REVCOMP_03: If revcomp lists given, change names to revcomp_name -- THIS STEP NEEDS TO BE AFTER ABOVE 2 STEPS : REVCOMP_01, REVCOMP_02
  querylengths <- revcomptignames(querylengths, rc.querylist)
  targetlengths <- revcomptignames(targetlengths, rc.targetlist)
  
  ## DEFINE PLOT PARAMTERS
  x.max <- max(querylengths$cumsum, targetlengths$cumsum)
  if(sum(is.na(xlim))>0){xlim <- c(0,x.max)}
  
  ## CAUTION: These params below should only be not-NA when looking at a zoom-in of a single-aln. Advanced params that will usually cause chaos if touched.
  qstartisna <- is.na(qstart)
  qendisna <- is.na(qend)
  tstartisna <- is.na(tstart)
  tendisna <- is.na(tend)
  
  ##CAUTION: perhaps only good to use these params below if there is one thing aligned to many -- or one small piece aligned to a bigger piece
  qoffsetisna <- is.na(qoffset) ## default is now 0 -- which should have same effect as NA
  toffsetisna <- is.na(toffset) ## default is now 0 -- which should have same effect as NA
  
  ## 2019-02-25 - am I now using xlim as two different variables? 
  ## -- earlier I was using xlim as what part of tig to plot; 
  ## later I introduced xlim as its interpreted for plotting... I'm not sure if I'm even using the original way anymore or not...
  if(!tendisna){xlim <- c(xlim[1], tend)}
  if(!tstartisna){xlim <- c(tstart, xlim[2])}
  
  ################################### PLOTTING ############################################
  ## INITIALIZE PLOT
  INITIALPLOT=left.as.NA(add.new.level)
  if(INITIALPLOT){
    
    assembalign.initialize_plot(x.max, xlim, ylim, font, xlab, xtext.line, xtext.cex, noyaxisticks, y.tgons, y.qgons, ylabels, tstep, 
                                targetlengthnorm, tspecificticks, xticknorm.target, ttickline, t.padj, force.int, 
                                plotqueryticks, qspecificticks, qoffset, querylengths, qstep, querylengthnorm,
                                xticknorm.query, reverse.qlabels, qtickline, q.padj, ...)
  }
  
  if( INITIALPLOT | (!INITIALPLOT & add.new.level == "query") ){
    ## Add query contig/scaffold boxes
    assembalign.addqgons(qstart=qstart, qoffset=qoffset,qstartisna=qstartisna,querylengths=querylengths,qendisna=qendisna,qtigcol=qtigcol,cadd=cadd,qgonlimits=qgonlimits,querytiglabels=querytiglabels, dotplot=FALSE, useQueryLengths.start = useQueryLengths.start, border=qbtigcol)
    if(!INITIALPLOT){
      ## ADD Y-AXIS TEXT
      if(!(noyaxisticks)){
        # axis(side=2, c(y.tgons, y.qgons), ylabels, las=1, ...)
        axis(side=2, c(y.tgons, y.qgons), c("",ylabels[2]), las=1, font=font, font.lab=font, font.axis=font, ...)
        ## TODO: if give optional tick approach that excludes the connecting line, use following line as option
        # axis(side=2, y.qgons, ylabels[2], las=1, ...)
        }
    }  
  }
  
  if( INITIALPLOT | (!INITIALPLOT & add.new.level == "target") ){
    ## Add target contig/scaffold boxes
    assembalign.addtgons(tstart=tstart, toffset=toffset, tstartisna=tstartisna, targetlengths=targetlengths, tendisna=tendisna, ttigcol=ttigcol, cadd=cadd, tgonlimits=tgonlimits, targettiglabels=targettiglabels, useTargetLengths.start = useTargetLengths.start, border=tbtigcol)
    if(!INITIALPLOT){
      ## ADD Y-AXIS TEXT
      if(!(noyaxisticks)){
        axis(side=2, c(y.tgons, y.qgons), c(ylabels[1],""), las=1, font=font, font.lab=font, font.axis=font, ...)
        ## TODO: if give optional tick approach that excludes the connecting line, use following line as option
        # axis(side=2, y.tgons, ylabels[1], las=1, ...)
        }
    }    
  }
  
  ## Draw polygon connections between query and target assemblies
  assembalign.alngons(paf=paf, ygonpars=ygonpars, querylengths=querylengths, targetlengths=targetlengths, positive.strand=positive.strand, cadd=cadd, pos.goncol=pos.goncol, neg.goncol=neg.goncol, qoffset=qoffset, toffset=toffset, bgoncol=bgoncol, pairwiseAln=pairwiseAlnGons)
  
  ## Add Scaffold Gap info if provided - THERE CAN BE NO NA VALUES IN GAPS OBJECT
  ## query gaps 
  if(left.as.NA(gappy)){gappy <- c(qgon.mid+qgon.span*gap.scale, qgon.mid-qgon.span*gap.scale, tgon.mid+tgon.span*gap.scale, tgon.mid-tgon.span*gap.scale)}
  print(gappy)
  if(sum(is.na(querygaps)) == 0){gapaln(tigs=querylengths, gaps=querygaps, top=gappy[1], bottom=gappy[2], col=gapcol, border=gapbcol, scale = querylengthnorm, borderwiddth=gapborderwidth)}
  if(sum(is.na(targetgaps)) == 0){gapaln(tigs=targetlengths, gaps=targetgaps, top=gappy[3], bottom=gappy[4], col=gapcol, border=gapbcol, scale=targetlengthnorm, borderwiddth=gapborderwidth)}
}





assembalign.dotplot <- function(querylengths=NA, targetlengths=NA, paf=NA, targetgaps=NA, querygaps=NA, 
                                rc.querylist=NA, rc.targetlist=NA, cadd=c(1,2,3), step=NA, goncol=NA, bgoncol=NA, 
                                pos.goncol=NA, pos.bgoncol="black", neg.goncol=NA, neg.bgoncol="black", tigcol=NA, 
                                btigcol="black", qtigcol=NA, qbtigcol="black", ttigcol=NA, tbtigcol="black", 
                                querylengthnorm=1, targetlengthnorm=1, querytiglabels=TRUE, targettiglabels=TRUE,
                                ylabels=c("target","query"), gapcol="white", gapbcol="white", xticknorm.target=1e6, 
                                xticknorm.query=1e6, plotqueryticks=FALSE, xlab="Pos (Mb)", tstep=50e6, qstep=50e6, 
                                ttickline=0, qtickline=0, t.padj=0, q.padj=0, xlabline=1.5, font=1, ygonpars=NA, 
                                ylim=NA, qgonlimits="small",tgonlimits="small", qstart=NA, qend=NA, tstart=NA, 
                                tend=NA, qoffset=0, toffset=0, longest2shortest=TRUE, qspecificticks=FALSE, 
                                tspecificticks=FALSE, xlim=NA, force.int=TRUE, noyaxisticks=FALSE, inner.qlabels=FALSE, 
                                inner.tlabels=FALSE, reverse.qlabels=FALSE, xtext.line=1.5, xtext.cex=0.75, 
                                orderOfAppearance=FALSE, targetOrder=NA, queryOrder=NA, segwd=2, changestrand=TRUE, 
                                gappy.scale=c(0.75, 0.25), gapborderwiddth=0.1, ...){
  ## paf need not be PAF object, but should be dataframe w/ following indexes: query, qstart, qend, target, tstart, tend, strand, and mapq --- mapq not currently used, but will have option soon
  ##  -- if do not have genome files (querylengths, targetlengths) -- then they are learned from PAF in which case you need the DataFrame to be set up as a PAF: q,qlen,qend,strand,t,tlen,tstart,tend,...
  ## querylengths and targetlengths should be renamed -- they specify data frames w/ 2 columns: seqnames, seqlengths
  ## targetgaps and querygaps would be better described as BED-like dataframes for the query and target -- i.e. where the BED ticks are put in the plt (top or bottom)
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
  ##    This will color both query and target tigs the same color
  ## qtigcol=NA, qbtigcol="black", ttigcol=NA, tbtigcol="black", 
  ##    assign different colors to query and target
  ##    If these are left NA they default to what is set for tigcol and btigcol
  ##querylengthnorm=1, targetlengthnorm=1,
  ##    These will multiply the sequence lengths by given number
  ##    TODO:::: THIS DOES NOT YET AFFECT ADD-ONS SUCH AS BED LOCATIONS
  ##querytiglabels=TRUE,targettiglabels=TRUE,ylabels=c("query","target")
  ##  Whether or not to use tiglabels for query or target
  ## ylabels=c("query","target")
  ##    set to whatever you want (e.g. c("","") for no ylabels)
  ## gapcol="white", gapbcol="white" -- colors of gaps (BED features) and their borders
  ## xticknorm.target=1e6, xticknorm.query=1e6, plotqueryticks=FALSE
  ## orderOfAppearance -- the contigs will appear in the order they appear in the PAF, not shortToLong, LongToShort
  ##    Overrides longestToShortest TRUE/FALSE.
  ## targetOrder -- given a c() vector of the names of target sequences in the order you want them to appear.
  ##              -- This fails to execute if orderByAppearance is flagged.
  ## queryOrder -- similar to target order, but for query seqs. Can give same vector, but it is kept separate to allow flexibility in orders.
  ## changestrand -- when revcomping, this changes the strand of the output so colors match accordingly. If something mapped to the neg strand of fwd, then it maps to pos of rev.
  
  ## If a paf is not given, use this as default
  if(sum(is.na(paf))>0){paf <- read.table("aln.filtered-pysorted.rnd2.scaf.forR", as.is = TRUE, col.names = c("query", "qlen", "qstart", "qend", "strand", "target", "tlen", "tstart", "tend", "nmatch","alnlen","mapq", "group"))}
  
  ## COLOR CODE: tigcol=NA, btigcol=NA, qtigcol=NA, qbtigcol="black", ttigcol=NA, tbtigcol="black"
  ## Phasing out tigcol and btigcol in favor of sep cols for target and query
  ## The next 2 lines will help older analyses continue to work in a legacy fashion
  if(!(is.na(tigcol))){qtigcol <- tigcol; ttigcol <- tigcol}
  if(!(is.na(btigcol))){qbtigcol <- btigcol; tbtigcol <- btigcol}
  ## Phasing out step in favor of qstep and tstep
  if(!(is.na(step))){tstep <- step; qstep <- step}
  if(!(is.na(goncol))){pos.goncol <- goncol; neg.goncol <- goncol}
  if(!(is.na(bgoncol))){pos.bgoncol <- bgoncol; neg.bgoncol <- bgoncol}
  
  
  
  
  ## If querylengths not given, use all from longest to shortest
  if(sum(is.na(querylengths))>0){
    querylengths <- get_query_lengths_from_paf(paf, querylengthnorm=1, 
                                               longest2shortest = longest2shortest, 
                                               orderOfAppearance = orderOfAppearance,
                                               orderByGiven = queryOrder)
  }
  
  ## If targetlengths not given, use all from longest to shortest
  if(sum(is.na(targetlengths))>0){
    #targetlengthnorm=1 on purpose, since that was already applied to PAF above in experimental section
    targetlengths <- get_target_lengths_from_paf(paf, targetlengthnorm = 1, 
                                                 longest2shortest = longest2shortest, 
                                                 orderOfAppearance = orderOfAppearance,
                                                 orderByGiven = targetOrder)
  }
  
  ## If revcomp lists are given, revcomp the PAF entries -- will return input PAF if both rc.lists are NA
  paf <- revcomppaf(paf, querylengths, targetlengths, rc.querylist, rc.targetlist, changestrand = changestrand)
  
  ## If gaps and revcomp lists are given, revcomp the BED entries -- will return input BED is rclist is NA
  if(sum(is.na(querygaps)) == 0){querygaps <- revcompbed(querygaps, querylengths, rc.querylist)}
  if(sum(is.na(targetgaps)) == 0){targetgaps <- revcompbed(targetgaps, targetlengths, rc.targetlist)}
  
  ## If revcomp lists given, change names to revcomp_name -- THIS STEP NEEDS TO BE AFTER ABOVE 2
  querylengths <- revcomptignames(querylengths, rc.querylist)
  targetlengths <- revcomptignames(targetlengths, rc.targetlist)
  
  ## Define plotting parameters
  x.max <- max(targetlengths$cumsum)
  y.max <- max(querylengths$cumsum)
  xy.max <- max(x.max, y.max)
  # x.max <- max(querylengths$cumsum, targetlengths$cumsum)
  
  # if(sum(is.na(xlim))>0){xlim <- c(0,x.max)}
  if(sum(is.na(xlim))>0){xlim <- c(-0.05*x.max,x.max)} ## need to bring boxes outside plotting area
  if(sum(is.na(ylim))>0){ylim <- c(-0.05*x.max,y.max)} ## xmax in pos1 on purpose so boxes same size
  
  
  ## qgonlimit setting
  qgonlimits <- c(xlim[1],0,0,xlim[1],xlim[1])
  tgonlimits <- qgonlimits
  
  
  
  ## These params below should only be not-NA when looking at a zoom-in of a single-aln. Advanced params that will usually cause chaos if touched.
  qstartisna <- is.na(qstart)
  qendisna <- is.na(qend)
  tstartisna <- is.na(tstart)
  tendisna <- is.na(tend)
  ## perhaps only good to use these params below if there is one thing aligned to many -- or one small piece aligned to a bigger piece
  qoffsetisna <- is.na(qoffset) ## default is now 0 -- which should have same effect as NA
  toffsetisna <- is.na(toffset) ## default is now 0 -- which should have same effect as NA
  
  ## 2019-02-25 - am I now using xlim as two different variables? 
  ## -- earlier I was using xlim as what part of tig to plot; 
  ## later I introduced xlim as its interpreted for plotting... I'm not sure if I'm even using the original way anymore or not...
  if(!tendisna){xlim <- c(xlim[1], tend)}
  if(!tstartisna){xlim <- c(tstart, xlim[2])}
  
  ## Begin plotting
  
  # plot(seq(0,x.max, length.out = 16), seq(0,x.max, length.out = 16), xlim=xlim, ylim=xlim, type="n", yaxt="n", xaxt="n", xlab="", ylab="", font=font, ...)
  plot(seq(0,x.max, length.out = 16), seq(0,y.max, length.out = 16), xlim=xlim, ylim=ylim, type="n", 
       yaxt="n", xaxt="n", xlab="", ylab="", font=font, ...)
  
  mtext(side=1, text = xlab, line = xtext.line, cex=xtext.cex, font=font)
  mtext(side=2, text = xlab, line = xtext.line, cex=xtext.cex, font=font)
  
  
  ans.tt <- assembalign.tticks(xlim=xlim,tstep=tstep, targetlengthnorm=targetlengthnorm, tspecificticks=tspecificticks, 
                               xticknorm.target=xticknorm.target,force.int=force.int)
  axis(side=1, ans.tt$x.at, ans.tt$labs, line = ttickline, padj=t.padj, font=font)
  
  ans.qt <- assembalign.plotqueryticks(xlim=xlim, plotqueryticks=plotqueryticks, qspecificticks=qspecificticks,
                                       qoffset=qoffset,querylengths=querylengths, qstep=qstep, 
                                       querylengthnorm=querylengthnorm, xticknorm.query=xticknorm.query, 
                                       force.int=force.int,reverse.qlabels=reverse.qlabels)
  if(plotqueryticks){axis(side=2, at = ans.qt$x.at, labels = ans.qt$labs, line = qtickline, padj=q.padj, font=font, las=1)}
  
  
  ## Add query contig/scaffold boxes
  assembalign.addqgons(qoffset=qoffset, qstartisna=qstartisna, querylengths=querylengths, qendisna=qendisna, 
                       qtigcol=qtigcol, cadd=cadd, qgonlimits=qgonlimits, querytiglabels=querytiglabels, 
                       dotplot=TRUE, border=qbtigcol)
  
  ## Add target contig/scaffold boxes
  assembalign.addtgons(toffset=toffset, tstartisna=tstartisna, targetlengths=targetlengths, tendisna=tendisna, 
                       ttigcol=ttigcol, cadd=cadd, tgonlimits=tgonlimits, targettiglabels=targettiglabels, border=tbtigcol)
  
  
  ## Draw dotplot line segments: code was basically transferred from pairwise approach, and "polygons" replaces with "segments". Can probably be optimized...
  for (i in 1:dim(paf)[1]){
    query <- paf$query[i]
    target <- paf$target[i]
    qstart <- querylengths$starts[querylengths$chr == query] + qoffset
    tstart <- targetlengths$starts[targetlengths$chr == target] + toffset
    positive.strand <- paf$strand[i] == "+"

    ## x0 and x1 (target): start/end always stay same :: query changes if negstrand
    x0 = paf$tstart[i]+tstart
    x1 = paf$tend[i]+tstart
    if(positive.strand){
      y0 = paf$qstart[i]+qstart
      y1 = paf$qend[i]+qstart
      if(sum(is.na(pos.goncol))>0){
        pgoncol <- i+cadd[3]
      }
      else{
        pgoncol<-pos.goncol
      }
    }
    else{ # negative strand
      y0 = paf$qend[i]+qstart
      y1 = paf$qstart[i]+qstart
      if(sum(is.na(neg.goncol))>0){
        pgoncol <- i+cadd[3]
      }
      else{
        pgoncol<-neg.goncol
      }
    } 
    #segments(x0 = paf$tstart[i]+tstart, x1 = paf$tend[i]+tstart, y0 = paf$qstart[i]+qstart, y1 = paf$qend[i]+qstart, col=pgoncol, lwd=segwd)
    segments(x0 = x0, x1 = x1, y0 = y0, y1 = y1, col=pgoncol, lwd=segwd)
    
  }
  
  # ## Add Scaffold Gap info if provided - THERE CAN BE NO NA VALUES IN GAPS OBJECT
  # ## query gaps
  # gappy.scale=c(0.75, 0.25)
  if(sum(is.na(querygaps)) == 0){gapaln(tigs=querylengths, gaps=querygaps, top=xlim[1]*gappy.scale[1], 
                                        bottom=xlim[1]*gappy.scale[2], col=gapcol, border=gapbcol, 
                                        scale = querylengthnorm, yaxisAln = TRUE) #, polygon.border.width=gapborderwidth)}
  }
  if(sum(is.na(targetgaps)) == 0){gapaln(tigs=targetlengths, gaps=targetgaps, top=xlim[1]*gappy.scale[1], 
                                         bottom=xlim[1]*gappy.scale[2], col=gapcol, border=gapbcol, 
                                         scale=targetlengthnorm) #, polygon.border.width==gapborderwidth)}
  }
}
  








####### HELPERS #####################################################

left.as.NA <- function(x){
  # E.g.: sum(is.na(querylengths))>0
  return( sum(is.na(x)) > 0 )
}

get.gonlimits <- function(bottom.and.top=NA, mid.and.up.down=NA){
  ## ASSERTION: Can only use one of the options.
  OPT1 <- !left.as.NA(bottom.and.top)
  OPT2 <- !left.as.NA(mid.and.up.down)
  both.used <- OPT1 & OPT2
  none.used <- !(OPT1 | OPT2)
  ## GET GON LIMITS
  if(both.used | none.used){
    return(NA)
  } else if(OPT1){
    b <- bottom.and.top[1]
    t <- bottom.and.top[2]
  } else if(OPT2){
    b <- mid.and.up.down[1] - mid.and.up.down[3]
    t <- mid.and.up.down[1] + mid.and.up.down[2]
  }
  return(c(b,t,t,b,b))
}

assembalign.get.gonlimits <- function(gonlimits, y.gons, gon.halfwidth=NA, method="query", k=0.05){
  ## method = query or target
  ## GON LIMITS: SETTING UPPER AND LOWER BOUNDARIES ON CONTIG POLYGONS
  ## IF HALFWIDTH PROVIDED, IT WILL DEFAULT TO USING Y.GON +/- HALFWIDTH
  ## IF NOT, IT WILL EXPLORE [Q/T]GONLIMITS VARABLE BY DEFAULT
  ## IF GONLIMITS VAR IS LEN=1 AND IN [SMALL, MEDIUM, LARGE, XLARGE], IT WILL USE Y.GON +/- PREDIFINED UP/DOWN SIZES
  ## IF GONLIMITS VAR IS LEN=2, THEN IT WILL USE THOSE AS BOTTOM AND TOP VALUES IN VECTOR: c(b,t,t,b,b)
  ## IF GONLIMITS VAR IS LEN=5, IT IS ASSUMED THAT IT IS THE VECTOR: c(b,t,t,b,b)
  
  if(method=="query"){d.params <- data.frame(small=c(k,k), medium=c(k,0.15), large=c(k,0.25), xlarge=c(k,0.5))}
  else if (method=="target"){d.params <- data.frame(small=c(k,k), medium=c(0.15,k), large=c(0.25,k), xlarge=c(0.5,k))}
  len.gonvar <- length(gonlimits)
  
  if (!left.as.NA(gon.halfwidth)){
    gonlimits <- get.gonlimits(mid.and.up.down = c(y.gons, gon.halfwidth, gon.halfwidth))
  } else if(len.gonvar == 1 & gonlimits %in% c("small","medium","large","xlarge")){
    gonlimits <- get.gonlimits(mid.and.up.down = c(y.gons, d.params[[gonlimits]][1], d.params[[gonlimits]][2]))
  } else if (len.gonvar==2){
    gonlimits <- get.gonlimits(bottom.and.top = gonlimits)
  } ##else if len5, stays same, so no need for code
  
  return(gonlimits)
}

update.paf.given.qnorm <- function(paf, querylengthnorm){
  paf$qlen <- paf$qlen*querylengthnorm
  paf$qstart <- paf$qstart*querylengthnorm
  paf$qend <- paf$qend*querylengthnorm
  return(paf)
}
update.paf.given.tnorm <- function(paf, targetlengthnorm){
  paf$tlen <- paf$tlen*targetlengthnorm
  paf$tstart <- paf$tstart*targetlengthnorm
  paf$tend <- paf$tend*targetlengthnorm
  return(paf)
}
update.paf.given.norms <- function(paf, querylengthnorm, targetlengthnorm){
  paf <- update.paf.given.qnorm(paf=paf, querylengthnorm=querylengthnorm)
  paf <- update.paf.given.tnorm(paf=paf, targetlengthnorm = targetlengthnorm)
  return(paf)
}


assembalign.initialize_plot <- function(x.max, xlim, ylim, font, xlab, xtext.line, xtext.cex, noyaxisticks, y.tgons, y.qgons, ylabels, tstep, 
                                        targetlengthnorm, tspecificticks, xticknorm.target, ttickline, t.padj, force.int, 
                                        plotqueryticks, qspecificticks, qoffset, querylengths, qstep, querylengthnorm,
                                        xticknorm.query, reverse.qlabels, qtickline, q.padj, ...){
  plot(seq(0,x.max, length.out = 16), rep(0, 16), xlim=xlim, ylim=ylim, type="n", yaxt="n", xaxt="n", xlab="", ylab="", font=font, ...)
  
  ## ADD X-AXIS TEXT
  mtext(side=1, text = xlab, line = xtext.line, cex=xtext.cex, font=font)
  
  ## ADD Y-AXIS TEXT
  if(!(noyaxisticks)){
    axis(side=2, c(y.tgons, y.qgons), ylabels, las=1, ...)
    ## TODO: give optional tick approach that excludes the connecting line, use following lines as option
    # axis(side=2, c(y.tgons, y.qgons)[1], ylabels[1], las=1, ...)
    # axis(side=2, c(y.tgons, y.qgons)[2], ylabels[2], las=1, ...)
  }
  
  ## ADD TARGET TICKS (SIDE=1, X-AXIS)
  ans.tt <- assembalign.tticks(xlim=xlim, tstep=tstep, targetlengthnorm=targetlengthnorm, tspecificticks=tspecificticks, xticknorm.target=xticknorm.target, force.int=force.int)
  axis(side=1, ans.tt$x.at, ans.tt$labs, line = ttickline, padj=t.padj, font=font)
  
  ## ADD QUERY TICKS (SIDE=1, X-AXIS)
  ans.qt <- assembalign.plotqueryticks(xlim=xlim, plotqueryticks=plotqueryticks, qspecificticks=qspecificticks, qoffset=qoffset, querylengths=querylengths, qstep=qstep, querylengthnorm=querylengthnorm, xticknorm.query=xticknorm.query, force.int=force.int, reverse.qlabels=reverse.qlabels)
  if(plotqueryticks){axis(side=3, at = ans.qt$x.at, labels = ans.qt$labs, line = qtickline, padj=q.padj, font=font, las=1)}
}



assembalign.addqgons <- function(qstart, qoffset,qstartisna,querylengths,qendisna,qtigcol,cadd,qgonlimits, querytiglabels, dotplot=FALSE, useQueryLengths.start=FALSE, border="black"){
  start <- qoffset
  if(!(qstartisna)){start<-qstart}
  if(useQueryLengths.start){start <- querylengths$starts[1]}
  ## no longer needed - qoffset defaults to 1 (perhaps should be 0) -- if(!(qoffsetisna)){start <- qoffset}
  for (i in 1:dim(querylengths)[1]){
    if(qendisna){end <- querylengths$cumsum[i]+qoffset}else{end <- qend}
    if(sum(is.na(qtigcol))>0){ptigcol <- i+cadd[1]}else{ptigcol<-qtigcol}
    gon.x <- c(start, start, end, end, start)
    gon.y <- qgonlimits
    if(dotplot){
      gon.x <- qgonlimits
      gon.y <- c(start, start, end, end, start)
    }
    polygon(x  = gon.x, y = gon.y, col = ptigcol, border=border)
    if(querytiglabels){axis(side=3, tick=FALSE, at = mean(c(querylengths$cumsum[i]+qoffset,querylengths$starts[i]+qoffset)), labels = querylengths$chr[i], las=2, cex.axis=0.6)}
    if(qstartisna){start <- end+1}else{start<-qstart}
  }
}

assembalign.addtgons <- function(tstart, toffset, tstartisna, targetlengths, tendisna, ttigcol, cadd, tgonlimits, targettiglabels, useTargetLengths.start=FALSE, border="black"){
  start <- toffset
  if(!(tstartisna)){start<-tstart}
  if(useTargetLengths.start){start <- targetlengths$starts[1]}
  #if(!(toffsetisna)){start <- toffset}
  for (i in 1:dim(targetlengths)[1]){
    if(tendisna){end <- targetlengths$cumsum[i]+toffset}else{end <- tend}
    if(sum(is.na(ttigcol))>0){ptigcol <- i+cadd[2]}else{ptigcol<-ttigcol}
    gon.x <- c(start, start, end, end, start)
    gon.y <- tgonlimits
    polygon(x = gon.x, y = gon.y, col = ptigcol, border=border)
    if(targettiglabels){axis(side=1, tick=FALSE, at = mean(c(targetlengths$cumsum[i],targetlengths$starts[i])), labels = targetlengths$chr[i], las=2, cex.axis=0.6)}
    if(tstartisna){start <- end+1}else{start<-tstart}
  }
}




assembalign.alngons <- function(paf, ygonpars, querylengths, targetlengths, positive.strand, cadd, pos.goncol, neg.goncol, qoffset, toffset, bgoncol, pairwiseAln=TRUE){
  if(!pairwiseAln){
    qdim <- dim(querylengths)
    tdim <- dim(targetlengths)
    norm <- mean( c(querylengths$cumsum[qdim[1]], targetlengths$cumsum[tdim[1]]) )
  }
  
  for (i in 1:dim(paf)[1]){
    query <- paf$query[i]
    target <- paf$target[i]
    qstart <- querylengths$starts[querylengths$chr == query] + qoffset
    tstart <- targetlengths$starts[targetlengths$chr == target] + toffset
    
    positive.strand <- paf$strand[i] == "+"
    #if(sum(is.na(goncol))>0){pgoncol <- i+cadd[3]}else{pgoncol<-goncol} ## from old way
    if(positive.strand){
      if(sum(is.na(pos.goncol))>0){
        pgoncol <- i+cadd[3]
      }
      else{
        pgoncol<-pos.goncol
      }
    }
    else{
      if(sum(is.na(neg.goncol))>0){
        pgoncol <- i+cadd[3]
      }
      else{
        pgoncol<-neg.goncol
      }
    } ## newer way
    if(pairwiseAln){
      # y = bottom, top, top, bottom, bottom
      ygon <- c(ygonpars[2], ygonpars[1], ygonpars[1], ygonpars[2], ygonpars[2])
      # x = qstart, tstart, tend, qend, qstart
      x <- c(paf$qstart[i]+qstart, paf$tstart[i]+tstart, paf$tend[i]+tstart, paf$qend[i]+qstart, paf$qstart[i]+qstart)
    } else {
      # y = bottom, midhigh, bottom, bottom, midlow, bottom, bottom
      y.mid <- mean(c(ygonpars[2], ygonpars[1]))
      meanwidth <- mean(c(paf$qend[i]-paf$qstart[i], paf$tend[i]-paf$tstart[i]))
      yfudge <- ((0.5 * meanwidth)/norm) * abs(ygonpars[2]-ygonpars[1])
      midhigh <- y.mid + yfudge
      midlow <- y.mid - yfudge
      bottom <- ygonpars[1]
      ygon <- c(bottom, midhigh, bottom, bottom, midlow, bottom, bottom)
      #print(querylengths)
      # x = x.lowstart, x.mid, x.highend, x.highstart, x.mid, x.lowend, x.lowstart
      if(paf$qstart[i]+qstart < paf$tstart[i]+tstart){
        x.lowstart <- paf$qstart[i]+qstart
        x.lowend <- paf$qend[i]+qstart
        x.highstart <- paf$tstart[i]+tstart
        x.highend <- paf$tend[i]+tstart
      } else {
        x.lowstart <- paf$tstart[i]+tstart
        x.lowend <- paf$tend[i]+tstart
        x.highstart <- paf$qstart[i]+qstart
        x.highend <- paf$qend[i]+qstart
      }
      x.mid <- mean(c(mean(c(paf$qstart[i]+qstart,paf$qend[i]+qstart)), mean(c(paf$tstart[i]+tstart,paf$tend[i]+tstart))))
      #x <- c(paf$qstart[i]+qstart, x.mid, paf$tend[i]+tstart, paf$tstart[i]+tstart, x.mid, paf$qend[i]+qstart, paf$qstart[i]+qstart)
      x <- c(x.lowstart, x.mid, x.highend, x.highstart, x.mid, x.lowend, x.lowstart)
    }
    polygon(x = x, y = ygon, col = pgoncol, border=bgoncol)
  }
}





get_opposite_strand <- function(strand){
  return( as.vector(sapply(strand, function(x){if(x=='+'){return('-')}else if(x=='-'){return('+')}})))
}

revcomppaf <- function(paf=NA, querylengths=NA, targetlengths=NA, rc.querylist=NA, rc.targetlist=NA, rename=TRUE, changestrand=FALSE){
  ## Require all arguments be provided explicitly to this fxn
  #print(1)
  if(sum(is.na(c(paf))) > 0){return("Need to provide PAF.")}
  newpaf <- paf
  if(sum(is.na(rc.querylist))==0){
    #print(2)
    for (tig in rc.querylist){
      #print(3)
      len <- querylengths$len[querylengths$chr == tig]
      newpaf$qstart[paf$query == tig] <- len - paf$qend[paf$query == tig] + 1 ## adjust pos and make the end = start
      newpaf$qend[paf$query == tig] <- len - paf$qstart[paf$query == tig] + 1 ## adjust pos and make the start = end
      if(rename){newpaf$query[paf$query == tig] <- paste0('rc_',tig)}
      if(changestrand){newpaf$strand[paf$query == tig] <- get_opposite_strand(newpaf$strand[paf$query == tig])}
    }
  }
  if(sum(is.na(rc.targetlist))==0){
    #print(4)
    for (tig in rc.targetlist){
      len <- targetlengths$len[targetlengths$chr == tig]
      newpaf$tstart[paf$target == tig] <- len - paf$tend[paf$target == tig] + 1 ## adjust pos and make the end = start
      newpaf$tend[paf$target == tig] <- len - paf$tstart[paf$target == tig] + 1 ## adjust pos and make the start = end
      if(rename){newpaf$target[paf$target == tig] <- paste0('rc_',tig)}
      if(changestrand){newpaf$strand[paf$target == tig] <- get_opposite_strand(newpaf$strand[paf$target == tig])}
    }
  }
  #print(5)
  #RETURN
  newpaf
}


finddenom <- function(x){
  try <- c(1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10)
  for(i in 1:length(try)){
    y <- x/try[i]
    if(y>=1 & y<10){return(try[i])}
  }
}


revcompbed <- function(bed=NA, tigs=NA, rc.list=NA, rename=TRUE, changestrand=TRUE){
  ## Require all arguments be provided explicitly to this fxn
  if(sum(is.na(c(bed))) > 0){return("Need to provide BED.")}
  newbed <- bed
  if(sum(is.na(rc.list))==0){
    for (tig in rc.list){
      len <- tigs$len[tigs$chr == tig]
      newbed$start[bed$chr == tig] <- len - bed$end[bed$chr == tig] ## adjust pos and make the end = start ## no +1 b/c BED 0-based
      newbed$end[bed$chr == tig] <- len - bed$start[bed$chr == tig] ## adjust pos and make the start = end ## no +1 b/c BED 0-based
      if(rename){newbed$chr[bed$chr == tig] <- paste0('rc_',tig)}
      if(changestrand & "strand" %in% colnames(bed)){
        newbed$strand[bed$chr == tig] <- get_opposite_strand(newbed$strand[bed$chr == tig])
      }
    }
  }
  #RETURN
  newbed
}

gapaln <- function(tigs=NA, gaps=NA, top=NA, bottom=NA, col="white", border="white", scale=1, offset=0, 
                   addlabels=FALSE, labelcex=1, labelpos=1, breakstr=" |-|_|,|\t|\n", laby=NA, font=2, 
                   intensity.color=NA, intensity.column="score", intensity.invert=FALSE, addstrand=FALSE, 
                   addstrandlabel=FALSE, yaxisAln=FALSE, intensity.fxn=identity, as.segments=FALSE, polygon.border.width=0.1){
  ## tigs (e.g. faltigs object), gaps (e.g. falcgaps object), top (e.g.1), bottom (e.g. 0.9)
  ## Require all arguments be provided explicitly to this fxn
  ## if addlabels set to true, labels need to be in "name" column of dataframe
  print("NEW GAPALN - lave-assembalign-dev.R")
  if(sum(is.na(c(tigs, gaps, top, bottom))) > 0){return("Need to provide all arguments.")}
  ## Ensure tigs has 4 colums
  if(!(dim(tigs)[2] == 4)){return("Tigs should be a 4-column dataframe")}
  ## Ensure tigs named correctly -- ASSUMES COLS IN CORRECT ORDER
  colnames(tigs) <- c("chr", "len", "cumsum", "starts")
  #print(gaps)
  ## ADD GAPS
  #pos <- 3
  #getpos <- c(3,0,1)
  #pos <- getpos[pos]
  #print(dim(gaps)[1])
  
  ## transform if need be
  gaps[[intensity.column]] <- intensity.fxn( gaps[[intensity.column]] )
  
  if(yaxisAln){left <- bottom; right <- top}
  ## May 2, 2019: yaxisAln is just added now to do specific task. i.e. experimental and may not work properly with all options yet.
  for (i in 1:dim(gaps)[1]){
    tig <- gaps$chr[i]
    if (tig %in% tigs$chr){
      start <- tigs$starts[tigs$chr == tig]
      if(yaxisAln){
        top <- scale*gaps$start[i]+start + offset
        bottom <- scale*gaps$end[i]+start + offset
      } else {
        left <- scale*gaps$start[i]+start + offset
        right <- scale*gaps$end[i]+start + offset
      }
      if(sum(is.na(intensity.color)==0)){
        intensity <- gaps[[intensity.column]][i]/max(gaps[[intensity.column]])
        if(intensity.invert){intensity <- 1 - intensity}
        pcol <- rgb(intensity.color[1], intensity.color[2], intensity.color[3], intensity)
      } else {pcol <- col}
      # polygon(x = c(left, left, right, right, left), y = c(bottom, top, top, bottom, bottom), col = pcol, border = border)
      # if(arrowgon){
      if(addstrand){
        # print(c("left", left))
        # print(c("right", right))
        # print(c("top", top))
        # print(c("bottom", bottom))
        draw_arrowgon(left, right, top, bottom, gaps$strand[i], pcol, border)
      }else if (as.segments){
        if(sum(is.na(intensity.color)==0)){
          y <- bottom + intensity*(top-bottom) }
        else {
          y <- (bottom+top)/2
        }
          segments(x0 = left, x1 = right, y0 = y, y1 = y, col = rgb(intensity.color[1], intensity.color[2], intensity.color[3]))
          
      } else {
        polygon(x = c(left, left, right, right, left), y = c(bottom, top, top, bottom, bottom), col = pcol, border = border, lwd=polygon.border.width)
      }
      if(addlabels){if(is.na(laby)){laby<-bottom}; text(x = (left+right)/2, y = laby, labels=paste(strsplit(gaps$name[i], split=breakstr)[[1]], collapse="\n"), cex=labelcex, pos=labelpos, srt=0, font=font)} #; pos<-getpos[pos]} 
      if(addstrandlabel){text(x = (left+right)/2, y = (top+bottom)/2, labels=gaps$strand[i], cex=labelcex, font=font)} #pos=labelpos, srt=0, }
    }
  }  
}

draw_arrowgon <- function(left, right, top, bottom, strand, col, border, xscale=0.35, yscale=0){
  stemtop <- top - yscale*abs(top-bottom)
  stembottom <- bottom + yscale*abs(top-bottom)
  midy <- top - 0.5*abs(top-bottom)
  if(strand=="+"){
    midx <- right - xscale*abs(right-left)
    polygon(x = c(left, midx, midx, right, midx, midx, left, left), y = c(stembottom, stembottom, bottom, midy, top, stemtop, stemtop, stembottom), col = col, border = border)
  } else if(strand=="-"){
    midx <- left + xscale*abs(right-left)
    polygon(x = c(left, midx, midx, right, right, midx, midx, left), y = c(midy, bottom, stembottom, stembottom, stemtop, stemtop, top, midy), col = col, border = border)
  }
}

get_query_lengths_from_paf <- function(paf,querylengthnorm=1, longest2shortest=TRUE, orderOfAppearance=FALSE, orderByGiven=NA){
  if(orderOfAppearance){
    querylengths <- unique(paf[,1:2])
  } else if ( sum(is.na(orderByGiven)) == 0 ) {
    querylengths <- unique(paf[ order(factor(paf$query, levels=orderByGiven)), 1:2])
  } else {
    querylengths <- unique(paf[order(paf$qlen, decreasing = longest2shortest),1:2])
  }
  colnames(querylengths) <- c("chr", "len")
  querylengths$len <- querylengths$len*querylengthnorm
  querylengths$cumsum <- cumsum(querylengths$len)
  if(length(querylengths$cumsum) == 1){querylengths$starts <- 0}
  else {querylengths$starts <- cumsum(c(0,querylengths$len[1:(length(querylengths$len)-1)]))}
  return(querylengths)
}

get_target_lengths_from_paf <- function(paf, targetlengthnorm=1, longest2shortest=TRUE, orderOfAppearance=FALSE, orderByGiven=NA){
  if(orderOfAppearance){
    targetlengths <- unique(paf[,6:7])
  } else if ( sum(is.na(orderByGiven)) == 0 ) {
    targetlengths <- unique(paf[ order(factor(paf$target, levels=orderByGiven)), 6:7])
  } else {
    targetlengths <- unique(paf[order(paf$tlen, decreasing = longest2shortest),6:7])
  }

  colnames(targetlengths) <- c("chr", "len")
  targetlengths$len <- targetlengths$len*targetlengthnorm
  targetlengths$cumsum <- cumsum(targetlengths$len)
  if(length(targetlengths$cumsum)==1){targetlengths$starts <- 0}
  else{targetlengths$starts <- cumsum(c(0,targetlengths$len[1:(length(targetlengths$len)-1)]))}
  return(targetlengths)
}



paf2bed <- function(paf, grab=c(6,8,9,1,12,5), names=c("chr","start","end","name", "score", "strand")){
  bed <- paf[,grab]
  colnames(bed) <- names
  return(bed)
}

draw_fwd_triangle <- function(x0, x1, y0, y1, ...){
  segments(x0=x0, x1=x0, y0=y0,y1=y1, ...)
  segments(x0=x0, x1=x1, y0=y0,y1=(y1+y0)/2, ...)
  segments(x0=x0, x1=x1, y0=y1,y1=(y1+y0)/2, ...)
}
draw_rev_triangle <- function(x0, x1, y0, y1, ...){
  segments(x0=x1, x1=x1, y0=y0,y1=y1, ...)
  segments(x0=x0, x1=x1, y0=(y1+y0)/2,y1=y0, ...)
  segments(x0=x0, x1=x1, y0=(y1+y0)/2,y1=y1, ...)
}


draw_fwd_arrow_manual <- function(x0, x1, y0, y1, headwidth, headheight, ...){
  segments(x0=x0, x1=x1, y0=y0,y1=y1, ...)
  hx0 <- x1-headwidth
  hx1 <- x1
  hy0 <- y1 - headheight/2
  hy1 <- y1 + headheight/2
  draw_fwd_triangle(hx0, hx1, hy0, hy1, ...)
}

draw_rev_arrow_manual <- function(x0, x1, y0, y1, headwidth, headheight, ...){
  segments(x0=x0, x1=x1, y0=y0,y1=y1, ...)
  hx0 <- x0
  hx1 <- x0 + headwidth
  hy0 <- y1 - headheight/2
  hy1 <- y1 + headheight/2
  draw_rev_triangle(hx0, hx1, hy0, hy1, ...)
}

draw_filled_fwd_triangle <- function(x0, x1, y0, y1, ...){
  if(length(y0) == 1){y0 <- rep(y0, length(x0))}
  if(length(y1) == 1){y1 <- rep(y1, length(x0))}
  for(i in 1:length(x0)){
    polygon(x = c(x0[i], x1[i], x0[i], x0[i]), y = c(y0[i], (y0[i]+y1[i])/2, y1[i], y0[i]), border=NA, ...)
  }
}
draw_filled_rev_triangle <- function(x0, x1, y0, y1, ...){
  if(length(y0) == 1){y0 <- rep(y0, length(x0))}
  if(length(y1) == 1){y1 <- rep(y1, length(x0))}
  for(i in 1:length(x0)){
    polygon(x = c(x0[i], x1[i], x1[i], x0[i]), y = c((y0[i]+y1[i])/2, y0[i], y1[i], (y0[i]+y1[i])/2), border=NA, ...)
  }
}

get_fwd_hx0 <- function(x0, x1, headwidth){
  n <- length(x0)
  if(n>0){
    h <- vector(length = n)
    for(i in 1:n){
      h[i] <- max(x0[i], x1[i]-headwidth)
    } 
  } else {h <- x0}
  return(h)
}

get_rev_hx1 <- function(x0, x1, headwidth){
  n <- length(x0)
  if(n>0){
    h <- vector(length = n)
    for(i in 1:n){
      h[i] <- min(x1[i], x0[i]+headwidth)
    }
  } else {h <- x1}
  return(h)
}

draw_filled_fwd_arrow_manual <- function(x0, x1, y0, y1, headwidth, headheight, trimline=TRUE, ...){
  hx0 <- get_fwd_hx0(x0,x1,headwidth)
  hx1 <- x1
  hy0 <- y0 - headheight/2
  hy1 <- y1 + headheight/2
  if(trimline){x1 <- hx0}
  segments(x0=x0, x1=x1, y0=y0,y1=y1, ...)
  draw_filled_fwd_triangle(hx0, hx1, hy0, hy1, ...)
}

draw_filled_rev_arrow_manual <- function(x0, x1, y0, y1, headwidth, headheight, trimline=TRUE, ...){
  hx0 <- x0
  hx1 <- get_rev_hx1(x0,x1,headwidth)
  hy0 <- y0 - headheight/2
  hy1 <- y1 + headheight/2
  if(trimline){x0<-hx1}
  segments(x0=x0, x1=x1, y0=y0,y1=y1, ...)
  draw_filled_rev_triangle(hx0, hx1, hy0, hy1, ...)
}

draw_fwd_arrow <- function(bed, y, headwidth, headheight, jiggle=0, ...){
  pos <- bed$strand=="+"
  x0 <- bed$start[pos]
  x1 <- bed$end[pos]
  n <- length(x0)
  y <- rep(y, n) + rnorm(n, 0, jiggle)
  draw_fwd_arrow_manual(x0,x1,y,y,headwidth, headheight, ...)
}

draw_rev_arrow <- function(bed, y, headwidth, headheight, jiggle=0, ...){
  pos <- bed$strand=="-"
  x0 <- bed$start[pos]
  x1 <- bed$end[pos]
  n <- length(x0)
  y <- rep(y, n) + rnorm(n, 0, jiggle)
  draw_rev_arrow_manual(x0,x1,y,y,headwidth, headheight, ...)
}

draw_filled_fwd_arrow <- function(bed, y, headwidth, headheight, jiggle=0, scale=1, offset=0, segmentsOnly=FALSE, ...){
  pos <- bed$strand=="+"
  x0 <- bed$start[pos]*scale + offset
  x1 <- bed$end[pos]*scale + offset
  n <- length(x0)
  y <- rep(y, n) + rnorm(n, 0, jiggle)
  if(segmentsOnly){
    segments(x0=x0, x1=x1, y0=y,y1=y, ...)
  } else {
    draw_filled_fwd_arrow_manual(x0,x1,y,y,headwidth, headheight, ...) 
  }
}

draw_filled_rev_arrow <- function(bed, y, headwidth, headheight, jiggle=0, scale=1, offset=0, segmentsOnly=FALSE, ...){
  pos <- bed$strand=="-"
  x0 <- bed$start[pos]*scale + offset
  x1 <- bed$end[pos]*scale + offset
  n <- length(x0)
  y <- rep(y, n) + rnorm(n, 0, jiggle)
  if(segmentsOnly){
    segments(x0=x0, x1=x1, y0=y,y1=y, ...)
  } else {
    draw_filled_rev_arrow_manual(x0,x1,y,y,headwidth, headheight, ...)
  }
}

draw_filled_arrows_from_paf <- function(paf, y, headwidth, headheight, jiggle=0, scale=1, offset=0, targetcoords=TRUE, draw.fwd=TRUE, draw.rev=TRUE, segmentsOnly=FALSE, ...){
  #targetcoords=FALSE to get paf from query coords
  if(targetcoords){grab<-c(6,8,9,1,12,5)} else {grab<-c(1,3,4,6,12,5)}
  bed <- paf2bed(paf, grab = grab)
  if(draw.fwd){draw_filled_fwd_arrow(bed, y, headwidth, headheight, jiggle=jiggle, scale=scale, offset=offset, segmentsOnly=segmentsOnly, ...)}
  if(draw.rev){draw_filled_rev_arrow(bed, y, headwidth, headheight, jiggle=jiggle, scale=scale, offset=offset, segmentsOnly=segmentsOnly, ...)}
}

assembalign.plotqueryticks <- function(xlim, plotqueryticks, qspecificticks, qoffset, querylengths, qstep, querylengthnorm, xticknorm.query, force.int, reverse.qlabels, ...){
  denom <- finddenom(xlim[2])
  x.ceil <- ceiling(xlim[2]/denom)
  if(plotqueryticks){
    if(qspecificticks){
      ticks.xlim <- c(qoffset,max(querylengths$cumsum+qoffset))
      denom <- finddenom(ticks.xlim[2])
      x.ceil <- ceiling(ticks.xlim[2]/denom)
      # xseq <- seq(qoffset,x.ceil*denom, qstep*querylengthnorm)
      xseq <- seq(ticks.xlim[1],x.ceil*denom, qstep*querylengthnorm)
      xseq <- xseq[xseq <= ticks.xlim[2]]
      labs <- (xseq-qoffset)/xticknorm.query/querylengthnorm
    } else {
      xseq <- seq(qoffset,ceiling(xlim[2]/denom)*denom, qstep*querylengthnorm)
      labs <- (xseq)/xticknorm.query/querylengthnorm
    }
    if(force.int){labs<-as.integer(round(labs))}
    if(reverse.qlabels){
      # x.at <- abs( (max(querylengths$cumsum)+qoffset) - xseq + qoffset ) ## need +2*qoffset here b/c you could also do (max(querylengths$cumsum)-(xseq-qoffset))+qoffset
      ## better way to express it is this:
      x.at <- abs( ticks.xlim[2] - xseq + qoffset )
    } else {
      x.at <- xseq
    }
    # print(labs)
    # print(x.at)
    # print(xseq)
    # print(xticknorm.query)
    # print(querylengthnorm)
    return(list(x.at=x.at, labs=labs))
  } else {
    return(NA)
  }
}

assembalign.tticks <- function(xlim,tstep,targetlengthnorm,tspecificticks,xticknorm.target,force.int){
  denom <- finddenom(xlim[2])
  x.ceil <- ceiling(xlim[2]/denom)
  xseq <- seq(0,x.ceil*denom, tstep*targetlengthnorm)
  if(tspecificticks){
    xseq <- xseq[xseq <= xlim[2]]
  }
  labs <-  xseq/xticknorm.target/targetlengthnorm
  if(force.int){labs<-as.integer(round(labs))}
  return(list(x.at=xseq, labs=labs))
}




















############## OLDER VERSIONS; DEPRECATED ################
# 
# assembalign.old.v02 <- function(querylengths=NA, targetlengths=NA, paf=NA, targetgaps=NA, querygaps=NA, rc.querylist=NA, rc.targetlist=NA, cadd=c(1,2,3), step=NA, goncol=NA, bgoncol=NA, pos.goncol=NA, pos.bgoncol="black", neg.goncol=NA, neg.bgoncol="black", tigcol=NA, btigcol="black", qtigcol=NA, qbtigcol="black", ttigcol=NA, tbtigcol="black", querylengthnorm=1, targetlengthnorm=1,querytiglabels=TRUE,targettiglabels=TRUE,ylabels=c("target","query"), gapcol="white", gapbcol="white", xticknorm.target=1e6, xticknorm.query=1e6, plotqueryticks=FALSE, xlab="Pos (Mb)", tstep=50e6, qstep=50e6, ttickline=0, qtickline=0, t.padj=0, q.padj=0, xlabline=1.5, font=1, ygonpars=NA, ylim=c(-0.01,1.01), qgonlimits="small",tgonlimits="small", qstart=NA, qend=NA, tstart=NA, tend=NA, qoffset=0, toffset=0, longest2shortest=TRUE, qspecificticks=FALSE, tspecificticks=FALSE, xlim=NA, force.int=TRUE, noyaxisticks=FALSE, inner.qlabels=FALSE, inner.tlabels=FALSE, reverse.qlabels=FALSE, xtext.line=1.5, xtext.cex=0.75, orderOfAppearance=FALSE, changestrand=TRUE,  gappy=c(1.01, 0.99, 0.01,-0.01), ...){
#   ## paf need not be PAF object, but should be dataframe w/ following indexes: query, qstart, qend, target, tstart, tend, strand, and mapq --- mapq not currently used, but will have option soon
#   ##  -- if do not have genome files (querylengths, targetlengths) -- then they are learned from PAF in which case you need the DataFrame to be set up as a PAF: q,qlen,qend,strand,t,tlen,tstart,tend,...
#   ## querylengths and targetlengths should be renamed -- they specify data frames w/ 2 columns: seqnames, seqlengths
#   ## targetgaps and querygaps would be better described as BED-like dataframes for the query and target -- i.e. where the BED ticks are put in the plt (top or bottom)
#   ##  In some cases when the assemblies are the same and the links are from a dataset (e.g. discordant PE or splitreads), 
#   ##    the same DataFrame can be given so the top and bottom mirror each other
#   ##    Or one can highlight different features on the top v bottom.
#   ##    TODO: In the future, perhaps diff features could be identified by colors
#   ##      -- actually I can already do that modularly with the "gapaln" function which would be better called bedaln or featurealn
#   ## the rc.lists are lists of sequences to orient as the revcomp (which specifies operations to do on paired alignments)
#   ##    This is useful to untangle knots, twists, and turns in assembly comparisons
#   ## cadd is just an integer offset to how to color polygons -- it specifies the offset for the top, bottom, and connecting polygons
#   ## step is a tick mark spacer basically 
#   ## goncol and bgoncol -- can tell it what to color the polygons and borders (overriding the cadd approach)
#   ## tigcol and btigcol -- same but for contig/seq polygons
#   ##    This will color both query and target tigs the same color
#   ## qtigcol=NA, qbtigcol="black", ttigcol=NA, tbtigcol="black", 
#   ##    assign different colors to query and target
#   ##    If these are left NA they default to what is set for tigcol and btigcol
#   ##querylengthnorm=1, targetlengthnorm=1,
#   ##    These will multiply the sequence lengths by given number
#   ##    TODO:::: THIS DOES NOT YET AFFECT ADD-ONS SUCH AS BED LOCATIONS
#   ##querytiglabels=TRUE,targettiglabels=TRUE,ylabels=c("query","target")
#   ##  Whether or not to use tiglabels for query or target
#   ## ylabels=c("query","target")
#   ##    set to whatever you want (e.g. c("","") for no ylabels)
#   ## gapcol="white", gapbcol="white" -- colors of gaps (BED features) and their borders
#   ## xticknorm.target=1e6, xticknorm.query=1e6, plotqueryticks=FALSE
#   ## orderOfAppearance -- the contigs will appear in the order they appear in the PAF, not shortToLong, LongToShort
#   ##    Overrides longestToShortest TRUE/FALSE.
#   ## changestrand -- when revcomping, this changes the strand of the output so colors match accordingly. If something mapped to the neg strand of fwd, then it maps to pos of rev.
#   
#   ## If a paf is not given, use this as default
#   if(sum(is.na(paf))>0){paf <- read.table("aln.filtered-pysorted.rnd2.scaf.forR", as.is = TRUE, col.names = c("query", "qlen", "qstart", "qend", "strand", "target", "tlen", "tstart", "tend", "nmatch","alnlen","mapq", "group"))}
#   
#   ## COLOR CODE: tigcol=NA, btigcol=NA, qtigcol=NA, qbtigcol="black", ttigcol=NA, tbtigcol="black"
#   ## Phasing out tigcol and btigcol in favor of sep cols for target and query
#   ## The next 2 lines will help older analyses continue to work in a legacy fashion
#   if(!(is.na(tigcol))){qtigcol <- tigcol; ttigcol <- tigcol}
#   if(!(is.na(btigcol))){qbtigcol <- btigcol; tbtigcol <- btigcol}
#   ## Phasing out step in favor of qstep and tstep
#   if(!(is.na(step))){tstep <- step; qstep <- step}
#   if(!(is.na(goncol))){pos.goncol <- goncol; neg.goncol <- goncol}
#   if(!(is.na(bgoncol))){pos.bgoncol <- bgoncol; neg.bgoncol <- bgoncol}
#   
#   ## qgonlimit setting
#   if(qgonlimits == "small"){qgonlimits=c(0.95,1.05,1.05,0.95,0.95)}
#   else if (qgonlimits == "medium"){qgonlimits=c(0.85,1.05,1.05,0.85,0.85)}
#   else if (qgonlimits == "large"){qgonlimits=c(0.75,1.05,1.05,0.75,0.75)}
#   else if (qgonlimits == "xlarge"){qgonlimits=c(0.5,1.05,1.05,0.5,0.5)}
#   if (tgonlimits=="small"){tgonlimits=c(-0.05,0.05,0.05,-0.05,-0.05)}
#   else if (tgonlimits=="medium"){tgonlimits=c(-0.05,0.15,0.15,-0.05,-0.05)}
#   else if (tgonlimits=="large"){tgonlimits=c(-0.05,0.25,0.25,-0.05,-0.05)}
#   else if (tgonlimits=="xlarge"){tgonlimits=c(-0.05,0.5,0.5,-0.05,-0.05)}
#   if(is.na(ygonpars)){ygonpars=c(tgonlimits[2],qgonlimits[1])}
#   
#   
#   ## EXPERIMENTAL
#   paf$qlen <- paf$qlen*querylengthnorm
#   paf$qstart <- paf$qstart*querylengthnorm
#   paf$qend <- paf$qend*querylengthnorm
#   paf$tlen <- paf$tlen*targetlengthnorm
#   paf$tstart <- paf$tstart*targetlengthnorm
#   paf$tend <- paf$tend*targetlengthnorm
#   
#   
#   ## If querylengths not given, use all from longest to shortest
#   if(sum(is.na(querylengths))>0){
#     #querylengthnorm=1 on purpose, since that was already applied to PAF above in experimental section
#     querylengths <- get_query_lengths_from_paf(paf, querylengthnorm=1, longest2shortest = longest2shortest, orderOfAppearance = orderOfAppearance)
#     
#     # old/dont need
#     # querylengths <- unique(paf[order(paf$qlen, decreasing = longest2shortest),1:2])
#     # colnames(querylengths) <- c("chr", "len")
#     # #######querylengths$len <- querylengths$len*querylengthnorm
#     # querylengths$cumsum <- cumsum(querylengths$len)
#     # if(length(querylengths$cumsum) == 1){querylengths$starts <- 0}
#     # else {querylengths$starts <- cumsum(c(0,querylengths$len[1:(length(querylengths$len)-1)]))}
#   }
#   
#   ## If targetlengths not given, use all from longest to shortest
#   if(sum(is.na(targetlengths))>0){
#     #targetlengthnorm=1 on purpose, since that was already applied to PAF above in experimental section
#     targetlengths <- get_target_lengths_from_paf(paf, targetlengthnorm = 1, longest2shortest = longest2shortest, orderOfAppearance = orderOfAppearance)
#     # targetlengths <- unique(paf[order(paf$tlen, decreasing = longest2shortest),6:7])
#     # colnames(targetlengths) <- c("chr", "len")
#     # ######targetlengths$len <- targetlengths$len*targetlengthnorm
#     # targetlengths$cumsum <- cumsum(targetlengths$len)
#     # if(length(targetlengths$cumsum)==1){targetlengths$starts <- 0}
#     # else{targetlengths$starts <- cumsum(c(0,targetlengths$len[1:(length(targetlengths$len)-1)]))}
#   }
#   
#   ## If revcomp lists are given, revcomp the PAF entries -- will return input PAF if both rc.lists are NA
#   paf <- revcomppaf(paf, querylengths, targetlengths, rc.querylist, rc.targetlist, changestrand = changestrand)
#   
#   ## If gaps and revcomp lists are given, revcomp the BED entries -- will return input BED is rclist is NA
#   if(sum(is.na(querygaps)) == 0){querygaps <- revcompbed(querygaps, querylengths, rc.querylist)}
#   if(sum(is.na(targetgaps)) == 0){targetgaps <- revcompbed(targetgaps, targetlengths, rc.targetlist)}
#   
#   ## If revcomp lists given, change names to revcomp_name -- THIS STEP NEEDS TO BE AFTER ABOVE 2
#   querylengths <- revcomptignames(querylengths, rc.querylist)
#   targetlengths <- revcomptignames(targetlengths, rc.targetlist)
#   
#   ## Define plotting parameters
#   x.max <- max(querylengths$cumsum, targetlengths$cumsum)
#   if(sum(is.na(xlim))>0){xlim <- c(0,x.max)}
#   ## These params below should only be not-NA when looking at a zoom-in of a single-aln. Advanced params that will usually cause chaos if touched.
#   qstartisna <- is.na(qstart)
#   qendisna <- is.na(qend)
#   tstartisna <- is.na(tstart)
#   tendisna <- is.na(tend)
#   ## perhaps only good to use these params below if there is one thing aligned to many -- or one small piece aligned to a bigger piece
#   qoffsetisna <- is.na(qoffset) ## default is now 0 -- which should have same effect as NA
#   toffsetisna <- is.na(toffset) ## default is now 0 -- which should have same effect as NA
#   
#   ## 2019-02-25 - am I now using xlim as two different variables? 
#   ## -- earlier I was using xlim as what part of tig to plot; 
#   ## later I introduced xlim as its interpreted for plotting... I'm not sure if I'm even using the original way anymore or not...
#   if(!tendisna){xlim <- c(xlim[1], tend)}
#   if(!tstartisna){xlim <- c(tstart, xlim[2])}
#   
#   ## Begin plotting
#   # plot(seq(0,xlim[2], length.out = 16), rep(0, 16), xlim=xlim, ylim=ylim, type="n", yaxt="n", xaxt="n", xlab="", ylab="", font=font, ...)
#   ## 2019-02-25 - it seems to produce better plots more reproducibly if I plot everything but use xlim to zoom -- rather than also using xlim to set what to plot up to
#   plot(seq(0,x.max, length.out = 16), rep(0, 16), xlim=xlim, ylim=ylim, type="n", yaxt="n", xaxt="n", xlab="", ylab="", font=font, ...)
#   
#   ################################### NEW ############################################
#   mtext(side=1, text = xlab, line = xtext.line, cex=xtext.cex, font=font)
#   
#   if(!(noyaxisticks)){axis(side=2, c(0,1), ylabels, las=1, ...)}
#   
#   ans.tt <- assembalign.tticks(xlim=xlim,tstep=tstep,targetlengthnorm=targetlengthnorm,tspecificticks=tspecificticks,xticknorm.target=xticknorm.target,force.int=force.int)
#   axis(side=1, ans.tt$x.at, ans.tt$labs, line = ttickline, padj=t.padj, font=font)
#   
#   ans.qt <- assembalign.plotqueryticks(xlim=xlim, plotqueryticks=plotqueryticks, qspecificticks=qspecificticks,qoffset=qoffset,querylengths=querylengths,qstep=qstep,querylengthnorm=querylengthnorm, xticknorm.query=xticknorm.query,force.int=force.int,reverse.qlabels=reverse.qlabels)
#   if(plotqueryticks){axis(side=3, at = ans.qt$x.at, labels = ans.qt$labs, line = qtickline, padj=q.padj, font=font, las=1)}
#   
#   ## Add query contig/scaffold boxes
#   assembalign.addqgons(qoffset=qoffset,qstartisna=qstartisna,querylengths=querylengths,qendisna=qendisna,qtigcol=qtigcol,cadd=cadd,qgonlimits=qgonlimits,querytiglabels=querytiglabels, dotplot=FALSE)
#   
#   ## Add target contig/scaffold boxes
#   assembalign.addtgons(toffset=toffset, tstartisna=tstartisna, targetlengths=targetlengths, tendisna=tendisna, ttigcol=ttigcol, cadd=cadd, tgonlimits=tgonlimits, targettiglabels=targettiglabels)
#   
#   ## Draw polygon connections between query and target assemblies
#   assembalign.alngons(paf=paf, ygonpars=ygonpars, querylengths=querylengths, targetlengths=targetlengths, positive.strand=positive.strand, cadd=cadd, pos.goncol=pos.goncol, neg.goncol=neg.goncol, qoffset=qoffset, toffset=toffset)
#   
#   ################################### NEW ############################################
#   
#   
#   ################################### OLD ############################################
#   mtext(side=1, text = xlab, line = xtext.line, cex=xtext.cex, font=font)
#   denom <- finddenom(xlim[2])
#   x.ceil <- ceiling(xlim[2]/denom)
#   xseq <- seq(0,x.ceil*denom, tstep*targetlengthnorm)
#   if(tspecificticks){
#     xseq <- xseq[xseq <= xlim[2]]
#   }
#   labs <-  xseq/xticknorm.target/targetlengthnorm
#   if(force.int){labs<-as.integer(round(labs))}
#   if(inner.tlabels){
#     axis(side=3, pos=ygonpars[1], xseq, labs, line = ttickline, padj=t.padj, font=font)
#   } else {
#     axis(side=1, xseq, labs, line = ttickline, padj=t.padj, font=font)
#   }
#   if(!(noyaxisticks)){axis(side=2, c(0,1), ylabels, las=1, ...)}
#   if(plotqueryticks){
#     if(qspecificticks){
#       ticks.xlim <- c(qoffset,max(querylengths$cumsum+qoffset))
#       denom <- finddenom(ticks.xlim[2])
#       x.ceil <- ceiling(ticks.xlim[2]/denom)
#       # xseq <- seq(qoffset,x.ceil*denom, qstep*querylengthnorm)
#       xseq <- seq(ticks.xlim[1],x.ceil*denom, qstep*querylengthnorm)
#       xseq <- xseq[xseq <= ticks.xlim[2]]
#       labs <- (xseq-qoffset)/xticknorm.query/querylengthnorm
#     } else {
#       xseq <- seq(qoffset,ceiling(xlim[2]/denom)*denom, qstep*querylengthnorm)
#       labs <- (xseq)/xticknorm.query/querylengthnorm
#     }
#     if(force.int){labs<-as.integer(round(labs))}
#     if(reverse.qlabels){
#       # x.at <- abs( (max(querylengths$cumsum)+qoffset) - xseq + qoffset ) ## need +2*qoffset here b/c you could also do (max(querylengths$cumsum)-(xseq-qoffset))+qoffset
#       ## better way to express it is this:
#       x.at <- abs( ticks.xlim[2] - xseq + qoffset )
#     } else {
#       x.at <- xseq
#     }
#     if(inner.qlabels){
#       axis(side=1, pos=ygonpars[2], x.at, labs, line = qtickline, padj=q.padj, font=font)
#     } else {
#       axis(side=3, x.at, labs, line = qtickline, padj=q.padj, font=font)
#     }
#     
#   }
#   
#   ## Add query contig/scaffold boxes
#   start <- qoffset
#   if(!(qstartisna)){start<-qstart}
#   ## no longer needed - qoffset defaults to 1 (perhaps should be 0) -- if(!(qoffsetisna)){start <- qoffset}
#   for (i in 1:dim(querylengths)[1]){
#     if(qendisna){end <- querylengths$cumsum[i]+qoffset}else{end <- qend}
#     if(sum(is.na(qtigcol))>0){ptigcol <- i+cadd[1]}else{ptigcol<-qtigcol}
#     polygon(x = c(start, start, end, end, start), y = qgonlimits, col = ptigcol)
#     if(querytiglabels){axis(side=3, tick=FALSE, at = mean(c(querylengths$cumsum[i]+qoffset,querylengths$starts[i]+qoffset)), labels = querylengths$chr[i], las=2, cex.axis=0.6)}
#     if(qstartisna){start <- end+1}else{start<-qstart}
#   }
#   ## Add target contig/scaffold boxes
#   start <- toffset
#   if(!(tstartisna)){start<-tstart}
#   #if(!(toffsetisna)){start <- toffset}
#   for (i in 1:dim(targetlengths)[1]){
#     if(tendisna){end <- targetlengths$cumsum[i]+toffset}else{end <- tend}
#     if(sum(is.na(ttigcol))>0){ptigcol <- i+cadd[2]}else{ptigcol<-ttigcol}
#     polygon(x = c(start, start, end, end, start), y = tgonlimits, col = ptigcol)
#     if(targettiglabels){axis(side=1, tick=FALSE, at = mean(c(targetlengths$cumsum[i],targetlengths$starts[i])), labels = targetlengths$chr[i], las=2, cex.axis=0.6)}
#     if(tstartisna){start <- end+1}else{start<-tstart}
#   }
#   
#   ## Draw polygon connections between query and target assemblies
#   ygon <- c(1,0,0,1,1)
#   ygon <- c(0.95,0.05,0.05,0.95,0.95)
#   ygon <- c(ygonpars[2], ygonpars[1], ygonpars[1], ygonpars[2], ygonpars[2])
#   for (i in 1:dim(paf)[1]){
#     query <- paf$query[i]
#     target <- paf$target[i]
#     qstart <- querylengths$starts[querylengths$chr == query] + qoffset
#     tstart <- targetlengths$starts[targetlengths$chr == target] + toffset
#     positive.strand <- paf$strand[i] == "+"
#     #if(sum(is.na(goncol))>0){pgoncol <- i+cadd[3]}else{pgoncol<-goncol} ## from old way
#     if(positive.strand){
#       if(sum(is.na(pos.goncol))>0){
#         pgoncol <- i+cadd[3]
#       }
#       else{
#         pgoncol<-pos.goncol
#       }
#     }
#     else{
#       if(sum(is.na(neg.goncol))>0){
#         pgoncol <- i+cadd[3]
#       }
#       else{
#         pgoncol<-neg.goncol
#       }
#     } ## newer way
#     polygon(x = c(paf$qstart[i]+qstart, paf$tstart[i]+tstart, paf$tend[i]+tstart, paf$qend[i]+qstart, paf$qstart[i]+qstart), y = ygon, col = pgoncol, border=bgoncol)
#   }
#   ################################### OLD ############################################
#   
#   
#   ## Add Scaffold Gap info if provided - THERE CAN BE NO NA VALUES IN GAPS OBJECT
#   ## query gaps 
#   if(sum(is.na(querygaps)) == 0){gapaln(tigs=querylengths, gaps=querygaps, top=gappy[1], bottom=gappy[2], col=gapcol, border=gapbcol, scale = querylengthnorm)}
#   if(sum(is.na(targetgaps)) == 0){gapaln(tigs=targetlengths, gaps=targetgaps, top=gappy[3], bottom=gappy[4], col=gapcol, border=gapbcol, scale=targetlengthnorm)}
# }
# 


