#!/usr/bin/Rscript
# methylation utils
library(Rsamtools)
library(GenomicRanges)

############################################################
#
# data reading
#
############################################################
read_data <- function(infp){
    if (infp == "stdin"){
        infp <- file(infp)
    }
    read_tsv(infp,col_names=F)
}

load_bed <- function(fpath,extracols=c("regstart","regend")){
    bedcnames=c("chrom","start","end","id","score","strand")
    cnames=c(bedcnames,extracols)
    clen=count_fields(fpath,tokenizer_tsv())[1]
    bed.tb=read_tsv(fpath,col_names=cnames[1:clen])%>%
        mutate(start=start+1)
    bed.gr=GRanges(bed.tb)
    bed.gr
}

GRangesTobed <- function(data){
    bedcols = c("seqnames","start","end","id","score","strand")
    dat.tb = as.tibble(data)
    extracols = colnames(dat.tb)[which( ! colnames(dat.tb) %in% c(bedcols,"width"))]
    dat.tb[,c(bedcols,extracols)]
}

############################################################
#
# tabix functions
#
############################################################
tabix <- function(querypath,dbpath,col_names=NULL,verbose=TRUE){
        if ("GRanges" %in% class(dbpath)){
            # input region is a GRanges object
            if (verbose) cat(paste0("reading regions defined by GRanges object",
                                    " in ",querypath,"\n"))
##            strand(dbpath) = "*"
##            regions = gsub(":.,","",toString(dbpath))
##            command = paste("tabix",querypath,regions)
            raw.list = scanTabix(querypath,param=dbpath)
            region.raw = do.call(c,raw.list)
        } else {
            if (verbose) cat(paste0("reading regions defined by ",
                                    dbpath," in ",querypath,"\n"))
            command=paste("tabix",querypath,"-R",dbpath)
            if (verbose) cat(paste0(substr(command,1,500),"...\n"))
            region.raw=system(command,intern=TRUE)
        }
        if (verbose) cat("converting to tibble\n")
        region=do.call(rbind,strsplit(region.raw,"\t"))
        region.tb=as.tibble(region)
        if (dim(region.tb)[1] == 0){ return(NA) }
        if (!is.null(col_names)) colnames(region.tb)=col_names
        region.tb %>% type_convert() %>%
            distinct()
}

mbedByCall <- function(mbed,region = NULL,pad = 2000, verbose=T){
    if (verbose) cat("parsing data into single call per line\n")
    calls = NA
    if (nrow(mbed) != 0){
      if (is.null(region)) {
        out.list=lapply(seq(dim(mbed)[1]),function(i){
            read=mbed[i,]
            mstring=read$mstring    
            call=str_extract_all(mstring,"[a-z]")[[1]]
            pos=cumsum(as.numeric(strsplit(mstring,"[a-z]")[[1]]))
            out=tibble(chrom=read$chrom,
                   start=read$start+pos+1,
                   end=start,
                   qname=read$readname,
                   mcall=call,
                   score=as.numeric(strsplit(read$scores,",")[[1]]),
                   context=strsplit(read$context,",")[[1]]
                   )
            out
        })
      } else {
        regstart <- region$start - pad
        regend <- region$end + pad
        out.list=lapply(seq(dim(mbed)[1]),function(i){
            read=mbed[i,]
            mstring=read$mstring    
            call=str_extract_all(mstring,"[a-z]")[[1]]
            pos=cumsum(as.numeric(strsplit(mstring,"[a-z]")[[1]]))
            out=tibble(chrom=read$chrom,
                   start=read$start+pos+1,
                   end=start,
                   qname=read$readname,
                   mcall=call,
                   score=as.numeric(strsplit(read$scores,",")[[1]]),
                   context=strsplit(read$context,",")[[1]]
                   )
            out[out$start >= regstart & out$start <= regend,] # region
        })
      }
      calls = do.call(rbind,out.list)
      calls$mcall[which(calls$mcall=="m")]=1
      calls$mcall[which(calls$mcall=="u")]=0
      calls$mcall[which(calls$mcall=="x")]=NA
      calls$mcall=as.numeric(calls$mcall)
    }
    calls
}

redo_mcall <- function(calls,thr) {
    mcall <- rep(NA,nrow(calls))
    mcall[calls$score > thr] <- 1
    mcall[calls$score < -thr] <- 0
    calls$mcall <- mcall
    calls
}

remove_fully_methylated <- function(reads,thr = 0.7){
  mcount <- str_count(reads$mstring,"m")
  ucount <- str_count(reads$mstring,"u")
  frac <- (mcount + 1)/(mcount + ucount + 1)
  keepi <- frac < thr
  reads[keepi,]
}
    
tabix_mbed <- function(querypath,dbpath=NULL,by=c("read","call"),extcol = NULL,verbose=TRUE){
    mbedcnames=c("chrom","start","end","readname","mstring","scores","context",extcol)
    if (!is.null(dbpath)) {
        region.tb=tabix(querypath,dbpath,mbedcnames,verbose=verbose)
    }else{
        if (verbose) cat(paste0("reading the entire data of ",querypath,"\n"))
        region.tb=read_tsv(querypath,col_names=mbedcnames)
    }
    if (verbose) cat("removing redundant loci\n")
    out.tb=unique(region.tb)
    parsembed=function(mbed.tb,what,verbose=TRUE){
        if (what == "call"){
            mbedByCall(mbed.tb,verbose=verbose)
        } else if (what == "read"){
            mbed.tb
        }
    }
    if (length(by)>1){
        out <- lapply(by,function(x){
          parsembed(out.tb,x,verbose)
        })
        names(out) = by
    }else{ out <- parsembed(out.tb,by) }
    return(out)
}
    
tabix_mfreq <- function(querypath,dbpath=NULL,cov=2,trinuc_exclude="GCG",verbose=TRUE){
    mfreqcnames=c("chrom","start","strand","meth","unmeth","dinuc","trinuc")
    if (!is.null(dbpath)) {
        region.tb=tabix(querypath,dbpath,mfreqcnames,verbose=verbose)
        if (is.na(region.tb)){
            if (verbose) cat("no data\n")
            return(tibble(freq=-1)) # if there is no data, output a tibble wiht just freq where freq=-1 for now
        }
        if (verbose) cat("removing redundant loci\n")
        out.tb=unique(region.tb)
    }else{
        if (verbose) cat(paste0("reading the entire data of ",querypath,"\n"))
        out.tb=read_tsv(querypath,col_names=mfreqcnames)
    }
    if (!is.null(trinuc_exclude)){
        if (verbose) cat(paste0("removing ",trinuc_exclude,"\n"))
        out.tb=out.tb[which(out.tb$trinuc!=trinuc_exclude),]
    }
    if (verbose) cat("calculating coverage and frequency\n")
    out.tb$cov=out.tb$meth+out.tb$unmeth
    if (!is.null(cov)){
        if (verbose) cat(paste0("removing sites with less than ",cov,"x coverage\n"))        
        out.tb=out.tb[which(out.tb$cov>=cov),]
    }
    out.tb$freq=out.tb$meth/out.tb$cov
    out.tb$end=out.tb$start
    na.omit(out.tb)
}

############################################################
#
# Bulk Methylation
#
############################################################
getRegionMeth <- function(query,subject,thr=2,verbose=TRUE){
    if (class(query)[1] != "GRanges") {
        if (verbose) cat("converting data to GRanges\n")
        query.gr=GRanges(query)
    }
    if (class(subject)[1] != "GRanges"){
        if (verbose) cat("converting subject to GRanges\n")
        subject=GRanges(subject)
    }
    if (verbose) cat("finding overlaps\n")
    ovl=findOverlaps(query.gr,subject)
    freq.tb=query[queryHits(ovl),]
    freq.tb$feature.index = subjectHits(ovl)
    if (verbose) cat("calculating methylation by region\n")
    freq.tb[which(freq.tb$cov>=thr),] %>% group_by(feature.index,add=TRUE)%>%
        summarize(totcov=sum(cov),
                  numsites=n(),
                  freq=mean(freq))
}

calculate_rollmean <- function(dat.tb,rollrange=NULL,win=50){
    rollmeans=numeric()
    tbout = FALSE
    if (! "dist" %in% names(dat.tb)) {
        if ("start" %in% names(dat.tb)) {
              dat.tb= dat.tb%>%
                mutate(dist=start)
        } else if ("pos" %in% names(dat.tb)){
              dat.tb= dat.dst %>%
                mutate(dist=pos)
        }
    }
    if (is.null(rollrange)) {
        rollrange=seq(from=min(dat.tb$dist)+win,
                       to=max(dat.tb$dist)-win,by=1)
        tbout = TRUE
    }
    for ( center in rollrange ){
        dat.win=dat.tb[which(dat.tb$dist >= center-win &
                             dat.tb$dist <= center+win),]
        rollmean=mean(dat.win$freq)
#        if (is.nan(rollmean)) rollmean=rollmeans[length(rollmeans)]
#        rollmean=sum(dat.win$totmeth)/
#            (sum(dat.win$totmeth)+sum(dat.win$totunmeth))
        rollmeans=c(rollmeans,rollmean)
    }
    if (tbout){
        tibble(start=rollrange,freq=rollmeans)
    }else{
        rollmeans
    }
}
aggregate_methylation <- function(dat.dist,win=50){
    cat("aggregating and calculating rolling mean\n")
    if (! "freq" %in% names(dat.dist)) {
        dat.dist = dat.dist %>% 
          mutate(freq=meth/(meth+unmeth)) %>%
          na.omit()
    }
    dat.ag=dat.dist%>%
        group_by(dist)%>%
            summarize(freq=mean(freq, na.rm = T)) %>%
        na.omit()
    roll.range=seq(from=min(dat.ag$dist)+win/2,
                   to=max(dat.ag$dist)-win/2,by=1)
    rollmeans=calculate_rollmean(dat.ag,rollrange=roll.range,win=win/2)
    tibble(dist=roll.range,freq=rollmeans)
}

bin_pairwise_methylation <- function(data,bin_number=75,saturation_quantile=0.9){
    n = 75
    breaks= seq(0,1,length.out=n)
    names(data) = c("x","y",names(data)[3:dim(data)[2]])
    data$x = cut(data$x,breaks)
    data$y = cut(data$y,breaks)
    levels(data$x)=levels(data$y)=breaks[1:length(breaks)-1]
    hist = data %>% group_by_all() %>%
        summarize(count=n()) %>% ungroup() %>%
        mutate(x=as.numeric(as.character(x)),
               y=as.numeric(as.character(y)))
    cutoff=round(quantile(hist$count,0.9))
    hist = hist%>%mutate(count=ifelse(count>cutoff,cutoff,count))
    hist
}
gpcPeakCaller <- function(bsobj,minwin = 50,a = 0.01,cutoff = NULL, qcutoff = 0.99){
  gpc.meth <- getMeth(bsobj,type="smooth",what="perBase")[,1]
  # remove NA
  keepi <- !is.na(gpc.meth)
  bsobj <- bsobj[keepi,]
  gpc.meth <- gpc.meth[keepi]
  gpc.cov <- getCoverage(bsobj,type="Cov",what="perBase")[,1]
  gpc.m <- getCoverage(bsobj,type="M",what="perBase")[,1]
  ## compare to baseline ----
  baseline <- median(gpc.meth,na.rm = T)
  gpc.diff <- gpc.meth - baseline
  if ( is.null(cutoff)) {
    cutoff <- quantile(gpc.diff,qcutoff,na.rm = T)
  }
  gpc.direction <- ifelse(gpc.diff > cutoff, 1, -1) # cut off by qcutoff
  gpc.gr <- granges(bsobj)
  chrs <- as.character(seqnames(gpc.gr))
  pos <- start(gpc.gr)
  ## find regions of + ----
  regions <- bsseq:::regionFinder3(gpc.direction,chrs,pos)$up
  regions <- as_tibble(regions)
  ## then add average and peak accesibility, along with coverage, then binomial test ---
  regions <- regions %>%
    rowwise() %>%
    mutate(
      coverage = sum(gpc.cov[idxStart:idxEnd]),
      methylated = sum(gpc.m[idxStart:idxEnd]),
      average = mean(gpc.meth[idxStart:idxEnd]),
      peak = max(gpc.meth[idxStart:idxEnd]),
      p.value = binom.test(methylated,coverage,baseline,alternative = "greater")$p.value
    ) %>%
    ungroup()
  ## multiple test adjusting using FDR
  regions$adjusted.pval  <- p.adjust(regions$p.value,"BH")

  ## significance based on :
  ## 1. width
  ## 2. alpha
  regions %>%
    mutate(
      width = end - start + 1,
      sig = ifelse(
        adjusted.pval <= a &
          width >= minwin,"sig","insig"))
}

findDMRs <- function(bsobj,onei,twoi,regions = NULL,cutoffs = NULL,qcutoffs = c(0.01,0.99),min_width = 100, min_n = 4, a = 0.01){
  bs.sub <- bsobj[,c(onei,twoi)]
  if ( ! is.null(regions)){
    message("subsetting data in regions ")
    keepi <- which(overlapsAny(bs.sub,regions))
    bs.sub <- bs.sub[keepi,]
    # if regions provided, I don't want cutoffs
    cutoffs = c(0,0)
  }
  meth <- as.matrix(getMeth(bs.sub,type="smooth",what="perBase"))
  # remove NA
  keepi <- rowSums(!is.na(meth)) == 2
  bs.sub <- bs.sub[keepi,]
  meth <- meth[keepi,]
  cov.mat <- as.matrix(getCoverage(bs.sub,type="Cov",what="perBase"))
  m.mat <- as.matrix(getCoverage(bs.sub,type="M",what="perBase"))
  ## compare two to one
  meth.diff <- meth[,2] - meth[,1]
  if ( is.null(cutoffs)) {
    cutoffs <- quantile(meth.diff,qcutoffs,na.rm = T)
  }
  directions <- rep(0L,length(meth.diff))
  directions[which(meth.diff < cutoffs[1])] <- -1
  directions[which(meth.diff > cutoffs[2])] <- 1
  gr <- granges(bs.sub)
  chrs <- as.character(seqnames(gr))
  pos <- start(gr)
  # find regions of significant differences  ----
  regions.list <- bsseq:::regionFinder3(directions,chrs,pos)
  regions <- as_tibble(bind_rows(regions.list,.id = "direction"))
  # then add average and peak difference, along with coverage, and fishers exact test ---
  regions <- regions %>%
    rowwise() %>%
    mutate(
      cov_one = sum(cov.mat[idxStart:idxEnd,1]),
      cov_two = sum(cov.mat[idxStart:idxEnd,2]),
      meth_one = sum(m.mat[idxStart:idxEnd,1]),
      meth_two = sum(m.mat[idxStart:idxEnd,2]),
      average_one = mean(meth[idxStart:idxEnd,1]),
      average_two = mean(meth[idxStart:idxEnd,2]),
      meandiff = mean(meth.diff[idxStart:idxEnd]),
      maxdiff = abs(max(meth.diff[idxStart:idxEnd])),
      p.value = ifelse(direction == "up",
        fisher.test(matrix(c(meth_two,cov_two-meth_two,meth_one,cov_one-meth_one),nrow = 2),alternative = "greater")$p.value,
        fisher.test(matrix(c(meth_two,cov_two-meth_two,meth_one,cov_one-meth_one),nrow = 2),alternative = "less")$p.value)
    ) %>%
    ungroup()
  # multiple test adjusting using FDR
  regions$adjusted.pval  <- p.adjust(regions$p.value,"BH")
  # significance test :
  # 1. width
  # 2. alpha
  # 3. minsites
  regions %>%
    mutate(
      width = end - start + 1,
      sig = ifelse(
        adjusted.pval <= a &
          width >= min_width &
          n >= min_n,"sig","insig"))
}

compareRegions <- function(bsobj,onei,twoi,regions, a = 0.01){
  # instead of finding dmrs, compare regions provided
  bs.sub <- bsobj[,c(onei,twoi)]
  # subset region
  regions.gr <- GRanges(regions)
  message("subsetting data in region")
  keepi <- which(overlapsAny(bs.sub,regions.gr))
  bs.sub <- bs.sub[keepi,]
  meth <- as.matrix(getMeth(bs.sub,type="smooth",what="perBase"))
  cov.mat <- as.matrix(getCoverage(bs.sub,type="Cov",what="perBase"))
  m.mat <- as.matrix(getCoverage(bs.sub,type="M",what="perBase"))
  meth.diff <- meth[,2]-meth[,1]
  chr <- as.character(seqnames(bs.sub))
  pos <- as.character(start(bs.sub))
  # get indexes
  message("getting indexes of data in region")
  ovl <- findOverlaps(bs.sub,regions.gr)
  regs.idx <- as_tibble(ovl) %>% 
    group_by(subjectHits) %>% 
    summarize(idxStart = min(queryHits),
          idxEnd = max(queryHits))
  message("comparing")
  # then do the comparison
  regions <- regions[regs.idx$subjectHits,] %>% 
    dplyr::select(chr,start,end) %>% 
    bind_cols(regs.idx) %>% 
    dplyr::select(-subjectHits) %>%
    rowwise() %>%
    mutate(
      cov_one = sum(cov.mat[idxStart:idxEnd,1]),
      cov_two = sum(cov.mat[idxStart:idxEnd,2]),
      meth_one = sum(m.mat[idxStart:idxEnd,1]),
      meth_two = sum(m.mat[idxStart:idxEnd,2]),
      average_one = mean(meth[idxStart:idxEnd,1]),
      average_two = mean(meth[idxStart:idxEnd,2]),
      meandiff = mean(meth.diff[idxStart:idxEnd]),
      maxdiff = sign(meandiff) * max(abs(meth.diff[idxStart:idxEnd])),
      direction = ifelse(meandiff > 0, "up","down"),
      p.value = ifelse(direction == "up",
        fisher.test(matrix(c(meth_two,cov_two-meth_two,meth_one,cov_one-meth_one),nrow = 2),alternative = "greater")$p.value,
        fisher.test(matrix(c(meth_two,cov_two-meth_two,meth_one,cov_one-meth_one),nrow = 2),alternative = "less")$p.value)
    ) %>%
    ungroup()
  # multiple test adjusting using FDR
  regions$adjusted.pval  <- p.adjust(regions$p.value,"BH")
  # significance test :
  # 1. alpha
  regions %>%
    mutate(
      width = end - start + 1,
      sig = ifelse(adjusted.pval <= a ,  "sig","insig"))
}


############################################################
#
# BSseq
#
############################################################
mfreqToBsseq <- function(data,pd,label="samp",fill=0){
    # rename to make compatible
    names(data)[which(names(data)==label)]="sample"
    # make matrices
    dat.meth = data %>%
        select(chrom,start,strand,meth,sample)%>%
        spread(key=sample,value=meth)
    dat.cov = data %>%
        select(chrom,start,strand,cov,sample)%>%
        spread(key=sample,value=cov)
    meth.mat = as.matrix(dat.meth[,pd$samp])
    cov.mat = as.matrix(dat.cov[,pd$samp])
    # replace NA
    naind = which(is.na(cov.mat))
    cov.mat[naind] = fill
    meth.mat[naind] = fill
    # assert same coords and make granges objecte
    stopifnot(all.equal(dat.meth[,c("chrom","start","strand")],
                        dat.cov[,c("chrom","start","strand")]))
    # load into bsseq
    bismark = BSseq(chr = dat.cov$chrom,
                    pos = dat.cov$start,
                    M = meth.mat,
                    Cov = cov.mat,
                    pData = pd)
    bismark
}

############################################################
#
# Smoothing and single-read
#
############################################################
# https://stackoverflow.com/questions/43875716/find-start-and-end-positions-indices-of-runs-consecutive-values
getRuns <- function(calls, maxGap = NULL, pad = 0){
  if (!is.null(maxGap)){
    pad = maxGap/4
    indices <- c(0,cumsum(diff(calls$start)> maxGap)) # based on difference to previous call
    calls$indices <- indices
    calls <- calls %>%
      filter(mcall != -1) %>%
      group_by(qname,indices)
  } else { 
    calls <- calls %>%
      filter(mcall != -1) %>%
      group_by(qname)
  }
  calls.list <- calls %>%
    group_split(keep = F)
  calls.keys <- calls %>%
    group_keys()
  runs.list <- lapply(calls.list,function(x){
    if (length(unique(x$mcall)) == 1){
      tibble(lengths = nrow(x),
        values = x$mcall[1],
        endi = lengths,
        starti = 1,
        start = min(x$start) - pad,
        end = max(x$start) + 1 + pad,
        width = end - start)
    } else {
      rle(x$mcall) %>%
      unclass() %>% as_tibble() %>%
      mutate( endi = cumsum(lengths),
              starti = c(1,dplyr::lag(endi)[-1]+1),
              start = x$start[starti] - pad,
              end = x$start[endi] + 1 + pad,
              width = end - start ) %>%
        filter( width >= 0) # remove negative widths (in case of dups, etc.)
    }
  })
  runs <- bind_rows(runs.list,.id = "run_index") 
  runs$qname = calls.keys$qname[as.numeric(runs$run_index)]
  runs[,-1]
}
getRuns_fast <- function(calls,mc.cores = 12,verbose = F){
  if (verbose) { message("splitting reads by qname and chrom") }
  calls <- calls[calls$mcall != -1,] %>%
    group_by(chrom,qname)
  calls.list <- calls %>%
    group_split(keep = T)
  calls.keys <- calls %>%
    group_keys()
  runs.list <- mclapply(mc.cores = mc.cores, calls.list,function(x){
    rle(x$mcall) %>%
    unclass() %>% as_tibble() %>%
    mutate( chrom = x$chrom[1],
      qname = x$qname[1],
      endi = cumsum(lengths),
      starti = c(1,dplyr::lag(endi)[-1]+1), 
      start = x$start[starti], 
      end = x$start[endi] + 1, 
      width = end - start ) %>%
      filter( width >= 0) # remove negative widths (in case of dups, etc.)
  })
  runs <- bind_rows(runs.list) 
}
order_reads <- function(x,offset = 0, bounds=NULL, method = "jaccard"){
  # get boundaries of reads if not provided
  if (is.null(bounds)){
    bounds <- x%>% group_by(qname) %>%
      summarize(start = min(start),
                end = max(end))
    # label y based on order of smallest start
    # label y based what's given
    bounds<- bounds %>% arrange(start, end) %>%
      mutate(
        readi = seq_len(length(unique(x$qname))),
        ymin = -readi - 0.8 - offset, 
        ymax = ymin + 0.6)
  }
  x <- x %>%
    mutate(ymin = bounds$ymin[match(qname,bounds$qname)],
           ymax = bounds$ymax[match(qname,bounds$qname)])
  bounds <- bounds %>%
    mutate(ymin = bounds$ymin[match(qname,bounds$qname)],
           ymax = bounds$ymax[match(qname,bounds$qname)])
  return(list(x = x,bounds = bounds))
}
smoothCalls <- function(calls,reg=NULL,bandwidth = 40,pad = 1000){
  calls <- calls %>%
    mutate(mcall = ifelse(abs(score)>1,sign(score),score)) # ceiling based on log-lik ratio 
  if (is.null(reg)) {
    xpoints <- seq(min(calls$start),max(calls$start))
  } else {
    reg <- as_tibble(reg)
    # pad by 1kb each side to include runs going outside the region
    xpoints <- seq(reg$start-pad,reg$end+pad)
  }
  ks <- ksmooth(calls$start,calls$mcall,bandwidth = bandwidth,kernel = "normal",x.points = xpoints)
  out <- tibble(
    chrom = calls$chrom[1],
    start = ks$x,
    end = start,
    qname = calls$qname[1],
    llr_smooth = ks$y, 
    mcall = case_when(
      llr_smooth > 0 ~ 1,
      llr_smooth < 0 ~ 0,
      TRUE ~ -1)
  )
  out[out$mcall != -1,]

}
smooth_reads_in_reg <- function(calls.reg,reg,bandwidth = 40){
  # separate by read
  calls.reg <- calls.reg %>%
    group_by(qname)
  group_names <- group_keys(calls.reg)$qname
  calls.list <- calls.reg %>% 
    group_split(keep = F)
  names(calls.list) <- group_names
  # smooth by read
  smooth.list <- lapply(calls.list,smoothCalls,reg, bandwidth)
  # return the combined tibble
  bind_rows(smooth.list,.id = "qname")
}
fix_protein_runs <- function(gcruns,mod){
  # first split by value
  gcruns_closed <- gcruns %>%
    filter(values == 0)
  gcruns_open <- gcruns %>%
    filter(values == 1)
  # predict using mod and get the probability of the run being in 1st cluster
  prot_prob <- predict(mod,gcruns_closed$width)$z[,1]
  # replace protein-bound calls with a open calls
  gcruns_closed$values[which(prot_prob > 0.5)] <- 1
  # recombine 
  bind_rows(gcruns_open,gcruns_closed)
}

############################################################
#
# Misc. functions
#
############################################################
getCenter <- function(db.gr){
    resize(shift(db.gr,shift=width(db.gr)/2),width=1,ignore.strand=T)
}
    
getDistance <- function(query,subject){
    # input subject to this must be center of features (done with getCenter)
    cat("getting distances\n")
    if (class(query)[1] != "GRanges") query=GRanges(query)
    if (class(subject)[1] != "GRanges") subject=Granges(subject)
    nearest.idx=nearest(query,subject,ignore.strand=T)
    sub.hit=subject[nearest.idx]
    dist.tb=tibble(index.feature=nearest.idx,
                   chrom.feature=as.character(seqnames(sub.hit)),
                   start.feature=start(sub.hit),
                   end.feature=end(sub.hit),
                   strand.feature=as.character(strand(sub.hit)),
                   dist=ifelse(strand.feature=="-",
                               end.feature-start(query),
                               start(query)-start.feature))
    dist.tb
}
