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

remove_fully_methylated <- function(reads,thr = 0.9){
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
getRuns <- function(calls, maxGap = NULL){
  if (!is.null(maxGap)){
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
        start = min(x$start),
        end = max(x$start) + 1,
        width = end - start)
    } else {
      rle(x$mcall) %>%
      unclass() %>% as_tibble() %>%
      mutate( endi = cumsum(lengths),
              starti = c(1,dplyr::lag(endi)[-1]+1),
              start = x$start[starti],
              end = x$start[endi] + 1,
              width = end - start ) %>%
        filter( width >= 0) # remove negative widths (in case of dups, etc.)
    }
  })
  runs <- bind_rows(runs.list,.id = "run_index") 
  runs$qname = calls.keys$qname[as.numeric(runs$run_index)]
  runs[,-1]
}
getRuns_fast <- function(calls){
  calls <- calls %>%
    filter(mcall != -1) %>%
    group_by(qname)
  calls.list <- calls %>%
    group_split(keep = F)
  calls.keys <- calls %>%
    group_keys()
  runs.list <- lapply(calls.list,function(x){
    rle(x$mcall) %>%
    unclass() %>% as_tibble() %>%
    mutate( endi = cumsum(lengths),
            starti = c(1,dplyr::lag(endi)[-1]+1),
            start = x$start[starti],
            end = x$start[endi] + 1,
            width = end - start ) %>%
      filter( width >= 0) # remove negative widths (in case of dups, etc.)
  })
  runs <- bind_rows(runs.list,.id = "run_index") 
  runs$qname = calls.keys$qname[as.numeric(runs$run_index)]
  runs[,-1]
}
order_reads <- function(x,offset = 0, bounds=NULL, method = "jaccard"){
  # get boundaries of reads if not provided
  if (is.null(bounds)){
    bounds <- x%>% group_by(qname) %>%
      summarize(start = min(start),
                end = max(end))
    # label y based on order of smallest start
    # label y based what's given
    bounds<- bounds %>% #arrange(start, end) %>%
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
smoothCalls <- function(calls,reg=NULL,bandwidth = 40){
  calls <- calls %>%
    mutate(mcall = ifelse(abs(score)>1,sign(score),score)) # ceiling based on log-lik ratio 
  if (is.null(reg)) {
    xpoints <- seq(min(calls$start),max(calls$start))
  } else {
    reg <- as_tibble(reg)
    # pad by 1kb each side to include runs going outside the region
    xpoints <- seq(reg$start-1000,reg$end+1000)
  }
  ks <- ksmooth(calls$start,calls$mcall,bandwidth = bandwidth,kernel = "normal",x.points = xpoints)
  tibble(
    start = ks$x,
    llr_smooth = ks$y, 
    mcall = case_when(
      llr_smooth > 0 ~ 1,
      llr_smooth < 0 ~ 0,
      TRUE ~ -1)) 
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
