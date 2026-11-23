#### 1. Load the observed (not imputed) final panels ####
args <- commandArgs(TRUE)
stopifnot(length(args)==2, args[1] %in% c("pilot","final"))
pilot <- args[1]=="pilot"
Lout <- normalizePath(args[2],mustWork=TRUE)
input <- Sys.getenv("LAMPRO_DIVERSITY_INPUTS","data/local/diversity")
sources <- c(Lgeno="ref5.rds",Lbg="lampro_landscape_background.rds",
             Lq="snmf_selected_qmatrix_k4.csv",Lmeta="TableS1_Individual_metadata.csv",
             Lvcf="ref5_exact_demographic.vcf.gz",Lcross="TreeMix_population_label_crosswalk.csv")
sources[] <- file.path(input,sources)
stopifnot(all(file.exists(sources)))
set.seed(112326)

read <- function(key) read.csv(sources[[key]],stringsAsFactors=FALSE)
order_ids <- function(key,ids){
  stopifnot(!anyDuplicated(key),!anyDuplicated(ids),all(ids %in% key))
  match(ids,key)
}
# Decode serialized adegenet SNPbin bit arrays; check against original VCF below.
decode <- function(z) {
  ids <- attr(z,"ind.names")
  loci <- attr(z,"loc.names")
  g <- vapply(attr(z,"gen"),function(x) {
    v <- attr(x,"snp")
    stopifnot(length(v)==2L,attr(x,"ploidy")==2L)
    y <- as.integer(rawToBits(v[[1]]))[seq_along(loci)]+as.integer(rawToBits(v[[2]]))[seq_along(loci)]
    y[attr(x,"NA.posi")] <- NA_integer_
    y
  },integer(length(loci)))
  dimnames(g) <- list(loci,ids)
  t(g)
}
lg <- readRDS(sources[["Lgeno"]]); L <- decode(lg)
Lbg <- decode(readRDS(sources[["Lbg"]]))
stopifnot(identical(dim(L),c(170L,19352L)),identical(rownames(Lbg),rownames(L)),ncol(Lbg)==19131L,
          identical(L[,order_ids(colnames(L),colnames(Lbg))],Lbg))
linfo <- data.frame(id=colnames(L),scaffold=as.character(attr(lg,"chromosome")),
                    position=attr(lg,"position"),tag=sub(":.*","",colnames(L)))
stopifnot(!anyDuplicated(linfo$tag),all(attr(lg,"ploidy")==2))
lq <- read("Lq"); lq <- lq[order_ids(lq$id,rownames(L)),paste0("K",1:4)]
lm <- read("Lmeta"); lm <- lm[order_ids(lm$Sample_ID,rownames(L)),]
lc <- read("Lcross")
s <- match(rownames(L),lc$sample_id)
gpop <- as.character(attr(lg,"pop"))
site <- attr(lg,"other")$latlon
stopifnot(!anyDuplicated(lc$sample_id),sum(!is.na(s))==163L,nrow(lc)==163L,
          all(gpop==ifelse(lm$Analytical_locality=="","unknown",lm$Analytical_locality)),
          all(gpop[!is.na(s)]==lc$population[s[!is.na(s)]]),
          all(abs(site$lat[!is.na(s)]-lc$latitude[s[!is.na(s)]])<.0001),
          all(abs(site$lon[!is.na(s)]-lc$longitude[s[!is.na(s)]])<.0001),
          all(lc$population_n==as.integer(table(lc$population)[lc$population])),
          all(gpop[match(as.character(411828:411834),rownames(L))]=="colossus"),
          gpop[match("411716",rownames(L))]=="deli3b")
Lmeta <- data.frame(id=rownames(L),population=gpop,geographic_label=lc$geographic_label[s],
                    locality=lm$Locality_description,region=lc$state_from_coordinates[s],
                    lat=site$lat,lon=site$lon,included=!is.na(s),
                    reason=ifelse(!is.na(s),"",ifelse(gpop=="unknown","Unknown population","Original group n < 3")),
                    metadata_POP=lm$Analytical_locality,lq)
stopifnot(all(table(gpop)[gpop[is.na(s) & gpop!="unknown"]]<3))

#### 2. Independent genotype extraction and hand-worked estimator checks ####
check_vcf <- function(path,g,snp,can_flip=FALSE){
  # Five fixed, spread-out individuals. Read the VCF directly without changing it.
  ids <- rownames(g)[round(seq(1,nrow(g),length.out=5))]
  f <- tempfile()
  on.exit(unlink(f))
  ret <- system2("bcftools",c("query","-s",paste(ids,collapse=","),
                  "-f",shQuote("%CHROM\\t%POS\\t%ID\\t%REF\\t%ALT[\\t%GT]\\n"),shQuote(path)),stdout=f)
  stopifnot(ret==0)
  b <- read.delim(f,header=FALSE,check.names=FALSE,colClasses="character")
  b <- b[order_ids(b[[3]],snp$id),]
  stopifnot(all(b[[1]]==snp$scaffold),all(as.numeric(b[[2]])==snp$position),
            all(nchar(b[[4]])==1),all(nchar(b[[5]])==1),all(b[[4]]!=b[[5]]))
  if(!can_flip)stopifnot(all(b[[4]]==snp$ref),all(b[[5]]==snp$alt))
  d <- matrix(match(gsub("|","/",as.matrix(b[6:10]),fixed=TRUE),c("0/0","0/1","1/0","1/1")),
              nrow=nrow(b))
  d <- matrix(c(0,1,1,2)[d],nrow=nrow(b))
  original <- t(g[ids,snp$id])
  stopifnot(identical(is.na(original),is.na(d)) || all(is.na(original)==is.na(d)))
  direct <- rowSums(d!=original,na.rm=TRUE)==0
  reverse <- rowSums(2-d!=original,na.rm=TRUE)==0
  stopifnot(all(direct | (can_flip & reverse)))
  cat("Independent VCF:",basename(path),"5 samples across",nrow(b),"positions,",
      sum(!direct & reverse),"oppositely encoded loci in subsample\n")
  invisible(b[4:5])
}
lalleles <- check_vcf(sources[["Lvcf"]],L,linfo,TRUE)
linfo$ref <- lalleles[[1]]; linfo$alt <- lalleles[[2]]
stopifnot(all(L[!is.na(L)] %in% 0:2),!anyDuplicated(colnames(L)))

counts <- function(g) {
  n <- colSums(!is.na(g)); k <- colSums(g,na.rm=TRUE); h <- colSums(g==1,na.rm=TRUE)
  list(n=n,k=k,h=h,p=k/(2*n),o=h/n)
}
wc <- function(u,v,i) {
  n1 <- u$n[i]; n2 <- v$n[i]; p1 <- u$p[i]; p2 <- v$p[i]
  n <- (n1+n2)/2; nc <- 2*n-(n1^2+n2^2)/(2*n)
  bar <- (n1*p1+n2*p2)/(2*n)
  ss <- (n1*(p1-bar)^2+n2*(p2-bar)^2)/n
  h <- (u$h[i]+v$h[i])/(2*n)
  a <- (n/nc)*(ss-(bar*(1-bar)-ss/2-h/4)/(n-1))
  b <- (n/(n-1))*(bar*(1-bar)-ss/2-((2*n-1)/(4*n))*h)
  c <- h/2
  list(a=a,total=a+b+c)
}
toy1 <- counts(matrix(c(0,0,1,1),ncol=1))
toy2 <- counts(matrix(c(1,1,2,2),ncol=1))
toy <- wc(toy1,toy2,1L)
stopifnot(toy1$o==.5,toy1$p==.25,
          abs(2*toy1$p*(1-toy1$p)*8/7-3/7)<1e-12,
          abs(toy$a/toy$total-1/3)<1e-12,
          abs(choose(6,6)/choose(8,6)-1/28)<1e-12)

#### 3. Population diversity and pairwise WC84 theta ####
pg <- function(g,meta,groups,region,info,panel,B=500L) {
  stopifnot(identical(rownames(g),meta$id),ncol(g)==nrow(info))
  units <- sort(unique(groups[!is.na(groups)]))
  use <- lapply(units,function(k) which(groups==k))
  N <- lengths(use)
  z <- lapply(use,function(i) counts(g[i,,drop=FALSE]))
  elig <- vapply(seq_along(z),function(i) z[[i]]$n>=max(3,ceiling(.7*N[i])),logical(ncol(g)))
  keep <- which(N>=3)
  shared <- which(rowSums(elig[,keep,drop=FALSE])==length(keep))
  stopifnot(length(shared)>100L)
  calls <- do.call(cbind,lapply(z,function(e)e$n))
  alleles <- do.call(cbind,lapply(z,function(e)e$k))
  miss_alt <- choose(2*calls[shared,keep,drop=FALSE]-alleles[shared,keep,drop=FALSE],6)/
    choose(2*calls[shared,keep,drop=FALSE],6)
  miss_ref <- choose(alleles[shared,keep,drop=FALSE],6)/choose(2*calls[shared,keep,drop=FALSE],6)
  stat <- data.frame(Population=units,Geographic_region=unname(region[units]),n=N,
          n_genotyped=vapply(use,function(i)sum(rowSums(!is.na(g[i,,drop=FALSE]))>0),integer(1)),
          n_panel_snps=ncol(g),n_population_eligible_snps=colSums(elig),
          n_loci=ifelse(N>=3,length(shared),NA_integer_),
          Genotype_call_rate=colSums(calls)/(N*ncol(g)),Missing_fraction=1-colSums(calls)/(N*ncol(g)),
          Usable_diploid_genotypes=colSums(calls),Ho=NA_real_,He=NA_real_,FIS=NA_real_,
          Polymorphic_SNPs=NA_integer_,Polymorphic_SNP_percent=NA_real_,
          Rarefied_allelic_richness=NA_real_,Private_allelic_richness=NA_real_,
          Genomic_panel=panel,Status=ifelse(N<3,"insufficient n",ifelse(N<5,"exploratory n=3–4","descriptive n>=5")))
  for (j in keep) {
    u <- z[[j]]; k <- match(j,keep)
    freq <- u$p[shared]
    h_e <- 2*freq*(1-freq)*2*u$n[shared]/(2*u$n[shared]-1)
    stat$Ho[j] <- mean(u$o[shared]); stat$He[j] <- mean(h_e)
    stat$FIS[j] <- if(sum(h_e)>0) 1-sum(u$o[shared])/sum(h_e) else NA_real_
    stat$Polymorphic_SNPs[j] <- sum(freq>0 & freq<1)
    stat$Polymorphic_SNP_percent[j] <- 100*stat$Polymorphic_SNPs[j]/length(shared)
    stat$Rarefied_allelic_richness[j] <- mean(2-miss_alt[,k]-miss_ref[,k])
    others <- setdiff(seq_along(keep),k)
    stat$Private_allelic_richness[j] <- mean((1-miss_alt[,k])*apply(miss_alt[,others,drop=FALSE],1,prod)+
                                      (1-miss_ref[,k])*apply(miss_ref[,others,drop=FALSE],1,prod))
    raw <- as.matrix(g[use[[j]],shared[1:10],drop=FALSE])
    manual <- mean(vapply(seq_len(ncol(raw)),function(i)mean(raw[!is.na(raw[,i]),i]==1),numeric(1)))
    stopifnot(abs(manual-mean(u$o[shared[1:10]]))<1e-12)
  }
  stopifnot(all(stat$Ho[keep]>=0 & stat$Ho[keep]<=1),all(stat$He[keep]<=1),
            all(stat$Rarefied_allelic_richness[keep]>=1-1e-12 & stat$Rarefied_allelic_richness[keep]<=2+1e-12),
            all(stat$Private_allelic_richness[keep]>=0 & stat$Private_allelic_richness[keep]<=2))
  tags <- as.character(info$tag[shared]); u_tag <- unique(tags)
  boot <- if(B>0) rmultinom(B,length(u_tag),rep(1,length(u_tag))) else NULL
  mat <- matrix(NA_real_,nrow=length(units),ncol=length(units),dimnames=list(units,units))
  tab <- list()
  for (j in seq_len(length(units)-1L)) for(k in (j+1L):length(units)) {
    ok <- all(N[c(j,k)]>=3)
    value <- lo <- hi <- NA_real_; variable <- NA_integer_
    if(ok) {
      v <- wc(z[[j]],z[[k]],shared)
      variable <- sum(v$total>0)
      value <- sum(v$a)/sum(v$total)
      if(B>0 && length(u_tag)>=50 && variable>=50) {
        a0 <- rowsum(v$a,tags,reorder=FALSE)[,1]
        t0 <- rowsum(v$total,tags,reorder=FALSE)[,1]
        x <- drop(crossprod(a0,boot))/drop(crossprod(t0,boot))
        if(all(is.finite(x))) { lim <- quantile(x,c(.025,.975)); lo <- lim[1]; hi <- lim[2] }
      }
    }
    mat[j,k] <- mat[k,j] <- value
    tab[[length(tab)+1L]] <- data.frame(Population_A=units[j],Population_B=units[k],nA=N[j],nB=N[k],
       shared_eligible_SNPs=if(ok)length(shared) else NA_integer_,informative_SNPs=variable,
       RAD_blocks=if(ok)length(u_tag) else NA_integer_,FST=value,CI_low=lo,CI_high=hi,
       bootstrap_replicates=if(is.finite(lo)) B else 0L,estimator="Weir & Cockerham 1984 theta",
       Genomic_panel=panel,Status=if(!ok)"insufficient n" else if(any(N[c(j,k)]<5))"exploratory n=3–4" else "descriptive n>=5")
  }
  stopifnot(!any(is.infinite(mat)),isSymmetric(mat),all(mat[is.finite(mat)]<=1))
  qsum <- do.call(rbind,lapply(seq_along(use),function(j) {
    qq <- as.matrix(meta[use[[j]],paste0("K",1:4)])
    data.frame(Population=units[j],n=N[j],as.list(colMeans(qq)),n_below_Q070=sum(apply(qq,1,max)<.7))
  }))
  list(diversity=stat,pairwise=do.call(rbind,tab),matrix=mat,shared=info[shared,],composition=qsum)
}

regions <- function(meta,grp) sapply(split(meta$region,grp),function(s)paste(sort(unique(s[!is.na(s) & nzchar(s)])),collapse="; "))
run <- function(g,meta,info,species,panel,out,background) {
  code <- paste0("K",1:4)
  qm <- as.matrix(meta[code]); stopifnot(all(is.finite(qm)),max(abs(rowSums(qm)-1))<.003)
  meta$dominant_K <- code[max.col(qm,ties.method="first")]
  meta$max_Q <- apply(qm,1,max)
  meta$high_ancestry_group <- ifelse(meta$max_Q>=.70,meta$dominant_K,"mixed / below Q 0.70")
  original <- ifelse(meta$included,meta$population,NA_character_)
  region <- regions(meta,original)
  if(pilot) {
    use <- seq(1,ncol(g),by=23)
    pg(g[,use],meta,original,region,info[use,],panel,25L)
    cat(species,"pilot passed on",length(use),"SNPs\n")
    return(invisible(NULL))
  }
  main <- pg(g,meta,original,region,info,panel)
  K <- ifelse(meta$max_Q>=.7,meta$dominant_K,NA_character_)
  cluster <- pg(g,meta,K,setNames(rep("Ancestry group",4),code),info,panel)
  bg <- pg(g[,background],meta,original,region,info[background,],paste(panel,"sensitivity removing named candidate set"),0L)
  bg_comparison <- merge(main$diversity[,c("Population","n_loci","Ho","He","FIS")],
                        bg$diversity[,c("Population","n_loci","Ho","He","FIS")],by="Population",suffixes=c("_full","_background"))
  # Hold SNPs fixed when testing the effect of removing mixed ancestry from each geographic group.
  methods <- lapply(seq_len(nrow(main$diversity)),function(j) {
    pop <- main$diversity$Population[j]
    both <- which(original==pop); strict <- both[meta$max_Q[both]>=.7]
    loci <- match(main$shared$id,colnames(g))
    if(length(strict)>=3) {
      us <- counts(g[strict,loci,drop=FALSE])
      ok <- us$n>=max(3,ceiling(.70*length(strict)))
      ua <- counts(g[both,loci,drop=FALSE])
      val <- function(u) c(mean(u$o[ok]),mean((2*u$p*(1-u$p)*2*u$n/(2*u$n-1))[ok]))
      n_snps <- sum(ok)
      a <- if(n_snps) val(ua) else c(NA,NA); b <- if(n_snps) val(us) else c(NA,NA)
    } else {n_snps <- NA_integer_; a <- b <- c(NA,NA)}
    data.frame(Population=pop,n_geo=length(both),n_Q070=length(strict),n_same_SNPs=n_snps,
       Ho_geo=a[1],He_geo=a[2],Ho_Q070=b[1],He_Q070=b[2])
  })
  fcols <- c("Population_A","Population_B","shared_eligible_SNPs","FST")
  fcomp <- merge(main$pairwise[fcols],bg$pairwise[fcols],by=fcols[1:2],suffixes=c("_full","_background"))
  stem <- paste0("TableS_",species)
  tab <- list(main$diversity,main$pairwise,data.frame(Population=rownames(main$matrix),main$matrix,check.names=FALSE),
       cluster$diversity,cluster$pairwise,data.frame(Cluster=rownames(cluster$matrix),cluster$matrix,check.names=FALSE),
       main$composition,bg_comparison,fcomp,do.call(rbind,methods),main$shared,cluster$shared)
  suffix <- c("_population_genetic_diversity","_pairwise_FST","_pairwise_FST_matrix",
       "_K4_diversity_Q070","_K4_pairwise_FST_Q070","_K4_FST_matrix_Q070","_geographic_ancestry_composition",
       "_candidate_exclusion_diversity_sensitivity","_candidate_exclusion_FST_sensitivity","_ancestry_filter_diversity_sensitivity",
       "_shared_SNPs","_K4_shared_SNPs")
  for(i in seq_along(tab))write.csv(tab[[i]],file.path(out,paste0(stem,suffix[i],".csv")),row.names=FALSE,na="")
  write.csv(meta,file.path(out,paste0(species,"_population_definitions.csv")),row.names=FALSE,na="")
  if(species=="Lampropholis")stopifnot(nrow(main$pairwise)==91L,sum(main$diversity$n)==163L,
                nrow(main$diversity)==14L,all(is.finite(main$pairwise$FST)),sum(!meta$included)==7L)
  cat(species,nrow(main$diversity),"groups,",sum(main$diversity$n),"individuals,",ncol(g),
      "SNPs, shared",nrow(main$shared),"SNPs and",length(unique(main$shared$tag)),
      "RAD blocks,",sum(is.finite(main$pairwise$FST)),"estimable pairs; K4 counts:",cluster$diversity$n,"\n")
}
run(L,Lmeta,linfo,"Lampropholis","Lampropholis ref5 19,352 SNPs",Lout,colnames(L)%in%colnames(Lbg))
print(sessionInfo())
