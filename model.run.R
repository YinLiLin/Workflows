
P_plus_T <- function(bfile = "", pheno.col = 6, val.num = 5, r2 = 0.2, rfile = ""){
  
  fam <- read.table(paste0(bfile, ".fam"), head = FALSE)
  pheno <- fam[, pheno.col] + 3
  index <- which(!is.na(pheno))
  if(length(index) == 0)  stop("No valid observations.")
  val.fold <- cut(1 : length(index), val.num, labels=FALSE)
  
  n <- nrow(read.table(rfile, head=FALSE))
  Acc <- matrix(NA, val.num, n)
  myrfile <- rfile

  # start cross-validation
  for(i in 1 : (val.num + 1)){

    # phenotype assigning
    val.pheno <- pheno
    if(i <= val.num){
      val.pheno[index[val.fold == i]] <- NA
    }
    if(file.exists(paste0(basename(bfile), "_PT_col", pheno.col, "_val", i, ".pheno"))) stop("file exists.")
    write.table(cbind(fam[, 1 : 2], val.pheno), paste0(basename(bfile), "_PT_col", pheno.col, "_val", i, ".pheno"), col.names=FALSE, row.names=FALSE, sep=" ", quote=FALSE)
    
    # association test
    system(
      paste0(
        "plink --bfile ", bfile,
        " --pheno ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i, ".pheno"),
        " --all-pheno ",
        " --allow-no-sex",
        " --keep-allele-order",
        " --freq",
        " --assoc",
        " --out ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i)
      )
    )

    sink(paste0(basename(bfile), "_PT_col", pheno.col, "_val", i, ".sh"))
    cat("#! /bin/bash", "\n")
    cat(paste0("paste -d' ' <(awk '{print $2,$5,$6}' ", paste0(bfile, ".bim"), " | sed '1iSNP A1 A2') <(awk '{print $5}' ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i), ".frq) <(awk '{print $5,$6,$9,$4}' ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i),".P1.qassoc) > ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i), ".ma"), "\n")
    sink()
    system(paste0("bash ", basename(bfile), "_PT_col", pheno.col, "_val", i, ".sh"))

    # clumping
    system(
      paste0(
        "plink --bfile ", bfile,
        " --clump ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i, ".ma"),
        " --clump-field P",
        " --clump-p1 1",
        " --clump-p2 1",
        " --clump-r2 ", r2,
        " --clump-kb 1000",
        " --allow-no-sex --out ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i)
      )
    )

    system(paste("awk 'NR!=1{print $3}' ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i), ".clumped > ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i), ".snplist", sep=""))

    # predicting for thresholds
    system(
      paste0(
        "plink --bfile ", bfile,
        " --score ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i, ".ma 1 2 5 header sum"),
        " --q-score-range ", myrfile,
        " ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i, ".ma 1 7 header"),
        " --extract ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i, ".snplist"),
        " --allow-no-sex --out ",
        paste0(basename(bfile), "_PT_col", pheno.col, "_val", i)
      )
    )

    if(i <= val.num){
      for(j in 1 : n){
        infebv <- try(read.table(paste0(paste0(basename(bfile), "_PT_col", pheno.col, "_val", i),".s",j,".profile"), head=T)[, 6], silent=TRUE)
        if(class(infebv) != "try-error"){
          Acc[i, j] <- pROC::auc(pheno[index[val.fold == i]] - 3, infebv[index[val.fold == i]])
        }
      }
    }else{
      predict <- read.table(paste0(paste0(basename(bfile), "_PT_col", pheno.col, "_val", i),".s1.profile"), head=T)[, 6]
    }

    if(i == val.num){
      n1 <- which.max(apply(Acc, 2, mean))
      system(paste0("sed -n ", n1, "p ", rfile, " > ", basename(bfile), "_PT_col", pheno.col, "_val", i + 1, ".opt.range"))
      myrfile <- paste0(basename(bfile), "_PT_col", pheno.col, "_val", i + 1, ".opt.range")
    }

    system(paste0("rm -f ", paste0(basename(bfile), "_PT_col", pheno.col, "_val", i, "*")))
  }

  return(predict)
}


BayesR <- 
function(
  bfile="", n=1, seed=123456, file.remove=TRUE, memo='BayesR'
)
{
  cat("------------BayesR started------------\n")
  phe <- read.delim(paste(bfile, ".fam", sep=""), head=FALSE, sep=" ")
  map <- read.delim(paste(bfile, ".bim", sep=""), head=FALSE, sep="\t")
  indexNA <- which(is.na(phe[, 5+n]))
  cat(paste("Number of Individuals:", nrow(phe), "\n"))
  cat(paste("Number of References:", nrow(phe)-length(indexNA), "\n"))
  cat(paste("Number of Candidates:", length(indexNA), "\n"))
  cat(paste("Number of SNPs:", nrow(map), "\n"))
  rm(phe);rm(map);gc()
  
  cat("Running...\n")
  
  t1=as.numeric(Sys.time())
  
  system(
    paste(
      "bayesR -bfile", bfile,
      "-n", n,
      "-out", memo,
      "-seed", seed
    )
  )
  
  system(
    paste(
      "bayesR -bfile", bfile,
      "-n", n,
      "-out", paste(memo, "pred", sep="_"),
      "-predict -model", paste(memo, ".model", sep=""),
      "-freq", paste(memo, ".frq", sep=""),
      "-param", paste(memo, ".param", sep="")
    )
  )
  
  t2=as.numeric(Sys.time())
  time <- t2 - t1
  
  refebv <- as.matrix(read.table(paste(memo, ".gv", sep=""), head=FALSE)[, 1])
  infebv <- as.matrix(read.delim(paste(paste(memo, "pred", sep="_"), ".gv", sep=""), head=FALSE)[, 1])
  if(length(indexNA) == 0){
    refebv <- refebv
    infebv <- NULL
  }else{
    refebv <- refebv[-indexNA, ]
    infebv <- infebv[indexNA, ]
  }
  if(file.remove){
    files <- paste(memo, c(".model", ".frq", ".param", ".hyp", ".gv", ".log", "_pred.gv", "_pred.log"), sep="")
    unlink(files)
  }
  cat("-----------------DONE-----------------\n")
  return(list(refebv=refebv, infebv=infebv, time=time))
}


BSLMM <- 
function(
  bfile="", n=1, seed=123456, file.remove=TRUE, memo='BSLMM'
)
{
  cat("-------------BSLMM started------------\n")
  phe <- read.delim(paste(bfile, ".fam", sep=""), head=FALSE, sep=" ")
  map <- read.delim(paste(bfile, ".bim", sep=""), head=FALSE, sep="\t")
  indexNA <- which(is.na(phe[, 5+n]))
  cat(paste("Number of Individuals:", nrow(phe), "\n"))
  cat(paste("Number of References:", nrow(phe)-length(indexNA), "\n"))
  cat(paste("Number of Candidates:", length(indexNA), "\n"))
  cat(paste("Number of SNPs:", nrow(map), "\n"))
  rm(phe);rm(map);gc()
  
  cat("Running...\n")
  
  t1=as.numeric(Sys.time())
  
  system(
    paste(
      "gemma -bfile", bfile,
      "-n", n,
      "-notsnp -bslmm 1 ",
      "-o", memo,
      "-seed", seed,
      ">", paste(memo, ".LOG", sep="")
    )
  )
  
  File <- read.delim(paste("./output/", memo, ".log.txt", sep=""), head=FALSE)
  cat(paste("Number of analyzed SNPs:", as.numeric(unlist(strsplit(as.character(File[14,]), " "))[7]), "\n"))
  
  if(length(indexNA) != 0){
    system(
      paste(
        "gemma -bfile ", bfile,
        " -n ", n,
        " -epm ./output/", memo, ".param.txt",
        " -emu ./output/", memo, ".log.txt",
        " -seed ", seed,
        " -predict 1 -o ",paste(memo ,"pred",sep="_"),
        " >>", paste(memo, ".LOG", sep="")
      ,sep="")
    )
  }
  t2=as.numeric(Sys.time())
  time <- t2 - t1
  
  refebv <- as.matrix(read.table(paste("./output/", memo, ".bv.txt", sep=""), head=FALSE)[, 1])

  if(length(indexNA) == 0){
    refebv <- refebv
    infebv <- NULL
  }else{
    infebv <- as.matrix(read.table(paste("./output/", paste(memo, "pred", sep="_"), ".prdt.txt", sep=""), head=FALSE)[, 1])
    refebv <- refebv[-indexNA, ]
    infebv <- infebv[indexNA, ]
  }
  if(file.remove){
    unlink(
      c(paste(memo, ".LOG", sep=""),
      paste("./output/",memo,
      c(".bv.txt",".gamma.txt",".hyp.txt",".log.txt",".param.txt",
      "_pred.log.txt","_pred.prdt.txt"),sep="")
      )
    )
  }
  cat("-----------------DONE-----------------\n")
  return(list(refebv=refebv, infebv=infebv, time=time))
}


DPR <- 
function(
  bfile="", n=1, seed=123456, file.remove=TRUE, memo='DPR'
)
{
  cat("-------------DPR started------------\n")
  phe <- read.delim(paste(bfile, ".fam", sep=""), head=FALSE, sep=" ")
  map <- read.delim(paste(bfile, ".bim", sep=""), head=FALSE, sep="\t")
  indexNA <- which(is.na(phe[, 5+n]))
  cat(paste("Number of Individuals:", nrow(phe), "\n"))
  cat(paste("Number of References:", nrow(phe)-length(indexNA), "\n"))
  cat(paste("Number of Candidates:", length(indexNA), "\n"))
  cat(paste("Number of SNPs:", nrow(map), "\n"))
  rm(phe);rm(map);gc()
  
  cat("Running...\n")
  
  t1=as.numeric(Sys.time())
  
  system(
    paste(
      "DPR -bfile", bfile,
      "-n", n,
      "-notsnp -dpr 3 ",
      "-o", memo,
      "-seed", seed,
      ">", paste(memo, ".LOG", sep="")
    )
  )
  
  File <- read.delim(paste("./output/", memo, ".log.txt", sep=""), head=FALSE)
  cat(paste("Number of analyzed SNPs:", as.numeric(unlist(strsplit(as.character(File[13,]), " "))[7]), "\n"))
  
  if(length(indexNA) != 0){
    system(
      paste(
        "DPR -bfile ", bfile,
        " -n ", n,
        " -epm ./output/", memo, ".param.txt",
        " -emu ./output/", memo, ".log.txt",
        " -seed ", seed,
        " -predict -o ", memo,
        " >>", paste(memo, ".LOG", sep="")
      ,sep="")
    )
  }
  t2=as.numeric(Sys.time())
  time <- t2 - t1
  
  #refebv <- as.matrix(read.table(paste("./output/", memo, ".prdt.txt", sep=""), head=FALSE)[, 1])

  if(length(indexNA) == 0){
    #refebv <- refebv
    infebv <- NULL
  }else{
    infebv <- as.matrix(read.table(paste("./output/", memo, ".prdt.txt", sep=""), head=FALSE)[, 1])
    #refebv <- refebv[-indexNA, ]
    infebv <- infebv[indexNA, ]
  }
  if(file.remove){
    unlink(
      c(paste(memo, ".LOG", sep=""),
      paste("./output/",memo,
      c(".log.txt",".param.txt",
      ".prdt.txt"),sep="")
      )
    )
  }
  cat("-----------------DONE-----------------\n")
  return(list(infebv=infebv, time=time))
}


MultiBLUPv5 <- 
function(
  bfile='', n=1, wind=75000, sig1="bonf", sig2=0.01, background=TRUE, file.remove=TRUE, memo='multiBLUP'
)
{

  cat("----------------------MultiBLUP started----------------------\n")
  phex <- read.delim(paste(bfile, ".fam", sep=""), head=FALSE, sep=" ")
  map <- read.delim(paste(bfile, ".bim", sep=""), head=FALSE, sep="\t")
  indexNA <- which(is.na(phex[, 5+n]))
  phe <- phex[!is.na(phex[, 5+n]), ]
  cat(paste("Number of Individuals:", nrow(phex), "\n"))
  cat(paste("Number of References:", nrow(phex)-length(indexNA), "\n"))
  cat(paste("Number of Candidates:", length(indexNA), "\n"))
  cat(paste("Number of SNPs:", nrow(map), "\n"))

  dir.create(memo)
  write.table(phe[,c(1,2,5+n)], paste("./", memo, "/", memo, ".pheno", sep=""), row.names=FALSE, col.names=FALSE, quote=FALSE, sep=" ")
  rm(phe);rm(map);gc()

  t1 <- as.numeric(Sys.time())
  
  #divide the genome into chunks. 
  cat(paste("Divide the genome into chunks(", wind/1000, "kb)...", sep=""))
  multiBLUP.try <- try(
    system(
      paste(
        "ldak --cut-genes ", memo, 
        " --chunks-bp ", wind, 
        " --bfile ", bfile, 
        " --ignore-weights YES",
        " >", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
      sep="")
    ), silent=TRUE
  )
  cat("Done!\n")

  #test the chunks for association.
  cat("Test the chunks for association...")
  multiBLUP.try <- try(
    system(
        paste(
        "ldak --calc-genes-reml ", memo,
        " --bfile ", bfile,
        " --pheno ", paste(memo, "/", memo, ".pheno", sep=""),
        " --ignore-weights YES --partition 1 --power -0.25",
        " >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
      sep="")
    ), silent=TRUE
  )
  cat("Done!\n")

  #identify significant chunks.
  if(is.null(sig1)){
    chunks.details <- read.delim(paste(memo, "/remls.1", sep=""), head=TRUE, sep=" ")
    sig1 <- min(c(0.05 / nrow(chunks.details), sort(chunks.details[, 8])[20]))
  }else{
    if(sig1=="bonf"){
      chunks.details <- read.delim(paste(memo, "/remls.1", sep=""), head=TRUE, sep=" ")
      sig1 <- 0.05 / nrow(chunks.details)
      rm(chunks.details); gc()
    }
  }
  cat(paste("Identify significant chunks(", sig1, " ", sig2, ")...", sep=""))
  multiBLUP.try <- try(
    system(
      paste(
        "ldak --join-genes-reml ", memo,
        " --bfile ", bfile,
        " --sig1 " , sig1,
        " --sig2 ", sig2,
        " >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
      sep="")
    ), silent=TRUE
  )
  cat("Done!\n")

  regionN <- read.delim(paste("./", memo, "/region.number", sep=""), head=FALSE)
  cat(paste("Number of regions:", as.numeric(regionN), "\n"))
  
  if(regionN == 0 | background){
    #calculate the background kinship matrix.
    cat("Calculate the background kinship matrix...")
    multiBLUP.try <- try(
      system(
        paste(
          "ldak --calc-kins-direct ", paste(memo, "/region.0", sep=""),
          " --bfile ", bfile,
          " --extract " , paste(memo, "/region.0", sep=""),
          " --ignore-weights YES --power -0.25",
          " >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
        sep="")
      ), silent=TRUE
    )
    cat("Done!\n")
  }
  
  #estimate variance components for the MultiBLUP Model.
  cat("Estimate variance components for the MultiBLUP Model...")
  if(as.numeric(regionN) != 0 & !background){
    multiBLUP.try <- try(
      system(
        paste(
          "ldak --reml ", paste(memo, "/REML", sep=""),
          " --bfile ", bfile,
          " --pheno ", paste(memo, "/", memo, ".pheno", sep=""),
          " --region-number " , as.numeric(regionN),
          " --region-prefix " , paste(memo, "/region.", sep=""),
          " --ignore-weights YES --power -0.25",
          " >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
        sep="")
      ), silent=TRUE
    )
  }else{
    multiBLUP.try <- try(
      system(
        paste(
          "ldak --reml ", paste(memo, "/REML", sep=""),
          " --bfile ", bfile,
          " --grm " , paste(memo, "/region.0", sep=""),
          " --pheno ", paste(memo, "/", memo, ".pheno", sep=""),
          " --region-number " , as.numeric(regionN),
          " --region-prefix " , paste(memo, "/region.", sep=""),
          " --ignore-weights YES --power -0.25",
          " >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
        sep="")
      ), silent=TRUE
    )
  }
  cat("Done!\n")

  #estimate SNP effect sizes and calculate predictions.
  cat("Estimate SNP effect sizes and calculate predictions...")
  if(as.numeric(regionN) != 0 & !background){
    multiBLUP.try <- try(
      system(
        paste(
          "ldak --calc-blups ", paste(memo, "/PRED", sep=""),
          " --bfile ", bfile,
          " --remlfile ", paste(memo, "/REML.reml", sep=""),
          " >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
        sep="")
      ), silent=TRUE
    )
  }else{
    multiBLUP.try <- try(
      system(
        paste(
          "ldak --calc-blups ", paste(memo, "/PRED", sep=""),
          " --bfile ", bfile,
          " --grm " , paste(memo, "/region.0", sep=""),
          " --remlfile ", paste(memo, "/REML.reml", sep=""),
          " >>", paste(memo, "/", memo, "_multiBLUP.log", sep=""),
        sep="")
      ), silent=TRUE
    )
  }
  cat("Done!\n")
  
  t2 <- as.numeric(Sys.time())
  time <- t2 - t1
  
  cat("----------------------------DONE-----------------------------\n")
  if(class(multiBLUP.try) == "try-error"){
    return("error!")
  }else{
    gebv <- read.delim(paste(memo,"/PRED.pred",sep=""), head=TRUE)
    gebv <- gebv[, ncol(gebv), drop=FALSE]
    if(length(indexNA) == 0){
      refebv <- gebv
      infebv <- NULL
    }else{
      refebv <- gebv[-indexNA, ]
      infebv <- gebv[indexNA, ]
    }
    if(file.remove){
      unlink(memo, recursive=TRUE)
    }
    return(list(refebv=refebv, infebv=infebv, indexNA=indexNA, time=time))
  }
}


