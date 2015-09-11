steponeR <- function(files=NULL, target.ratios=NULL, fluor.norm=NULL,
                     copy.number=NULL, ploidy=NULL, extract=NULL) {
  require(plyr); require(reshape2)
  if(is.null(files)) stop("No data files specified")
#   # Import data from files
#   data <- do.call("rbind", lapply(files, function(file) {
#     data.frame(Filename=file, read.csv(file, skip=7, na.strings="Undetermined",
#                                        blank.lines.skip=T))
#   }))
  # New import strategy
  data <- lapply(files, function(x) {
    temp <- readLines(x)
    linesToSkip <- grep("^Well", temp)-1
    read.csv(text = temp, skip = linesToSkip)
  })
  data <- do.call("rbind", data)
  # Change C_ to CT
  colnames(data) <- sub(x=colnames(data), pattern="C_", replacement="CT")
  # Subset CT and sample metadata
  data <- data[, c("Filename", "Well", "Sample.Name", "Target.Name", "Task", "CT")]
  # Check and remove NTC wells
  ntc <- data[which(data$Task=="NTC"), ]
  if(any(!is.na(ntc$CT))) warning("Template detected in NTC: interpret data with caution")
  if(!empty(ntc)) data <- data[!rownames(data) %in% rownames(ntc), ]
  # Remove wells with no target
#   nosample <- data[which(data$Sample.Name==""), ]
#   if(!empty(nosample)) {
#     apply(nosample, 1, function(x) message(paste("Well", x["Well"], "in", x["Filename"], "discarded: no sample")))
#     data <- data[!rownames(data) %in% rownames(nosample), ]
#   }
  notarget <- data[which(data$Target.Name==""), ]
  if(!empty(notarget)) {
    apply(notarget, 1, function(x) message(paste("Well", x["Well"], "in", x["Filename"], "discarded: no target")))
    data <- data[!rownames(data) %in% rownames(notarget), ]
  }
  # Drop unused levels from data
  data <- droplevels(data)
  # Create unique sample-plate IDs to distinguish samples run on multiple plates
  data$Sample.Plate <- interaction(data$Sample.Name, data$Filename, sep="~")
  

  # Calculate mean and sd of technical replicates for each target for each sample run
  ctmeans <- dcast(data, Sample.Plate ~ Target.Name, mean, na.rm=F, value.var="CT")  # na.rm=F
  colnames(ctmeans) <- c(colnames(ctmeans)[1], paste(colnames(ctmeans)[-1], "CT.mean", sep="."))
  ctsds <- dcast(data, Sample.Plate ~ Target.Name, sd, na.rm=F, value.var="CT")
  colnames(ctsds) <- c(colnames(ctsds)[1], paste(colnames(ctsds)[-1], "CT.sd", sep="."))
  # Combine CT means, SDs
  result <- merge(ctmeans, ctsds)
  # Split Sample.Plate column into Plate and Sample.Name columns
  result <- cbind(colsplit(as.character(result$Sample.Plate), pattern="~", names=c("Sample.Name", "File.Name")),
                  result[, -1])
  #return(result)
  

  # List targets present in data
  targets <- levels(data$Target.Name)
  
  # Calculations
  # Fluorescence normalization
  if(!is.null(fluor.norm)) {
    if(is.list(fluor.norm)) {
      if(any(!names(fluor.norm) %in% targets)) {
        warning(paste(names(fluor.norm)[which(!names(fluor.norm) %in% targets)], "not a valid Target"))
      }
      for (fluor in names(fluor.norm)) {
        result[, paste(fluor, "CT.mean", sep=".")] <- result[, paste(fluor, "CT.mean", sep=".")] - fluor.norm[[fluor]]
      }
    } else {
      stop("fluor.norm must be a list")
    }
  }
  # Target ratios
  if(!is.null(target.ratios)) {
    # Separate numerator and denominator of desired ratios by period
    ratios <- strsplit(target.ratios, split=".", fixed=T)
    # Check that all ratios have length two
    if(any(unlist(lapply(ratios, length))!=2)) {
      warning(paste(target.ratios[which(unlist(lapply(ratios, length))!=2)], 
                    "is not a valid target ratio and will be ignored\n"))
      target.ratios <- target.ratios[-which(unlist(lapply(ratios, length))!=2)]
    }
    # Calculate ratios
    for(ratio in target.ratios) {
      num <- strsplit(ratio, split=".", fixed=T)[[1]][1]
      denom <- strsplit(ratio, split=".", fixed=T)[[1]][2]
      result[, ratio] <- 2^(result[, paste(denom, "CT.mean", sep=".")] - result[, paste(num, "CT.mean", sep=".")])
    }
  }
  if(!is.null(copy.number)) {
    for(ratio in target.ratios) {
      num <- strsplit(ratio, split=".", fixed=T)[[1]][1]
      denom <- strsplit(ratio, split=".", fixed=T)[[1]][2]
      cnratio <- copy.number[[num]] / copy.number[[denom]]
      result[, ratio] <- result[, ratio] / cnratio
    }
  }
  if(!is.null(ploidy)) {
    for(ratio in target.ratios) {
      num <- strsplit(ratio, split=".", fixed=T)[[1]][1]
      denom <- strsplit(ratio, split=".", fixed=T)[[1]][2]
      pratio <- ploidy[[num]] / ploidy[[denom]]
      result[, ratio] <- result[, ratio] / pratio
    }
  }
  if(!is.null(extract)) {
    for(ratio in target.ratios) {
      num <- strsplit(ratio, split=".", fixed=T)[[1]][1]
      denom <- strsplit(ratio, split=".", fixed=T)[[1]][2]
      eeratio <- extract[[num]] / extract[[denom]]
      result[, ratio] <- result[, ratio] / eeratio
    }
  }
  return(result)
}



files=list("20150807_KBayRecov_Mcap_2_data.csv", "20150808_KBayRecov_Mcap_1_data.csv")

df <- steponeR(files=files, target.ratios=c("C.Mcap", "D.Mcap"), 
               fluor.norm=list(C=2.26827, D=0, Mcap=0.84815),
               copy.number=list(C=10, D=2, Mcap=1),
               ploidy=list(C=1, D=1, Mcap=2),
               extract=list(C=0.813, D=0.813, Mcap=0.982))
