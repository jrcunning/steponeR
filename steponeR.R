steponeR <- function(files=NULL, target.ratios=NULL, fluor.norm=NULL,
                     copy.number=NULL, ploidy=NULL, extract=NULL) {
  require(plyr); require(reshape2)
  if(is.null(files)) stop("No data files specified")
  # Import data from files
  data0 <- lapply(files, function(x) {
    temp <- suppressWarnings(readLines(x))
    linesToSkip <- grep("^Well", temp) - 1
    dat <- data.frame(Filename=basename(x),
                      read.csv(text=temp, skip=linesToSkip, na.strings="Undetermined"))
    dat <- dat[which(dat$Target.Name!=""), ]  # Omit empty wells
    dat$Sample.Name <- as.character(dat$Sample.Name)  # Convert sample names to character
    dat
  })
  data0 <- rbind.fill(data0)
  # Change C_ to CT
  colnames(data0) <- sub(x=colnames(data0), pattern="C_", replacement="CT")
  # Check and remove NTC wells
  ntc <- data0[which(data0$Task=="NTC"), ]
  if(any(!is.na(ntc$CT))) warning("Template detected in NTC: interpret data with caution")
  data <- droplevels(data0[!rownames(data0) %in% rownames(ntc), ])
  # Check remaining tasks
  tasks <- levels(data$Task)
  # Subset CT and sample metadata
  if("STANDARD" %in% tasks) {
    data <- data[, c("Filename", "Well", "Sample.Name", "Target.Name", "Task", "CT", "Quantity")]
  } else {
    data <- data[, c("Filename", "Well", "Sample.Name", "Target.Name", "Task", "CT")]
  }
  # Remove wells with no target
  notarget <- data[which(data$Target.Name==""), ]
  if(!empty(notarget)) {
    apply(notarget, 1, function(x) message(paste("Well", x["Well"], "in", x["Filename"], "discarded: no target")))
    data <- data[!rownames(data) %in% rownames(notarget), ]
  }
  # Drop unused levels from data
  data <- droplevels(data)
  # Create unique sample-plate IDs to distinguish samples run on multiple plates
  data$Sample.Plate <- interaction(data$Sample.Name, data$Filename, sep="~")
  # Separate UNKNOWNS from STANDARDS
  unk <- data[which(data$Task=="UNKNOWN"), ]
  # Process STANDARDS
  if("STANDARD" %in% tasks) {
    # Process STANDARDS
    std <- data[which(data$Task=="STANDARD"), ]
    std <- std[, c("Quantity", "CT", names(which(sapply(std[,c("Target.Name", "Filename")], function(x) nlevels(x) > 1))))]
    std.lm <- lm(log10(Quantity) ~ ., data=std)
    # Use standard curves to calculate quantities for unknowns
    unk$copies <- 10^predict(std.lm, newdata = unk[, c("CT", "Target.Name", "Filename")])
    copymeans <- dcast(unk, Sample.Plate ~ Target.Name, mean, na.rm=T, value.var="copies")
    colnames(copymeans) <- c(colnames(copymeans)[1], paste(colnames(copymeans)[-1], "copies", sep="."))
  }
  # Calculate mean and sd of technical replicates for each target for each sample run
  ctmeans <- dcast(unk, Sample.Plate ~ Target.Name, mean, na.rm=T, value.var="CT")  # na.rm=F
  colnames(ctmeans) <- c(colnames(ctmeans)[1], paste(colnames(ctmeans)[-1], "CT.mean", sep="."))
  ctsds <- dcast(unk, Sample.Plate ~ Target.Name, sd, na.rm=T, value.var="CT")
  colnames(ctsds) <- c(colnames(ctsds)[1], paste(colnames(ctsds)[-1], "CT.sd", sep="."))
  techreps <- dcast(unk, Sample.Plate ~ Target.Name, value.var="CT",
                    fun.aggregate = function(x) length(which(!is.na(x))))
  colnames(techreps) <- c(colnames(techreps)[1], paste(colnames(techreps)[-1], "reps", sep="."))
  
  # Combine CT means, SDs, techreps, quantities
  if(exists("copymeans")) {
    result <- join_all(list(ctmeans, ctsds, techreps, copymeans), by="Sample.Plate")
  } else {
    result <- join_all(list(ctmeans, ctsds, techreps), by="Sample.Plate")
  }
  # Split Sample.Plate column into Plate and Sample.Name columns
  result <- cbind(colsplit(as.character(result$Sample.Plate), pattern="~", names=c("Sample.Name", "File.Name")),
                  result[, -1])
  # List targets present in data
  targets <- levels(data$Target.Name)
  
  # Optional data adjustments
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
  if("STANDARD" %in% tasks) return(list(standards=list(data=std, lm=std.lm), unknowns=unk, result=result))
  else return(list(unknowns=unk, result=result))
}


# # # EXAMPLE USAGE ----------
# # 
# files=list("20150807_KBayRecov_Mcap_2_data.csv", "20150808_KBayRecov_Mcap_1_data.csv")
# # 
# df <- steponeR(files=files, target.ratios=c("C.Mcap", "D.Mcap"), 
#                fluor.norm=list(C=2.26827, D=0, Mcap=0.84815),
#                copy.number=list(C=10, D=2, Mcap=1),
#                ploidy=list(C=1, D=1, Mcap=2),
#                extract=list(C=0.813, D=0.813, Mcap=0.982))
# #
# # Raw data for unknowns
# df$unknowns
# # Target ratio results
# df$result


