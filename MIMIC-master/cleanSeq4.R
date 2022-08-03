## Function to read in MassArray ooutput .csv file and return Beta values

## Testing only (MB)
#filename <- "~/MIMIC/www/test_samples.csv"
#threshold <- 0.2

cleanSeq4 <- function(threshold=0.2, filename)
{
  ## Modified by MB 17th Feb 2015
  ## Check for comment row

  skip <- 0
  if (length(grep("iPLEX", readLines(filename)[1])) != 0) {
    # Comment row present set skip to 1
    skip <- 2
  }

  data <- read.csv(filename, skip=skip, header=T)

  ## Get and clean assay IDs if necessary
  assayID <- as.character(data$Assay.Id)
  assayID2 <- sapply(1:length(assayID), function(x) strsplit(assayID[x],split="\\_")[[1]][1])

  # Replace existing Assay.Id
  data$Assay.Id <- assayID2

  ## MB detect presence of col called Sample.Descripton and change to Sample
  # Define a funtion for use later
  CheckAndFixCol <- function(x, find.string = NULL, replace.string = NULL) {
    col_names <- colnames(x)
    grep_out <- grep(find.string, col_names)
    if (length(grep_out) != 0) {
      col_names[grep_out] <- replace.string
    }
    colnames(x) <- col_names
    return(x)
  }

  data <- CheckAndFixCol(data, "Sample.Description", "Sample")

  ## Check to see if Sample col is empty is so replace with content from Sample.Id col
  if (!sum(is.na(data[,2])) == 0) {
    cat("lacking Sample info in Sample col copying from Sample.Id")
    data[,2] <- data[,1]
  }


  ## Remove NTC (None Targeting Control) samples.
  # All lines which have NTC at the start of their value in the Sample col will be removed
  data <- data[grep("^NTC", data$Sample, invert = T), ]
  # Also remove PTC samples
  data <- data[grep("^PTC", data$Sample, invert = T), ]


  ## Skip this for production code since could have only 17, or 17 + BS efficiency control
  ## Remove any incomplete samples without 28 probes
  #x <- table(data$Sample)
  #data <- data[!data$Sample %in% names(x)[x<28 | x>28],]
  data$Sample <- factor(data$Sample)

  ## T formats - these are calculated as 1 - avgFreq1. (i.e. avgFreq2)
  Tform <- c("cg01561259","cg01986767","cg07262395",
             "cg08307469","cg17185060","cg18302652","cg24280645",
             "cg04541368","cg09923107","cg20912770","cg25923609")

  ## Known failures
  known <- c("cg00147216","cg15597109","cg19333614","cg08307469","cg01462184")


  ## For each sample, calculate yield
  ## MB detect presence of different col names and correct each one
  data <- CheckAndFixCol(data, "HEIGHT.P", "height_1")
  data <- CheckAndFixCol(data, "HEIGHT.Allele1", "height_2")
  data <- CheckAndFixCol(data, "HEIGHT.Allele2", "height_3")

  yield <- 1 - ((data$height_1) / (data$height_1 + data$height_2 + data$height_3))

  ## Set estimates below filter to NA
  yield[yield < threshold] <- NA


  output <- matrix(ncol=length(unique(data$Sample)), nrow=length(unique(assayID2)))
  colnames(output) <- 1:ncol(output)

  tmpNames <- unique(assayID2)

  rownames(output) <- tmpNames[order(tmpNames)]

  for(i in 1:length(unique(data$Sample)))
  {
    tmp <- cbind(data[data$Sample == unique(data$Sample)[i],], yield[data$Sample==unique(data$Sample)[i]])

    ## Order rownames in alphabetical order
    tmp <- tmp[order(tmp$Assay.Id),]

    ## Get peak height 1 over sum of heights 1 and 2 and replace with 1-avgFreq1 where appropriate
    output[,i] <- tmp$height_2 / ( tmp$height_2 + tmp$height_3)

    ## Replace with 1-avgFreq1
    to_replace <- which(tmp$Assay.Id %in% Tform)
    output[to_replace,i] <- (1-(tmp$height_2 / ( tmp$height_2 + tmp$height_3)))[to_replace]

    ## Check yield, and if below threshold, set to NA
    fail <- which(is.na(tmp[,ncol(tmp)]))
    if(length(fail) >=1) output[fail,i] <- NA

    ## Set known failed probes to NA
    to_replace <- which(tmp$Assay.Id %in% known)
    output[to_replace,i] <- NA
    colnames(output)[i] <- unique(as.character(data$Sample))[i]

  }

  ## How many failures are there beyond the 5 we know about
  cat("\n\n")
  cat("Number of failures:\n")
  print(apply(output,2, function(x) length(x[is.na(x)])-5))

  ## Remove Conv row, if present
  x <- which(rownames(output) == "Conv")
  if(length(x)==1) output <- output[-x,]

  ## Old File based output removed, we now return
  output2 <- output[order(rownames(output)),order(colnames(output))]
  #write.table(output2, file=outputName,sep=",")
  #return(output2)

  ## MB 6th of Jan 2016, code to extract BS conversion efficiency and return a
  ## list with existing data and new BS conversion efficiency data

  # Extract data for each sample on the Conv probe from df data
  BS_Eff <- data[data$Assay.Id == "Conv",c("Sample","Assay.Id","height_2","height_3")]
  BS_Eff_vect <- 100 - (BS_Eff$height_2/BS_Eff$height_3 * 100)
  BS_Eff[,"BS_Eff"] <- BS_Eff_vect
  # Set rownames
  rownames(BS_Eff) <- BS_Eff[,1]
  # Remove unwanted cols
  BS_Eff <- BS_Eff[,c(1,5)]

  # Make a unified list to return which has the old output (output2) and the new
  # BS_Eff data frame
  output3 <- list(output2, BS_Eff)
  return(output3)

}
