#!/usr/bin/env Rscript

# just some defaults when loading tables
load_csv <- function(fp, ...) {
  read.csv(fp, check.names = FALSE, stringsAsFactors = FALSE, ...)
}

# make one table of all sequences and parse the IDs into separate fields
aggr_fastas <- function(fp_ighv, fp_ighj) {
  ighv <- dnar::read.fa(fp_ighv)
  ighj <- dnar::read.fa(fp_ighj)
  ighv$Segment <- "IGHV"
  ighj$Segment <- "IGHJ"
  combo <- rbind(ighv, ighj)
  combo$Family <- with(
    combo,
    ifelse(
      grepl("IGHV[0-9]-", name),
      sub("-.*$", "", name),
      ifelse(
        grepl("VH[0-9]\\.", name),
        sub("VH([0-9]).*", "IGHV\\1", name),
        ifelse(
          grepl("IGHJ", name),
          sub("[-\\*].*", "", name),
          ""))))
  combo$Allele <- sub("(\\*.*[0-9]+)[abcd]+$", "\\1", combo$name)
  combo$Gene <- sub("\\*.*$", "", combo$name)
  combo$GeneNumber <- sub(".*[-.]([0-9]*)[A-Z]*", "\\1", combo$Gene)
  combo$GeneLetter <- sub(".*[^A-Z]([A-Z]*)$", "\\1", combo$Gene)
  combo$GeneNumber[combo$Segment == "IGHJ"] <- ""
  combo$GeneLetter[combo$Segment == "IGHJ"] <- ""
  combo$GeneNumber <- as.integer(combo$GeneNumber)
  combo$Suffix <- sub(".*\\*[0-9]+", "", combo$name)
  combo$SeqInIMGT <- c("F", "T")[grepl("a", combo$Suffix) + 1]
  combo$SeqInKIGenome <- c("F", "T")[grepl("b", combo$Suffix) + 1]
  combo$SeqInThisPaper <- c("F", "T")[grepl("c", combo$Suffix) + 1]
  combo$SeqInKIIgM <- c("F", "T")[grepl("d", combo$Suffix) + 1]
  combo
}

# aggregate one of the paper's tables into the existing combo table of all
# sequences and attributes
aggr_table <- function(combo, tbl, group) {
  idx <- match(tbl$`Sequence ID`, combo$name)
  if (! "Group" %in% colnames(combo)) {
    combo$Group <- ""
  }
  combo$Group[idx] <- group
  for (colname in colnames(tbl)[-1]) {
    if (! colname %in% colnames(combo)) {
      combo[[colname]] <- ""
    }
    # don't overwrite any entries
    if (! identical(unique(combo[[colname]][idx]), "")) {
      stop(paste0("Values already defined for column ", colname))
    }
    combo[[colname]][idx] <- tbl[[colname]]
  }
  combo
}

# gather all the info from the CSV/FASTA files into one big data frame
make_combo_table <- function(dp_paper) {
  table_paths <- list.files(dp_paper, pattern = "*.csv", full.names = TRUE)
  names(table_paths) <- sub("\\.csv$", "", basename(table_paths))
  tables <- lapply(table_paths, load_csv)
  combo <- aggr_fastas(
    file.path(dp_paper, "tableS3_V.fasta"),
    file.path(dp_paper, "tableS3_J.fasta"))
  combo <- aggr_table(combo, tables$table1, "table1")
  combo <- aggr_table(combo, tables$table2, "table2")
  combo <- aggr_table(combo, tables$table3, "table3")
  combo <- aggr_table(combo, tables$tableS2, "tableS2")
  combo$Flag <- sub("novel_allele", "Novel allele", combo$Flag)
  combo$Flag <- sub("novel_gene", "Novel gene", combo$Flag)
  combo$Flag <- sub("known", "Known", combo$Flag)
  combo <- combo[
    with(combo,
         order(Segment != "IGHV", Family, GeneNumber, GeneLetter, name)), ]
  colnames(combo)[1:2] <- c("SeqID", "Seq")
  columns <- c(
    "SeqID", "Allele", "Gene", "Family", "Segment",
    "GeneNumber", "GeneLetter", "Suffix",
    "SeqInIMGT", "SeqInKIGenome", "SeqInThisPaper", "SeqInKIIgM",
    "Group", "Flag", "Nearest Corresponding Alleles", "Functionality",
    "Mismatch/Gap", "Supporting Individuals", "Supporting Samples",
    "Chromosome", "Position", "Nearest corresponding alleles", "Seq")
  if (! identical(sort(columns), sort(colnames(combo)))) {
    stop("unexpected column names")
  }
  combo <- combo[, columns]
  combo
}

main <- function() {
  combo <- make_combo_table("from-paper")
  write.csv(
    combo, "output/alleles.csv", quote = FALSE, na = "", row.names = FALSE)
}

if (any(grepl("^--file", commandArgs()))) {
  main()
}