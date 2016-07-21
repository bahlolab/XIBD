#' XIBD BOT Genotype Switching
#'
#' The HapMap allele frequencies in XIBDs HapMap allele frequency files are calculated for the A allele only,
#' where the A allele is determined by the following rules:
#' \enumerate{
#' \item When one of the possible variations of the SNP is adenine (A), then adenine is labeled the A allele
#' and the remaining variation is labeled the B allele, regardless of what this might be.
#' \item If adenine (A) is not a variation of the SNP but cytosine (C) is, then cytosine is labeled the A allele
#' and the remaining variation is labeled the B allele.
#' \item If neither adenine (A) or cytosine (C) are variants of the SNP then thymine (T) is labeled the A allele.
#' }
#' Illuminas convention for the naming of A and B alleles differs to that of the HapMap data
#' (http://www.illumina.com/documents/products/technotes/technote_topbot.pdf). Rather, the classification
#' of A and B alleles depend on the top (TOP) and bottom (BOT) designations of the SNP. This
#' means that the A allele in the HapMap data is not always the same as the A allele in the Illumina data. In
#' fact, alleles that have been named according to the BOT designation actually correspond the the B allele
#' in the HapMap data. To correct for this, \code{switchBOTgenotypes()} switchs the A and B alleles in
#' the input genotypes for all SNPs corresponding to BOT designations. This mean a homozygous genotype, 0, will be
#' changed to a homozygous alternative genotype, 2, and vis versa. Heterozygous genotypes will be unchanged.
#' NOTE: this function should only be implemented with Illumina SNPchip data when XIBD's HapMap reference data is used
#' and if there is a noticeable discrepancy between population allele frequencies calculated from the HapMap reference data
#' and those calculated from the input dataset.
#' @param ped.genotypes a named list containing \code{pedigree}, \code{genotypes} and \code{model}.
#' See \code{Value} description in \code{\link{getGenotypes}} for more details.
#' The family IDs and individual IDs in \code{pedigree} must match the family IDs and individual IDs in the header of \code{genotypes}.
#' @param topbot.annotation the annotation file provided by XIBD containing information of TOPBOT status for each SNP
#' and population allele frequencies of the 11 HapMap populations.
#' This file can be downloaded from \url{http://bioinf.wehi.edu.au/software/XIBD}.
#' @return A named list of the same format as the input \code{ped.genotypes} with A and B alleles switched for BOT SNPs.
#' @export
switchBOTgenotypes <- function(ped.genotypes, topbot.annotation) {

  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 3)
  stopifnot(c("pedigree","genotypes","model") %in% names(ped.genotypes))
  pedigree  <- ped.genotypes[["pedigree"]]
  genotypes <- ped.genotypes[["genotypes"]]
  model     <- ped.genotypes[["model"]]
  stopifnot(is.numeric(model))
  stopifnot(model == 1 | model == 2)

  # check the pedigree has 6 coloumns
  if (ncol(pedigree) != 6)
    stop ("'ped.genotypes' has incorrect format")
  colnames(pedigree) <- c("fid", "iid", "pid", "mid", "sex", "aff")

  # check there are ped.genotypes and pairs to perform analysis
  if (model == 1){
    if(ncol(genotypes) < 8 & nrow(genotypes) <= 1)
      stop ("'ped.genotypes' has incorrect format")
    if(!all(colnames(genotypes)[1:5] %in% c("chr", "snp_id", "pos_M","pos_bp", "freq")))
      stop ("'ped.genotypes' has incorrect format")
    samples <- colnames(genotypes)[!(colnames(genotypes) %in% c("chr", "snp_id", "pos_M","pos_bp", "freq"))]
  }
  if (model == 2) {
    if(ncol(genotypes) < 14 & nrow(genotypes) <= 1)
      stop ("'ped.genotypes' has incorrect format")
    if(!all(colnames(genotypes)[1:11] %in% c("chr","snp_id","pos_M","pos_bp","freq","condition_snp", "pba", "pbA", "pBa", "pBA", "freq_condition_snp")))
      stop ("'ped.genotypes' has incorrect format")
    samples <- colnames(genotypes)[!(colnames(genotypes) %in% c("chr","snp_id","pos_M","pos_bp","freq","condition_snp", "pba", "pbA", "pBa", "pBA", "freq_condition_snp"))]
  }

  # check topbot annotation
  stopifnot(is.data.frame(topbot.annotation))
  if(ncol(topbot.annotation) != 18)
    stop("'topbot.annotation' has incorrect format")
  if(!all(colnames(topbot.annotation) %in% c("chr","snp_id","pos_bp","pos_M","TOPBOT","A","B","CEU","ASW","CHB","CHD","GIH","JPT","LWK","MEX","MKK","TSI","YRI")))
    stop ("'topbot.annotation' has incorrect format")

  # merge annotation file with genotype data
  genotype.topbot    <- merge(genotypes, topbot.annotation, by="snp_id")
  genotype.topbot.v1 <- genotype.topbot[order(genotype.topbot[,"chr.x"],genotype.topbot[,"pos_bp.x"]),]
  if (nrow(genotype.topbot.v1) != nrow(genotypes))
    stop("missing TOPBOT information for some SNPs")

  # get TOP/BOT SNPs
  topbot       <- as.character(genotype.topbot.v1[,"TOPBOT"])
  topbotMatrix <- as.matrix(genotype.topbot.v1[,samples])
  genotypesNew <- topbotChange(topbotMatrix, topbot)
  if (model == 1) {
    called.genotypes.v2 <- cbind(genotype.topbot.v1[,c("chr.x", "snp_id", "pos_M.x","pos_bp.x", "freq")],genotypesNew)
    colnames(called.genotypes.v2)[1:5] <- c("chr","snp_id","pos_M","pos_bp","freq")
  }
  if (model == 2) {
    called.genotypes.v2 <- cbind(genotype.topbot.v1[,c("chr.x", "snp_id", "pos_M.x", "pos_bp.x", "freq", "condition_snp", "pba", "pbA", "pBa", "pBA", "freq_condition_snp")],genotypesNew)
    colnames(called.genotypes.v2)[1:11] <- c("chr", "snp_id", "pos_M", "pos_bp", "freq", "condition_snp", "pba", "pbA", "pBa", "pBA", "freq_condition_snp")
  }

  return_list <- list(pedigree, called.genotypes.v2, model)
  names(return_list) <- c("pedigree","genotypes","model")
  return(return_list)
}
