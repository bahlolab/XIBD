#' XIBD Plot IBD Segments
#'
#' Plot IBD segments for pairs across the genome.
#' Annotation genes can be added to the plot and specific genes on interest can be highlighted.
#'
#' @param ibd.segments a named list containing the \code{ibd_segments} inferred from pairs of samples.
#' See \code{value} description in \code{\link{getIBDsegments}} for more details.
#' @param ped.genotypes a named list containing \code{pedigree}, \code{genotypes} and \code{model}.
#' See \code{Value} description in \code{\link{getGenotypes}} for more details.
#' The family IDs and individual IDs in \code{pedigree} must match the family IDs and individual IDs in the header of \code{genotypes}.
#' @param interval a vector of length 3 containing a chromosome, a start position (bp) and an end position (bp)
#' of an interval to plot. The default is \code{interval=NULL} which will plot the IBD segments over all chromosomes
#' in \code{ped.genotypes}.
#' @param annotation.genes a data frame with at least 5 columns of information:
#' \enumerate{
#' \item Chromosome (type \code{numeric} or \code{integer})
#' \item Gene name (type \code{character})
#' \item Start location of the gene in base-pairs (type \code{numeric} or \code{integer})
#' \item End location of the gene in base-pairs (type \code{numeric} or \code{integer})
#' \item Gene strand (+ or -)
#' }
#' \code{annotation.genes} should contain the following headers \code{chr, name, start, end} and \code{strand}.
#' This data frame does not have to be in a specific order, however it must contain all of the above information
#' with respective labels.
#' @param highlight.genes a data frame with at least 4 columns of information:
#' \enumerate{
#' \item Chromosome (type \code{numeric} or \code{integer})
#' \item Gene name (type \code{character})
#' \item Start location of the gene in base-pairs (type \code{numeric} or \code{integer})
#' \item End location of the gene in base-pairs (type \code{numeric} or \code{integer})
#' }
#' \code{highlight.genes} should contain the following headers \code{chr, name, start} and \code{end}.
#' This data frame does not have to be in a specific order, however it must contain all of the above information
#' with respective labels.
#' @param segment.hight the hight of IBD segment blocks, such that 0 < \code{segment.hight} <= 1. The default is \code{segment.hight=0.5}.
#' @param number.per.page the maximum number of IBD pairs to plot in a single praphics window. The default is
#' \code{number.per.page=NULL} which will plot all IBD pairs in a single window. This may not be ideal when there are many IBD pairs. If
#' \code{number.per.page} is set, it is recommended to plot the output to a file as opposed to the R console.
#' @param add.fid.name a logical value indicating if family IDs should be included in the y-axis labels.
#' The default is \code{add.fid.name=TRUE}.
#' @param add.iid.name a logical value indicating if individual IDs should be included in the y-axis labels.
#' The default is \code{add.iid.name=TRUE}.
#' @param add.rug a logical value indicating whether SNP positions should be added to the plot. The default is \code{add.rug=TRUE}.
#' @param plot.title a character string of a title to be added to the plot. The default is \code{plot.title=NULL} which does not add a title to the plot.
#' @param add.legend a logical value indicating whether the IBD status legend should be plotted. The default is \code{add.legend=TRUE}.
#' @importFrom ggplot2 ggplot
#' @export
plotIBDsegments <- function (ibd.segments, ped.genotypes, interval = NULL,
                             annotation.genes = NULL, highlight.genes = NULL, segment.hight = 0.5,
                             number.per.page = NULL, add.fid.name = TRUE, add.iid.name = TRUE, add.rug = FALSE,
                             plot.title = NULL, add.legend = TRUE) {

  # check ibd.segments file is a dataframe with correct fields
  stopifnot(is.list(ibd.segments))
  stopifnot("ibd_segments" %in% names(ibd.segments))
  stopifnot(is.data.frame(ibd.segments[["ibd_segments"]]))
  ibd.segments <- ibd.segments[["ibd_segments"]]
  if (ncol(ibd.segments) != 15)
    stop ("ibd.segments has incorrect format")
  colnames(ibd.segments) <- c("fid1", "ind1", "fid2", "ind2", "chr", "start.snp", "end.snp", "start.position.bp",
                              "end.position.bp", "start.position.M", "end.position.M", "number.snps", "length.bp",
                              "length.M", "ibd.status")

  # check format of input data
  stopifnot(is.list(ped.genotypes) | length(ped.genotypes) == 3)
  stopifnot(c("pedigree","genotypes","model") %in% names(ped.genotypes))
  pedigree  <- ped.genotypes[["pedigree"]]
  genotypes <- ped.genotypes[["genotypes"]]
  if(!all(colnames(genotypes)[1:5] %in% c("chr", "snp_id", "pos_M","pos_bp", "freq")))
    stop ("ped.genotypes has incorrect format")

  # check groups
  #if (!is.null(groups)) {
  #  stopifnot(is.data.frame(groups))
  #  stopifnot(ncol(groups) > 2)
  #  if (ncol(groups) > 4){
  #    cat("using first 4 columns of groups")
  #    groups <- groups[,1:4]
  #  }
  #  colnames(groups)[1:2] <- c("fid","iid")

    # check isolates belong to a group
  #  group.names <- paste(groups[,"fid"],groups[,"iid"],sep="/")
  #  isolate.names <- paste(pedigree[,"fid"],pedigree[,"iid"],sep="/")
  #  if (!all(isolate.names %in% group.names))
  #    stop("'groups' is missing information for some isoaltes")
  #}

  # if interval specified
  if (!is.null(interval) & length(interval) == 3) {
    chr   <- as.character(interval[1])
    start <- as.numeric(interval[2])
    stop  <- as.numeric(interval[3])
    stopifnot(chr %in% ibd.segments[,"chr"])
    if(!(chr %in% genotypes[,"chr"]))
      stop(paste0("chromosome ",chr," not in 'ped.genotypes'\n"))
    if(start > stop)
      stop(paste("interval start=",interval[2]," is greater than interval end=",interval[3],sep=""))

    # subset ibd by interval
    ibd.segments.0 <- ibd.segments[ibd.segments[,"chr"] == chr,]
    ibd.segments.1 <- NULL
    if(nrow(ibd.segments.0) > 0) {
      # function to find genes that overlap interval
      fun.overlap <- function(region.1, region.2) {
        a <- max(region.1[1], region.2[1])
        b <- min(region.1[2], region.2[2])
        return(c(a,b))
      }
      for(g in 1:nrow(ibd.segments.0)){
        gene.overlap <- fun.overlap(ibd.segments.0[g,c("start.position.bp","end.position.bp")],c(start,stop))
        if (gene.overlap[2] - gene.overlap[1] > 0)
          ibd.segments.1 <- rbind(ibd.segments.1, ibd.segments.0[g,])
      }
    } else {
      stop("no ibd segments in interval")
    }
    ibd.segments.1 <- data.frame(ibd.segments.1)
    if(nrow(ibd.segments.1) == 0) {
      stop("no ibd segments in interval\n")
    } else {
      ibd.segments.1[,"start.position.bp"] <- as.numeric(ibd.segments.1[,"start.position.bp"])
      ibd.segments.1[,"end.position.bp"] <- as.numeric(ibd.segments.1[,"end.position.bp"])
    }
  } else {
    ibd.segments.1 <- ibd.segments
  }

  # check annotation.genes
  if (!is.null(annotation.genes)) {
    stopifnot(is.data.frame(annotation.genes))
    stopifnot(c("name","strand","chr","start","end") %in% colnames(annotation.genes))

    # subset genes if interval specified
    if (!is.null(interval)) {
      annotation.genes.0 <- annotation.genes[annotation.genes[,"chr"] == chr,]
      annotation.genes.1 <- NULL
      if(nrow(annotation.genes.0) > 0) {
        # function to find genes that overlap interval
        fun.overlap <- function(region.1, region.2) {
          a <- max(region.1[1], region.2[1])
          b <- min(region.1[2], region.2[2])
          return(c(a,b))
        }
        for(g in 1:nrow(annotation.genes.0)){
          gene.overlap <- fun.overlap(annotation.genes.0[g,c("start","end")],c(start,stop))
          if (gene.overlap[2] - gene.overlap[1] > 0)
            annotation.genes.1 <- rbind(annotation.genes.1, annotation.genes.0[g,])
        }
      }
      annotation.genes.1 <- data.frame(annotation.genes.1)
    } else {
      annotation.genes.1 <- annotation.genes
    }
    if(nrow(annotation.genes.1) == 0) {
      cat("no annotation.genes in interval\n")
    } else {
      annotation.genes.1[,"start"] <- as.numeric(annotation.genes.1[,"start"])
      annotation.genes.1[,"end"] <- as.numeric(annotation.genes.1[,"end"])
    }
  }

  # check highlight.genes
  if (!is.null(highlight.genes)) {
    stopifnot(is.data.frame(highlight.genes))
    stopifnot(c("name","chr","start","end") %in% colnames(highlight.genes))

    # subset genes if interval specified
    if (!is.null(interval)) {
      highlight.genes.0 <- highlight.genes[highlight.genes[,"chr"] == chr,]
      highlight.genes.1 <- NULL
      if(nrow(highlight.genes.0) > 0) {
        # function to find genes that overlap interval
        fun.overlap <- function(region.1, region.2) {
          a <- max(region.1[1], region.2[1])
          b <- min(region.1[2], region.2[2])
          return(c(a,b))
        }
        for(g in 1:nrow(highlight.genes.0)){
          gene.overlap <- fun.overlap(highlight.genes.0[g,c("start","end")],c(start,stop))
          if (gene.overlap[2] - gene.overlap[1] > 0)
            highlight.genes.1 <- rbind(highlight.genes.1, highlight.genes.0[g,])
        }
      }
      highlight.genes.1 <- data.frame(highlight.genes.1)
    } else {
      highlight.genes.1 <- highlight.genes
    }
    if(nrow(highlight.genes.1) == 0) {
      cat("no highlight.genes in interval\n")
    } else {
      highlight.genes.1[,"start"] <- as.numeric(highlight.genes.1[,"start"])
      highlight.genes.1[,"end"] <- as.numeric(highlight.genes.1[,"end"])
    }
  }

  # check segment hight
  stopifnot(is.numeric(segment.hight))
  if (length(segment.hight) > 1)
    segment.hight <- segment.hight[1]

  # check segment hight
  if (!is.null(number.per.page)) {
    stopifnot(is.numeric(number.per.page))
    if (length(number.per.page) > 1)
      number.per.page <- number.per.page[1]
  }

  # check logical
  stopifnot(is.logical(add.fid.name))
  stopifnot(is.logical(add.iid.name))
  stopifnot(is.logical(add.rug))
  stopifnot(is.logical(add.legend))

  # check title
  if (!is.null(plot.title)) {
    stopifnot(is.vector(plot.title))
    plot.title <- as.character(plot.title)
  }

  # get chromosomes
  chromosomes <- unique(genotypes[,"chr"])

  # define continious plot positions if interval not specified
  if (is.null(interval)) {
    chrstart <- 0
    chradd   <- 0
    labelpos <- NULL
    newpos   <- NULL
    for (i in 1:length(chromosomes)) {
      maxpos   <- max(genotypes[genotypes[,"chr"] == chromosomes[i],"pos_bp"])
      minpos   <- min(genotypes[genotypes[,"chr"] == chromosomes[i],"pos_bp"])
      newpos   <- c(newpos, genotypes[genotypes[,"chr"] == chromosomes[i],"pos_bp"] + chrstart)
      ibd.segments.1[ibd.segments.1[,"chr"] == chromosomes[i],"start.position.bp"] <- ibd.segments.1[ibd.segments.1[,"chr"] == chromosomes[i],"start.position.bp"] + chrstart
      ibd.segments.1[ibd.segments.1[,"chr"] == chromosomes[i],"end.position.bp"] <- ibd.segments.1[ibd.segments.1[,"chr"] == chromosomes[i],"end.position.bp"] + chrstart
      labelpos[i] <- (maxpos - minpos + 2*chrstart)/2
      chradd[i]   <- chrstart
      if (!is.null(annotation.genes)) {
        if (nrow(annotation.genes.1) != 0){
          annotation.genes.1[annotation.genes.1[,"chr"] == chromosomes[i],"start"] <- annotation.genes.1[annotation.genes.1[,"chr"] == chromosomes[i],"start"] + chrstart
          annotation.genes.1[annotation.genes.1[,"chr"] == chromosomes[i],"end"] <- annotation.genes.1[annotation.genes.1[,"chr"] == chromosomes[i],"end"] + chrstart
        }
      }
      if (!is.null(highlight.genes)) {
        if (nrow(highlight.genes.1) != 0){
          highlight.genes.1[highlight.genes.1[,"chr"] == chromosomes[i],"start"] <- highlight.genes.1[highlight.genes.1[,"chr"] == chromosomes[i],"start"] + chrstart
          highlight.genes.1[highlight.genes.1[,"chr"] == chromosomes[i],"end"] <- highlight.genes.1[highlight.genes.1[,"chr"] == chromosomes[i],"end"] + chrstart
        }
      }
      chrstart <- chrstart + maxpos
    }
    chradd <- chradd[-1]
  }

  # get group pairs for each ibd segment
  #if (!is.null(groups)) {
  #  for (i in 1:nrow(groups)) {
  #    ibd.segments.1[ibd.segments.1[,"fid1"] %in% groups[i,1] & ibd.segments.1[,"ind1"] %in% groups[i,2],"group1"] <- groups[i,3]
  #    ibd.segments.1[ibd.segments.1[,"fid2"] %in% groups[i,1] & ibd.segments.1[,"ind2"] %in% groups[i,2],"group2"] <- groups[i,3]
  #  }
  #  ibd.segments.1$group_pairs <- paste(ibd.segments.1[,"group1"],ibd.segments.1[,"group2"],sep="/")
  #} else {
  #  ibd.segments.1$group_pairs <- "all"
  #}
  #number.groups <- length(unique(ibd.segments.1$group_pairs))

  # for each unique group pair, get numeric values for pairs
  ibd.segments.1$unique.pairs.1 <- paste(ibd.segments.1[,"fid1"],ibd.segments.1[,"ind1"],ibd.segments.1[,"fid2"],ibd.segments.1[,"ind2"],sep="/")
  ibd.segments.2 <- NULL
  #for (i in sort(unique(ibd.segments.1[,"group_pairs"]),decreasing=TRUE)) {
  #  ibd.segments.i <- ibd.segments.1[ibd.segments.1[,"group_pairs"] == i,]
  #  unique.pairs.i <- unique(paste(ibd.segments.i[,"fid1"],ibd.segments.i[,"ind1"],ibd.segments.i[,"fid2"],ibd.segments.i[,"ind2"],sep="/"))
  unique.pairs.i <- unique(ibd.segments.1$unique.pairs.1)
  num.ID <- 1:length(unique.pairs.i)
  numID  <- data.frame(unique.pairs.i, num.ID)
  ibd.segments.2 <- merge(ibd.segments.1,numID,by.x="unique.pairs.1",by.y="unique.pairs.i")
  #ibd.segments.2 <- rbind(ibd.segments.2, merge(ibd.segments.1,numID,by.x="unique.pairs.1",by.y="unique.pairs.i"))
  #}

  # set number per page
  if (!is.null(number.per.page)) {
    page.start <- 1
    for (i in seq(1,length(unique(ibd.segments.2[,"unique.pairs.1"])),number.per.page)) {
      if (i == max(seq(1,length(unique(ibd.segments.2[,"unique.pairs.1"])),number.per.page))) {
        page.pairs <- unique(ibd.segments.2[,"unique.pairs.1"])[page.start:length(unique(ibd.segments.2[,"unique.pairs.1"]))]
      } else
        page.pairs <- unique(ibd.segments.2[,"unique.pairs.1"])[page.start:(page.start+number.per.page-1)]
      ibd.segments.2[ibd.segments.2[,"unique.pairs.1"] %in% page.pairs,"page.num"] <- i
      for (j in 1:length(page.pairs)) {
        ibd.segments.2[ibd.segments.2[,"unique.pairs.1"] == page.pairs[j],"num.ID"] <- j
      }
      page.start <- page.start + number.per.page
    }
  } else {
    ibd.segments.2[,"page.num"] <- 1
  }


  # plotting segments:

  for (i in unique(ibd.segments.2[,"page.num"])) {
    # base plot
    p <- ggplot2::ggplot()
    if (add.legend) {
      p <- p + ggplot2::geom_rect(data=ibd.segments.2[ibd.segments.2[,"page.num"] == i,],
                                  ggplot2::aes(xmin = start.position.bp, xmax = end.position.bp, ymin = num.ID, ymax = num.ID+segment.hight,
                                               fill = ifelse(ibd.status == 1,"#69B4FF","#99DD55")), alpha=0.8)
      p <- p + ggplot2::scale_fill_manual("", values = c("#69B4FF","#99DD55"), labels=c("IBD = 1", "IBD = 2"))
    } else {
      p <- p + ggplot2::geom_rect(data=ibd.segments.2[ibd.segments.2[,"page.num"] == i,],
                                  ggplot2::aes(xmin = start.position.bp, xmax = end.position.bp, ymin = num.ID, ymax = num.ID+segment.hight),
                                  fill = ifelse(ibd.segments.2[ibd.segments.2[,"page.num"] == i,"ibd.status"] == 1,"#69B4FF","#99DD55"), alpha=0.8)
    }
    p <- p + ggplot2::theme_bw()
    p <- p + ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                            panel.grid.major = ggplot2::element_blank(),
                            strip.background = ggplot2::element_rect(fill = "white",color = "white"),
                            legend.title = ggplot2::element_blank(),
                            strip.text.y = ggplot2::element_text(angle = 360),
                            axis.title.y = ggplot2::element_blank())

    # add facet
    #p <- p + ggplot2::facet_grid(group_pairs~.)

    # remove legend
    #if (number.groups == 1)
    #p <- p + ggplot2::theme(legend.position = "none")

    # add rug
    if (add.rug) {
      if (is.null(interval)) {
        p <- p + ggplot2::geom_rug(ggplot2::aes(x = newpos), size = 0.1, colour = "gray30")
      } else
        p <- p + ggplot2::geom_rug(data=genotypes[genotypes[,"chr"] == chr,], ggplot2::aes(x = pos_bp), size = 0.1, colour = "gray30")
    }

    # add title
    if (!is.null(plot.title))
      p <- p + ggplot2::ggtitle(plot.title)

    # add highlighted genes
    if (!is.null(highlight.genes)) {
      if (nrow(highlight.genes.1) > 0) {
        p <- p + ggplot2::geom_rect(data=highlight.genes.1, ggplot2::aes(xmin = start, xmax = end), ymin = -Inf, ymax = Inf, fill = "gray40",alpha = 0.1)
        p <- p + ggplot2::geom_vline(data=highlight.genes.1, ggplot2::aes(xintercept = start), colour = "gray40", linetype = "solid", alpha = 0.1)
        x.pos <- rowMeans(highlight.genes.1[,c("start","end")])
        y.pos <- max(ibd.segments.2[ibd.segments.2[,"page.num"] == i,"num.ID"]) + segment.hight + 0.1
        p <- p + ggplot2::geom_text(data=highlight.genes.1, ggplot2::aes(x = x.pos, y = y.pos, label = name),
                                    colour = "gray20", angle = 0, hjust = 0.5, vjust = 0, size = 3, alpha = 0.8)
      }
    }

    # add annotation genes
    if (!is.null(annotation.genes)) {
      if (nrow(annotation.genes.1) > 0) {
        min.y <- 0.1*max(ibd.segments.2[ibd.segments.2[,"page.num"] == i,"num.ID"])
        y.limits <- c((0.5 - min.y), max(ibd.segments.2[ibd.segments.2[,"page.num"] == i,"num.ID"]) + segment.hight + 0.1)
        p <- p + ggplot2::geom_rect(data=annotation.genes.1, ggplot2::aes(xmin = start, xmax = end), ymin = 0.5, ymax = (0.5 - min.y), alpha = 0.9, fill = ifelse(annotation.genes.1[,"strand"]=="+","#FFE455","#FF7575"))
      }
    }

    # change x axis labels
    if (is.null(interval)) {
      p <- p + ggplot2::xlab("Chromosome")
      p <- p + ggplot2::geom_vline(xintercept = chradd, colour = "gray87", linetype = "longdash")
      p <- p + ggplot2::scale_x_continuous(breaks = labelpos, labels = chromosomes)
      p <- p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5))
    } else {
      p <- p + ggplot2::xlab(paste("Chromosome", chr))
      p <- p + ggplot2::coord_cartesian(xlim = c(start, stop))
    }

    # change y axis labels
    if (add.fid.name & add.iid.name){
      y.labs <- ibd.segments.2[ibd.segments.2[,"page.num"] == i,]
      y.labs <- y.labs[!duplicated(y.labs[,"num.ID"]),]
      y.breaks <- y.labs[,"num.ID"] + segment.hight/2
      y.labels <- paste(y.labs[,"fid1"], y.labs[,"ind1"], y.labs[,"fid2"], y.labs[,"ind2"],sep="/")
      if (!is.null(annotation.genes)) {
        if (nrow(annotation.genes.1) > 0) {
          p <- p + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels, limits=y.limits)
        } else
          p <- p + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels)
      } else
        p <- p + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels)
    }
    if (add.fid.name & !add.iid.name){
      y.labs <- ibd.segments.2[ibd.segments.2[,"page.num"] == i,]
      y.labs <- y.labs[!duplicated(y.labs[,"num.ID"]),]
      y.breaks <- y.labs[,"num.ID"] + segment.hight/2
      y.labels <- paste(y.labs[,"fid1"], y.labs[,"fid2"],sep="/")
      if (!is.null(annotation.genes)) {
        if (nrow(annotation.genes.1) > 0) {
          p <- p + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels, limits=y.limits)
        } else
          p <- p + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels)
      } else
        p <- p + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels)
    }
    if (!add.fid.name & add.iid.name) {
      y.labs <- ibd.segments.2[ibd.segments.2[,"page.num"] == i,]
      y.labs <- y.labs[!duplicated(y.labs[,"num.ID"]),]
      y.breaks <- y.labs[,"num.ID"] + segment.hight/2
      y.labels <- paste(y.labs[,"ind1"], y.labs[,"ind2"],sep="/")
      if (!is.null(annotation.genes)) {
        if (nrow(annotation.genes.1) > 0) {
          p <- p + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels, limits=y.limits)
        } else
          p <- p + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels)
      } else
        p <- p + ggplot2::scale_y_continuous(breaks = y.breaks, labels = y.labels)
    }
    if (!add.fid.name & !add.iid.name) {
      if (!is.null(annotation.genes)) {
        if (nrow(annotation.genes.1) > 0)
          p <- p + ggplot2::ylim(y.limits[1], y.limits[2])
        p <- p + ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                                axis.ticks=ggplot2::element_blank())
      } else
        p <- p + ggplot2::theme(axis.text.y=ggplot2::element_blank(),
                                axis.ticks=ggplot2::element_blank())
    }

    print(p)
  }

}
