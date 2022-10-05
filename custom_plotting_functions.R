gseaplot2_custom <- function(x, order,geneSetID, title = "", color="green", base_size = 11,
                             rel_heights=c(1.5, .5, 1), subplots = 1:3,
                             pvalue_table = FALSE, ES_geom="line") {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  
  geneList <- position <- NULL ## to satisfy codetool
  
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
    
  }
  gsdata$Description=factor(gsdata$Description,levels=order)
  gsdata=gsdata[order(gsdata$Description),]
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))
  
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                          size=1)
  } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description),
                           size=1.2, data = subset(gsdata, position == 1),show.legend = FALSE)
  }
  
  p.res <- p + es_layer +
    theme(legend.position = c(1,1), legend.title = element_blank(),
          legend.background = element_rect(fill = "transparent"),legend.text=element_text(size=12))+guides(color = guide_legend(byrow = TRUE))+
    guides(color = guide_legend(override.aes = list(size = 6)))
  
  p.res <- p.res + ylab("Running Enrichment Score") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),axis.title.y=element_text(size=15))
  
  i <- 0
  for (term in unique(rev(gsdata$Description))) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(legend.position = "none",
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))+geom_vline(xintercept = sum(geneList2>0), size=0.5)
  
  if (length(geneSetID) == 1) {
    ## geneList <- gsdata$geneList
    ## j <- which.min(abs(geneList))
    ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
    ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]
    
    ## v <- sort(c(v1, v2))
    ## inv <- findInterval(geneList, v)
    
    v <- seq(1, sum(gsdata$position), length.out=9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) inv <- inv + 1
    
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = col[unique(inv)])
    p2 <- p2 + geom_rect(
      aes_(xmin=~xmin,
           xmax=~xmax,
           ymin=~ymin,
           ymax=~ymax,
           fill=~I(col)),
      data=d,
      inherit.aes=FALSE)
  }
  
  ## p2 <- p2 +
  ## geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
  ##           ymin=ymin, ymax = ymin + yy, alpha=.5) +
  ## theme(legend.position="none") +
  ## scale_fill_gradientn(colors=color_palette(c("blue", "red")))
  
  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                            color="grey")
  p.pos <- p.pos + ylab("Ranked List Metric") +
    xlab(expression("Gene rank in gene list ordered by log"[2]*"FC")) +
    theme(axis.title.x=element_text(size=15))
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values=color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    # pd <- pd[order(pd[,1], decreasing=FALSE),]
    rownames(pd) <- pd$Description
    
    pd <- pd[,-1]
    pd <- round(pd, 4)
    
    tp <- tableGrob2(pd, p.res)
    
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, .5),
                        xmax = quantile(p.res$data$x, .95),
                        ymin = quantile(p.res$data$runningScore, .75),
                        ymax = quantile(p.res$data$runningScore, .9))
  }
  
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x=element_line(),
          axis.text.x = element_text())+xlab(expression("Gene rank in gene list ordered by log"[2]*"FC"))+theme(axis.title.x=element_text(size=15))
  
  if (length(subplots) == 1)
    return(plotlist[[1]] )
  
  
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  
  plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights=rel_heights)
}

##' @rdname gseaplot
##' @exportMethod gseaplot
setMethod("gseaplot", signature(x = "gseaResult"),
          function (x, geneSetID, by = "all", title = "", color='black',
                    color.line="green", color.vline="#FA5860", ...){
            gseaplot.gseaResult(x, geneSetID = geneSetID,
                                by = by, title = title,
                                color = color, color.line = color.line,
                                color.vline = color.vline, ...)
          })

##' @rdname gseaplot
##' @param color color of line segments
##' @param color.line color of running enrichment score line
##' @param color.vline color of vertical line which indicating the
##' maximum/minimal running enrichment score
##' @return ggplot2 object
##' @importFrom ggplot2 ggplot
##' @importFrom ggplot2 geom_linerange
##' @importFrom ggplot2 geom_line
##' @importFrom ggplot2 geom_vline
##' @importFrom ggplot2 geom_hline
##' @importFrom ggplot2 xlab
##' @importFrom ggplot2 ylab
##' @importFrom ggplot2 xlim
##' @importFrom ggplot2 aes
##' @importFrom ggplot2 ggplotGrob
##' @importFrom ggplot2 geom_segment
##' @importFrom ggplot2 ggplot_gtable
##' @importFrom ggplot2 ggplot_build
##' @importFrom ggplot2 ggtitle
##' @importFrom ggplot2 element_text
##' @importFrom ggplot2 rel
##' @importFrom aplot plot_list
##' @author Guangchuang Yu
gseaplot.gseaResult <- function (x, geneSetID, by = "all", title = "",
                                 color='black', color.line="green",
                                 color.vline="#FA5860", ...){
  by <- match.arg(by, c("runningScore", "preranked", "all"))
  gsdata <- gsInfo(x, geneSetID)
  p <- ggplot(gsdata, aes_(x = ~x)) +
    theme_dose() + xlab("Position in the Ranked List of Genes")
  if (by == "runningScore" || by == "all") {
    p.res <- p + geom_linerange(aes_(ymin=~ymin, ymax=~ymax), color=color)
    p.res <- p.res + geom_line(aes_(y = ~runningScore), color=color.line,
                               size=1)
    enrichmentScore <- x@result[geneSetID, "enrichmentScore"]
    es.df <- data.frame(es = which.min(abs(p$data$runningScore - enrichmentScore)))
    p.res <- p.res + geom_vline(data = es.df, aes_(xintercept = ~es),
                                colour = color.vline, linetype = "dashed")
    p.res <- p.res + ylab("Running Enrichment Score")
    p.res <- p.res + geom_hline(yintercept = 0)
  }
  if (by == "preranked" || by == "all") {
    df2 <- data.frame(x = which(p$data$position == 1))
    df2$y <- p$data$geneList[df2$x]
    p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                              color=color)
    p.pos <- p.pos + ylab("Ranked List Metric") +
      xlim(0, length(p$data$geneList))
  }
  if (by == "runningScore")
    return(p.res + ggtitle(title))
  if (by == "preranked")
    return(p.pos + ggtitle(title))
  
  p.pos <- p.pos + xlab(NULL) + theme(axis.text.x = element_blank(),
                                      axis.ticks.x = element_blank())
  p.pos <- p.pos + ggtitle(title) +
    theme(plot.title=element_text(hjust=0.5, size=rel(2)))
  #plot_list(gglist =  list(p.pos, p.res), ncol=1)
  
  aplot::gglist(gglist = list(p.pos, p.res), ncol=1)
}


##' extract gsea result of selected geneSet
##'
##'
##' @title gsInfo
##' @param object gseaResult object
##' @param geneSetID gene set ID
##' @return data.frame
##' @author Guangchuang Yu
## @export
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}


gseaplot2_custom_legend <- function(x, order,geneSetID, title = "", color="green", base_size = 11,
                                    rel_heights=c(1.5, .5, 1), subplots = 1:3,
                                    pvalue_table = FALSE, ES_geom="line") {
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  
  geneList <- position <- NULL ## to satisfy codetool
  
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  } else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
    
  }
  gsdata$Description=factor(gsdata$Description,levels=order)
  gsdata=gsdata[order(gsdata$Description),]
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) +
    theme_classic(base_size) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank()) +
    scale_x_continuous(expand=c(0,0))
  
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color= ~Description),
                          size=0)
  } else {
    es_layer <- geom_point(aes_(y = ~runningScore, color= ~Description),
                           size=0, data = subset(gsdata, position == 1))
  }
  
  p.res <- p + es_layer +
    theme( legend.title = element_blank(),
           legend.background = element_rect(fill = "transparent"),legend.text=element_text(size=12))+guides(color = guide_legend(byrow = TRUE))+
    guides(color = guide_legend(override.aes = list(size = 6)))
  
  p.res <- p.res + ylab("Running Enrichment Score") +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.line.x=element_blank(),axis.title.y=element_text(size=15))
  
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) +
    geom_linerange(aes_(ymin=~ymin, ymax=~ymax, color=~Description)) +
    xlab(NULL) + ylab(NULL) + theme_classic(base_size) +
    theme(legend.position = "none",
          axis.ticks = element_blank(),
          axis.text = element_blank(),
          axis.line.x = element_blank()) +
    scale_x_continuous(expand=c(0,0)) +
    scale_y_continuous(expand=c(0,0))
  
  if (length(geneSetID) == 1) {
    ## geneList <- gsdata$geneList
    ## j <- which.min(abs(geneList))
    ## v1 <- quantile(geneList[1:j], seq(0,1, length.out=6))[1:5]
    ## v2 <- quantile(geneList[j:length(geneList)], seq(0,1, length.out=6))[1:5]
    
    ## v <- sort(c(v1, v2))
    ## inv <- findInterval(geneList, v)
    
    v <- seq(1, sum(gsdata$position), length.out=9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) inv <- inv + 1
    
    col <- c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * .3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[as.character(unique(inv))])
    d <- data.frame(ymin = ymin, ymax = yy,
                    xmin = xmin,
                    xmax = xmax,
                    col = col[unique(inv)])
    p2 <- p2 + geom_rect(
      aes_(xmin=~xmin,
           xmax=~xmax,
           ymin=~ymin,
           ymax=~ymax,
           fill=~I(col)),
      data=d,
      inherit.aes=FALSE)
  }
  
  ## p2 <- p2 +
  ## geom_rect(aes(xmin=x-.5, xmax=x+.5, fill=geneList),
  ##           ymin=ymin, ymax = ymin + yy, alpha=.5) +
  ## theme(legend.position="none") +
  ## scale_fill_gradientn(colors=color_palette(c("blue", "red")))
  
  df2 <- p$data #data.frame(x = which(p$data$position == 1))
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data=df2, aes_(x=~x, xend=~x, y=~y, yend=0),
                            color="grey")
  p.pos <- p.pos + ylab("Ranked List Metric") +
    xlab(expression("Gene rank in gene list ordered by log"[2]*"FC")) +
    theme(axis.title.x=element_text(size=15))
  
  if (!is.null(title) && !is.na(title) && title != "")
    p.res <- p.res + ggtitle(title)
  
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values=color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    } else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    # pd <- pd[order(pd[,1], decreasing=FALSE),]
    rownames(pd) <- pd$Description
    
    pd <- pd[,-1]
    pd <- round(pd, 4)
    
    tp <- tableGrob2(pd, p.res)
    
    p.res <- p.res + theme(legend.position = "none") +
      annotation_custom(tp,
                        xmin = quantile(p.res$data$x, .5),
                        xmax = quantile(p.res$data$x, .95),
                        ymin = quantile(p.res$data$runningScore, .75),
                        ymax = quantile(p.res$data$runningScore, .9))
  }
  
  
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] +
    theme(axis.line.x = element_line(),
          axis.ticks.x=element_line(),
          axis.text.x = element_text())+xlab(expression("Gene rank in gene list ordered by log"[2]*"FC"))+theme(axis.title.x=element_text(size=15))
  
  if (length(subplots) == 1)
    return(plotlist[[1]] )
  
  
  if (length(rel_heights) > length(subplots))
    rel_heights <- rel_heights[subplots]
  
  plot_grid(plotlist = plotlist, ncol=1, align="v", rel_heights=rel_heights)
}


plotGeneCountpseudo <- function(curve, counts = NULL, gene = NULL, clusters = NULL,
                                models = NULL, title = NULL){
  rd <- slingshot::slingReducedDim(curve)
  if (!is.null(gene)) {
    logcounts <- log1p(counts[gene, ])
    cols <- logcounts
    scales <- scale_color_brewer(palette="Reds")
    if (is.null(title)) title <- paste0("logged count of gene ", gene)
  } else {
    cols <- slingPseudotime(crv, na = FALSE)[,1]
    scales <- NULL
    rescale <- function(x) (x-min(x))/(max(x) - min(x)) * 100
    cols=rescale(cols)
    if (is.null(title)) title <- "Pseudotime"
  }
  # Getting the main plot
  df <- data.frame(dim1 = rd[, 1], dim2 = rd[, 2], col = as.numeric(cols))
  p <- ggplot(df, aes(x = dim1, y = dim2, col = col)) +
    geom_point(size = 3,aes(color=col))+scale_color_distiller(palette="Reds",trans = "reverse")+
    theme_classic() +
    labs(col = title) +
    theme(axis.title.x=element_text(size=20))+ 
    theme(axis.title.y=element_text(size=20,vjust=0))+ 
    theme(axis.text.x=element_text(size=10,margin=margin(0,0,10,0)))+ 
    theme(axis.text.y=element_text(size=10,margin=margin(0,0,0,15)))+ 
    theme(legend.title=element_text(size=20))+ 
    theme(legend.text=element_text(size=15))+xlab("Dimension 1") +
    ylab("Dimension 2")
  
  # Adding the curves
  #  for (i in seq_along(slingCurves(curve))) {
  #    curve_i <- slingCurves(curve)[[i]]
  #    curve_i <- curve_i$s[curve_i$ord, seq_len(2)]
  #    colnames(curve_i) <- c("dim1", "dim2")
  #    p <- p + geom_path(data = as.data.frame(curve_i), col = "black", size = 1)
  #  }
  
  # Adding the knots
  nCurves <- length(slingCurves(curve))
  if (!is.null(models)) {
    if (is(models, "list")) {
      sce <- FALSE
    } else if(is(models, "SingleCellExperiment")){
      sce <- TRUE
    }
    if(!sce){
      m <- .getModelReference(models)
      knots <- m$smooth[[1]]$xp
    } else if(sce){
      knots <- S4Vectors::metadata(models)$tradeSeq$knots
    }
    # times <- slingPseudotime(curve, na = FALSE)
    knots_dim <- matrix(ncol = 2, nrow = nCurves * length(knots))
    for (ii in seq_along(slingCurves(curve))) {
      S <- project_to_curve(x = slingCurves(curve)[[ii]]$s,
                            s = slingCurves(curve)[[ii]]$s[slingCurves(curve)[[ii]]$ord, ],
                            stretch = 0)
      for (jj in seq_along(knots)) {
        kn <- knots[jj]
        times <- S$lambda
        knot <- which.min(abs(times - kn))
        knots_dim[jj + (ii-1)*length(knots), ] <- S$s[knot, seq_len(2)]
      }
    }
    knots_dim <- as.data.frame(knots_dim)
    colnames(knots_dim) <- c("dim1", "dim2")
    p <- p +
      geom_point(data = knots_dim, col = "black", size = 2)
  }
  return(p)
}


