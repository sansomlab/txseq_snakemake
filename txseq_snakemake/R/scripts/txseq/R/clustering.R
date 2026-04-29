#' Function to return a clustered dendrogram
#' with optional optimised leaf order.
#'
#' @param optimise Should the dendrogram leaf order be optimised (TRUE|FALSE) 
#' @param dist_method The distance method, passed to the "dist" function
#' @param clust_method The clustering method, passed to the "hclust" function
#' @param p The "p" parameter for the "minkowski" distance metric
#'
#' @import dendsort
#'
#' @export
get_den <- function(m,
                    optimize=F,
                    dist_method="manhattan",
                    clust_method="complete",
                    p=1.5)
    {

    require(fastcluster)

    if(dist_method=="cor" | dist_method=="correlation")
        {
        dissym <- (1 - cor(m))/2
        d <- as.dist(dissym)
        }
    else {if(dist_method=="minkowski") {
            d <- dist(m, method=dist_method, p=p)
            }
          else {
              d <- dist(m[,], method=dist_method)
            }
         }

    h <- hclust(d, method=clust_method)

    if(optimize)
    {
      require(dendsort)
      h <- dendsort(h)
          #co <- order.optimal(d, h$merge)

    # overwrite the hclust object with the optimal leaf order
    #h$merge <- co$merge
    #h$order <- co$order
    }

    as.dendrogram(h)
   }


#' function to draw heatmap from RNA-seq counts or FPKMS
#' @export
rnaseq_heatmap <- function(m,
                           rowden=F,
                           colden=F,
                           dendrogram="column",
                           labRow=NULL,
                           ramp_colors=c("black","blue","white","yellow","orange","red","brown"),
                           key.xticks=c(1:4),
                           key.xlab="expression level",
                           log_factor=10,
                           mar=c(14,12),
                           title=NULL,
                           lhei=c(1.15,4),
                           lwid=c(1.3,4),
                           cexRow=1,
                           cexCol=1,
                           ...
                           )
{

  nbreaks=200 # number of graduations
  rm <- range(m)
  breaks=seq(rm[1],rm[2],diff(rm)/nbreaks) #range (-2 -> 2) and breakpoint size = range_size / no. graduations
  colors=colorRampPalette(ramp_colors)(nbreaks) # get the color pallete

  heatmap.2(m,
            col=colors,
            breaks=breaks,
            scale="none",
            Rowv=rowden,
            Colv=colden,
            mar=mar,
            trace="none",
            dendrogram=dendrogram,
            colsep=c(seq(1,dim(m)[2],1)),
            sepwidth=c(0.02,0.02),
            key.title = "",
            #lmat = rbind(c(0,3),c(2,1),c(0,4)),
            key.par=list(mar=c(5.5,.3,5.5,.3)),
            #key.ytickfun = "",
            key.xtickfun = function()
                 {
                    int=key.xticks
                    at = int / max(m)
                    labels = as.character(log_factor^int)
                    return(list(labels=labels, at=at))
                 },
            density.info=c("none"),
            lwid = lwid,
            lhei = lhei,
            key.xlab = key.xlab,
            key.ylab = "",
            labRow = labRow,
            main=title,
            cexRow = cexRow,
            cexCol = cexCol
    )
    }