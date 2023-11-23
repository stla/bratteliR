ids <- function(I, prefix) {
  ndigits <- floor(log10(I)) + 1
  paste0(prefix, sprintf(paste0("%0", ndigits, "d"), 1L:I))
}

if(getRversion() >= "2.15.1") {
  utils::globalVariables(c(".N", "node1", "node2"))
}

#' Generate TikZ code of a Bratteli graph
#'
#' @param edgelabels \code{"default"}, \code{"default_letters"},
#'   \code{"order"}, \code{"kernels"}, \code{NA}, or a VECTORIZED function
#' @param vertexlabels \code{"colnames"} (default) to use the column names of
#'   the matrices, \code{"dims"} to use the dimensions of the vertices,
#'   or a VECTORIZED(?) function
#' @param bending curvature when there are multiple edges
#' @param northsouth node connections
#'
#' @export
#' @importFrom gmp numerator denominator
#' @importFrom data.table `:=` data.table
#' @importFrom diagram coordinates
#'
#' @examples
#' Pascal_Mn <- function(n){
#'  M <- matrix(0, nrow=n+1, ncol=n+2)
#'  for(i in 1:(n+1)){
#'    M[i, ][c(i, i+1)] <- 1
#'  }
#'  return(M)
#' }
#' BgraphTikZ("/tmp/PascalGraph.tex", Pascal_Mn, 3)
#'
bratteliGraph <- function(
    outfile, Mn, N,
    edgelabels = "default",
    vertexlabels = "colnames",
    fvertexlabels = NULL,
    colorpath = NULL,
    ROOTLABEL = "\\varnothing", LATEXIFY=TRUE,
    packages = NULL,
    scale = c(50,50), bending = 1,
    hor = FALSE, mirror = FALSE,
    northsouth = FALSE
) {
  Ms <- lapply(0:(N-1), Mn)
  for(i in 1L:N){
    if(is.null(colnames(Ms[[i]]))) {
      colnames(Ms[[i]]) <- seq_len(ncol(Ms[[i]]))
    }
  }
  nvertices <- vapply(1L:N, function(n) { # number of vertices per level
    nrow(Ms[[n]])
  }, integer(1L))
  left <- c(0L, cumsum(nvertices))
  vertex <- function(n, k){ # n: level ; k: vertex at this level
    left[n] + k
  }
  nvertices <- c(nvertices, ncol(Ms[[N]]))
  elpos <- coordinates(nvertices, relsize = 1, hor = !hor) # positions of vertices
  elpos <- data.table(elpos)
  names(elpos) <- c("x", "y")
  if(mirror){
    if(!hor) elpos[, y := max(y)-y] else elpos[, x := max(x)-x]
  }
  elpos[["level"]] <- rep(seq_along(nvertices), times = nvertices) - 1L
  # scale
  elpos[, `:=`(x = scale[1]*x, y = scale[2]*y)]
  # node id's
  elpos[, `:=`(node = ids(.N, LETTERS[level[1L]+1L])), by = "level"]
  if(is.null(fvertexlabels) && is.character(vertexlabels)){
    if(vertexlabels == "colnames") {
      fvertexlabels <- function(n) colnames(Ms[[n]])
    }
    if(vertexlabels == "dims"){
      dims <- bratteliDimensions(Mn, N)
      fvertexlabels <- function(n) dims[[n]]
    }
  }
  elpos[["nodelabel"]] <- c(ROOTLABEL, lapply(1L:N, fvertexlabels))
  if(LATEXIFY) elpos[, nodelabel:=paste0("$", nodelabel, "$")]
  # code for nodes
  elpos[, code:=sprintf(
    "\\node[VertexStyle](%s) at (%s, %s) {%s};", node, x, y, nodelabel
  )]
  # code for edges
  connections <- data.frame(
    level = integer(), from = integer(), to = integer(),
    multiplicity = integer(), node1 = character(), node2 = character(),
    stringsAsFactors = FALSE
  )
  counter <- 1L
  for(n in 0L:(N-1)){
    M <- Ms[[n+1L]]
    for(i in 1L:nrow(M)){
      from <- vertex(n+1L, i)
      for(k in which(M[i, ]>0)){
        to <- vertex(n+2L, k)
        for(m in 1L:M[i, k]){
          connections[counter, ] <- data.frame(
            n, i, k, M[i, k], elpos[from, ][["node"]], elpos[to, ][["node"]],
            stringsAsFactors = FALSE
          )
          counter <- counter + 1L
        }
      }
    }
  }
  connections <- data.table(connections)
  connections[, id := paste0(level, from, to), by = 1L:nrow(connections)] # connections id's
  ## multiplicity index
  connections[, `:=`(mindex = seq_len(.N)), by = "id"]
  # edge labels
  if(is.character(edgelabels)){
    labels_on_edge <- TRUE
    if(edgelabels == "default") {
      connections[, edgelabel := seq_along(to)-1L, by = node1]
    } else if(edgelabels == "default_letters") {
      connections[, edgelabel:=letters[seq_along(to)], by = node1]
    } else if(edgelabels == "order") {
      connections[, edgelabel := seq_len(.N)-1L, by = node2]
    } else if(edgelabels == "kernels"){
      if(!is.element("nicefrac", packages)) packages <- c(packages, "nicefrac")
      ckernels <- bratteliKernels(Mn, N)
      ckernels_numer <- lapply(ckernels, function(x) {
        as.character(numerator(x))
      })
      ckernels_denom <- lapply(ckernels, function(x) {
        as.character(denominator(x))
      })
      f <- Vectorize(function(level, from, to){
        if(ckernels_denom[[level+1L]][to, from] == "1") return("1")
        sprintf(
          "\\nicefrac{%s}{%s}",
          ckernels_numer[[level+1L]][to, from],
          ckernels_denom[[level+1L]][to, from]
        )
      })
      connections[, edgelabel := f(level, from, to)]
    }
  }
  if(is.function(edgelabels)){
    labels_on_edge <- TRUE
    connections[, edgelabel := edgelabels(level, from, to, mindex)]
  }
  if(is.atomic(edgelabels) && is.na(edgelabels)){
    labels_on_edge <- FALSE
  }
  # curvatures
  fbend <- function(m) {
    if(m == 1L) return(NA_real_)
    bend <- seq(0, 10*bending, length.out = m) * (2*m - 2)
    bend - mean(bend)
  }
  connections[, bend := fbend(.N), by = "id"]
  # paths
  if(!is.null(colorpath)){
    dd <- as.data.frame(connections)
    ff <- function(e) {
      e1 <- which(dd[["node1"]] == dd[["node2"]][e]) # edges connected to edge i
      paths <- lapply(e1, function(j) c(e, j))
      if(N == 2){
        return(paths)
      }
      for(m in 2L:(N-1)){
        e2 <- lapply(paths, function(p) {
          which(dd[["node1"]] == dd[["node2"]][tail(p, 1L)])
        })
        paths <-
          do.call(c, lapply(seq_along(paths), function(p) {
            lapply(e2[[p]], function(ee) c(paths[[p]], ee))}
          ))
      }
      paths
    }
    paths <- do.call(c, lapply(which(dd$level == 0L), ff))
    path <- logical(nrow(dd))
    if(colorpath > length(paths)) {
      warning(sprintf(
        "There are %d paths and you supplied `colorpath=%d`",
        length(paths), colorpath
      ))
    } else {
      path[paths[[colorpath]]] <- TRUE
    }
  } else {
    path <- logical(nrow(connections))
  }
  connections[, path := path]
  #
  if(labels_on_edge) {
    if(LATEXIFY) connections[, edgelabel := paste0("$", edgelabel, "$")]
    drawcode <- Vectorize(function(bend) {
      if(is.na(bend)) {
        return(ifelse(
          northsouth,
          "\\draw[%s](%s.south) to node[EdgeLabelStyle]{%s} (%s.north);",
          "\\draw[%s](%s) to node[EdgeLabelStyle]{%s} (%s);"
        ))
      }
      paste0(
        sprintf("\\draw[%s, bend left=%s]", "%s", bend),
        ifelse(northsouth,
               "(%s.south) to node[EdgeLabelStyle]{%s} (%s.north);",
               "(%s) to node[EdgeLabelStyle]{%s} (%s);")
      )
    })
    connections[, code := sprintf(
      drawcode(bend),
      ifelse(path, "EdgeStylePath", "EdgeStyle"), node1, edgelabel, node2
    )]
  } else {
    drawcode <- Vectorize(function(bend) {
      if(is.na(bend)) return(ifelse(
        northsouth,
        "\\draw[%s](%s.south) to (%s.north);",
        "\\draw[%s](%s) to (%s);"
      ))
      paste0(
        sprintf("\\draw[%s, bend left=%s]", "%s", bend),
        ifelse(northsouth, "(%s.south) to (%s.north);", "(%s) to (%s);")
      )
    })
    connections[, code := sprintf(
      drawcode(bend), ifelse(path, "EdgeStylePath", "EdgeStyle"), node1, node2
    )]
  }
  # TikZ code
  Code <- paste0("\t", c(elpos[["code"]], connections[["code"]]), collapse="\n")
  # add packages
  if(!is.null(packages)){
    packages <- paste0(sapply(packages, function(x) sprintf("\\usepackage{%s}\n", x)), collapse="")
  }else{
    packages <- ""
  }
  # write code to template
  template <- system.file("templates", "template_BratteliTikZ3.RDS", package="bratteli")
  texfile <- sprintf(readRDS(template), packages, Code)
  writeLines(texfile, outfile)
  #return(connections)
  #return(paths)
  if(!is.null(colorpath)){
    return(invisible(paths))
  }
  return(invisible())
}
