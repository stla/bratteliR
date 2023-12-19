ids <- function(I, prefix) {
  ndigits <- floor(log10(I)) + 1
  paste0(prefix, sprintf(paste0("%0", ndigits, "d"), 1L:I))
}

if(getRversion() >= "2.15.1") {
  utils::globalVariables(
    c(
      ".N", "node1", "node2", "bend", "code", "edgelabel", "id", "level",
      "mindex", "node", "nodelabel", "x", "y"
    )
  )
}

#' @title LaTeX code drawing a Bratteli graph
#' @description Generate a LaTeX file containing TikZ code that renders a
#'   picture of a Bratteli graph.
#'
#' @param outfile path to the output file
#' @param Mn a function returning for each integer \code{n} the incidence
#' matrix between levels \code{n} and \code{n+1}; the matrix \code{Mn(0)}
#' must have one and only one row
#' @param N the level up to which the graph is wanted
#' @param edgelabels \code{"default"}, \code{"letters"}, \code{"order"},
#'   \code{"kernels"}, \code{NA}, or a \emph{vectorized} function with
#'   four arguments: the level of the graph, the index of the "from" vertex,
#'   the index of the "to" vertex, and the index of the edge among the
#'   multiple edges, if there are multiple edges
#' @param vertexlabels \code{"colnames"} to use the column names of the
#'   matrices, \code{"dims"} to use the dimensions of the vertices, \code{NA},
#'   or a function with one argument, the level of the graph, returning
#'   for level \code{n} the vector of labels at the \code{n}-th level
#' @param colorpath an index of a path to be colored, or \code{NA}
#' @param rootlabel a label for the root vertex of the graph
#' @param latex Boolean, whether to enclose all labels between dollars
#' @param xscale,yscale scaling factors for the graph
#' @param bending curvature when there are multiple edges
#' @param hor Boolean, whether to render a horizontal graph
#' @param mirror Boolean, whether to "reverse" the graph
#' @param northsouth Boolean, whether to draw the edges with endpoints at the
#'   bottom and the top of the vertex labels
#'
#' @return No return value; called to generate the LaTeX file.
#'
#' @export
#' @importFrom gmp numerator denominator as.bigq
#' @importFrom data.table `:=` data.table setnames
#' @importFrom diagram coordinates
#' @importFrom utils tail
bratteliGraph <- function(
    outfile, Mn, N,
    edgelabels = NA,
    vertexlabels = "colnames",
    colorpath = NA,
    rootlabel = "\\varnothing", latex = TRUE,
    xscale = 50, yscale = 50, bending = 1,
    hor = FALSE, mirror = FALSE,
    northsouth = FALSE
) {
  packages <- NULL
  stopifnot(
    is.function(edgelabels) || is.na(edgelabels) ||
      edgelabels %in% c("default", "letters", "order", "kernels")
  )
  stopifnot(
    is.function(vertexlabels) || is.na(vertexlabels) ||
      vertexlabels %in% c("colnames", "dims")
  )
  Ms <- lapply(0L:(N-1), Mn)
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
  setnames(elpos, c("x", "y"))
  if(mirror) {
    if(!hor) elpos[, y := max(y)-y] else elpos[, x := max(x)-x]
  }
  elpos[, level := rep(seq_along(nvertices), times = nvertices) - 1L]
  # scale
  elpos[, `:=`(x = xscale*x, y = yscale*y)]
  # node id's
  elpos[, `:=`(node = ids(.N, LETTERS[level[1L]+1L])), by = "level"]
  if(is.na(vertexlabels)) {
    fvertexlabels <- function(n) rep("", ncol(Ms[[n]]))
  } else if(is.character(vertexlabels)){
    if(vertexlabels == "colnames") {
      fvertexlabels <- function(n) colnames(Ms[[n]])
    } else if(vertexlabels == "dims"){
      dims <- bratteliDimensions(Mn, N)
      fvertexlabels <- function(n) dims[[n]]
    }
  }
  elpos[, nodelabel := c(rootlabel, unlist(lapply(1L:N, fvertexlabels)))]
  if(latex) elpos[, nodelabel := paste0("$", nodelabel, "$")]
  # code for nodes
  elpos[, code := sprintf(
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
    } else if(edgelabels == "letters") {
      connections[, edgelabel:=letters[seq_along(to)], by = node1]
    } else if(edgelabels == "order") {
      connections[, edgelabel := seq_len(.N)-1L, by = node2]
    } else if(edgelabels == "kernels"){
      if(!is.element("nicefrac", packages)) packages <- c(packages, "nicefrac")
      ckernels <- bratteliKernels(Mn, N)
      ckernels_numer <- lapply(ckernels, function(x) {
        as.character(numerator(as.bigq(x)))
      })
      ckernels_denom <- lapply(ckernels, function(x) {
        as.character(denominator(as.bigq(x)))
      })
      f <- Vectorize(function(level, from, to){
        if(ckernels_numer[[level+1L]][to, from] == "0") return("0")
        if(ckernels_denom[[level+1L]][to, from] == "1") return("1")
        sprintf(
          "\\nicefrac{%s}{%s}",
          ckernels_numer[[level+1L]][to, from],
          ckernels_denom[[level+1L]][to, from]
        )
      })
      connections[, edgelabel := f(level, from, to)]
    }
  } else if(is.function(edgelabels)){
    labels_on_edge <- TRUE
    connections[, edgelabel := edgelabels(level, from, to, mindex)]
  } else {
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
  if(!is.na(colorpath)){
    colorpath <- as.integer(colorpath)
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
    if(latex) connections[, edgelabel := paste0("$", edgelabel, "$")]
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
  templateFile <-
    system.file("templateGraph", "templateGraph.tex", package = "bratteli")
  template <- paste0(readLines(templateFile), collapse = "\n")
  tex <- sprintf(template, packages, Code)
  writeLines(tex, outfile)
  #return(connections)
  #return(paths)
  if(!is.na(colorpath)){
    return(invisible(paths))
  }
  return(invisible())
}
