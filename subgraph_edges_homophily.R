subgraph_edges_homophily <- function(graph, vattr_name, heterophily = FALSE,
                                     drop_isolates = FALSE) {
  stopifnot( # arg checks
    igraph::is.igraph(graph) || is.character(vattr_name) || 
      length(vattr_name) == 1L || !is.na(vattr_name) || 
      vattr %in% igraph::vertex_attr_names(vattr_name)
  )
  
  vattrs <- igraph::vertex_attr(graph, name = vattr_name)
  total_el <- igraph::as_edgelist(graph, names = FALSE)
  
  # rows from total_el where the attribute of the edge source == attribute of edge target
  edges_to_keep <- vattrs[total_el[, 1L]] == vattrs[total_el[, 2L]]
  
  # for heterophilous ties, just negate the "in_group" version
  if (heterophily) edges_to_keep <- !edges_to_keep
  
  igraph::subgraph.edges(graph, 
                         eids = which(edges_to_keep), 
                         delete.vertices = drop_isolates)
}