#takes a query rd and a subject rd and returns
#a list with one entry per space containing a data.frame
#of the nearest ranges in subject to each of the ranges in
#query. Rows are named by query range name
nearest.rd <- function(query, subject ){
  spaces.query <- names(query)
  spaces.subject <- names(subject)
  common.spaces <- intersect(spaces.query, spaces.subject)
  if(length(common.spaces)<1){stop("No common spaces between subject and query")}

  res <- list()
  for(sp in common.spaces){

    q.iranges <- unlist(ranges(query[sp]))
    s.iranges <- unlist(ranges(subject[sp]))

    inds <- nearest(q.iranges, s.iranges)
    nr <-  as.data.frame(subject[sp][inds,])
    rownames(nr) <- names(q.iranges)
    res[[sp]] <- nr 
  }
  return(res)  
}

# A note on the nearest function:
#     1. Find the ranges in ‘subject’ that overlap ‘xi’. If a single
#          range ‘si’ in ‘subject’ overlaps ‘xi’, ‘si’ is returned as
#          the nearest neighbor of ‘xi’. If there are multiple overlaps,
#          one of the overlapping ranges is chosen arbitrarily.
#
#       2. If no ranges in ‘subject’ overlap with ‘xi’, then the range
#          in ‘subject’ with the shortest distance from its end to the
#          start ‘xi’ or its start to the end of ‘xi’ is returned.
#
