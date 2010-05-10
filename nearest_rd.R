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

    nr <-  subject[nearest(q.iranges, s.iranges ),]
    rownames(nr) <- names(q.iranges)
    res[[sp]] <- nr 
  }
  return(res)  
}

