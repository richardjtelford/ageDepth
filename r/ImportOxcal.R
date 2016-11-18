#Import OxCal functions

get.vectors <- function(js, target) {
  x <- js[grep(target, js)]
  lab <- gsub("^(ocd\\[\\d+\\]).*", "\\1", x )
  val <- gsub(".+=\\[(.+)\\];", "\\1", x)
  val <- strsplit(val, ", ")
  res <- sapply(val, as.numeric)
  setNames(res, lab)
}

get.values <- function(js, target, name = "value", numeric = TRUE) {
  x <- js[grep(target, js)]
  val <- gsub(".+=(.+);", "\\1", x)
  lab <- gsub("^(ocd\\[\\d+\\]).*", "\\1", x )
  if(numeric) {
    val <- as.numeric(val)
  } 
  res <- data_frame(lab, val)
  setNames(res, c("lab", name))
}


get.OxCal <- function(js, posterior = TRUE){
  if(posterior){
    what <- "posterior"
  }else{
    what <- "likelihood"
  }
  depths <- get.values(js, "\\.data\\.z", "depth")
  comment <- get.values(js, paste0(what,"\\.comment\\[0\\]"), "comment", FALSE)
  depth_comment <- left_join(depths, comment)
  
  start <- get.values(js, paste0(what,"\\.start"), "start")
  resolution <- get.values(js, paste0(what,"\\.resolution"), "resolution")
  meta <- full_join(start, resolution)
  meta <- inner_join(depth_comment, meta)
  probs <- get.vectors(js, paste0(what,"\\.prob="))
  
  out <- plyr::adply(meta, 1, function(r){
    prob <- probs[[r$lab]]
    date <- r$start + 0:(length(prob) - 1) * r$resolution
    data_frame(lab = r$lab, depth = r$depth, date = date, prob = prob, comment  = r$comment, posterior = posterior)
  })
  out$date <- 1950 - out$date#convert to BP
  out
}
