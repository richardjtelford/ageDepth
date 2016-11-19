#Import OxCal functions
## ---- getVectors
get.vectors <- function(js, target) {
  x <- js[grep(target, js)]
  lab <- gsub("^(ocd\\[\\d+\\]).*", "\\1", x )
  val <- gsub(".+=\\[(.+)\\];", "\\1", x)
  val <- strsplit(val, ", ")
  res <- sapply(val, as.numeric)
  setNames(res, lab)
}

## ---- getValues
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

## ---- getOxcal
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
  out$start <- NULL
  out$resolution <- NULL
  out$date <- 1950 - out$date#convert to BP
  out
}

## ---- getAgemodel
get.ageModel <- function(x, ci = c(0.023, 0.977)){
  x %>% 
    filter(posterior == TRUE) %>% 
    group_by(depth, comment) %>% 
    mutate(sumProb = cumsum(prob)) %>% 
    mutate(sumProb = sumProb/max(sumProb)) %>% 
    summarise(
      medianDate = approx(x = sumProb, y = date, xout =  0.5)$y,
      lowCI = approx(x = sumProb, y = date, xout = ci[1])$y,
      highCI = approx(x = sumProb, y = date, xout = ci[2])$y
    ) 
}

## ---- plotdatesmodel
plot.dates_model <- function(dates, model, title){
  dates %>% filter(grepl("Dep", comment)) %>% #just the dates
    ggplot(aes(x = depth, y = date, size = prob, colour = posterior, group = comment)) +   
    geom_ribbon(data = model, aes(x = depth, ymax = highCI, ymin = lowCI), inherit.aes = FALSE, alpha  = 0.2) +
    geom_line(data = mod, aes(x = depth, y = medianDate), inherit.aes = FALSE, colour = "grey40") +
    geom_line(alpha = 0.6) + 
    scale_x_reverse() +
    scale_size(range = c(0, 6)) +
    coord_flip() +
    labs(x = "Depth cm", y = "Date yr BP") + 
    theme_bw() + 
    guides(size = "none", colour = guide_legend(title = "Posterior", override.aes = list(size = 4))) + 
    ggtitle(title)
}