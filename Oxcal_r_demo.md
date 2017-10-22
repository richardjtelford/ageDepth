# OxCal & R
Richard J. Telford  
October 2, 2016  

[OxCal](https://c14.arch.ox.ac.uk/oxcal.html) is perhaps the most powerful of the three main Bayesian age-depth modelling procedures, with many many options and the potential for building far more complicated than your typical palaeolimnologist needs to use. Unlike [Bchron](https://cran.r-project.org/web/packages/Bchron/index.html) and [Bacon](http://www.chrono.qub.ac.uk/blaauw/bacon.html), OxCal is not R based. Much as I want to use OxCal, I also want to use R. 

In this post, I show how OxCal models can be built in R, the OxCal command line program called from R, and the output imported into R and plotted. The code does not deal will all possible outputs from OxCal (it may need altering if your model is different from mine) and absolutely does not absolve you from checking the log files for any problems. 

You need to start by [downloading OxCal](https://c14.arch.ox.ac.uk/oxcalhelp/readme.html#local) and installing it somewhere convenient on your computer.



We need a couple of packages for data manipulation and plotting.


```r
library("dplyr")
library("ggplot2")
```

And we need some radiocarbon dates. Here are a few of the dates from a model I'm currently working on (the full set of dates takes about six hours to process).


```r
dates <- read.table(header = TRUE, text = 
"depth_corr	Lab_Code	raw_14C_age	error
312.5	X-1	16010	100
352   X-2	16670	 70
406.5	X-3	17730	130
440.5	X-4	18690	120
490.5	X-5	20150	220
535.5	X-6	20930	170
")
```

OxCal was written by archaeologists so you need to think like an archaeologist to use it: oldest dates need to come first in the model. (I suspect many palaeolimnologists get this wrong and then wait days before the model refuses to converge. At least I've done this a couple of times.)


```r
##sort so oldest first
dates <- dates[order(dates$depth_corr, decreasing = TRUE), ]
```

I want to use a P sequence model to fit an age-depth model to the dates. I'm using the Marine13 calibration curve with a local reservoir offset of 40 +/-30 years. I'm using the P sequence model with a _k_ selected from a prior distribution by OxCal and with an outlier model. All these commands and their arguments are explained in the on-line help for OxCal. `writeLines` prints the formatted model to the screen. I've used lots of white space (don't worry, it's free) to make the model human readable.


```r
#Calibration information
curve <- ' Curve("marine13");'
DR <- ' Delta_R("Local Marine",40,30);'

#The dates
ageDepths <- with(dates, sprintf('  R_date("Dep%s",%d,%d){\n    z=%s;\n    Outlier(0.05);\n  };', depth_corr, raw_14C_age, error, depth_corr))
ageDepths <- paste(ageDepths, collapse = "\n")

#Putting it all together
commandp <- paste(
  curve, 
  DR,
  ' Outlier_Model("General",T(5),U(0,4),"t");',
  ' P_Sequence("k",1,1,U(-2,2)){',
  '  Boundary("Bottom");',
     ageDepths,
  '  Boundary("Top");',
  ' };', 
  sep = "\n")

#Wrapping in Plot()
commandp <- paste('Plot(){', commandp, '};', sep = "\n")
writeLines(commandp)
```

```
## Plot(){
##  Curve("marine13");
##  Delta_R("Local Marine",40,30);
##  Outlier_Model("General",T(5),U(0,4),"t");
##  P_Sequence("k",1,1,U(-2,2)){
##   Boundary("Bottom");
##   R_date("Dep535.5",20930,170){
##     z=535.5;
##     Outlier(0.05);
##   };
##   R_date("Dep490.5",20150,220){
##     z=490.5;
##     Outlier(0.05);
##   };
##   R_date("Dep440.5",18690,120){
##     z=440.5;
##     Outlier(0.05);
##   };
##   R_date("Dep406.5",17730,130){
##     z=406.5;
##     Outlier(0.05);
##   };
##   R_date("Dep352",16670,70){
##     z=352;
##     Outlier(0.05);
##   };
##   R_date("Dep312.5",16010,100){
##     z=312.5;
##     Outlier(0.05);
##   };
##   Boundary("Top");
##  };
## };
```

Once we're happy with the model, we can save the file and call the OxCal command line utility. There are different versions of this for Linux, Mac and Windows.


```r
writeLines(commandp, con = "oxcal_demo.input")
#system(paste("~/oxcal/OxCal/bin/OxCalLinux", "oxcal_demo.input"))
```

Depending on how complicated your model is, you now have time for a cup of tea or a round the world cruise. During this time, R will sit there apparently doing nothing, but in the background, OxCal is busy shunting electrons round.

OxCal will produce three files, `oxcal_demo.log`, `oxcal_demo.txt`, and `oxcal_demo.js`. You need to look in the `.log` file and check for any problems. The `.txt` file contains a summary of the model output and the `.js` file contains all the details. It is this file we need to read into R. Unfortunately it is a large and complicated file which takes a bit of effort to read. I've made some functions to help import the data in to a `data.frame` that I can plot with `ggplot`. The functions extract what I need - they will need modifying if you need other information or if you have used a different naming convention for the dates.


```r
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
```

```r
get.vectors <- function(js, target) {
  x <- js[grep(target, js)]
  lab <- gsub("^(ocd\\[\\d+\\]).*", "\\1", x )
  val <- gsub(".+=\\[(.+)\\];", "\\1", x)
  val <- strsplit(val, ", ")
  res <- sapply(val, as.numeric)
  setNames(res, lab)
}
```

```r
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
```
These functions use regular expressions to find and extract the required information. If you don't know how to use regular expressions, it is well worth spending a morning learning how to use them as they are very useful for processing data.

With these functions, we need to first read the `.js` file into memory and then extract the posterior and likelihood information separately.


```r
js <- readLines("oxcal_demo.js")
posterior <- get.OxCal(js, posterior = TRUE)
likelihood <- get.OxCal(js, posterior = FALSE)
allprob <- rbind(posterior, likelihood)

head(allprob)
```

```
##      lab depth    date     prob               comment posterior
## 1 ocd[9] 535.5 25544.5 0.000910 "Dep535.5 Posterior "      TRUE
## 2 ocd[9] 535.5 25539.5 0.000455 "Dep535.5 Posterior "      TRUE
## 3 ocd[9] 535.5 25534.5 0.000455 "Dep535.5 Posterior "      TRUE
## 4 ocd[9] 535.5 25529.5 0.003095 "Dep535.5 Posterior "      TRUE
## 5 ocd[9] 535.5 25524.5 0.001729 "Dep535.5 Posterior "      TRUE
## 6 ocd[9] 535.5 25519.5 0.002367 "Dep535.5 Posterior "      TRUE
```

Another function will summarise the age-depth model, extracting, by default, the 95.4% credibility interval and the median.


```r
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
```

```r
mod <- get.ageModel(x = allprob)
```

A final function plots the model


```r
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
```


```r
g <- plot.dates_model(allprob, mod, "My Core")
print(g)
```

![](Oxcal_r_demo_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

The red lozenges show the probability distribution function of the each date when calibrated individually. The blue lozenges show the posterior probability distribution functions of each date when modelled with the p-sequence model.

This code does not extract all the important information from OxCal, but it is a start. 


