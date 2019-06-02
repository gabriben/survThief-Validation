# source: https://stackoverflow.com/questions/46877451/generating-tabs-in-r-markdown-with-a-loop


library(survThief)
functions <- lsf.str("package:survThief")
# for (i in functions) print(eval(parse(text=paste0("survThief:::", i))))
setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/March Presentation")
# Function to create multiple tabs
make.tabs <- function(functions){
  res <- NULL
  for(i in functions){
    res <- c(res, '#### ', i, ' {-}', '\n',
             "```{r, comment=NA}", '\n',
             "eval(parse(text=paste0('survThief:::', '", i, "')))", '\n',
             '```', '\n\n')
  }
  return(res)
}

# print the code from package
make.tabs <- function(functions){
  res <- NULL
  for(i in functions){
    f <- eval(parse(text=paste0('survThief:::', i)))
    res <- c(res, '#### ', i, ' {-}', '\n',
             "```{r, echo = T, eval=FALSE}", '\n',
            as.character(body(f)), '\n',
             '```', '\n\n')
  }
  return(res)
}

# print the code from package dir
setwd("/Users/gabrielbenedict/Google_Drive/docs/UNIS/KU Leuven/Courses/Master Thesis/Digitize Survival Curves in sections/Code/survThief/R/")
functions <- list.files(pattern = ".R")

make.tabs <- function(functions){
  res <- NULL
  for(i in functions){
    res <- c(res, '#### ', i, ' {-}', '\n',
             "```{r, echo = T, eval=FALSE}", '\n',
             paste(paste0(capture.output(print(source(functions[1])$value)), sep="\n"), collapse = ", "), '\n',
             '```', '\n\n')
  }
  return(res)
}

# Create the Rmd to knit
cat(
  '

',
  make.tabs(functions),
  sep = "",
  file = "survThiefFunctions.Rmd")

# Render the Rmd created into html here
# rmarkdown::render("survThiefFunctions.Rmd")
