library(ggplot2)
setwd("~/ctlab/EA17/")
load("./work/GA1.rda")

vectorToInt <- function(x) {
  return(sum(2 ^ (0:(length(x) - 1)) * x))
}

chromToReal <- function(x){
  return(globalMin + vectorToInt(x) * globalDelta)
}

crossingOver <- function(x, y, probCross){
  if(runif(1) > probCross)
    return(list(x, y))
  crossPoint <- sample(1:(chromosomeSize-1), 1)
  retX <- c(x[1:crossPoint], y[(crossPoint+1):chromosomeSize])
  retY <- c(y[1:crossPoint], x[(crossPoint+1):chromosomeSize])
  return(list(retX, retY))
}

randomMutation <- function(x, probMutation){
  if(runif(1) > probMutation)
    return(x)
  position <- sample(1:chromosomeSize, 1)
  x[position] <- (1 + x[position]) %% 2
  return(x)
}

predicate <- function(iter, iterWOChange, mode){
  if(mode == "boptim"){
    return(iterWOChange>5)
  }else if(mode == "generation"){
    return(iter>100)
  }else{
    stop("mode must be \"boptim\" or \"generation\"")
  }
}

drawMe <- function(x, y, iter){
  df <- data.frame(x, y)
  optimumPoint <- df[df$y==(max(df$y)),]
  ggplot(df, aes(x, y)) + geom_point(colour = "blue") + xlim(globalMin, globalMax)+ 
    stat_function(fun=fitnessFunc) +
    labs(title=paste0("Iteration: ", iter))+ 
    theme(plot.title = element_text(hjust = 0.5))+
    geom_point(data=optimumPoint, aes(x, y), colour="green")
  ggsave(file=paste0("work/GA1/plot", iter, ".png"), width = 6, height = 4)
}

fitnessFunc <-  function(x) sin(x) / (1 + exp(-x)) + 1

chromosomeSize <- 15
globalMin <- 0.5
globalMax <- 10
globalDelta <- (globalMax - globalMin) / 2 ^ chromosomeSize

run <- function(sizeOfPopulation, probCross, probMutation, mode, save){
  start.time <- Sys.time()
  file.remove(list.files("work/GA1", full.names = TRUE))
  
  chroms <- lapply(1:sizeOfPopulation, 
                   function(x) sample.int(2, size=chromosomeSize, replace=TRUE) -1)
  xcords <- unlist(lapply(chroms, chromToReal))
  fits <- unlist(lapply(xcords, fitnessFunc))
  optimumFit <- max(fits)
  iter <- 0
  iterWOChange <- 1
  if(save){
    drawMe(xcords, fits, iter)
  }
  while(!predicate(iter, iterWOChange, mode)){
    reproduction <- sample(chroms, size=sizeOfPopulation, prob = fits, replace = TRUE)
    # crossing over
    hey <- lapply(1:(sizeOfPopulation / 2), function(x) 
      crossingOver(reproduction[[x * 2 - 1]], reproduction[[x * 2]], probCross))
    reproduction <- unlist(hey, recursive = FALSE)
    # mutation
    chroms <- lapply(reproduction, function(x) randomMutation(x, probMutation))
    xcords <- unlist(lapply(chroms, chromToReal))
    fits <- unlist(lapply(xcords, fitnessFunc))
    if(save){
      drawMe(xcords, fits, iter)
    }
    iter <- iter + 1
    iterWOChange <- iterWOChange + 1
    if(optimumFit != min(fits)){
      iterWOChange <- 1
      optimumFit <- min(fits)
    }
  }
  argopt <- 7.8543697
  err <- 1 - abs(argopt - xcords[which(fits == max(fits))[1]]) / (argopt - 0.5)
  return(c(err,
         iter,
         Sys.time()-start.time))
}

argt <- seq(10, 200, by=10)
argpc <- seq(0.25, 0.75, by=0.05)
argpm <- seq(0.0005, 0.004, by = 0.0005)
args <- expand.grid(argt, argpc, argpm)
results <- list()
for(i in 1:nrow(args)){
  results[[i]] <- rowMeans(replicate(1, 
                    run(args[1][i,], args[2][i,], args[3][i,], "boptim", FALSE)))
}

rrow <- lapply(1:length(results[[1]]), function(i) sapply(results, function(x) x[i]))
args$accuracy <- rrow[[1]]
args$iter <- rrow[[2]]
args$time <- rrow[[3]]

colnames(args) <- c("number", "pC", "pM", "accuracy", "iteration", "time")

tmp <- aggregate(x=args, by=list(args$number, args$pM), FUN = mean)
ggplot(tmp, aes(x = number, y=pM))+geom_point(aes(size=iteration, color=accuracy))+
  # geom_text(aes(label=format(accuracy, digits = 2)), hjust=0, nudge_x = 0, vjust=0, nudge_y = 0)+
  #  scale_size_continuous(range = c(1, 6)) + 
  labs(y = "probability of mutation", x = "size of population")+
  scale_colour_gradientn(colours = terrain.colors(11))
  # scale_color_gradient2(low="blue", mid="red", high="green") 
ggsave("./work/1popVSmut.png", width = 6, height = 4)  

save(args, file = "./work/GA1.rda")

