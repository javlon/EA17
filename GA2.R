library(plotly)
setwd("~/ctlab/EA17/")
load("./work/GA2.rda")

arithmeticCrossover <- function(x, y){
  if(runif(1) > probCross)
    return(list(x, y))
  retX <- x * w + y * (1 - w)
  retY <- x * (1 - w) + y * w
  return(list(retX, retY))
}
randVar <- function(){
  return(runif(1) * (maxX - minX) + minX)
}
randomMutation <- function(x){
  if(runif(1) <= probMutation)
    x[sample(1:length(x), 1)] <- randVar()
  return(x)
}

fitnessFunc <- function(x){
  -a * exp(-b * sqrt(sum(x^2) / length(x))) -
    exp(sum(cos(c * x)) / length(x)) + a + exp(1)
}

getAsDF <- function(chroms, fits){
  x1 <- sapply(chroms, function(x) x[1])
  x2 <- sapply(chroms, function(x) x[2])
  return(data.frame(x=x1,y=x2,z=fits))
}

create3DPlot <- function(){
  x <- y <- seq(minX, maxX, by=(maxX-minX)/50)
  z <- sapply(y, function(y){ sapply(x, function(z) fitnessFunc(c(y,z)))})
  return(plot_ly() %>% add_trace(x=~x, y=~y, z=~z, opacity=1, type = "surface"))
}

a <- 20 ; b <- .2; c <- 2 * pi;
minX <- -1
maxX <- 1
predicate <- function(iter, iterWOChange, mode){
  if(mode == "boptim"){
    return(iterWOChange>5)
  }else if(mode == "generation"){
    return(iter>100)
  }else{
    stop("mode must be \"boptim\" or \"generation\"")
  }
}

run <- function(sizeOfPopulation, probCross, probMutation, n, mode){
  start.time <- Sys.time()
  assign("w", .1, envir = .GlobalEnv) 
  assign("probCross", probCross, envir = .GlobalEnv) 
  assign("probMutation", probMutation, envir = .GlobalEnv) 
  # init
  chroms <- lapply(1:sizeOfPopulation, function(x) replicate(n, randVar()))
  fits <- unlist(lapply(chroms, fitnessFunc))
  if(n == 2){
    p <- create3DPlot()
    df <- getAsDF(chroms, fits)
    p %>% add_trace(df, x=~df$x, y=~df$y, z=~df$z, type="scatter3d", color=~df, colors = "green")
  }

  iter <- 0
  iterWOChange <- 1
  minimumFit <- min(fits)
  
  while(!predicate(iter, iterWOChange, mode)){
    prob <- max(fits) - fits
    if(max(prob) == 0)
      prob <- numeric(sizeOfPopulation) + 1
    # selection
    reproduction <- sample(chroms, size = sizeOfPopulation, replace = TRUE, prob = prob )
    # crossover
    hey <- lapply(1:(sizeOfPopulation/2), 
           function(x) arithmeticCrossover(reproduction[[2*x-1]], reproduction[[2*x]]))
    reproduction <- unlist(hey, recursive = FALSE)
    # mutation
    chroms <- lapply(reproduction, randomMutation)
    fits <- unlist(lapply(chroms, fitnessFunc))
    iter <- iter + 1
    iterWOChange <- iterWOChange + 1
    if(minimumFit != min(fits)){
      iterWOChange <- 1
      minimumFit <- min(fits)
    }
  }
  if(n == 2){
    df <- getAsDF(chroms, fits)
    p %>% add_trace(df, x=~df$x, y=~df$y, z=~df$z, type="scatter3d", color=~df, colors = "green")
  }
  #ans <- c(chroms[which(fits == minimumFit)[1]][[1]], minimumFit)
  err <- 1 -sqrt(sum(chroms[which(fits == minimumFit)[1]][[1]]^2))/sqrt(2)
  return(c(err, iter, (Sys.time()-start.time)))
}

argt <- seq(10, 190, by=20)
argpc <- seq(0.35, 0.65, by=0.05)
argpm <- seq(0.0005, 0.004, by = 0.0005)
args <- expand.grid(argt, argpc, argpm)
results <- list()
for(i in 1:nrow(args)){
  results[[i]] <- rowMeans(replicate(1, run(args[1][i,], args[2][i,], args[3][i,], 2, "boptim")))
}
rrow <- lapply(1:length(results[[1]]), function(i) sapply(results, function(x) x[i]))
hey <- cbind(rrow[[1]], rrow[[2]], rrow[[3]])
args$accuracy <- rrow[[1]]
args$iter <- rrow[[2]]
args$time <- rrow[[3]]

colnames(args) <- c("number", "pC", "pM", "accuracy", "iteration", "time")
tmp <- aggregate(x=args, by=list(args$number, args$pC), FUN = mean)
ggplot(tmp, aes(x = number, y=pC))+geom_point(aes(size=iteration, color=accuracy))+
  # scale_size_continuous(range = c(1, 6)) +
  labs(y = "probability of crossing over", x = "size of population")+
  scale_colour_gradientn(colours = terrain.colors(11))
  # scale_color_gradient(low="blue", high="red")
ggsave("./work/2popVScross.png", width = 6, height = 4)  


save(args, file = "./work/GA2.rda")

