library(foreach)
library(doSNOW)

initPositions <- function(x, numParticles, numStars, starRadius, centers)
{
  x = data.frame(star=rep(0, sum(numParticles)),
                 xPos=rep(0, sum(numParticles)),
                 yPos=rep(0, sum(numParticles)))
  
  totalPartsPlaced = 1;
  
  for(i in 1:numStars)
  {
    placedParticles = 1;
    
    while(placedParticles <= numParticles[i])
    {
      tempX = runif(1, -(starRadius[i]), starRadius[i])
      tempY = runif(1, -(starRadius[i]), starRadius[i])
      
      tempPos = c(tempX, tempY)
      
      if(sqrt(sum(tempPos^2)) <= starRadius[i])
      {
        x[totalPartsPlaced, ] <- c(i, tempPos[1] + centers[i, 1], tempPos[2] + centers[i, 2])
        placedParticles = placedParticles + 1
        totalPartsPlaced = totalPartsPlaced + 1
      }
    }
  }
  return(x)
}

initMasses <- function(m, numParticles, numStars, starMasses)
{
  m = data.frame(star=rep(0, sum(numParticles)),
                 mass=rep(0, sum(numParticles)))
  
  totalPartsPlaced = 1;
  
  for(i in 1:numStars)
  {
    placedParticles = 1;
    
    while(placedParticles <= numParticles[i])
    {
      partMass = starMasses[i]/numParticles[i]
      
      m[totalPartsPlaced, ] <- c(i, partMass)
      placedParticles = placedParticles + 1
      totalPartsPlaced = totalPartsPlaced + 1
    }
  }
  return(m)
}

initLambda <- function(lambda, numStars, starMass, starRadius, presureConstant, PolyIndex)
{
  lambda = c()
  
  for(i in 1:numStars)
  {
    lambda <- c(lambda,
                ((2*presureConstant * (pi^(-1/PolyIndex))) *
                   ((starMass[i])*(1+PolyIndex)/((starRadius[i])^2))^(1+(1/PolyIndex)))/
                  starMass[i])
  }
  
  return(lambda)
}

#kernal functions
kernel <- function(position, smoothingLength, dimensions)
{
  ch <- switch(
    dimensions,
    1/6* smoothingLength,
    5/(14 * pi * smoothingLength^2),
    1/(4 * pi * smoothingLength^3)
  )
  
  q<-(sqrt(sum(position^2)))/smoothingLength
  
  #Stops program if ch is null (ie not assigned by switch)
  stopifnot(!is.null(ch))
  
  if(is.na(q))
  {
    return(0)
  }
  
  if(q >= 0 && q < 1)
  {
    return (ch *((2-q)^3 - 4*(1-q)^3))
  } 
  else if(q >= 1 && q < 2)
  {
    return (ch * (2 -q)^3)
  }
  else if(q >= 2)
  {
    return(0)
  }
}

gradKernel <- function(position, smoothingLength,dimensions)
{
  ch <- switch(
    dimensions,
    1/(6* smoothingLength),
    5/(14 * pi * smoothingLength^2),
    1/(4 * pi * smoothingLength^3)
  )
  
  unitR <- position / (sqrt(sum(position^2)))
  q<-(sqrt(sum(position^2)))/smoothingLength
  
  #Stops program if ch is null (ie not assigned by switch)
  stopifnot(!is.null(ch))
  
  if(is.na(q))
  {
    return(0)
  }
  
  if(q >= 0 && q < 1)
  {
    return (ch * (1/smoothingLength) * (-3*(2-q)^2 + 12*(1-q)^2) * unitR)
  } 
  else if(q >= 1 && q < 2)
  {
    return (ch * (1/smoothingLength) * (-3*(2 - q)^2) * unitR)
  }
  else if(q >= 2)
  {
    return(0)
  }
}

#simulation Paramters
numParticles = c(200, 50)
totalParticles = sum(numParticles)
dimensions= 2
numStars = 2
starMass = c(1.6, .4)
starRadius = c(0.75,0.75)
smoothingLength = .04/sqrt(totalParticles/1000) #orginal .04/sqrt(numParticles/1000)
timeStep = .04
damping = 2
presureConstant = 0.1
PolyIndex = 1
maxTimeSetps <- 450

centers = data.frame(x = c(0, 2),
                     y = c(0, 0))

rho = data.frame(rep(0, totalParticles))

#placeholders which will be set with init methods
x <- 0
m <- 0
lambda = 0

v = data.frame(xVel=rep(0, totalParticles),
               yVel=rep(0, totalParticles))

accel = data.frame(xAccel=rep(0, totalParticles),
                   yAccel=rep(0, totalParticles))

#I think this needs to be calcuated, but wasnt included in the code
#for now ill just assume that in our problem all particles are at rest for t < 0
v_mhalf = data.frame(xVel=rep(0, totalParticles),
                     yVel=rep(0, totalParticles))

v_phalf = data.frame(xVel=rep(0, totalParticles),
                     yVel=rep(0, totalParticles))

x <- initPositions(x, numParticles,numStars, starRadius, centers)
m <- initMasses(m, numParticles, numStars, starMass)
lambda <- initLambda(lambda, numStars, starMass,  starRadius, presureConstant , PolyIndex)

mainNonParForeach <- function()
{
  
  foreach(i = 1:(totalParticles-1)) %do%
  {
    #"initialize density with i = j contribution"
    rho[i, 1] <- m[i, 2] * kernel(0, smoothingLength, dimensions)
    
    foreach(j=(i+1):totalParticles) %do%
    {
      #"calculate vector between two particles"
      uij = x[i,2:(dimensions+1)] - x[j,2:(dimensions+1)]
      rho_ij = m[i, 2] * kernel(uij, smoothingLength, dimensions)
      
      #"add contribution to density"
      rho[i, 1] <- rho[i, 1] + rho_ij 
      rho[j, 1] <- rho[j, 1] + rho_ij
    }
  }
  
  print(rho)
}

mainParForeach <- function()
{
  cl<-makeCluster(5) #change the 2 to your number of CPU cores
  registerDoSNOW(cl)
  
  #data frames are like reference types in java. (pointers)
  #print(tracemem(m))
  #mpar=m
  #print(tracemem(mpar))
  
  mpar = m[,]
  
  #foreach(i=1:(totalParticles-1), .verbose = TRUE) %dopar%
 # {
    #"initialize density with i = j contribution"
   # print(mpar)
    
 #   rho[i, 1] <- mpar[i, 2] * kernel(0, .04/sqrt(250/1000), 2)
  #}
  
  foreach(i=1:(totalParticles-1)) %do%
  {
    foreach(j=(i+1):totalParticles, .export = "x") %dopar%
    {
      #"calculate vector between two particles"
      uij = x[i,2:(dimensions+1)] - x[j,2:(dimensions+1)]
      rho_ij = mpar[i, 2] * kernel(uij, .04/sqrt(250/1000), 2)
      
      #"add contribution to density"
      rho[i, 1] <- rho[i, 1] + rho_ij 
      rho[j, 1] <- rho[j, 1] + rho_ij
    }
  }
  
  stopCluster(cl)
  
  print(rho)
}

mainNonParStand <- function()
{
  for(i in 1:(totalParticles-1))
  {
    #"initialize density with i = j contribution"
    rho[i, 1] <- m[i, 2] * kernel(0, smoothingLength, dimensions)
    
    for(j in (i+1):totalParticles)
    {
      #"calculate vector between two particles"
      uij = x[i,2:(dimensions+1)] - x[j,2:(dimensions+1)]
      rho_ij = m[i, 2] * kernel(uij, smoothingLength, dimensions)
      
      #"add contribution to density"
      rho[i, 1] <- rho[i, 1] + rho_ij 
      rho[j, 1] <- rho[j, 1] + rho_ij
    }
  }
  
  print(rho)
}
