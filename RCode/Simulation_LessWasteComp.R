setwd("C:/Users/Johnny/Dropbox/Hood/DeptHonors/CommonEnvelopeSPH")

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
    print(((2*presureConstant * (pi^(-1/PolyIndex))) *
             ((starMass[i])*(1+PolyIndex)/((starRadius[i])^2))^(1+(1/PolyIndex)))/
      starMass[i])
    
    lambda <- c(lambda,
                ((2*presureConstant * (pi^(-1/PolyIndex))) *
                   ((starMass[i])*(1+PolyIndex)/((starRadius[i])^2))^(1+(1/PolyIndex)))/
                  starMass[i])
  }
  
  return(lambda)
}

calcDistanceMatrix <- function(x, totalParticles, dimensions)
{
  distVectorMatrix = array(dim=c(totalParticles, totalParticles, dimensions), rep(0,totalParticles *totalParticles* dimensions ))
  
  print(dim(distVectorMatrix))
  for(i in 1:(totalParticles-1))
  {
    for(j in (i+1):totalParticles)
    {
      #print(distVectorMatrix[i, j, 1])
      ix = x[i,2:(dimensions+1)]
      jx= x[j,2:(dimensions+1)]
      print(typeof(ix))
      distVectorMatrix[i, j, 1:2 ] = as.vector(ix-jx)
      #print(distVectorMatrix[i,j])
    }
  }
  
  return(distVectorMatrix)
}

#kernal functions
kernel <- function(position, smoothingLength, ch)
{
 # ch <- switch(
 #   dimensions,
 #   1/(6* smoothingLength),
 #   5/(14 * pi * smoothingLength^2),
 #   1/(4 * pi * smoothingLength^3)
 # )
  q<-(sqrt(sum((position)^2)))/smoothingLength
  
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

gradKernel <- function(position, smoothingLength,ch)
{
 #ch <- switch(
 #   dimensions,
 #   1/(6* smoothingLength),
 ##   5/(14 * pi * smoothingLength^2),
 #   1/(4 * pi * smoothingLength^3)
 # )
  
  unitR <- position / (sqrt(sum(position^2)))
  q<-(sqrt(sum(position^2)))/smoothingLength
  
  #Stops program if ch is null (ie not assigned by switch)
  if(is.null(ch))
  {
    print("null")
  }
  stopifnot(!is.null(ch)) #need to look up a better way to stop execution
  
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

#sudocode implemtation
calculate_density <- function(x, m, h, rho, numParticles, totalParticles, dimensions, ch, distVectorMatrix)
{
  for(i in 1:(totalParticles-1)){
    #"initialize density with i = j contribution"
    rho[i, 1] <- m[i, 2] * kernel(0, h, dimensions)
    
    for(j in (i+1):totalParticles)
    {
      #"calculate vector between two particles"
      uij = x[i,2:(dimensions+1)] - x[j,2:(dimensions+1)]
      
      #uij = distVectorMatrix[i, j, ]
      print(uij)
      rho_ij = m[i, 2] * kernel(uij, h, ch)
      
      #"add contribution to density"
      rho[i, 1] <- rho[i, 1] + rho_ij 
      rho[j, 1] <- rho[j, 1] + rho_ij
    }
  }
  
  return (rho)
}

calculate_Acceleration <- function(x, v, m, rho, P, nu, lambda, h, accel, numParticles, totalParticles, dimensions, ch, distVectorMatrix)
{
  #"add damping and gravity"
  accel <- data.frame(xAccel=rep(0, totalParticles),
                      yAccel=rep(0, totalParticles)) 
  if(dimensions == 3)
  {
    accel$zAccel <- rep(0, totalParticles)
  }

  #will need to make this into a looped statment to be cleaner
  accel[c(1:numParticles[1]), 1:dimensions] <- -nu * v[c(1:numParticles[1]), 1:dimensions] - lambda[1]* x[c(1:numParticles[1]), 1:dimensions]
  
  accel[c((numParticles[1] + 1): numParticles[2]), 1:dimensions] <- -nu * v[c((numParticles[1] + 1): numParticles[2]), 1:dimensions]
                                                                  - lambda[2]* x[c((numParticles[1] + 1): numParticles[2]), 1:dimensions]
  
  
  #not completly sure why it needs to be -1 but it ends up as a 101 row matrix and caues issues
  #for(i in 1:numParticles)
  #{
  #  accel[i, 1:dimensions] <-  -nu * v[i, 1:dimensions] - lambda[x[i, 1]] * x[i, 2:(dimensions+1)]
  #}
  
  #"add pressure"
  for(i in 1:(totalParticles-1))
  {
    pOverRSquare = P[i, 1]/(rho[i, 1])^2
    for(j in (i+1):totalParticles)
    {
      #"calculate vector between two particles"
      uij = x[i,2:(dimensions+1)] - x[j,2:(dimensions+1)]
      #uij = distVectorMatrix[1, j]
      #"calculate acceleration due to pressure"
      p_a = (-m[j, 2])*(pOverRSquare + P[j, 1]/(rho[j, 1])^2)*gradKernel(uij, h, ch)

      accel[i,] <- accel[i,] + p_a
      accel[j,] <- accel[j,] - p_a
      
    }
  }

  return(accel)
}

main <- function(){
  
  #simulation Paramters
  numParticles = c(120, 30)
  totalParticles = sum(numParticles)
  dimensions= 2
  numStars = 2
  starMass = c(1.5, .5)
  starRadius = c(0.75,0.75)
  smoothingLength = .04/sqrt(totalParticles/1000) #orginal .04/sqrt(numParticles/1000)
  timeStep = .04
  damping = 2.0
  presureConstant = 0.1
  PolyIndex = 1
  maxTimeSetps <- 250
  profilingTimeSteps <- 100
  
  
  ch = switch(
    dimensions,
    1/(6* smoothingLength),
    5/(14 * pi * smoothingLength^2),
    1/(4 * pi * smoothingLength^3)
  )
  
  centers = data.frame(x = c(0, 2),
                       y = c(0, 0))
  
  rho = data.frame(rep(0, totalParticles))
  
  #placeholders which will be set with init methods
  x = 0
  m = 0
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
  
  if(dimensions == 3)
  {
    x$zPos <- runif(totalParticles, -starRadius, starRadius)
    v$zVel <- runif(totalParticles, -0.25, .25)
    accel$zAccel <- rep(0, totalParticles)
    v_mhalf$zVel <- zVel=runif(totalParticles, -0.25, .25)
    v_phalf <- rep(0, totalParticles)
  }
  
  
  print(paste("Start ",Sys.time()))

  x = initPositions(x, numParticles,numStars, starRadius, centers)
  m = initMasses(m, numParticles, numStars, starMass)
  lambda = initLambda(lambda, numStars, starMass,  starRadius, presureConstant , PolyIndex)
  
  print(paste("Done initlizing particle positions at", Sys.time()))
  
  png(file = './OutputPlots/_Start.png')
  plot(x$xPos, x$yPos, xlim = c(-1, 3), ylim = c(-2, 2))
  dev.off()

  print("Starting main loop")

  #for(i in 1:maxTimeSetps)
  for(i in 1:1)
  {
    v_phalf = v_mhalf + (accel * timeStep)
    x[c(2,3)] = x[c(2,3)] + v_phalf * timeStep
    v = .5 * (v_mhalf + v_phalf)
    v_mhalf = v_phalf
  
    #"update densities, pressures, accelerations"
    #distVectorMatrix = calcDistanceMatrix(x, 10, dimensions)

    rho = calculate_density(x, m, smoothingLength, rho, numParticles, totalParticles, dimensions, ch, distVectorMatrix)

    print(rho)
    print(rho^(1+1/PolyIndex))
    
        P = presureConstant * rho^(1+1/PolyIndex)
    accel = calculate_Acceleration(x, v, m, rho, P, damping, lambda, smoothingLength, accel, numParticles, totalParticles, dimensions, ch, distVectorMatrix)
    
    png(file = paste("./OutputPlots/After", i,"loops.png", sep = ""))
    plot(x$xPos, x$yPos, xlim = c(-1, 3), ylim = c(-2, 2))
    dev.off()
    
    print(paste("Done loop", i, "at",Sys.time(), sep=" "))
  }

  print(paste("Done at", Sys.time()))

}
