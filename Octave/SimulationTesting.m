function SimulationTesting
  graphics_toolkit('gnuplot')
  page_screen_output(0);
  page_output_immediately(1);
  
  numStars = 2;
  numParticles = [50, 50];
  totalParticles = sum(numParticles(1:numStars));
  
  centers = [5, 0; 0, 0];
  starRadius = [1, 1];
  #starRadius = [0.00464913034, 0.00464913034] #sun radius = 0.00464913034 AU 
  #starMasses = [1.989*(10**12), 1.989*(10**12)];
  starMasses = [5,5];
  smoothingLength = 0.04/sqrt(totalParticles/1000);
  timeStep = 0.04;
  damping = 2;
  dimensions = 2;
  
  presureConstant = 0.1;
  polyIndex = 1;
  
  testingTimeSteps = 100;
  
  dataFile = "testFile.txt"
  
  switch(dimensions)
    case 1
      ch = (1/(6*smoothingLength));
    case 2
      ch = (5/(14 * pi * smoothingLength^2));
    case 3
      ch = (1/(4 * pi * smoothingLength^3));
    otherwise
        error("Dimensions must be 1, 2, or 3");  
  endswitch

  rho = zeros(1, totalParticles);
  P = zeros(1, totalParticles);
  
  v=zeros(totalParticles, dimensions);
  accel = zeros(totalParticles, dimensions);
  
  %for now ill just assume that in our problem all particles are at rest for t < 0
  %in the furture, we would need to 1 time step back which Eulers method
  vMHalf = zeros(totalParticles, dimensions);
  vPHalf = zeros(totalParticles, dimensions);
  
  %initilize the position, mass, and lambda matrices
  x = initPosition(numParticles, numStars, starRadius, dimensions, centers, false, dataFile);
  mass = initMass(numParticles, numStars, starMasses);
  lambda = initLambda(numStars, starMasses, starRadius, presureConstant, polyIndex)
  #lambda = ((2* presureConstant * (pi ^(-1/polyIndex))) * ((sum(starMasses(1:numStars)))* (1+polyIndex)/((0.04/sqrt(totalParticles/1000))^2)^(1+(1/polyIndex)))/sum(starMasses(1:numStars)));
  
  #savePositions(x, dataFile)
  disp(strcat("Done initlizing particle positions at-", strftime ("%T %x", localtime (time()))))
 
  #outputFolder = strcat('.\OutputPlots-', strftime("%T %x", localtime (time())));
  mkdir(strcat('OutputPlots\OutputPlots-', datestr(now(), 30)));
  outputFolder = strcat('OutputPlots\OutputPlots-', datestr(now(), 30));

  plot(x(1:numParticles(1),2), x(1:numParticles(1),3), '.b',
         x((numParticles(1)+1):numParticles(2)+numParticles(1),2), x((numParticles(1)+1):numParticles(2)+numParticles(1),3), '.r');
  #axis([-55, 5050, -150, 150]);
  #axis([-.025, .025, -.025, .025]);
  axis([-2, 7, -5, 5]);
  print(strcat(outputFolder , '\0Start.png'))
  
  disp("Starting main loop")
  for i = 1:testingTimeSteps
    #disp(x)
    #disp(v)
    #disp(accel)    
    
    vPHalf = vMHalf + (accel * timeStep);
    x(:, 2:(dimensions+1)) = x(:, 2:(dimensions+1)) + vPHalf * timeStep;
    v = 0.5 * (vMHalf+vPHalf);
    vMHalf = vPHalf;
    
    rho = calcDensity(x, mass, rho, smoothingLength,numParticles, totalParticles, dimensions, ch);
    P = presureConstant * (rho.^(1+1/polyIndex));
    accel = calcAccel(x, v, mass, rho, P, damping, lambda, smoothingLength, numParticles, totalParticles, dimensions, ch);
    #disp(accel)
    
    plot(x(1:numParticles(1),2), x(1:numParticles(1),3), '.b',
         x((numParticles(1)+1):numParticles(2)+numParticles(1),2), x((numParticles(1)+1):numParticles(2)+numParticles(1),3), '.r');
    #axis([-55, 5050, -150, 150]);
    #axis([-.025, .025, -.025, .025]);
    axis([-2, 7, -5, 5]);
    print(strcat(outputFolder , '\AfterLoop', num2str(i), '.png'))
    
    disp(strcat("Done loop-", num2str(i), " at ", strftime(" %T %x", localtime (time())))) 
  end
  
  #plot(x(1:numParticles(1),2), x(1:numParticles(1),3), '.b',
       #x((numParticles(1)+1):numParticles(2)+numParticles(1),2), x((numParticles(1)+1):numParticles(2)+numParticles(1),3), '.r');
  #axis([2490, 2505, -10, 10]); 
  #print(strcat(outputFolder , '\End.png'))  
  
  
endfunction

function x = initPosition(numParticles, numStars, starRadius, dimensions, centers, loadFile, fileString)
  
  x = zeros(sum(numParticles(1:numStars)), 1+dimensions);
  
  if(loadFile)
    x = dlmread(fileString,",");
  else
    totalPlacedParts = 1;
    for i=1:numStars
      
      placedParticiles = 1;
      while(placedParticiles <= numParticles(i))
      
        tempPos = unifrnd(-(starRadius(i)), starRadius(i), 1, 2);
        
        if(norm(tempPos, 2) <= starRadius(i))
          x(totalPlacedParts,:) = [i, centers(i,:) + tempPos];
          placedParticiles++;
          totalPlacedParts++; 
        end 
      end
    end
  end
endfunction

function savePositions(x, fileString)
  dlmwrite(fileString, x, ",");
endfunction 

function vel = initVel(numParticles, numStars, dimensions)
  
  vel = zeros(sum(numParticles), dimensions);
  
  totalPlacedParts = 1;
  for i=1:numStars
    for j= 1:(numParticles(i))
    
      massFrac = starMasses(i)/numParticles(i);
      m(totalPlacedParts, :) = [i, massFrac];
      
      totalPlacedParts++;
    end
  end
endfunction

function m = initMass(numParticles, numStars, starMasses)
  
  m = zeros(sum(numParticles), 2);
  
  totalPlacedParts = 1;
  for i=1:numStars
    for j= 1:(numParticles(i))
    
      massFrac = starMasses(i)/numParticles(i);
      m(totalPlacedParts, :) = [i, massFrac];
      
      totalPlacedParts++;
    end
  end
endfunction

function lambda = initLambda(numStars, starMass, starRadius, presureConstant, polyIndex)
  lambda = zeros(1, numStars);
  
  for i = 1:numStars
    lambda(i) =((2* presureConstant * (pi ^(-1/polyIndex))) * ((starMass(i))* (1+polyIndex)/((starRadius(i))^2))^(1+(1/polyIndex)))/starMass(i);
  end  
endfunction

function distVectorMatrix = calcDistanceMatrix(x, totalParticles, dimensions)
  
endfunction

function kernelOut = kernel(position, smoothingLength, ch)
  
  q = norm(position, 2)/smoothingLength;
  
  if(q >= 0 && q < 1)
    kernelOut = (ch * ((2-q)^3 - 4*(1-q)^3));
    
  elseif(q >= 1 && q < 2)
    kernelOut = (ch * (2-q)^3);
    
  elseif(q >=2)
    kernelOut = 0;
    
  endif  
endfunction

function gradKernelOut = gradKernel(position, smoothingLength, ch)
  posNorm = norm(position);
  
  unitR = position/posNorm;
  q = posNorm/smoothingLength;
  
  if(q >= 0 && q < 1)
    gradKernelOut = (ch * (1/smoothingLength) * (-3*(2-q)^2 + 12*(1-q)^2) * unitR);
    
  elseif(q >= 1 && q < 2)
    gradKernelOut = (ch * (1/smoothingLength) * (-3*(2-q)^2)* unitR);
    
  elseif (q >=2)
    gradKernelOut = 0;
    
  endif  
endfunction

function rho = calcDensity(x, mass, rho, smoothingLength, numParticles, totalParticles, dimensions, ch)
  for i = 1:totalParticles
    rho(i) = mass(i, 2) * kernel(0, smoothingLength, ch);
    
    for j = (i+1):totalParticles
      uij = x(i, 2:(dimensions+1)) - x(j, 2:(dimensions+1));
      %rho_ij = mass(i, 2) * kernel(uij, smoothingLength, ch);
      kernelOut = 0;
      q = norm(uij, 2)/smoothingLength;
  
      if(q >= 0 && q < 1)
        kernelOut = (ch * ((2-q)^3 - 4*(1-q)^3));
        
      elseif(q >= 1 && q < 2)
        kernelOut = (ch * (2-q)^3);
        
      elseif(q >=2)
        kernelOut = 0;
      endif  
       
      rho_ij = mass(i, 2) * kernelOut;   
           
      rho(i) = rho(i) + rho_ij;
      rho(j) = rho(j) + rho_ij;      
    end 
  end  
endfunction

function accel = calcAccel(x, v, mass, rho, P, nu, lambda, smoothingLength, numParticles, totalParticles, dimensions, ch)
  accel = zeros(totalParticles, dimensions);
  
  %need to figure out how to loop this.
  #accel(1:numParticles(1),:) = -nu * v(1:numParticles(1),:); - lambda(1) * x(1:numParticles(1), 2:(dimensions+1));
  #accel((numParticles(1) + 1):totalParticles,:) = -nu * v((numParticles(1) + 1):totalParticles,:); - lambda(2) * x((numParticles(1) + 1):totalParticles, 2:(dimensions+1));

   for i = 1:totalParticles
    accel(i, :) = accel(i, :) + -nu * v(i,:) - lambda(x(i,1)) * x(i, 2:(dimensions+1)); 
    #accel(i, :) = -nu * v(i,:) - lambda * x(i, 2:(dimensions+1)); 
  end
  
  #bigG = 6.674 * (10**-011);
  #bigG = 1.3268 * (10**20);
  #bigG = 8.825741331*(10**-25);
  #for i = 1:totalParticles
  #   for j = (i+1):totalParticles
  #     uij = x(i, 2:(dimensions+1)) - x(j, 2:(dimensions+1));
  #     
  #     #my gravity
  #     #forceOfGrav = ((bigG*mass(i, 2)*mass(j, 2)) / (norm(uij) * norm(uij)));
  #            
  #     #sph paper grav
  #     forceOfGrav = (bigG*mass(i, 2)*mass(j, 2))*norm(uij); #* norm(uij));
  #    
  #     accel(i, :) = accel(i, :) - ((forceOfGrav/ mass(i, 2)) * ((uij)/ norm(uij)));
  #     accel(j, :) = accel(j, :) + ((forceOfGrav/ mass(j, 2)) * ((uij)/ norm(uij)));
  #   end
  #end
  
  for i = 1:15
    #disp(accel(i, :))
  end  
 
  if(true)
    for i = 1:totalParticles
       pOverRSquare = P(i)/(rho(i))^2;
       
       for j = (i+1):totalParticles
        uij = x(i, 2:(dimensions+1)) - x(j, 2:(dimensions+1));
        %p_a = -mass(j, 2)*(pOverRSquare + P(j)/(rho(j))^2)*gradKernel(uij, smoothingLength, ch);
        
        gradKernelOut = 0;
        posNorm = norm(uij);
    
        unitR = uij/posNorm;
        q = posNorm/smoothingLength;
        
        if(q >= 0 && q < 1)
          gradKernelOut = (ch * (1/smoothingLength) * (-3*(2-q)^2 + 12*(1-q)^2) * unitR);
          
        elseif(q >= 1 && q < 2)
          gradKernelOut = (ch * (1/smoothingLength) * (-3*(2-q)^2)* unitR);
          
        elseif (q >=2)
          gradKernelOut = 0;
          
        endif 
       
        p_a = -mass(j, 2)*(pOverRSquare + P(j)/(rho(j))^2)* gradKernelOut;
        #p_a = -(mass(j, 2)/(3*(10**26)))*(pOverRSquare + P(j)/(rho(j))^2)* gradKernelOut;
        accel(i,:) = accel(i,:) + p_a;
        accel(j,:) = accel(j,:) - p_a;
      end 
    end
  end  
  
   for i = 190:210
    #disp(accel(i, :))
   end
  
endfunction


