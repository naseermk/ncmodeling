function [coords, L] = initCubicGrid(nPart,density)

% Initialize with zeroes
coords = zeros(3,nPart);

% Get the cooresponding box size
L = 100;

c1=10;
c2=45;

% Find the lowest perfect cube greater than or equal to the number of
% particles
%nCube = 2;

%         while (nCube^3 < nPart)
%             nCube = nCube + 10;
%         end


% Start positioning - use a 3D index for counting the spots
index = [0,0,0]';

% Assign particle positions
%need to change this to take into account birth position + distance
parameter=1;
 for part=1:nPart
parameter=parameter+1;
% Set coordinate
%scale=10;
%coords(:,part) = (index+[2*abs(randn),2*abs(randn),2*abs(randn)]');
   coords(:,part) = (index+[(c2-c1)*rand+c1,(c2-(c1+5))*rand + c1+5,10]');
 %coords(:,part) = (index+[2+parameter,2+parameter,2+parameter]');
                                                        
                                               
end
