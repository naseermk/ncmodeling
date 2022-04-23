function [forces,gamma3,pressure] = Forcepara1(coords,rad,poisson,modulus,...
    nPart,receptor,ligand)
 
rng('shuffle'); 
etab=0.00001; %changing from 0.0001
 
% Initialize all forces to 0
forces = zeros(size(coords));
ftrac = zeros(size(coords));
fbias = zeros(size(coords));
gamma3 = zeros(size(rad));
pressure = zeros(size(rad));
% SRijHat = zeros(size(coords));
zmax = 6;
subpoisson =  0.5;
submodulus =  1.5*10^(-4);%150Pa the value is in MPa
receptorsub=0.9;
ligandsub=0.9;
% coordsLW(1,1)=1;
% coordsRW(1,1)=100;
% coordsBW(1,1)=10;

 a1=0;
 b1=1;
% Get the number of particles
%nPart = size(coords,2);
 
f = 0.000005;  %adhesion factor typical run f=0.000005
fsub = 0.00000925;  %substrate adhesion factor
 
for part=1:nPart
    %instead of looping over all pairs in the force calculation, find the distance between all
    %pairs
 
    %For a given cell, we create array dlist which has its distancefrom all other cells. We then
    %sort h= R_{part}+R_i-dlist. The neighbours are only those cells for
    %which h is >0. For calculating forces, we only take contact forces.
    %Thus this is an efficient method. 
    
    dlist= sqrt(sum(bsxfun(@minus, coords(:,part), coords).^2, 1));
    
    
    [d, ind] = sort(rad(part,1)+rad'-dlist);
    
    %begin_index=find(d==2*rad(part,1));
    begin_index=find(d>0.0);
    
    %ind_closest = ind(2:(min(nPart,26))); %find the n nearest neighbors
    ind_closest=ind((begin_index):(end));
    
    coords_closest = coords(:,ind_closest);
    %rad_closest=rad(ind_closest,1);
    %poission_closest=poisson(partA,1)
    for partA=1:(size(ind_closest,2))
        
        dr =  coords(:,part) - coords_closest(:,partA);
        
        if norm(dr) > 0.0
        
        Rij = norm(dr);
        
        RijHat = dr/Rij;
        
        hij = (rad(ind_closest(partA),1) + rad(part,1) - Rij);
        
        Eij = ((1 - poisson(ind_closest(partA),1)^2)/(modulus(ind_closest(partA),1)) ...
            + (1 - poisson(part,1)^2)/(modulus(part,1)));
        
        Rijf = (1/rad(part,1) + 1/rad(ind_closest(partA),1));
        
        %hij = max(0,h0);
        
        
        areaint = pi*(1/Rijf)*hij; %overlap area between 2 cells
        
        
        invDr2 = (hij)^(3/2); % 1/r^2
        
        forceFact = (invDr2/(0.75*(Eij)*sqrt(Rijf)))-areaint*f...
            *0.5*(receptor(ind_closest(partA),1)*ligand(part,1)+ ...
            receptor(part,1)*ligand(ind_closest(partA),1));
        
        % if (areaint > 0)
        pressure(part,1) = pressure(part,1)+ ...
            abs(forceFact/areaint);
        % end
        
%         dsub = rad(part,1)-(coords(3,part)-zmax);
%         if dsub >0 
%             ftrac(1:2,part)=2*((submodulus/modulus(part,1))^1.75)...
%                 *dsub*RijHat(1:2,1);
%             ftrac(3,part)=0;
%             
%         else
%             ftrac(:,part)=zeros(3,1);
%         end
        
        forces(:,part) = forces(:,part) + (forceFact*RijHat);%+ftrac(:,part);
        
         end 
        
    end
    
    %calculating gamma3 parameter.  In the paper by gernot

    
%     dxL =  coords(1,part) - coordsLW(1,1);
%     dxR =  coords(1,part) - coordsRW(1,1);
%     dyB =  coords(2,part) - coordsBW(1,1);
%     
%     %include repulsion due to boundaries 
%     forces(1,part) = forces(1,part) + (dxL*10^(-2))/(abs(dxL))^4 + ...
%        (dxR*10^-2)/(abs(dxR))^4;
%     
%     forces(2,part) = forces(2,part) + (dyB*10^-2)/(abs(dyB))^4 ;
    %boundaries at x=1, x=100, y=1
    
    for partA=1:(size(ind_closest,2))
        
        dr =  coords(:,part) - coords_closest(:,partA);
        
        if norm(dr)>0.0
        
        Rij = norm(dr);
        
        RijHat = dr/Rij;
        
        hij = (rad(ind_closest(partA),1) + rad(part,1) - Rij);
        
        %Eij = ((1 - poisson(partA,1)^2)/(modulus(partA,1)) ...
        %        + (1 - poisson(part,1)^2)/(modulus(part,1)));
        
        Rijf = (1/rad(part,1) + 1/rad(ind_closest(partA),1));
        
        %hij = max(0,h0);
        
        %  if (h0 > 0)
        
        mag=norm(forces(:,part));
        
        areaint = pi*(1/Rijf)*hij;
        
         gamma2 = etab*areaint*0.5*(1+(sum(bsxfun(@times, forces(:,part),RijHat)))/mag)...
            *0.5*(receptor(ind_closest(partA),1)*ligand(part,1)+ ...
            receptor(part,1)*ligand(ind_closest(partA),1));

        
        gamma3(part,1) = gamma3(part,1) + gamma2;
        
        end
      
        
    end
    
       dsub = rad(part,1)-(coords(3,part)-zmax);
    
        if dsub > 0
        
        Eijs = ((1 - subpoisson^2)/submodulus ...
            + (1 - poisson(part,1)^2)/(modulus(part,1)));
        
        Fsubrep = (4/3)*(1/Eijs)*(rad(part,1)^0.5)*dsub^(1.5);
        
        areasub = pi*rad(part,1)*dsub;
        
        Fsubatt = areasub*fsub*0.5*(receptorsub*ligand(part,1)+ ...
            receptor(part,1)*ligandsub);
         
         Fsubtot = Fsubrep - Fsubatt;    
         
         %generating random numbers between 0,1
               
                r31=(b1-a1).*rand(1) + a1;
                r41=(pi/2)*(b1-a1).*rand(1) + (pi/2)*a1;
                r51=2*pi*(b1-a1).*rand(1) + 2*pi*a1;
                           
%                 SRijHat(1,part) = sin(r41)*cos(r51);
%                 SRijHat(2,part) = sin(r41)*sin(r51);
%                 SRijHat(3,part) = cos(r41);
                SRijHat=[sin(r41)*cos(r51);sin(r41)*sin(r51);cos(r41)];%removing z-component cos(r41)
                
                %fbias(1,part)=1e-5*dsub; %1e-5*dsub before
                fbias(1,part)=0; %changing the bias force to zero. 
                fbias(2:3,part)=0;
                
        %    dsub = rad(part,1)-(coords(3,part)-zmax);
        %if dsub >0 
            ftrac(1:2,part)=((submodulus/modulus(part,1))^1.75)...
                *dsub*SRijHat(1:2,1);
            ftrac(3,part)=0;
            
        %else
            
        %end
        
        
            forces(:,part) = forces(:,part)+Fsubtot*SRijHat+ftrac(:,part);
        else
            fbias(:,part)=zeros(3,1);
            ftrac(:,part)=zeros(3,1);
        end 
    
    
end
 
 
end
