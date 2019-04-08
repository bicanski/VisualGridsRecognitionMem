

% Code developed by Andrej Bicanski
% andrej.bicanski@gmail.com
% www.andrejbicanski.com
% 
% model published in Current Biology
%
% Bicanski A, Burgess N. - A computational model of recognition 
% memory via grid cells. Current Biology, 2019, 29, 1–12. 
% DOI: 10.1016/j.cub.2019.01.077
%
% This script creates generic grid cell firing rate maps that are used as a
% lookup table


N_mod   = 9;                                                                          % 9 modules
Fmods   = [0.0014 0.002 0.0028 0.004 0.0055 0.0077 0.0108 0.0151 0.0211]*2*pi;    % frequencies separated by factor sqrt(2)
FAC     = [0.0638 0.0911  0.128   0.182  0.26   0.365   0.515  0.73 1.0429];          % scaling offset steps to cover the complete rombus
ORIs    = 0*[0  pi/3  pi/4  pi/2  pi/6 1.2*pi  1.7*pi  1.9*pi  pi/5];                 % pick random orientation for the modules
NperMod = 100;                                                                        % 100 offsets per module

res = 1;
Xmax = 440;
Ymax = 440;

X = [0:res:Xmax];   Nx = length(X);
Y = [0:res:Ymax];   Ny = length(Y);

GC_FRmaps    = zeros(length(X),length(Y),NperMod,N_mod);
GC_FRmaps_SD = zeros(round(length(X)/10),round(length(Y)/10),NperMod,N_mod);

[Xg,Yg] = meshgrid(X,Y);

Nx = length(Xg(:));

XY = [Xg(:), Yg(:)];

b0  = [ cos(0), sin(0) ];
b1  = [ cos(pi/3), sin(pi/3) ];
b2  = [ cos(2*pi/3), sin(2*pi/3) ];

for i = 1:length(ORIs)
    
    F   = Fmods(i);
    fac = FAC(i);
    ORI = ORIs(i);
    R   = [cos(ORI) -sin(ORI); sin(ORI) cos(ORI)];

    x_off_base1 = 0;
    y_off_base1 = fac*1/F;

    x_off_base2 = fac*(1/F)*cos(pi/6);
    y_off_base2 = fac*(1/F)*sin(pi/6);
    
    Off_vec1 = R*[x_off_base1, y_off_base1]';
    Off_vec2 = R*[x_off_base2, y_off_base2]';
    
    for w = 1:sqrt(NperMod)
        for j = 1:sqrt(NperMod)
            
            Off = (j-1)/10*Off_vec1' + (w-1)/10*Off_vec2';
            
            z0  = sum( repmat((R*b0')',Nx,1).*(F*XY+repmat(Off,Nx,1)) ,2);
            z1  = sum( repmat((R*b1')',Nx,1).*(F*XY+repmat(Off,Nx,1)) ,2);
            z2  = sum( repmat((R*b2')',Nx,1).*(F*XY+repmat(Off,Nx,1)) ,2);
            
            FRmap = cos(z0)+cos(z1)+cos(z2);
            FRmap = FRmap/(max(max(FRmap)));
            FRmap(FRmap<0)=0;
            reshFRmap = reshape(FRmap,[length(X),length(Y)]);
            reshFRmap = reshFRmap/max(max(reshFRmap));  % normalize to 1 becasue these are rates
            
            GC_FRmaps(:,:,(j-1)*sqrt(NperMod)+w,i) = reshFRmap;
            
            tmp = zeros(round(length(X)/10),round(length(Y)/10));
            for k = 1:Xmax/10
                for l = 1:Ymax/10
                    tmp(k,l) = mean(mean(reshFRmap( 5+(k-1)*10-4:5+(k-1)*10+4 , 5+(l-1)*10-4:5+(l-1)*10+4 )));
                    % downsampling after peaks have been realiably
                    % identified in high res version, but didn't end up
                    % using this
                end
            end
            GC_FRmaps_SD(:,:,(j-1)*sqrt(NperMod)+w,i) = tmp/max(max(tmp));  % normalize to 1 becasue these are rates
        end
    end
end

% reduced resolution version not actually used
save DOE_GC_FRmaps GC_FRmaps_SD GC_FRmaps
