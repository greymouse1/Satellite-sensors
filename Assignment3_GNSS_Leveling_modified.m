
% Geoid modelling using Earth gravitational models (e.g. EGM2008)
% By: Mohammad Bagherbandi  2018-02-03
% Email:modbai@hig.se
% University of Gävle, Sweden
% Modified by N.G., Universitz of Gävle, Sweden
% Modification added in order for code run for all the values
% of nmax and save resulting tables in Excel file
clc;
clear all

% Define array of nmax values
nmax_values = [45,90,120,180,250,360,450,720,900,1200,2160];
% Loop through each nmax value
for nmax_idx = 1:length(nmax_values)
    nmax = nmax_values(nmax_idx);
    db = 180/nmax;

%nmax= 2160;  db = 180/nmax; % data resolution

phi_south = 52+db/2;
phi_north  = 72-db/2;
lambda_west = 2+db/2;
lambda_east = 34 -db/2;


phi_step = db; % block size in degree
lambda_step = db;

FEE= phi_south : phi_step : phi_north ;
LON= lambda_west : lambda_step : lambda_east;

     GM=0.3986004418E+15; 
     R=6378137;
 
      
phi_north=max(FEE);
phi_south=min(FEE);
lambda_west=min(LON);
lambda_east=max(LON);

% READING EGM08
% model='EGM08_2160x2160.txt';                            
% [C,S]=ReadEGM(model);     
load S_EGM08;  S = S_EGM08(1:nmax+1,1:nmax+1);
load C_EGM08;  C = C_EGM08(1:nmax+1,1:nmax+1);


%__________________________________________________________________________     
 % Subtract normal GRS80 coefficients from zonal GGMs coefficients 
  C(1,1)=0;
  C(3,1)=C(3,1)+4.841667749835220E-04;
  C(5,1)=C(5,1)-7.903037335100000E-07;
  C(7,1)=C(7,1)+1.687249611513883E-09;
  C(9,1)=C(9,1)-3.460524683954118E-12;
  C(11,1)=C(11,1)+2.650022257475667E-15; 
%__________________________________________________________________________


n=[0:nmax];
%Initial value of derived quantities
  T=0; 
for m=0:nmax
     %Call associate Legendre function
    P=plm(n,m,90-FEE);
    
    A=C(1:nmax+1,m+1)';
    B=S(1:nmax+1,m+1)';
   
        FUNC  =  ((P)*A')*cosd(m*LON)+((P)*B')*sind(m*LON);   T = T + FUNC;     % gravity anomaly
  end
 
    
   %Geoid
   N_EGM2008 = GM/R*T / 9.7;
                              
figure
h=worldmap([min(FEE),max(FEE)],[min(LON),max(LON)]);
title('Geoid Height using EGM08 model')
                setm(gca,'frame','on','ffill',600,...
                       'fedgecolor',[0 0 0],'flinewidth',3,'fontsize',24,...
                       'MeridianLabel','on','LabelFormat','none')
                gridm('GLineStyle','-','Gcolor','black')
                p = findobj(h,'type','patch'); % Find background
contourfm(FEE,LON,N_EGM2008,100,'LineStyle', 'none');
geoshow('landareas.shp','Facecolor','none'); 
tc=colorbar('Location','manual','Position',[0.25 0.68 0.012 0.25],'FontSize',22); 
set(get(tc,'xlabel'),'String','metre','FontSize',18); % colorbar label
colormap jet

%-----------------------------------------------------------------------------------------------


% Data format:
% No. 	Geodetic latitude	Longitude	
% GNSS height above the GRS 80 ellipsoid	GNSS-height	   
% Levelled normal height in RH 2000 (Official EVRS realisation) First order levelling


load gnsslev12_LM_SWE_140213.txt
fi_gnss=gnsslev12_LM_SWE_140213(:,2);    lambda_gnss = gnsslev12_LM_SWE_140213(:,3); 
h = gnsslev12_LM_SWE_140213(:,4);  % GNSS leveling
H = gnsslev12_LM_SWE_140213(:,5);  % Levelled_normal_height_RH2000

N_GNSS_Leveling = h - H;   % Geoid determination using h and H.
% N_accurate_geoid_LM = gnsslev12_LM_SWE_140213(:,6);

% determining the Geoid height for corresponding points in gnsslev12_LM_SWE_140213.txt
N_GNSS_Leveling_EGM2008 = interp2(LON,FEE,N_EGM2008,lambda_gnss,fi_gnss);

% Difference between the geoid from EGM08 and GNSS leveling
DN = N_GNSS_Leveling - N_GNSS_Leveling_EGM2008;


[r c]=size(DN);  RMS2 = (DN).^2; RMS3 = sum(RMS2(:))/(r*c); RMS = sqrt(RMS3);
fprintf('  Difference between the geoid from EGM08 and GNSS leveling\n  ');
 fprintf('Max =%6.5g         mean =%6.4g        Min =%6.5g         std =%6.4g    rms=%6.4g  metre\n',...
     max(DN(:)),mean(DN(:)),min(DN(:)),std(DN(:)),RMS);fprintf(' \n');


 figure
                h=worldmap([min(FEE),max(FEE)],[min(LON),max(LON)]);
                title('Geoid based on GNSS leveling points')
                setm(gca,'frame','on','ffill',600,...
                       'fedgecolor',[0 0 0],'flinewidth',3,'fontsize',24,...
                       'MeridianLabel','on','LabelFormat','none')
                gridm('GLineStyle','-','Gcolor','black')
                p = findobj(h,'type','patch'); % Find background
                set(p,'FaceColor',[1 1 1]); % Change background to white
                land = shaperead('landareas.shp', 'UseGeoCoords', true);
                geoshow(land, 'FaceColor', [0.7 0.7 0.7])   
                scatterm(fi_gnss,lambda_gnss,50, N_GNSS_Leveling,'filled','marker','o'); % plot scatter dots
                tc=colorbar('Location','manual','Position',[0.35 0.68 0.012 0.25],'FontSize',22); 
                set(get(tc,'xlabel'),'String','metre','FontSize',18); % colorbar label
                colormap jet
 
 
 figure
                h=worldmap([min(FEE),max(FEE)],[min(LON),max(LON)]);
                title('Difference between the geoid from EGM08 and GNSS leveling')
                setm(gca,'frame','on','ffill',600,...
                       'fedgecolor',[0 0 0],'flinewidth',3,'fontsize',24,...
                       'MeridianLabel','on','LabelFormat','none')
                gridm('GLineStyle','-','Gcolor','black')
                p = findobj(h,'type','patch'); % Find background
                set(p,'FaceColor',[1 1 1]); % Change background to white
                land = shaperead('landareas.shp', 'UseGeoCoords', true);
                geoshow(land, 'FaceColor', [0.7 0.7 0.7])   
                scatterm(fi_gnss,lambda_gnss,50, DN,'filled','marker','o'); % plot scatter dots
                tc=colorbar('Location','manual','Position',[0.35 0.68 0.012 0.25],'FontSize',22); 
                set(get(tc,'xlabel'),'String','metre','FontSize',18); % colorbar label
                colormap jet


           
% Save DN in Excel file
    excel_filename = sprintf('DN_current_nmax_%d.xlsx', nmax);
    writematrix(DN, excel_filename);

    % Display information
    fprintf('Results for nmax = %d saved in %s\n', nmax, excel_filename);
end
                


  


              


  




  



  


