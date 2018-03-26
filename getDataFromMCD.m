function [allValsOut,maxValsOut,minValsOut] = getDataFromMCD(targetVar,LSlon,LSlat,maxMinValsOnly) 
% getDataFromMCD Extracts data from the Mars Climate Database for a given
% MCD variable, longitude, and latitude

% John B McClean 
% Imperial College London 

% 26 March 2018: Initial version.

% data = array(month,lon,lat,time of day);
% L_s = 30*month - 15
data = zeros(12,64,49,12);

% Specify name of target variable
% ps    Surface pressure (Pa)
% tsurf Surface temperature (K)
% targetVar = 'ps';

% Read data in for target variable from MCD netCDF files
for i = 1:12
    s = num2str(i,'%02u');
    fileName = join(['clim_',s,'_me.nc']);
    disp(join(['Reading ',fileName,'...']));
    dataAppend = ncread(fileName,targetVar);
    switch targetVar
        case 'ps'
            data(i,:,:,:) = dataAppend;
        case 'tsurf'
            data(i,:,:,:) = dataAppend;
        case 'vmr_co2'
            data(i,:,:,:) = dataAppend(:,:,1,:); % use surface layer only
        case 'dustq'
            data(i,:,:,:) = dataAppend(:,:,1,:); % use surface layer only
        case 'rho'
            data(i,:,:,:) = dataAppend(:,:,1,:); % use surface layer only
        case 'reffdust'
            data(i,:,:,:) = dataAppend(:,:,1,:); % use surface layer only
    end
end
    disp('Done.');

clear dataAppend i s fileName


%% Set up interpolation

% Specify the longitude and latitude grid
lonVals = linspace(-180,180,65);
latVals = linspace(90,-90,49);
[X,Y] = meshgrid(lonVals,latVals);

% Set query points
Xq = LSlon;
Yq = LSlat;

% Get L_s values at midday on each sol from Excel file
L_s_midday = importfile('Sol_L_s_midday.xlsx','Sheet1',2,669);


%% Loop over sols

result = zeros(668,12);

if (maxMinValsOnly == 1)
   maxVals = zeros(668,1);
   minVals = zeros(668,1);
end

for thisSolNum = 1:668 % 668 sols in the year

    % Find out which data files to target using L_s_midday
    fileNumLo = floor((L_s_midday(thisSolNum)+15)/30);
    fileNumHi = fileNumLo + 1;
    
    L_s_lower = fileNumLo*30-15;
    L_s_upper = fileNumHi*30-15;
    
    if (fileNumLo == 0)
        fileNumLo = 12;
        L_s_lower = -15;
    end
    
    if (fileNumHi == 13)
       fileNumHi = 1; 
       L_s_upper = 375;
    end

    % Two columns, one for each bounding file
    valueOverThisSol = zeros(2,12);

    % Get spatially-interpolated value over one sol at lower L_s
    for hr = 1:12
        V(:,:) = data(fileNumLo,:,:,hr);
        Vextend = [V; V(1,:)]; % append an extra row so that we go from -180 to +180 deg
        valueOverThisSol(1,hr) = interp2(X,Y,Vextend',Xq,Yq);
    end

    % Get spatially-interpolated value over one sol at upper L_s
    for hr = 1:12
        V(:,:) = data(fileNumHi,:,:,hr);
        Vextend = [V; V(1,:)]; % append an extra row so that we go from -180 to +180 deg
        valueOverThisSol(2,hr) = interp2(X,Y,Vextend',Xq,Yq);
    end
    
    % Temporal interpolation
    
    % sample points - L_s values of the lower/upper files
    x = [L_s_lower,L_s_upper];
    
    for hr = 1:12
        % Sample values
        v = [valueOverThisSol(1,hr), valueOverThisSol(2,hr)];
        
        % Query point is the L_s at midday for this sol
        xq = L_s_midday(thisSolNum);
        
        result(thisSolNum,hr) = interp1(x,v,xq); 
    end
    
    if (maxMinValsOnly == 1)
        maxVals(thisSolNum,1) = max(result(thisSolNum,:));
        minVals(thisSolNum,1) = min(result(thisSolNum,:));
    end
    
end

result = result';
allValsOut = result(:);

maxValsOut = maxVals;
minValsOut = minVals;

end

