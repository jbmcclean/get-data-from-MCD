% matlabNetCDFreadTest Used with getDataFromMCD to plot values of various
% parameters for various landing sites

% John B McClean 
% Imperial College London

% 26 March 2018: Initial version.

%% Clear command window and workspace
% Aim: to generate a full Mars year plot of all atmospheric parameters
% relevant to MOXIE flow rate (P, T, mass fraction of CO2 in atmosphere,
% etc...)

%clear;
%clc;

%% Define constants
R = 8.31446; % gas constant, kg m^2 s^-2 K^-1 mol^-1
n_CO2 = 2; % number of moles of CO2 consumed
n_O2 = 1; % number of moles of O2 produced
M_O2 = 0.032; % molecular mass of O2, kg mol^-1
M_CO2 = 44.01; % molar mass of CO2, g mol^-1
M_atm = 43.34; % mean molar mass of Mars atmosphere, g mol^-1
rho_d = 2700; % density of dust, kg m^-3;

const = R*n_CO2/(n_O2*M_O2);

%% Define MOXIE parameters

eta = 0.5; % conversion efficiency
m_dot_O2_g_hr = 22; % O2 production rate, g hr^-1
A_F = 0.189; % filter area, m^2
T_op = 1; % operation time, hours

m_dot_O2 = m_dot_O2_g_hr/(1000*60*60);  % O2 production rate, kg s^-1

%% Define landing site coordinates

%Set target site coordinates
%        JEZERO   NE SYRT  COL HILLS
LSlon = [77.5031, 77.1599, 175.6255]; % deg E
LSlat = [18.4386, 17.8899, -14.5478]; % deg N

%% Options
maxMinValsOnly = 1; % Only get max and min values rather than every 2-hourly

%% Prepare for MCD data import

P           = zeros(3,8016);
T           = zeros(3,8016);
w_CO2       = zeros(3,8016);
dustq       = zeros(3,8016);
rho         = zeros(3,8016);
reffdust    = zeros(3,8016);
n_dust      = zeros(3,8016);

Pmin        = zeros(3,668);
Tmin        = zeros(3,668);
w_CO2min    = zeros(3,668);
dustqmin    = zeros(3,668);
rhomin      = zeros(3,668);
reffdustmin = zeros(3,668);
n_dustmin   = zeros(3,668);

Pmax        = zeros(3,668);
Tmax        = zeros(3,668);
w_CO2max    = zeros(3,668);
dustqmax    = zeros(3,668);
rhomax      = zeros(3,668);
reffdustmax = zeros(3,668); 
n_dustmax   = zeros(3,668);

%% Import data from MCD
for LSNum = 1:3
    
    lat = LSlat(1,LSNum);
    lon = LSlon(1,LSNum);
    
    disp(join(['Getting data for landing site at ', num2str(lat,'%2f'), ' N ', num2str(lon,'%2f'), ' E...', ]));
    
    disp('Getting pressure data...');
        [P(LSNum,:),Pmax(LSNum,:),Pmin(LSNum,:)] = ...
            getDataFromMCD('ps',lon,lat,maxMinValsOnly);
        
    disp('Getting temperature data...');
        [T(LSNum,:),Tmax(LSNum,:),Tmin(LSNum,:)] = ...
            getDataFromMCD('tsurf',lon,lat,maxMinValsOnly);
        
    disp('Getting CO2 volume mixing ratio data...');
        [w_CO2(LSNum,:),w_CO2max(LSNum,:),w_CO2min(LSNum,:)] = ...
            getDataFromMCD('vmr_co2',lon,lat,maxMinValsOnly);
        
    disp('Getting dust mass mixing ratio data...');
        [dustq(LSNum,:),dustqmax(LSNum,:),dustqmin(LSNum,:)] = ...
            getDataFromMCD('dustq',lon,lat,maxMinValsOnly);
        
    disp('Getting atmospheric density data...');
        [rho(LSNum,:),rhomax(LSNum,:),rhomin(LSNum,:)] = ...
            getDataFromMCD('rho',lon,lat,maxMinValsOnly); 
    
    disp('Getting dust effective radius data...');
        [reffdust(LSNum,:),reffdustmax(LSNum,:),reffdustmin(LSNum,:)] = ...
            getDataFromMCD('reffdust',lon,lat,maxMinValsOnly);

end

%% Carry out conversions
% Convert CO2 gas volume mixing ratio to mass mixing ratio
w_CO2 = (M_CO2/M_atm)*w_CO2;
w_CO2min = (M_CO2/M_atm)*w_CO2min;
w_CO2max = (M_CO2/M_atm)*w_CO2max;
    
% Calculate dust number density
n_dust    = (3./(4*pi*(reffdust.^3)))   .*(rho./rho_d)   .*(dustq);
n_dustmin = (3./(4*pi*(reffdustmin.^3))).*(rhomin./rho_d).*(dustqmin);
n_dustmax = (3./(4*pi*(reffdustmax.^3))).*(rhomax./rho_d).*(dustqmax);
    
%% Calculate the volumetric flow rate required, dV/dt (m^3 s^-1)

V_dot_g    = zeros(3,8016);
V_dot_gmin = zeros(3,668);
V_dot_gmax = zeros(3,668);

% 2-hourly values (original)
figure(1);
for LSNum = 1:3
    V_dot_g(LSNum,:) = (const*m_dot_O2/eta)*(T(LSNum,:)./(P(LSNum,:).*w_CO2(LSNum,:)));
    plot(V_dot_g(LSNum,:));
    title('2-hourly volumetric flow rate (m^3/s)')
    hold on;
end
hold off;

% Max/min values
figure(2)
x = 1:668;
for LSNum = 1:3
    % Reshape the 2-hourly matrix into one column per sol
    V_dot_g_sol(LSNum,:,:) = reshape(V_dot_g(LSNum,:),[12,668]);

    % Find the maximum and minimum value in each sol
    V_dot_gmax(LSNum,:) = max(V_dot_g_sol(LSNum,:,:));
    V_dot_gmin(LSNum,:) = min(V_dot_g_sol(LSNum,:,:));
    
    % Plot the values
    plot(x,V_dot_gmax(LSNum,:),'r-');
    title('Min/max volumetric flow rate (m^3/s)')
    hold on;
    plot(x,V_dot_gmin(LSNum,:),'k-');
    x2 = [x, fliplr(x)];
    inBetween = [V_dot_gmin(LSNum,:), fliplr(V_dot_gmax(LSNum,:))];
    fill(x2, inBetween, 'k', 'FaceAlpha', 0.333);
end

%% Calculate the rate of dust mass ingestion, dM/dt (kg s^-1)

dMdt = zeros(3,8016);
dMdt_min = zeros(3,668); 
dMdt_max = zeros(3,668); 

figure(3);
for LSNum = 1:3
    dMdt(LSNum,:) = (4/(3*pi))*rho_d.*V_dot_g(LSNum,:).*n_dust(LSNum,:).*(reffdust(LSNum,:).^3);
    plot(dMdt(LSNum,:));
    title('2-hourly dust mass ingestion rate (kg/s)')
    hold on;
end

figure(4)
x = 1:668;
for LSNum = 1:3
    % Reshape the 2-hourly matrix into one column per sol
    dMdt_sol(LSNum,:,:) = reshape(dMdt(LSNum,:),[12,668]);

    % Find the maximum and minimum value in each sol
    dMdt_max(LSNum,:) = max(dMdt_sol(LSNum,:,:));
    dMdt_min(LSNum,:) = min(dMdt_sol(LSNum,:,:));
    
    % Plot the values
    plot(x,dMdt_max(LSNum,:),'r-');
    title('Min/max dust mass ingestion rate (kg/s)')
    hold on;
    plot(x,dMdt_min(LSNum,:),'k-');
    x2 = [x, fliplr(x)];
    inBetween = [dMdt_min(LSNum,:), fliplr(dMdt_max(LSNum,:))];
    fill(x2, inBetween, 'k', 'FaceAlpha', 0.333);
end

%% Calculate dust loading rate, dm/dt (g m^-2 hr^-1)

dmdt = zeros(3,8016);
dmdt_min = zeros(3,668);
dmdt_max = zeros(3,668);

figure(5);
for LSNum = 1:3
    %           Convert from kg s^-1        Divide by area
    %                to mg hr^-1                of filter
    %                     |                       |
    %                     V                       V
    dmdt(LSNum,:) = (60*60*1e6*dMdt(LSNum,:))./A_F;
    plot(dmdt(LSNum,:));
    title('2-hourly dust loading rate (mg m^{-2} hr^{-1})')
    hold on;
end

figure(6)
lineMarkers = {'s','o','^'};
x = 1:668;
for LSNum = 1:3
    % Reshape the 2-hourly matrix into one column per sol
    dmdt_sol(LSNum,:,:) = reshape(dmdt(LSNum,:),[12,668]);

    % Find the maximum and minimum value in each sol
    dmdt_max(LSNum,:) = max(dmdt_sol(LSNum,:,:));
    dmdt_min(LSNum,:) = min(dmdt_sol(LSNum,:,:));
    
    % Plot the values
    plot(x(1:55:end),dmdt_max(LSNum,1:55:end),'color','r','marker',lineMarkers{LSNum});
    title('Min/max dust loading rate (mg m^{-2} hr^{-1})')
    hold on;
    plot(x(1:55:end),dmdt_min(LSNum,1:55:end),'color','b','marker',lineMarkers{LSNum});
    x2 = [x, fliplr(x)];
%     inBetween = [dmdt_min(LSNum,:), fliplr(dmdt_max(LSNum,:))];
%     fill(x2, inBetween, 'k', 'FaceAlpha', 0.333);
end

legend({'Jezero','NE Syrtis','Columbia Hills'});
xlim([1,668]);

%% Repeat mass loading rate plot but with L_s x-axis (0, 30, ..., 360)

L_s_midday = importfile('Sol_L_s_midday.xlsx','Sheet1',2,669);

% Add top and tail to L_s_midday since 0 < L_s < 360 deg
L_s_midday = [zeros(1); L_s_midday];
L_s_midday(1) = L_s_midday(end)-360;
L_s_midday = [L_s_midday; zeros(1)];
L_s_midday(end) = L_s_midday(2)+360;

%dmdt_max = dmdt_max';
%dmdt_min = dmdt_min';

% Add top and tail to values v since 0 < L_s < 360 deg
v = [dmdt_max dmdt_min];
v = [zeros(1,6); v];
v(1,:) = v(end,:);
v = [v; zeros(1,6)];
v(end,:) = v(2,:);

x = L_s_midday;
xq = (0:30:360);

dmdt_L_s_out = zeros(length(xq),6);

for i = 1:6
dmdt_L_s_out(:,i) = interp1(x,v(:,i),xq);
end






