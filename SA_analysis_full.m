%% Solar Panel Analysis - Full Data Set
% Lena Siskind, Instrument Systems 382A JPL
%%
close all
clear
clc

%% import data 
% import from two-week and year-long data set

% read in CSV from both data sets
M = readtable('SA_data_reduced_60x.csv'); % every minute for 1 year
M2 = readtable('iss_bad_data_joined_365.csv'); % every sec for 2 weeks

% position data (ECI frame)
posX = table2array(M(:,4)); % [km]
posY = table2array(M(:,5));
posZ = table2array(M(:,6));

posX2 = table2array(M2(:,5)); % [km]
posY2 = table2array(M2(:,6));
posZ2 = table2array(M2(:,7));

% define SARJs (solar alpha rotary joint)
SSARJ = table2array(M(:,7)); % [deg] 
PSARJ = table2array(M(:,8)); % port

SSARJ2 = table2array(M2(:,8)); % [deg]
PSARJ2 = table2array(M2(:,9));

% define BGAs angles per array (beta gimbal assembly)
P44A_BGA = table2array(M(:,9)); % [deg]
P42A_BGA = table2array(M(:,10));
P42B_BGA = table2array(M(:,11));
P64B_BGA = table2array(M(:,12));
S43A_BGA = table2array(M(:,13));
S41A_BGA = table2array(M(:,14));
S61B_BGA = table2array(M(:,15));
S63B_BGA = table2array(M(:,16));

P44A_BGA_2 = table2array(M2(:,10)); % [deg]
P42A_BGA_2 = table2array(M2(:,11));
P42B_BGA_2 = table2array(M2(:,12));
P64B_BGA_2 = table2array(M2(:,13));
S43A_BGA_2 = table2array(M2(:,14));
S41A_BGA_2 = table2array(M2(:,15));
S61B_BGA_2 = table2array(M2(:,16));
S63B_BGA_2 = table2array(M2(:,17));

%% convert timestamps

% grab time stamps
timestamp = M(:,1);
timestamp = table2array(timestamp);

timestamp2 = M2(:,2);
timestamp2 = table2array(timestamp2);

% convert timestamps to date-time format
DateTime = datetime(timestamp,'InputFormat','uuuu:DDD:HH:mm:ss');
DateTime2 = datetime(timestamp2,'InputFormat','uuuu:DDD:HH:mm:ss');

% convert date-time to number of days since Jan 0, 0000
DateNum = datenum(DateTime);
DateNum2 = datenum(DateTime2);

% set start time to 0 days, 0 seconds
timevec = DateNum - DateNum(1).*ones(size(DateNum));
timevec2 = DateNum2 - DateNum2(1).*ones(size(DateNum2));

%% Convert position to lat/long

% 365 day data
r = [posX posY posZ];
lla = eci2lla(r, utc);

lat = lla(:,1);
long = lla(:,2);
h = lla(:,3);

% two week data
r2 = [posX posY posZ];
lla2 = eci2lla(r, utc);

lat2 = lla(:,1);
long2 = lla(:,2);
h2 = lla(:,3);

%% interpolate data









%% interpolate data sets
% try interpolation over the first two weeks. 

% cut 365 array down to two weeks
marker = ismembertol(timevec, timevec2(end), .0002);
index = find(marker==1, 1);
timevec_short = timevec(1:index);
PSARJ_short = PSARJ(1:index);

% interpolate using matlab function

%% find intrusions into FOV and KOZ
% time instances and alpha angles of intrusions into EMIT FOV
t2A_fov_2 = timevec2(PSARJ2 >= 268 & PSARJ2 <= 276);
PSARJ_2Afov_2 = PSARJ2(PSARJ2 >= 268 & PSARJ2 <= 276);

t4A_fov_2 = timevec2(PSARJ2 >= 89 & PSARJ2 <= 96);
PSARJ_4Afov_2 = PSARJ2(PSARJ2 >= 89 & PSARJ2 <= 96);

% time instances and alpha angles of intrusions into stray light KOZ envelope
t2A_koz_2 = timevec2(PSARJ2>=255 & PSARJ2<=289); 
PSARJ_2Akoz_2 = PSARJ2(PSARJ2>=255 & PSARJ2<=289); 

t4A_koz_2 = timevec2(PSARJ2>=75 & PSARJ2<=110); 
PSARJ_4Akoz_2 = PSARJ2(PSARJ2>=75 & PSARJ2<=110); 

% time instances of beta intrusions into koz
t2A_koz_b_2 = timevec2((PSARJ2 == 270 & (P42A_BGA_2 >=24 & P42A_BGA_2 <= 150)) | (PSARJ2 == 270 & (P42A_BGA_2 >=229 & P42A_BGA_2 <= 337)));
BGA_2Akoz_2 = P42A_BGA_2((PSARJ2 == 270 & (P42A_BGA_2 >=24 & P42A_BGA_2 <= 150)) | (PSARJ2 == 270 & (P42A_BGA_2 >=229 & P42A_BGA_2 <= 337)));

t4A_koz_b_2 = timevec2((PSARJ2 == 90 & (P44A_BGA_2 >=41 & P44A_BGA_2<=149)) | (PSARJ2 == 90 & (P44A_BGA_2>=230 & P44A_BGA_2<=338)));
BGA_4Akoz_2 = P44A_BGA_2((PSARJ2 == 90 & (P44A_BGA_2 >=41 & P44A_BGA_2 <= 149)) | (PSARJ2 == 90 & (P44A_BGA_2 >=230 & P44A_BGA_2 <= 338)));

% Percent overlap
percFOV_2 = ((sum(t2A_fov_2)+sum(t4A_fov_2))/sum(timevec2))*100; % percentage of time
percKOZ_2 = ((sum(t2A_koz_2) + sum(t4A_koz_2))/sum(timevec2))*100;

% create patch for FOV intrusions
x_2 = [timevec2(1) timevec2(end) timevec2(end) timevec2(1)];
yf2A_a_2 = [268 268 276 276];
yf4A_a_2 = [89 89 96 96];

% create patches for KOZ intrusions
yk2A_a_2 = [255 255 289 289]; % relates to psarj
yk4A_a_2 = [75 75 110 110];

% yk2A_b1 = [24 24 150 150]; % relates to beta angle, commented out because
% only relevent based on sarj
% yk2A_b2 = [229 229 337 337];
% yk4A_b1 = [41 41 149 149]; 
% yk4A_b2 = [230 230 338 338];

%% Figure out position units

Re = 6378; % radius of earth km

% calculate altitude of ISS
r = [posX2 posY2 posZ2];
rmag = sqrt(sum(r.^2,2));
alt = rmag - Re; % alt=~410 km, true to iss orbit

%% Make plots

% plot orbit of ISS
figure()
plot3(posX2, posY2, posZ2);
hold on
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
surf(XS*Re, YS*Re, ZS*Re);
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title('Orbit of ISS Around Earth');
hold off

% plot PSARJ angles
figure()
plot(timevec2, PSARJ2)
hold on
patch(x_2,yf2A_a_2,'red','FaceAlpha',.2);
patch(x_2,yf4A_a_2,'green','FaceAlpha',.2);
legend('PSARJ Angles','Region of P4-2A Intrusion (268^o - 276^o)','Region of P4-4A Intrusion (89^o - 96^o)');
xlabel('Time [days]');
ylabel('\alpha [deg]');
title(['ISS SAWs intrude into EMIT FOV ' num2str(percFOV_2) '% of 2019']);
hold off

figure()
plot(timevec2, PSARJ2)
hold on
% patch(x,yk2A_a,'red','FaceAlpha',.2);
% patch(x, yk4A_a,'green','FaceAlpha',.2);
plot(t2A_koz_2, PSARJ_2Akoz_2)
plot(t4A_koz_2, PSARJ_4Akoz_2)
legend('All PSARJ Angles','P4-2A Intrusion into EMIT KOZ','P4-4A Intrusion into EMIT KOZ');
ylabel('\alpha [deg]');
xlabel('Time [days]');
title(['ISS SAWs intrude into EMIT KOZ ' num2str(percKOZ_2) '% of 2019']);
hold off

% beta angles
figure()
plot(timevec2, P42A_BGA_2);
hold on
plot(t2A_koz_b_2, BGA_2Akoz_2,'*');
legend('P4-2A BGA Angles','P4-2A intrusions into EMIT KOZ');
ylabel('\beta [deg]');
xlabel('Time [days]');
title('\alpha=270^o, Beta angles resulting in P4-2A SAW Intrusions into EMIT KOZ');
hold off

figure()
plot(timevec2, P44A_BGA_2);
hold on
plot(t4A_koz_b_2, BGA_4Akoz_2,'*');
legend('P4-4A BGA Angles','P4-4A intrusions into EMIT KOZ');
ylabel('\beta [deg]');
xlabel('Time [days]');
title('\alpha=90^o, Beta angles resulting in P4-4A SAW Intrusions into EMIT KOZ');
hold off