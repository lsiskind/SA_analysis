%% Solar Panel Analysis - Full Data Set
% Lena Siskind, Instrument Systems 382A JPL
%%

close all
clear
clc

% read in CSV
M = readtable('iss_bad_data_joined_365.csv');

% seperate time stamp from angle data
timestamp = M(:,2);
timestamp = table2array(timestamp);

% convert timestamps to date-time format
DateTime = datetime(timestamp,'InputFormat','uuuu:DDD:HH:mm:ss');

% convert date-time to number of days since Jan 0, 0000
DateNum = datenum(DateTime);

%start time at 0
timevec = DateNum - DateNum(1).*ones(size(DateNum));

% timetag
t_course = table2array(M(:,3)); % [s]
t_fine = table2array(M(:,4));

% attitude data
posX = table2array(M(:,5)); % [ft]
posY = table2array(M(:,6));
posZ = table2array(M(:,7));

% define SARJs
SSARJ = table2array(M(:,8));
PSARJ = table2array(M(:,9));

% define beta angles per array
P44A_BGA = table2array(M(:,10));
P42A_BGA = table2array(M(:,11));
P42B_BGA = table2array(M(:,12));
P64B_BGA = table2array(M(:,13));
S43A_BGA = table2array(M(:,14));
S41A_BGA = table2array(M(:,15));
S61B_BGA = table2array(M(:,16));
S63B_BGA = table2array(M(:,17));

% time instances and alpha angles of intrusions into EMIT FOV
t2A_fov = timevec(PSARJ >= 268 & PSARJ <= 276);
PSARJ_2Afov = PSARJ(PSARJ >= 268 & PSARJ <= 276);

t4A_fov = timevec(PSARJ >= 89 & PSARJ <= 96);
PSARJ_4Afov = PSARJ(PSARJ >= 89 & PSARJ <= 96);

% time instances and alpha angles of intrusions into stray light KOZ envelope
t2A_koz = timevec((PSARJ>=255 & PSARJ<=289) | (PSARJ == 270 & (P42A_BGA >=24 & P42A_BGA <= 150)) | (PSARJ == 270 & (P42A_BGA >=229 & P42A_BGA <= 337))); % FIX
PSARJ_2Akoz = PSARJ((PSARJ>=255 & PSARJ<=289) | (PSARJ == 270 & P42A_BGA >=24 & P42A_BGA <= 150) | (PSARJ == 270 & (P42A_BGA >=229 & P42A_BGA <= 337)));

t4A_koz = timevec((PSARJ>=75 & PSARJ<=110) | (PSARJ == 90 & P44A_BGA >=41 & P44A_BGA<=150) | (PSARJ == 90 & P44A_BGA>=229 & P44A_BGA<=337));
PSARJ_4Akoz = PSARJ((PSARJ>=75 & PSARJ<=110) |(PSARJ == 90 & P44A_BGA >=41 & P44A_BGA<=150) | (PSARJ == 90 & P44A_BGA>=229 & P44A_BGA<=337));

% Percent overlap
percFOV = ((sum(t2A_fov)+sum(t4A_fov))/sum(timevec))*100; % percentage of time
percKOZ = ((sum(t2A_koz) + sum(t4A_koz))/sum(timevec))*100;
perc4A = ((sum(t4A_fov) + sum(PSARJ_4Akoz)))/sum(timevec)*100;

% create patch for FOV intrusions
x = [timevec(1) timevec(end) timevec(end) timevec(1)];
yf2A_a = [268 268 276 276];
yf4A_a = [89 89 96 96];

% create patches for KOZ intrusions
yk2A_a = [255 255 289 289]; % relates to psarj
yk4A_a = [75 75 110 110];

% yk2A_b1 = [24 24 150 150]; % relates to beta angle, commented out because
% only relevent based on sarj
% yk2A_b2 = [229 229 337 337];
% yk4A_b1 = [41 41 149 149]; 
% yk4A_b2 = [230 230 338 338];

%% make plots

figure()
plot(timevec, PSARJ)
hold on
patch(x,yf2A_a,'red','FaceAlpha',.2);
patch(x,yf4A_a,'green','FaceAlpha',.2);
legend('PSARJ Angles','Region of P4-2A Intrusion (268^o - 276^o)','Region of P4-4A Intrusion (89^o - 96^o)');
xlabel('Time [days]');
ylabel('\alpha [deg]');
title(['ISS SAWs intrude into EMIT FOV ' num2str(percFOV) '% of 2019']);
hold off

figure()
plot(timevec, PSARJ)
hold on
% patch(x,yk2A_a,'red','FaceAlpha',.2);
% patch(x, yk4A_a,'green','FaceAlpha',.2);
plot(t2A_koz, PSARJ_2Akoz)
plot(t4A_koz, PSARJ_4Akoz)
legend('All PSARJ Angles','P4-2A Intrusion into EMIT KOZ','P4-4A Intrusion into EMIT KOZ');
ylabel('\alpha [deg]');
xlabel('Time [days]');
title(['ISS SAWs intrude into EMIT KOZ ' num2str(percKOZ) '% of 2019']);
hold off