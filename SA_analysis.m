%% Solar Panel Analysis
% Lena Siskind, Instrument Systems 382A JPL
%% import data
close all
clear
clc

% read in CSV
M = readtable('SA_data_reduced_60x.csv');

% seperate time stamp from angle data
timestamp = M(:,1);
timestamp = table2array(timestamp);

% attitude data
posX = table2array(M(:,4)); % [miles ??]
posY = table2array(M(:,5));
posZ = table2array(M(:,6));

% define SARJs
SSARJ = table2array(M(:,7));
PSARJ = table2array(M(:,8));

% define beta angles per array
P44A_BGA = table2array(M(:,9));
P42A_BGA = table2array(M(:,10));
P42B_BGA = table2array(M(:,11));
P64B_BGA = table2array(M(:,12));
S43A_BGA = table2array(M(:,13));
S41A_BGA = table2array(M(:,14));
S61B_BGA = table2array(M(:,15));
S63B_BGA = table2array(M(:,16));

%% convert time format

% convert timestamps to date-time format
DateTime = datetime(timestamp,'InputFormat','uuuu:DDD:HH:mm:ss');

% convert date-time to number of days since Jan 0, 0000
DateNum = datenum(DateTime);

%start time at 0
timevec = DateNum - DateNum(1).*ones(size(DateNum));

%% interpolate
% might work when we can start filtering out data
% 
% % decrease step size to every second
% [rows, ~] = size(timevec);
% timevec_int = timevec(1):1/(rows*6):timevec(end); % step size of 10 seconds
% 
% % interpolate PSARJ and BGA angles for 2A and 4A
% PSARJ_int = interp1(timevec, PSARJ, timevec_int);
% P44A_BGA_int = interp1(timevec,P44A_BGA,timevec_int);
% P42A_BGA_int = interp1(timevec,P42A_BGA,timevec_int);

%% find intrusions into kov and fov

% time instances and alpha angles of intrusions into EMIT FOV
t2A_fov = timevec(PSARJ >= 268 & PSARJ <= 276);
PSARJ_2Afov = PSARJ(PSARJ >= 268 & PSARJ <= 276);

t4A_fov = timevec(PSARJ >= 89 & PSARJ <= 96);
PSARJ_4Afov = PSARJ(PSARJ >= 89 & PSARJ <= 96);

% time instances and alpha angles of intrusions into stray light KOZ envelope
t2A_koz_a = timevec(PSARJ>=255 & PSARJ<=289);
PSARJ_2Akoz = PSARJ(PSARJ>=255 & PSARJ<=289);

t4A_koz_a = timevec(PSARJ>=75 & PSARJ<=110); 
PSARJ_4Akoz = PSARJ(PSARJ>=75 & PSARJ<=110);

% time instances and beta angles of intrusions into KOZ
t2A_koz_b = timevec((PSARJ == 270 & (P42A_BGA >=24 & P42A_BGA <= 150)) | (PSARJ == 270 & (P42A_BGA >=229 & P42A_BGA <= 337)));
BGA_2Akoz = P42A_BGA((PSARJ == 270 & (P42A_BGA >=24 & P42A_BGA <= 150)) | (PSARJ == 270 & (P42A_BGA >=229 & P42A_BGA <= 337)));

t4A_koz_b = timevec((PSARJ == 90 & (P44A_BGA >=41 & P44A_BGA<=149)) | (PSARJ == 90 & (P44A_BGA>=230 & P44A_BGA<=338)));
BGA_4Akoz = P44A_BGA((PSARJ == 90 & (P44A_BGA >=41 & P44A_BGA<=149)) | (PSARJ == 90 & (P44A_BGA>=230 & P44A_BGA<=338)));

% % repeat beta angles with interpolated data
% t2A_koz_b_int = timevec_int((PSARJ_int == 270 & (P42A_BGA_int >=24 & P42A_BGA_int <= 150)) | (PSARJ_int == 270 & (P42A_BGA_int >=229 & P42A_BGA_int <= 337)));
% BGA_2Akoz_int = P42A_BGA_int((PSARJ_int == 270 & (P42A_BGA_int >=24 & P42A_BGA_int <= 150)) | (PSARJ_int == 270 & (P42A_BGA_int >=229 & P42A_BGA_int <= 337)));
% 
% t4A_koz_b_int = timevec_int((PSARJ_int == 90 & (P44A_BGA_int >=41 & P44A_BGA_int <= 149)) | (PSARJ_int == 90 & (P44A_BGA_int >= 230 & P44A_BGA_int <=338)));
% BGA_4Akoz_int = P44A_BGA_int((PSARJ_int == 90 & (P44A_BGA_int >=41 & P44A_BGA_int <= 149)) | (PSARJ_int == 90 & (P44A_BGA_int >=230 & P44A_BGA_int <=338)));


%% Percent overlap

percFOV = ((sum(t2A_fov)+sum(t4A_fov))/sum(timevec))*100; % percentage of time
percKOZ = ((sum(t2A_koz_a) + sum(t4A_koz_a))/sum(timevec))*100;
perc4A = ((sum(t4A_fov) + sum(PSARJ_4Akoz)))/sum(timevec)*100;


%% make plots

% plot orbit of ISS
figure()
Re = 3958.8; % [miles]
plot3(posX, posY, posZ);
hold on
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
surf(XS*Re, YS*Re, ZS*Re);
xlabel('x [miles]');
ylabel('y [miles]');
zlabel('z [miles]');
title('Orbit of ISS Around Earth');
hold off

% alpha angles
figure()
plot(timevec, PSARJ)
hold on
plot(t2A_fov, PSARJ_2Afov);
plot(t4A_fov, PSARJ_4Afov);
legend('PSARJ Angles','Region of P4-2A Intrusion (268^o - 276^o)','Region of P4-4A Intrusion (89^o - 96^o)');
xlabel('Time [days]');
ylabel('\alpha [deg]');
title(['ISS SAW PSARJs intrude into EMIT FOV ' num2str(percFOV) '% of 2019']);
hold off

figure()
plot(timevec, PSARJ)
hold on
plot(t2A_koz_a, PSARJ_2Akoz)
plot(t4A_koz_a, PSARJ_4Akoz)
legend('All PSARJ Angles','P4-2A Intrusion into EMIT KOZ','P4-4A Intrusion into EMIT KOZ');
ylabel('\alpha [deg]');
xlabel('Time [days]');
title(['ISS SAW PSARJ intrude into EMIT KOZ ' num2str(percKOZ) '% of 2019']);
hold off

% beta angles
figure()
plot(timevec, P42A_BGA);
hold on
plot(t2A_koz_b, BGA_2Akoz);
legend('P4-2A BGA Angles','P4-2A intrusions into EMIT KOZ');
ylabel('\beta [deg]');
xlabel('Time [days]');
title('\alpha=270^o, Beta angles resulting in P4-2A SAW Intrusions into EMIT KOZ');
hold off

figure()
plot(timevec, P44A_BGA);
hold on
plot(t4A_koz_b, BGA_4Akoz);
legend('P4-4A BGA Angles','P4-4A intrusions into EMIT KOZ');
ylabel('\beta [deg]');
xlabel('Time [days]');
title('\alpha=90^o, Beta angles resulting in P4-4A SAW Intrusions into EMIT KOZ');
hold off

% beta anlges with interpolated data
% figure()
% plot(timevec_int, P42A_BGA_int);
% hold on
% plot(t2A_koz_b_int, BGA_2Akoz_int);
% legend('P4-2A BGA Angles','P4-2A intrusions into EMIT KOZ');
% ylabel('\beta [deg]');
% xlabel('Time [days]');
% title({'\alpha=270^o, Beta angles resulting in P4-2A SAW Intrusions into EMIT KOZ'},{'Interpolated'});
% hold off
% 
% figure()
% plot(timevec_int, P44A_BGA_int);
% hold on
% plot(t4A_koz_b_int, BGA_4Akoz_int);
% legend('P4-4A BGA Angles','P4-4A intrusions into EMIT KOZ');
% ylabel('\beta [deg]');
% xlabel('Time [days]');
% title({'\alpha=90^o, Beta angles resulting in P4-4A SAW Intrusions into EMIT KOZ'},{'Interpolated'});
% hold off


