%% Solar Panel Analysis
% Lena Siskind, Instrument Systems 382A JPL
%% import data
close all
clear
clc
tic

% read in CSV
M = readtable('SA_data_reduced_60x.csv');

% import timestamp
timestamp = M(:,1);
timestamp = table2array(timestamp);

% attitude data
posX = table2array(M(:,4)); % [km]
posY = table2array(M(:,5));
posZ = table2array(M(:,6));

% define SARJs (solar alpha rotary joints)
SSARJ = table2array(M(:,7)); % [deg]
PSARJ = table2array(M(:,8));

% define bgas per array (beta gimbly assembly)
P44A_BGA = table2array(M(:,9)); % [deg]
P42A_BGA = table2array(M(:,10));


%% convert time format

% convert timestamps to date-time format
DateTime = datetime(timestamp,'InputFormat','uuuu:DDD:HH:mm:ss');

% convert date-time to number of days since Jan 0, 0000
DateNum = datenum(DateTime);

%start time at 0 days, 0 seconds
timevec = DateNum - DateNum(1).*ones(size(DateNum));

% convert date-time to utc
[year, month, day, hour, mins, sec] = datevec(DateTime);
utc = datevec(DateTime);

%% convert position to lat/long

r = [posX posY posZ];
lla = eci2lla(r, utc);

lat = lla(:,1);
long = lla(:,2);
h = lla(:,3);

%% Calculate orbit beta angle

% ecliptic true solar longitude
GMTOffset = 0;
Gamma = solar_ecliptic(year, month, day, hour, mins, GMTOffset);

% calculate RAAN
% need velocity

%% interpolate

% % decrease step size to every 10 second
% [rows, ~] = size(timevec);
% timevec_int = timevec(1):1/(rows*6):timevec(end); % step size of 10 seconds
% timevec_int = timevec_int';
% 
% % interpolate PSARJ and BGA angles for 2A and 4A
% PSARJ_int = interp1(timevec, PSARJ, timevec_int);
% P44A_BGA_int = interp1(timevec,P44A_BGA,timevec_int);
% P42A_BGA_int = interp1(timevec,P42A_BGA,timevec_int);
% lat_int = interp1(timevec, lat, timevec_int);

%% filter out >|52| lat

timevec_flt = timevec(lat<=abs(52));
PSARJ_flt = PSARJ(lat<=abs(52));
P44A_BGA_flt = P44A_BGA(lat<=abs(52));
P42A_BGA_flt = P42A_BGA(lat<=abs(52));

% repeat with interpolated data
% timevec_int = timevec_int(lat_int<=abs(52));
% PSARJ_int = PSARJ_int(lat_int<=abs(52));
% P44A_BGA_int = P44A_BGA_int(lat_int<=abs(52));
% P42A_BGA_int = P42A_BGA_int(lat_int<=abs(52));

%% find intrusions into kov and fov

% time instances and alpha angles of intrusions into EMIT FOV
t2A_fov = timevec_flt(PSARJ_flt >= 268 & PSARJ_flt <= 276);
PSARJ_2Afov = PSARJ_flt(PSARJ_flt >= 268 & PSARJ_flt <= 276);

t4A_fov = timevec_flt(PSARJ_flt >= 89 & PSARJ_flt <= 96);
PSARJ_4Afov = PSARJ_flt(PSARJ_flt >= 89 & PSARJ_flt <= 96);

% time instances and alpha angles of intrusions into stray light KOZ envelope
t2A_koz_a = timevec_flt(PSARJ_flt>=255 & PSARJ_flt<=289);
PSARJ_2Akoz = PSARJ_flt(PSARJ_flt>=255 & PSARJ_flt<=289);

t4A_koz_a = timevec_flt(PSARJ_flt>=75 & PSARJ_flt<=110); 
PSARJ_4Akoz = PSARJ_flt(PSARJ_flt>=75 & PSARJ_flt<=110);

% time instances and beta angles of intrusions into KOZ (tolerances??)
t2A_koz_b = timevec_flt(((PSARJ_flt>=269.99 & PSARJ_flt <= 270.01) & (P42A_BGA_flt >=24 & P42A_BGA_flt <= 150)) | ((PSARJ_flt>=269.99 & PSARJ_flt <= 270.01)& (P42A_BGA_flt >=229 & P42A_BGA_flt <= 337)));
BGA_2Akoz = P42A_BGA_flt(((PSARJ_flt>=269.99 & PSARJ_flt <= 270.01) & (P42A_BGA_flt >=24 & P42A_BGA_flt <= 150)) | ((PSARJ_flt>=269.99 & PSARJ_flt <= 270.01) & (P42A_BGA_flt >=229 & P42A_BGA_flt <= 337)));

t4A_koz_b = timevec_flt(((PSARJ_flt>=89.99 & PSARJ_flt <= 90.01) & (P44A_BGA_flt >=41 & P44A_BGA_flt<=149)) | ((PSARJ_flt>=89.99 & PSARJ_flt <= 90.01) & (P44A_BGA_flt>=230 & P44A_BGA_flt<=338)));
BGA_4Akoz = P44A_BGA_flt(((PSARJ_flt>=89.99 & PSARJ_flt <= 90.01) & (P44A_BGA_flt >=41 & P44A_BGA_flt<=149)) | ((PSARJ_flt>=89.99 & PSARJ_flt <= 90.01) & (P44A_BGA_flt>=230 & P44A_BGA_flt<=338)));

% repeat beta angles with interpolated data (tolerances???)
% t2A_koz_b_int = timevec_int(((PSARJ_int>=269.99 & PSARJ_int <= 270.01) & (P42A_BGA_int >=24 & P42A_BGA_int <= 150)) | ((PSARJ_int>=269.99 & PSARJ_int <= 270.01) & (P42A_BGA_int >=229 & P42A_BGA_int <= 337)));
% BGA_2Akoz_int = P42A_BGA_int((P(PSARJ_int>=269.99 & PSARJ_int <= 270.01)& (P42A_BGA_int >=24 & P42A_BGA_int <= 150)) | ((PSARJ_int>=269.99 & PSARJ_int <= 270.01) & (P42A_BGA_int >=229 & P42A_BGA_int <= 337)));
% 
% t4A_koz_b_int = timevec_int(((PSARJ_int>=89.99 & PSARJ_int <= 90.01) & (P44A_BGA_int >=41 & P44A_BGA_int <= 149)) | ((PSARJ_int>=89.99 & PSARJ_int <= 90.01) & (P44A_BGA_int >= 230 & P44A_BGA_int <=338)));
% BGA_4Akoz_int = P44A_BGA_int(((PSARJ_int>=89.99 & PSARJ_int <= 90.01) & (P44A_BGA_int >=41 & P44A_BGA_int <= 149)) | ((PSARJ_int>=89.99 & PSARJ_int <= 90.01) & (P44A_BGA_int >=230 & P44A_BGA_int <=338)));


%% Percent overlap

percFOV = ((sum(t2A_fov)+sum(t4A_fov))/sum(timevec_flt))*100; % percentage of time
percKOZ = ((sum(t2A_koz_a) + sum(t4A_koz_a))/sum(timevec_flt))*100;

%% total time lost


%% make plots

% plot orbit of ISS
figure()
Re = 6397; % [km]
plot3(posX, posY, posZ);
hold on
[XS, YS, ZS] = sphere(30); % plot the Earth using Matlab sphere command
surf(XS*Re, YS*Re, ZS*Re);
xlabel('x [km]');
ylabel('y [km]');
zlabel('z [km]');
title('Orbit of ISS Around Earth');
hold off

% alpha angles
figure()
plot(timevec_flt, PSARJ_flt)
hold on
plot(t2A_fov, PSARJ_2Afov);
plot(t4A_fov, PSARJ_4Afov);
legend('PSARJ Angles','Region of P4-2A Intrusion (268^o - 276^o)','Region of P4-4A Intrusion (89^o - 96^o)');
xlabel('Time [days]');
ylabel('\alpha [deg]');
title(['ISS SAW PSARJs intrude into EMIT FOV ' num2str(percFOV) '% of 2019']);
hold off

figure()
plot(timevec_flt, PSARJ_flt)
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
plot(timevec_flt, P42A_BGA_flt);
hold on
plot(t2A_koz_b, BGA_2Akoz,'*');
legend('P4-2A BGA Angles','P4-2A intrusions into EMIT KOZ');
ylabel('\beta [deg]');
xlabel('Time [days]');
title('\alpha=270^o, Beta angles resulting in P4-2A SAW Intrusions into EMIT KOZ');
hold off

figure()
plot(timevec_flt, P44A_BGA_flt);
hold on
plot(t4A_koz_b, BGA_4Akoz,'*');
legend('P4-4A BGA Angles','P4-4A intrusions into EMIT KOZ');
ylabel('\beta [deg]');
xlabel('Time [days]');
title('\alpha=90^o, Beta angles resulting in P4-4A SAW Intrusions into EMIT KOZ');
hold off

% % beta anlges with interpolated data
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

% histograms

figure() % PSARJ into FOV
histogram(PSARJ_flt);
hold on
histogram(PSARJ_2Afov);
histogram(PSARJ_4Afov);
legend('PSARJ angles','P4-2A Intrusion (268^o - 276^o)','P4-4A Intrusion (89^o - 96^o)','Location','northwest')
xlabel('Angles [deg]');
ylabel('Frequency');
title(['ISS SAW PSARJs intrude into EMIT FOV ' num2str(percFOV) '% of 2019']);
hold off

figure() % PSARJ into KOZ
histogram(PSARJ_flt);
hold on
histogram(PSARJ_2Akoz);
histogram(PSARJ_4Akoz);
legend('PSARJ angles','P4-2A Intrusion (255^o - 289^o)','P4-4A Intrusion (75^o - 110^o)','Location','northwest')
xlabel('Angles [deg]');
ylabel('Frequency');
title(['ISS SAW PSARJs intrude into EMIT KOZ ' num2str(percKOZ) '% of 2019']);
hold off


toc
