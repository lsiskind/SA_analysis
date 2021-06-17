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

% position data
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
lla = eci2lla(r, [year month day hour mins sec]);

lat = lla(:,1);
long = lla(:,2);
h = lla(:,3);
toc
%% Calculate orbit beta angle
% tic
% ecliptic true solar longitude
% [~, ~, Gamma, ~ ] = earthEphemeris(DateTime);

% calculate RAAN - need velocity

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
mins_flt = mins(lat<=abs(52));

PSARJ_flt = PSARJ(lat<=abs(52));
P44A_BGA_flt = P44A_BGA(lat<=abs(52));
P42A_BGA_flt = P42A_BGA(lat<=abs(52));

% repeat with interpolated data
% timevec_int = timevec_int(lat_int<=abs(52));
% PSARJ_int = PSARJ_int(lat_int<=abs(52));
% P44A_BGA_int = P44A_BGA_int(lat_int<=abs(52));
% P42A_BGA_int = P42A_BGA_int(lat_int<=abs(52));

%% find intrusions into kov and fov

% time instances and alpha angles of intrusions into KOZ and FOV
[afov2, tafov2, akoz2, takoz2, afov4, tafov4, akoz4, takoz4] = alphaIntrusions(PSARJ_flt, timevec_flt);

% time instances and beta angles of intrusions into KOZ
[bkoz2, tbkoz2, bkoz4, tbkoz4] = betaIntrusions(PSARJ_flt, P42A_BGA_flt, P44A_BGA_flt, timevec_flt);

% alpha angles where bga intrusions occur
[abkoz2, abkoz4, tabkoz2, tabkoz4] = alphabetaIntrusions(PSARJ_flt, P42A_BGA_flt, P44A_BGA_flt, timevec_flt);

%% time lost

% amount of time missing from data set- unrelated to intrusions
% (pre-filter)
[totaldata, ~] = size(mins);
expecteddata = 365*24*60;
mins_missing = expecteddata - totaldata;
hours_missing = mins_missing/60;
days_missing = hours_missing/24;
perc_missing = (mins_missing/expecteddata)*100;

[totalmins, ~] = size(mins_flt); % number of minutes spent taking science

% time p4-2A intruded into koz due to dynamic psarj
[tint_akoz2,~] = size(takoz2); 
[tint_bkoz2,~] = size(tbkoz2);
plost_akoz2 = (tint_akoz2/totalmins)*100; % percent of time lost due to p4-2a intrusions
fprintf('PSARJ dynamic rotation resulted in P4-2A intrusion into EMIT KOZ for a total of %d minutes and %d%% of mission operation time\n', tint_akoz2, plost_akoz2); 

% time p4-4A intruded into koz due to dynamic psarj
[tint_akoz4,~] = size(takoz4); 
[tint_bkoz4,~] = size(tbkoz4);
plost_akoz4 = (tint_akoz4/totalmins)*100; % percent of time lost due to p4-4a intrusions
fprintf('PSARJ dynamic rotation resulted in P4-4A intrusion into EMIT KOZ for a total of %d minutes and %d%% of mission operation time\n', tint_akoz4, plost_akoz4); 

% time p4-2A intruded into fov due to dynamic psarj
[tint_afov2,~] = size(tafov2); 
plost_afov2 = (tint_afov2/totalmins)*100; % percent of time lost due to p4-4a intrusions
fprintf('PSARJ dynamic rotation resulted in P4-2A intrusion into EMIT FOV for a total of %d minutes and %d%% of mission operation time\n', tint_afov2, plost_afov2); 

% time p4-4A intruded into fov due to dynamic psarj
[tint_afov4,~] = size(tafov4); 
plost_afov4 = (tint_afov4/totalmins)*100; % percent of time lost due to p4-4a intrusions
fprintf('PSARJ dynamic rotation resulted in P4-44 intrusion into EMIT FOV for a total of %d minutes and %d%% of mission operation time\n', tint_afov4, plost_afov4); 

% total time KOZ intruded on due to dynamic psarj
ttotal_koz = tint_akoz2 + tint_akoz4; %[min]
perc_total_koz = ((ttotal_koz)/totalmins)*100;

% total time FOV intruded on due to dynamic psarj
ttotal_fov = tint_afov2 + tint_afov4;

%% make plots
close all

%----------------histograms----------------
nbins = .5;
figure()
t = tiledlayout(1,2);
title(t,'P4-2A SAW Interference','FontWeight','Bold');
nexttile
histogram(akoz2,'BinWidth',nbins);
hold on
histogram(afov2,'BinWidth',nbins,'FaceColor','red');
histogram(abkoz2, 'BinWidth', nbins,'FaceColor','green');
hold off
xlabel('Angle [deg]');
ylabel('Frequency [min]');
title(['Dynamic PSARJ Rotation Causes Interferences in EMIT KOZ for ',num2str(tint_akoz2),' min'],['Interference in EMIT FOV for ',num2str(tint_afov2),' min'],'FontWeight','normal');
legend('KOZ Interference 255^o-289^o','FOV Inteference 268^o-276^o','PSARJ Angles where BGA Interference Occurs');

nexttile
histogram(bkoz2);
xlabel('Angle [deg]');
ylabel('Frequency [min]');
legend('BGA Interference Range: 24^o-15^o & 229^o-227^o');
title(['In PSARJ Range (left), Static BGA Intrusion into KOZ for ', num2str(tint_bkoz2),' min'],'FontWeight','normal');
ylim([0 1200]);

figure()
t = tiledlayout(1,2);
title(t,'P4-4A SAW Interference','FontWeight','Bold');
nexttile
histogram(akoz4,'BinWidth',nbins);
hold on
histogram(afov4,'BinWidth',nbins,'FaceColor','red');
histogram(abkoz4, 'BinWidth', nbins,'FaceColor','green');
hold off
xlabel('Angle [deg]');
ylabel('Frequency [min]');
title(['Dynamic PSARJ Rotation Causes Interferences in EMIT KOZ for ',num2str(tint_akoz4),' min'],['Interference in EMIT FOV for ',num2str(tint_afov4),' min'],'FontWeight','normal');
legend('KOZ Interference 75^o-110^o','FOV Inteference 89^o-96^o','PSARJ Angles where BGA Interference Occurs');

nexttile
histogram(bkoz4);
xlabel('Angle [deg]');
ylabel('Frequency [min]');
title(['In PSARJ Range (left), Static BGA Intrusion into KOZ for ', num2str(tint_bkoz4),' min'],'FontWeight','normal');
legend('BGA Interference Range: 41^o-149^o & 230^o-338^o');

%----------- line plots--------------
figure()
plot(timevec_flt, PSARJ_flt,'.','LineWidth',.5);
hold on
plot(takoz2, akoz2,'.','LineWidth',.5);
plot(tabkoz2, abkoz2,'.y','LineWidth',.5);
plot(takoz4, akoz4,'.','LineWidth',.5);
plot(tabkoz4, abkoz4,'.y','LineWidth',.5);
xlabel('Time [min]');
ylabel('Angle [deg]');
legend('Filtered PSARJ Angles','P4-2A PSARJ Interference Angles','P4-2A PSARJ Angles of BGA Interference','P4-4A PSARJ Interference Angles','P4-4A PSARJ Angles of BGA Intereference');
title(['P4-2A and P4-4A Interfere in EMIT KOZ ', num2str(perc_total_koz),'% of Operation Time'],['P4-4A Interferes in EMIT KOZ ',num2str(plost_akoz4),'% of Operation Time']);
subtitle([num2str(mins_missing),' min of Missing Data Pre-Latitude Filter (', num2str(perc_missing),'% of year)']); 


toc

