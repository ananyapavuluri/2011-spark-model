% Plotting data from the 2011 Sticky Cluster Model
% PROVIDED BY DR. ERIC SOBIE, DEPARTMENT OF PHARMOCOLOGY, ICAHN SCHOOL OF MEDICINE AT MOUNT SINAI

% Read in csv files from C++
Nopen_all = csvread('N_open.csv');
Irel_all = csvread('Irel.csv');
CaJSR_all = csvread('CaJSR.csv');
Cass_all = csvread('Ca_ss.csv');
plottime = csvread('plottime.csv');

% Redefine CSQ buffering constants
CSQ = 30e3;
KCSQ = 630;

% Plots
figure
plot(plottime, Nopen_all)
title('Number of RyRs Open over Time');
xlabel('# of RyRs open');
ylabel('Time in ms');
print -depsc Nopen
    %undergoing a spark
sparks = find(max(Nopen_all) > 5);

figure
plot(plottime, CaJSR_all(:,sparks))
hold on
    %free CaJSR during a spark
CaJSR_avg = mean(CaJSR_all(:, sparks), 2);
plot(plottime, CaJSR_avg, 'b', 'LineWidth', 2.25)
title('Average Concentration of Free CaJSR Over Time');
xlabel('Time in ms');
ylabel('Average Concentration in uM');
print -depsc freeCaJSR

CaJSRtot_all = CaJSR_avg +  CaJSR_avg * CSQ./(KCSQ + CaJSR_avg);
CaJSRtot_all = CaJSRtot_all/max(max(CaJSRtot_all));

figure
hold on
plot(plottime, CaJSRtot_all)
plot(plottime, 1 - exp(-plottime/90), 'r', 'LineWidth', 2.5)
title('Total CaJSR')
print -depsc totalCaJSR

