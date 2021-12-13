% This file is part of OGNSS_scheduler.
%
% OGNSS_scheduler is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% OGNSS_scheduler is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Lesser General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with OGNSS_scheduler.  If not, see <http://www.gnu.org/licenses/>.
%
% Copyright (C) 2020 Andrea Nardin
% Navigation, Signal Analysis and Simulation (NavSAS) group,
% Politecnico di Torino


function mstData = showOutput(mstData,mstSequence,settings)
%SHOWOUTPUT Process output data and show data and plots.


% Init operations
Nmeo    = settings.Nmeo;
Nleo    = settings.Nleo;
NtimeDS = size(mstData.idleAntennas,2);

% No of epochs with disconnected MST
disconnectedEpochs  = sum(mstData.disconnected>1);

% Average no of idle antennas per epoch
meanIdleAntennas    = sum(mstData.idleAntennas)/NtimeDS;

% Perc. of epochs where at least 1 LEO is connected to 3 different MEO orbital planes
threePlanes = sum(mstData.leo3planes>0)/NtimeDS*100;

% max MEO revisiting time of all epochs
maxrevtime = mstData.maxRevTime;

% No of epochs with at least 1 node exceeding its max degree
dd = zeros(1,NtimeDS);
for ii= 1:NtimeDS
    dd(ii)=sum(degree(mstSequence{ii})>settings.maxNodeDeg)>0;
end
nonBoundedEpochs = sum(dd);

%--- Data delivery latency
A               = zeros(Nmeo);          % temp array
lp_id1          = zeros(NtimeDS,Nmeo*(Nmeo+1)/2);       % index of MEOs that are further from each other 
lp_id2          = zeros(NtimeDS,Nmeo*(Nmeo+1)/2);       % index of MEOs that are further from each other 
longestPaths    = zeros(NtimeDS,1);       % distance between these MEOs
avgdist         = zeros(NtimeDS,Nmeo);    % average distance from MEO i at epoch j (+Nleo to have index consistent with weight matrices)

for ii = 1:NtimeDS
    A(:,:) = mstData.distMat(ii,Nleo+1:end,Nleo+1:end);
    
    % longest possible path among MEOs for each epoch   
    B = triu(A);
    M = max(max(B));
    I = find(B==M);
    longestPaths(ii) = M;
    [I1,I2] = ind2sub(size(A),I);
    lp_id1(ii,1:length(I1)) = I1;
    lp_id2(ii,1:length(I2)) = I2;
    
    % average distance for each MEO and each epoch
    avgdist(ii,:) = sum(A)/(Nmeo-1);
end
[worstDist.dist,worstDist.epoch] = max(longestPaths);

%--- Assignation Fairness
% Amount o contact time received by each MEO (to any LEO)
% (tpoint is also accounted since it involves energy consumption)
A = mstData.adjMat;
A(~isnan(A)) = 1; % each selected contact in the plan is count as 1
B = sum(sum(A,2,'omitnan'),'omitnan'); % each non selected, non existing, or M2M link count is omitted as NaN

meoTotContactTime = zeros(Nmeo,1);
meoTotContactTime(:) = B(1,1,Nleo+1:end);

% Jain's fairness Index
jfi = sum(meoTotContactTime)^2/(Nmeo*sum(meoTotContactTime.^2));

% Output numeric data
T1 = table(disconnectedEpochs,meanIdleAntennas,threePlanes,maxrevtime,...
    nonBoundedEpochs,worstDist.dist,worstDist.epoch,jfi,'variableNames',...
    {'Disconnected ep.','Avg idle antennas','No. of epochs with at least 1 L23M',...
    'Max revisit time','Ep. exceeding max node degree', 'Longest distance',...
    'Longest distance ep.','JFI'})
T2 = table(settings.percentiles',mstData.percentiles,'variableNames',...
    {'Percentile', 'PDOP'})

%--- Update mstData
mstData.disconnectedEpochs  = disconnectedEpochs;
mstData.meanIdleAntennas    = meanIdleAntennas;
mstData.threePlanes         = threePlanes;
mstData.maxrevtime          = maxrevtime;
mstData.nonBoundedEpochs    = nonBoundedEpochs;
mstData.longestPaths_id1    = lp_id1;
mstData.longestPaths_id2    = lp_id2;
mstData.longestPaths        = longestPaths;
mstData.avgdist             = avgdist;
mstData.worstDist           = worstDist;
mstData.meoTotContactTime   = meoTotContactTime;
mstData.jfi                 = jfi;
mstData.table1              = T1;
mstData.table2              = T2;


%% export Data
AA = zeros(size(mstData.meoPdops)); 
AA(mstData.meoPdops~=0) = 1;

if ~exist('../results', 'dir')
   mkdir('../results')
end

writematrix(mstData.meoPdops, '../results/Pdops.csv');
writematrix(mstData.revTimes, '../results/revTimes.csv');
writematrix(AA, '../results/meoAssignations.csv');

writematrix(mstData.longestPaths_id1, '../results/longestPathId1.csv');
writematrix(mstData.longestPaths_id2, '../results/longestPathId2.csv');
writematrix(mstData.longestPaths, '../results/longestPaths.csv');
writematrix(mstData.avgdist, '../results/averageDistances.csv');

% writematrix(mstData.meoTotContactTime,'output/meoTotContactTime.csv');
writetable(mstData.table1,'../results/table1.csv');
writetable(mstData.table2,'../results/table2.csv');


%% Plots

defaultAxesFontsize = 16;
defaultLegendFontsize = 14;
defaultLegendInterpreter = 'latex';
defaultLinelinewidth = 2;
defaultAxesTickLabelInterpreter = 'latex';
defaultTextInterpreter = 'latex';
set(0,'defaultAxesFontsize',defaultAxesFontsize,'defaultLegendFontsize',defaultLegendFontsize,...
    'defaultLegendInterpreter',defaultLegendInterpreter,'defaultLinelinewidth',defaultLinelinewidth,...
    'defaultAxesTickLabelInterpreter',defaultAxesTickLabelInterpreter);
set(0,'defaultTextInterpreter',defaultTextInterpreter);



figure
subplot(3,1,1)
plot(mstData.avgPdop), hold on
xlabel('$t_{min}$ slot index')
ylabel('PDOP');
axis tight, grid on
legend 'mean'

subplot(3,1,2)
plot(mstData.maxPdop), hold on
xlabel('$t_{min}$ slot index')
ylabel('PDOP');
axis tight, grid on
legend 'max'

subplot(3,1,3)
plot(mstData.minPdop)
xlabel('$t_{min}$ slot index')
ylabel('PDOP');
axis tight, grid on
legend 'min'

set(gcf, 'Position',  [100, 100, 700, 500])
saveas(gcf,'../results/aggregatedPdop.png')


% figure
% surf(avgdist)
% ylabel('MEO id')
% xlabel('$t_{min}$ window index');
% zlabel('Average hop distance');
end
