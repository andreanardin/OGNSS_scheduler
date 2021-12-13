% *************************************************************************
% ************************ O-GNSS Scheduler ******************************* 
% * PDOP-driven Scheduler for Optical Inter-Satellite Links enabled GNSS **
% *************************************************************************

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


close all
clearvars

addpath 'library'
addpath 'STKdata'


%% To retrieve STK data:

% 1) Compute MEO2MEO visibility access report to get average H matrix
% values. Choose 3 adjacent MEO, compute relative elv and az pairwise,
% then compute average H matrix entries for the central MEO

% 2) Compute relative metrics for each visible M2L link at each instant
% First: select multiple items to set a common constrants
% For each MEO
%   -set cutoff elevation for each MEO!! (14.6° conical aperture from
%   nadir, i.e. max elev -75.4°) (max 0° REALLY illuminates the whole LEO
%   segment!!)
%   - set cutoff for each LEO (visibility cone aperture 90°, i.e. min
%   elev 0°
%   Go to MEO>Access> then select multiple LEO target >AER report 
% the report contains relative elev and azimuth (and range)
% For 24 MEO you get 24 reports. In each there are 4 access reports with all LEOs
% Set: scenario interval, step time
% ( You can then modify elevation constraint for all LEOs and directly update all reports)

tic

%% Settings

settings = initSettings();


%% Initial operations

Nleo        = settings.Nleo;
Nmeo        = settings.Nmeo;
orbitPlanes = settings.orbitPlanes;
tmin        = settings.tmin;
tpoint      = settings.tpoint;
accessTimeStep = settings.accessTimeStep;

accessMatrices = zeros(Nleo,settings.Ntime,Nmeo);


%% Compute PDOP

if ~settings.loadWeights
%--- Compute average incomplete Hmatrix for MEO satellites
Hfixed = computeMeoHfixed(settings.m2mAccess_filename);

%--- For each MEO: read access reports, compute PDOP, store data in a
% epochs-coherent matrix
N = Nmeo/orbitPlanes;
for ii = 1:Nmeo
    
    %--- Read Access (AER) report from STK (csv file)
    % build input filename for current LEO
    meo     = mod(ii-1,N);
    orbit   = ceil(ii/N)-1;
    C = textscan(settings.m2lAccess_filename,'%s','delimiter','*');
    filename = [C{1}{1} num2str(orbit) C{1}{2} num2str(meo) C{1}{3}];
    [m2lData, time] = readSTKaccessCsv(filename,Nleo,settings.scenarioStartTime);
    
    % --> additional cutoff filter??? here or new STK reports generated
    
    %--- Calculate PDOP (Nmeo cells)
    m2lPdopData = computePdopDataset(m2lData,Hfixed);
    
    
    %--- merge MEO PDOPs in a matrix coherent with the epochs
    % Thus [Nleo x Scenarioduration] matrix where each PDOP is in a position
    % coherent with its epoch
    
    for jj = 1: Nleo
        % convert time in integer epoch indexes
        epochs = floor(time{jj}/accessTimeStep)+1;
        % associate PDOPs with proper time index in the matrix
        accessMatrices(jj,epochs,ii) = m2lPdopData{jj};
        % Note for adjacent datetime which belong to the same step (e.g. final
        % instants of a visibility window) the most recent PDOP is considered
    end
end

%% Compute average PDOPs on accessMatrix
%--- Compute average slotted PDOP
[W1, W2] = computeSlottedWeights(accessMatrices,accessTimeStep,tmin,tpoint);
else
    load('output/weights','W1','W2');
end
if settings.saveWeights
    save('output/weights.mat','W1','W2','settings');
end


%% GRAPH & MST

%--- Compute sequence of MST given a dynamic matrix of weights (2 flavors)
[mstSequence, mstData] = computeMstSequence(W1,W2,settings);
toc


%% OUTPUT AND PLOTS

mstData = showOutput(mstData,mstSequence,settings);

