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

function [settings] = initSettings()
%INITSETTINGS initialize the general settings of the simulator
%   Specify your settings below

%--- Save & load
saveWeights = 0;
loadWeights = 0;

%--- Random CPD
fakeWeights = 0; % simulate a random scheduling (satisfy system constraints)

%--- MST % heuristic parameters
mstMethod           = 'dense'; % better for randomizing the root of the mst computation
% mstMethod           = 'sparse';

% Number of optical transponders for each LEO and MEO sat 
maxDegreeLeo = 3;
maxDegreeMeo = 3;
fixedLctsMeo = 2; % changing it does not update the adjacency matrix!

maxepsilon = realmin;

%--- Secondary optimization criteria
% secondaryOpt    = 'pdop';
% secondaryOpt    = 'revTime'; 
secondaryOpt    = 'none';
threePlanesPriority = 0;
randomties = 1; % Break ties randomly. Only for revTime (a tie with PDOP is very unlikely)

%--- Input files
% Complete access report for 1 LEO satellite to all MEO satellites.
% '*' means that a matching filename is used for each LEO satellite report
m2lAccess_filename = 'Satellite-MEO_Satellite10*0*-To-Satellite-LEO_Satellite20000 AER.csv';
% A Sample access report for 1 sample MEO satellite with its 2 neighbors
m2mAccess_filename = 'Satellite-MEO_Satellite1-To-Satellite-MEO_Satellite2 AER.csv';

%--- Constellation parameters
Nmeo        = 24;
Nleo        = 6;
orbitPlanes = 3;

%--- Optical link constraints
tmin    = 5*60; % (s)
tpoint  = 2*60; % (s)
tuseful = tmin - tpoint;

%--- Scenario date & time
scenarioStartTime   = '1 Jan 2020 00:00:00.000';
scenarioEndTime     = '11 Jan 2020 12:00:00.000';
accessTimeStep      = 30; % access report time step (seconds)

%--- Output options
percentiles = [25 50 75 95 99 99.9];


%% Init operations

scenarioStartTime   = datetime(scenarioStartTime,'InputFormat','dd MMM yyyy HH:mm:ss.SSS');
scenarioEndTime     = datetime(scenarioEndTime,'InputFormat','dd MMM yyyy HH:mm:ss.SSS');
scenarioDuration = scenarioEndTime - scenarioStartTime;
Ntime = floor(seconds(scenarioDuration)/accessTimeStep)+1;

settings.saveWeights = saveWeights;
settings.loadWeights = loadWeights;
settings.m2lAccess_filename = m2lAccess_filename;
settings.m2mAccess_filename = m2mAccess_filename;
settings.scenarioStartTime = scenarioStartTime;
settings.scenarioEndTime = scenarioEndTime;

settings.Nmeo = Nmeo;
settings.Nleo = Nleo;
settings.orbitPlanes = orbitPlanes;
settings.maxDegreeLeo = maxDegreeLeo;
settings.maxDegreeMeo = maxDegreeMeo;
settings.fixedLctsMeo = fixedLctsMeo;
settings.maxepsilon = maxepsilon;
settings.fakeWeights = fakeWeights;
settings.percentiles = percentiles;
settings.mstMethod = mstMethod;
settings.secondaryOpt = secondaryOpt;
settings.threePlanesPriority = threePlanesPriority;
% maxdeg could be even set differently for each single node
settings.maxNodeDeg = [maxDegreeLeo*ones(Nleo,1); maxDegreeMeo*ones(Nmeo,1)];
settings.Ntime = Ntime;
settings.accessTimeStep = accessTimeStep;
settings.tmin = tmin;
settings.tpoint = tpoint;
settings.randomties = randomties;

end


