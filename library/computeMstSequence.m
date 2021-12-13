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

function [mstSequence, mstData] = computeMstSequence(W1,W2,settings)
% Compute sequence of MST given a dynamic matrix of weights (2 flavors)
% Undirected, weighted, dynamic graph Gn(V,En) (n discrete time)
% The graph is represented by a Visibility matrix 3-D: V = Nleo x NtimeDs x Nmeo

[Nleo,NtimeDs,Nmeo] = size(W1);
Nplanes = settings.orbitPlanes;
N = Nmeo/Nplanes;
maxepsilon = settings.maxepsilon;

%--- Initial operations
% initial mst (empty)
mst = graph(zeros(Nleo+Nmeo));
% init output analysis variables
mstSequence = cell(1,NtimeDs);
mstAdjMat   = zeros(NtimeDs,Nleo+Nmeo,Nleo+Nmeo);

mstAdjMat_meoOnly   = zeros(NtimeDs,Nleo,Nmeo);

mstMeoPdops = zeros(NtimeDs,Nmeo);
mstRevTimes = zeros(NtimeDs,Nmeo);
mstDistMat     = zeros(NtimeDs,Nleo+Nmeo,Nleo+Nmeo);
mstDisconnected = zeros(1,NtimeDs);
mstAveragePdop = zeros(1,NtimeDs);
mstStdPdop = zeros(1,NtimeDs);
mstMaxPdop = zeros(1,NtimeDs);
mstMinPdop = zeros(1,NtimeDs);
leo3planes = zeros(1,NtimeDs);
idleAntennas = nan(1,NtimeDs);
n = zeros(1,Nleo + Nmeo);
n_L = zeros(1,Nleo);
maxRevTime = 0;

%% Build Graph G
%--- Create adjacency matrix
% (M2M intra-orbit fixed (epsilon), M2M = 0 elsewhere, L2L = 0 everywhere,
% L2M are computed at each time from W1 and W2)

%--- Create M2M links (fixed weights)
w1 = zeros(Nleo+Nmeo);
w2 = zeros(Nleo+Nmeo);


% For each tmin slot
for ii = 1:NtimeDs
    
    % You can generate random epsilon to let each MEO edge be part of a mst each
    % time (independtly of edge order, otheriwse epsilon ties are broken by node order)
    w1 = restoreM2M(w1,settings,'adjacency');
    w2 = restoreM2M(w2,settings,'adjacency');
    
    %--- Build the proper weight matrix according to past time slot
    % build symmetric weight matrices, coherent with graph adjacency matrix
    w1(1:Nleo,Nleo+1:end) = W1(:,ii,:);
    w1(Nleo+1:end,1:Nleo) = w1(1:Nleo,Nleo+1:end).';
    w2(1:Nleo,Nleo+1:end) = W2(:,ii,:);
    w2(Nleo+1:end,1:Nleo) = w2(1:Nleo,Nleo+1:end).';
    
    % se l'arco esisteva nel MST usa W1 altrimenti usa W2
    M = full(adjacency(mst));
    weights = w2;
    weights(M>0) = w1(M>0);
    
    %--- Generate graph from weight matrix for the current time slot
    if settings.fakeWeights
        % No criterion for mst. All l2m weights = 1 or rand.
        fakeweights = weights;
        %         fakeweights(weights>max(0,maxepsilon)) = 1;
        FK = 3+30*rand(size(weights,1),size(weights,2));
        fakeweights(weights>max(0,maxepsilon)) = FK(weights>max(0,maxepsilon));
        
        G = graph(fakeweights,'upper');
    else
        G = graph(weights,'upper');
    end
    
    
    %% Extract DCMST + secondary optimization
    
    %--- Compute Degree-Constrained MST
    mst = dcmstHeuristic(G,settings.maxNodeDeg,settings);
      
    % rebuild some M2M links (erased by mst construction)
    mst_cyclic = restoreM2M(mst,settings);
    
    %--- Secondary Optimization applied on idle antennas
    if ~strcmp(settings.secondaryOpt,'none')
        mst_cyclic = secondaryOptimization2(mst_cyclic,G,n,settings);
    end
    
    %--- Update Revisit time
    % degree
    D = degree(mst_cyclic);
    
    %--- Compute revisit time
    % Leo revisit time
    n_L(D(1:Nleo)>0) = 0;
    n_L(D(1:Nleo)<=0) = n_L(D(1:Nleo)<=0)+1;
    % Meo revisit time
    f = settings.fixedLctsMeo; % fixed LCTs in each MEO node
    n(D>f) = 0;
    n(D<=f) = n(D<=f)+1;
    % merge
    n(1:Nleo) = n_L;
    
    
    %%  Output data collection
    
    mstSequence{ii} = mst_cyclic;
    
    % Connected components
    [~,binsize] = conncomp(mst_cyclic);
    if length(binsize) > 1
        mstDisconnected(ii) = binsize;
    end
    
    %--- PDOP statistic
    A = adjacency(mst_cyclic,'weighted'); A(A<=max(0,maxepsilon))=NaN;
    if settings.fakeWeights
        % use real weights
        A(A>max(0,maxepsilon)) = weights(A>max(0,maxepsilon));
    end
    
    mstAdjMat(ii,:,:)   = A;
    
    % MEO PDOP statistic
    % Adjacency with MEO
    Ameo = A(1:Nleo,Nleo+1:end); 
    mstAdjMat_meoOnly(ii,:,:)   = Ameo;
    mstMeoPdops(ii,:) = mean(Ameo,'omitnan');    
    % Mean & Std dev PDOP of just L2M links
    mstAveragePdop(ii) = mean(mean(Ameo,'omitnan'),'omitnan');
    mstStdPdop(ii)      = std(std(Ameo,'omitnan'),'omitnan');
    mstMaxPdop(ii)      = max(max(Ameo));
    mstMinPdop(ii)      = min(min(Ameo));
    
    % LEO with at least 3 MEO planes connected x epoch
    for jj = 1:Nleo
        % compute connected & disconnected planes
        connections = find(~isnan(A(jj,:)));
        connectedPlanes = ceil((connections-Nleo)/N);
        leo3planes(ii) = leo3planes(ii) + (length(unique(connectedPlanes)) == settings.orbitPlanes);
    end
    
    % Idle LEO antennas
    idleAntennas(ii) = sum(settings.maxNodeDeg(1:Nleo) - D(1:Nleo));
    
    % Revisiting time statistic
    mstRevTimes(ii,:) = n(Nleo+1:end);
    if max(n(Nleo+1:end)) > maxRevTime
        maxRevTime = max(n(Nleo+1:end));
    end
    
    % Data delivery latency
    mstDistMat(ii,:,:) = distances(mst_cyclic,'method','unweighted');
    
    
end
mstData.adjMat          = mstAdjMat;
mstData.meoPdops        = mstMeoPdops;
mstData.revTimes        = mstRevTimes;
mstData.disconnected    = mstDisconnected;
mstData.avgPdop         = mstAveragePdop;
mstData.stdPdop         = mstStdPdop;
mstData.leo3planes      = leo3planes;
mstData.idleAntennas    = idleAntennas;
mstData.maxRevTime      = maxRevTime;
mstData.maxPdop         = mstMaxPdop;
mstData.minPdop         = mstMinPdop;
mstData.percentiles     = prctile(mstAdjMat,settings.percentiles,'all');
mstData.distMat         = mstDistMat;

end
