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

function [mst, idleAntennas] = secondaryOptimization2(mst,G,n,settings)
% Global best pdop or global best revisit time.
% The optimization is based on sorting the edges according to the chosen
% criteria and then add them one-by-one, satisfying the system constraints

% get settings
opt = lower(settings.secondaryOpt);
Nleo = settings.Nleo;
Nmeo = settings.Nmeo;
Nplanes = settings.orbitPlanes;
maxNodeDeg = settings.maxNodeDeg;
N = Nmeo/Nplanes;

% compute a list of edges + weights (Nedges x 3))
M = G.Edges{:,:};
priorityServed = 0;
switch opt
    case 'pdop'
        % sort by weight
        M = sortrows(M,3,'ascend');
    case 'revtime'
        % comprehensive revisit time for an edge: revtime s + revtime t
        n(1:Nleo)=0;
        edgeRevTime = n(M(:,1)) + n(M(:,2));       
        if settings.randomties
            [~,I] = sort(edgeRevTime + rand(size(edgeRevTime))-0.5,'descend');
        else
            [~,I] = sort(edgeRevTime,'descend');
        end
        M = M(I,:);
    otherwise
        error 'unknown secondary optimization criterion'
end

ii = 1;

while ii <= size(M,1)
    s = M(ii,1); t = M(ii,2);
    
    % if it's not in mst try to add it
    if findedge(mst,s,t) == 0
        if degree(mst,s)+1 <= maxNodeDeg(s) && degree(mst,t)+1 <= maxNodeDeg(t)
            % if L23M priority and priority was not yet given
            if settings.threePlanesPriority && ~priorityServed
                if connectNewPlane(mst,s,t,Nleo,N)
                    mst = addedge(mst,s,t,M(ii,3));
                end
                % otherwise add the edge directly
            else
                mst = addedge(mst,s,t,M(ii,3));
            end
        end
        
    end
    ii = ii+1;
    if ii>size(M,1) && ~priorityServed
        ii = 1;
        priorityServed = 1;
    end
end

% compute total LEO idleAntennas for this mst
dd = degree(mst);
idleAntennas = sum(dd(1:Nleo) - maxNodeDeg(1:Nleo));
end

function out = connectNewPlane(G,s,t,Nleo,N)
% True if s,t connect a new (and different) orbital plane to s or to t
out = 0;

% if is a L2M link
if (s <= Nleo) + (t<=Nleo) == 1
    % which is a LEO
    if s<=Nleo
        x = s;
        meo = t;
    else
        x = t;
        meo = s;
    end
    
    nn = neighbors(G,x);
    connectedPlanes = ceil((nn-Nleo)/N);
    isnotthere = isempty(find(connectedPlanes == ceil((meo-Nleo)/N),1));
    if isnotthere
        out = 1;
    else
        out = 0;
    end   
end
end

