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

function mst_cyclic = restoreM2M(mst,settings,format)
% Restores epsilon valus over connected orbital planes
% if input is graph, output is a graph. Vice versa if input is an adjacency
% matrix, output is an adjacency matrix.


if nargin<3 || strcmp(format,'graph')
    A = adjacency(mst,'weighted');
elseif strcmp(format,'adjacency')
    A = mst;
else
    error 'unknown I/O format'
end

Nleo = settings.Nleo;
Nplanes = settings.orbitPlanes;
N = settings.Nmeo/Nplanes;
for jj = 1:Nplanes
    % a = settings.maxepsilon/2*rand(1,N-1);
    a = settings.maxepsilon * ones(1,N-1);
    R = diag(a,1);
    %R(1,N) = settings.maxepsilon/2*rand;
    R(1,N) = settings.maxepsilon;
   
    win = Nleo+1+(jj-1)*N:Nleo+jj*N;
    A(win,win) = R + R.';
end

if nargin<3 || strcmp(format,'graph')
    mst_cyclic = graph(A,'upper');
else
    mst_cyclic = A;
end

end
