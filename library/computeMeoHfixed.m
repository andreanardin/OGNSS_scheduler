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

function Hfixed = computeMeoHfixed(m2mAccess_filename)
% Compute average incomplete H-matrix for a sample MEO satellite with its 2
% neighbors

NsampleMeo = 2;
% read azimuth and elevation data
m2mData = readSTKaccessCsv(m2mAccess_filename,NsampleMeo);
% compute mean azimuth and elevation per meo pairs
meanAzmt = zeros(1,NsampleMeo);
meanElev = zeros(1,NsampleMeo);
for ii = 1:NsampleMeo
   meanAzmt(ii) = mean(mod(m2mData{ii}{2}(:),360));
   meanElev(ii) = mean(m2mData{ii}{3}(:));
end

% build H-matrix
Hfixed = buildHmatrix(meanAzmt,meanElev);
end
