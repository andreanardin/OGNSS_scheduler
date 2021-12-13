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

function pdopDataset = computePdopDataset(accessData,Hfixed)
% Compute the PDOP values for each access epochs in accessData, for all
% access pairs, considering a partial H matrix fixed that characterizes
% each target satellite of the access pair.

Ntgt = length(accessData);
pdopDataset = cell(1,Ntgt);
for ii = 1:Ntgt
    az = accessData{ii}{2};
    el = accessData{ii}{3};
    pdopDataset{ii} = zeros(length(az),1);
    
    if length(az) ~= length(el)
        error 'Unconsistent data: elevation and azimuth sets have different lengths.'
    end
    for jj = 1:length(az)
        H = buildHmatrix(az(jj),el(jj),Hfixed);    
        pdopDataset{ii}(jj) = computeDOP(H,'pdop3');
    end
end


end
