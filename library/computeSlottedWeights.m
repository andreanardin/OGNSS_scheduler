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

function [W1, W2, NtimeDs] = computeSlottedWeights(accessMat,accessTimeStep,tmin,tpoint)
% Compute average slotted PDOP
% accessMat is a 3d (Nleo x Ntime x Nmeo) matrix with access information for all
% the satellites and epochs

[Nleo,Ntime,Nmeo] = size(accessMat);

% min integer number of epochs to respect tmin
epMin = ceil(tmin/accessTimeStep);
% min integer number of epochs to respect tpoint
epPoint = ceil(tpoint/accessTimeStep);
% downsampled time interval duration
NtimeDs = floor(Ntime/epMin);

%--- Downsample the accessMatrix to get a value for each tmin slot

% Downsampled matrices
W1 = zeros(Nleo,NtimeDs,Nmeo);
W2 = zeros(Nleo,NtimeDs,Nmeo);

% NaN where the access is not avaiable (essential for mean computation)
accessMat(accessMat==0) = nan;

for ii = 1:NtimeDs
    window1 = (ii-1)*epMin+1 : ii*epMin;
    window2 = (ii-1)*epMin+epPoint+1 : ii*epMin;

    % 1st mean: average over the tmin slot
    % if start after or finish before tmin there are NaN in the windows (so
    % mean is NaN)
    B = mean(accessMat(:,window1,:),2);     
    
    % 2nd mean: average over tuseful only
    C = mean(accessMat(:,window2,:),2); 
    % erase values where the tmin slot is not fully covered by availability
    C(isnan(B)) = NaN;
    
    W1(:,ii,:) = B;
    W2(:,ii,:) = C; 
end
clear B C
W1(isnan(W1)) = 0;
W2(isnan(W2)) = 0;

% Now i get 2 "downsampled" matrices where their PDOP values is the mean of
% the PDOP values for the upcoming tmin slot. If the visibility last less
% than the slot (or start after, or finish before), PDOP is set to zero 
% (link not considered for optimal choice) -> it's probably a waste of
% availability window

% W1 contains average PDOP for the whole tmin duration of the
% upcoming slot in case tpoint has already elapsed
% W2 contains average PDOP for the last tuseful time of the
% upcoming slot in case tpoint is still needed



end
