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

function dop = computeDOP(H,dopType)
% Compute Dilution Of Precision (DOP).
% computeDOP(H,dopType) computes DOP according to dopType, which can be
% 'gdop', 'pdop'/'pdop3', 'hdop'/'hdop3', 'vdop'/'vdop3', 'tdop'. If *dop3
% is selected the DOP value is computed even if only 3 navigation sources
% are available. This is done by ignoring the clock component of the 
% navigation solution.

lastwarn(''); % empty lastwarn
warning('off','MATLAB:singularMatrix');

if nargin < 2
    dopType = 'gdop';
end

dopType = lower(dopType);
if strcmp(dopType,'pdop3') || strcmp(dopType,'hdop3') || strcmp(dopType,'vdop3') 
    H(:,4) = [];
end

G = inv((H.'*H));

% exclude singular matrices from minimum search
[warnMsg, warnId] = lastwarn;
if strcmp(warnId,'MATLAB:singularMatrix')
    g = Inf*eye(size(H,2));
else
    g = diag(G);   
end


switch dopType
    case {'pdop','pdop3'}
        dop = sqrt(sum(g(1:3)));
    case {'hdop','hdop3'}
        dop = sqrt(sum(g(1:2)));
    case {'vdop','vdop3'}
        dop = sqrt(sum(g(3)));
    case 'tdop'
        dop = sqrt(sum(g(4)));
    case 'gdop'
        dop = sqrt(sum(g));
    otherwise
        error 'unknown DOP type'
end

end
