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

function H = buildHmatrix(azmt, elev, Hin)
% build H matrix from elevation and azimuth.
% H = BUILDHMATRIX(azmt, elev) build H matrix using a set of azimuth and 
% elevation values.
% H = BUILDHMATRIX(azmt, elev, Hin) build H matrix extending Hin matrix
% using an additional set of azimuth and elevation values.

if size(azmt) ~= size(elev)
    error 'Azimuth and elevation sets have different lengths.'
end
if size(azmt,1) < size(azmt,2)
    azmt = azmt.';
    elev = elev.';
end

if nargin < 3
    Hin = [];
end

H = [Hin;...
    cosd(elev).*cosd(azmt) cosd(elev).*sind(azmt) sind(elev) -ones(size(azmt,1),1)]; 

end
