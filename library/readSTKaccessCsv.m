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


function [accessData, time] = readSTKaccessCsv(filename,Ntgt,scenarioStartTime)
% leo2meo is a cell array containing:
% - Ntgt cells
% - in each one a cell for each measure (date,az,el,range)
% - each cell is vector with one value per epoch
% Ntgt: is the number of targets of the access report


% counter of target objects of the access report
tgtCnt = 0;

fin = fopen(filename);

accessData = cell(1,Ntgt); %{Ntgt}x{date,el,az,range}x(time)

while ~feof(fin)
    tgtCnt = tgtCnt+1;
    textscan(fin,'%*[^\n]',1,'delimiter',',');
    accessData{tgtCnt} = textscan(fin,'%{dd MMM yyyy HH:mm:ss.SSS}D%f%f%f%*[^\n]','delimiter',',');
    
    % when a L2M pair is over a NAT (not a time) is read at the beginning o
    % the next L2M section: erase it.
    if size(accessData{tgtCnt}{1},1) > size(accessData{tgtCnt}{4},1)
        accessData{tgtCnt}{1}(end) = [];
    end
end

if tgtCnt~=Ntgt
    error 'Data provided lacks at least 1 visibility window for each LEO-MEO pair (acording to Nleo). I cannot tell which.'
end

if nargin > 2
    time = cell(1,Ntgt);
    for ii = 1:Ntgt
        % convert datetime in duration from scenario beginning
        durationSinceStart = accessData{ii}{1} - scenarioStartTime;
        time{ii} = seconds(durationSinceStart);
    end
end

fclose(fin);
end
