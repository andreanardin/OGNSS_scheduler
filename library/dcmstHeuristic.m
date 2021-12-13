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

function mst = dcmstHeuristic(G,maxNodeDeg,settings)
% DCMSTHEURISTIC computes an approximate solution for the
% degree-constrained minimum spanning tree problem (NP-hard) using the
% heuristic BF1 or BF2 found in:
% Boldon et al., " Minimum-weight degree-constrained spanning tree problem:
% Heuristics and implementation on an SIMD parallel machine", Parallel 
% Computing 22, 1996.
% BF2 is implemented by finding the smallest edge for each node and avoid
% penalizing that edge. However, an edge (u,v) may be the smallest for u
% and for v or just for one of the two. So this information is mixed in the
% f value where the connection to an infeasible node is forgiven if the 
% edge is the smallest for that node.
% A generalization of Boldon et al. is implemented here, allowing to
% compute the NODE degree constrained minimum spanning tree (NDCMST) where
% the degree constraint is specified indipendently for each node. 


%--- Initialization
maxepsilon = settings.maxepsilon;
% initial weight matrix
Winit = adjacency(G,'weighted');
W = Winit;
% initial mst
mst = minspantree(G,'type','forest','method',settings.mstMethod,'Root',randi(size(Winit,1)));
%--- specific for this "cyclic" configuration -----------------------------
mst = restoreM2M(mst,settings);
% (epsilons don't need to be rand here, but remembering erased weight is
% pointless)
%--------------------------------------------------------------------------
D = degree(mst);

% while at least 1 node exceeds its maxNodeDeg
while sum(D > maxNodeDeg) > 0
    
    W_penalty = zeros(numnodes(G));
    
    M = adjacency(mst,'weighted');
    
    % weight range should be referred to PDOP range, thus 'epsilon'
    % weights are excluded (so an edge whose weight = wmin isn't penalized)
    wmin = min(min(M(M>max(0,maxepsilon))));
    wmax = max(max(M));
    
    % infeasible nodes
    iNodes = D > maxNodeDeg;
    
    % is the edge the smallest for that node?
    MM = M;
    MM(MM<=max(0,maxepsilon)) = nan;
    [~,smallest] = min(MM,[],'omitnan');
    
    % for each edge update the weights
    for ii = 1:size(M,1)
        for jj = 1:size(M,2)
            if M(ii,jj)>max(0,maxepsilon)
                % edge feasibility
                f = iNodes(ii)*(smallest(ii)~=jj) +iNodes(jj)*(smallest(jj)~=ii);
                if f~=0
                    % penalty
                    penalty = f*((M(ii,jj)-wmin)/(wmax-wmin))*wmax;
                    if settings.fakeWeights
                        penalty = f*10*randi([0 1]);
                    end
                    W_penalty(ii,jj) = M(ii,jj) + penalty;
                end
            end
        end      
    end
    
    % Update GRAPH weights
    W(W_penalty~=0) = W_penalty(W_penalty~=0);
    
    % Re-compute mst
    G = graph(W,'upper');
    mst = minspantree(G,'type','forest','method',settings.mstMethod,'Root',randi(size(Winit,1)));
    
    % ----- specific for this "cyclic" configuration 
    % rebuild some M2M links (erased by mst construction)
    mst_cyclic = restoreM2M(mst,settings);
    %---------------------------------------------------
    
    D = degree(mst_cyclic);
end

% Restore non-penalized weights
M = adjacency(mst,'weighted');
M(M>max(0,maxepsilon)) = Winit(M>max(0,maxepsilon));
mst = graph(M,'upper');

end




