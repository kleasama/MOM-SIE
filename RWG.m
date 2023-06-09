function [rwg] = RWG(t)
%RWG - Create RWG Common Edge Elements Topology from triangle-to-nodes 
%      connectivity (Supports up to 999999 triangles)
%--------------------------------------------------------------------------
% Usage:
%    rwg = RWG(t);
%--------------------------------------------------------------------------
% Input:
%    t   - (TrianglesTotal x 4) - Triangle-to-Nodes connectivity
%--------------------------------------------------------------------------
% Output:
%    rwg - (rwgTotal x 10)      - [Tp,Tm,Ep,Em,n1,n2,Vp,Vm,lbl,idx]
% where the following topology data are block-stored as follows:
%    Tp  - (rwgTotal x 1) - Positive Triangle Index
%    Tm  - (rwgTotal x 1) - Negative triangle Index
%    Ep  - (rwgTotal x 1) - Positive Common Edge Local Index
%    Em  - (rwgTotal x 1) - Negative Common Edge Local Index
%    n1  - (rwgTotal x 1) - Node_1 Index of Common Edge
%    n2  - (rwgTotal x 1) - Node_2 Index of Common Edge
%    Vp  - (rwgTotal x 1) - Positive Free Vertex Index
%    Vm  - (rwgTotal x 1) - Negative Free Vertex Index
%    lbl - (rwgTotal x 1) - Boundary Condition Selection (1-PEC, 2-PEC/Diel, 4-Diel)
%    idx - (rwgTotal x 1) - RWG Basis Function Index
%--------------------------------------------------------------------------
%         - RWG-Element -
%            n2_ _ _ _ _ _ _ _V-
%            /\              /
%           /  \            /
%          /    \     T-   /
%         /    L \        /
%        /        \      /
%       /  T+      \    /
%      /            \  /
%     /_ _ _ _ _ _ _ \/
%    V+               n1
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, December 2022
%--------------------------------------------------------------------------

% Store the triangulation size in TrianglesTotal
TrianglesTotal = length(t(:,1));
% Create the 'edges' block that concatenates the node indeces of edge1, 
% edge2 and edge3 of the triangles, the triangle indeces per edge, the 
% non-unique edge indeces and the free node indeces as follows:
%
% edges =
%     [ node2, node3, triangle_idx, edge1_idx, free_node1]
%     [ node3, node1, triangle_idx, edge2_idx, free_node2]
%     [ node1, node2, triangle_idx, edge3_idx, free_node3]
%
% the 'edges' block size is (3*Triangles Total x 5)
edges = [[t(:,2), t(:,3), (1:TrianglesTotal)',   ones(TrianglesTotal,1), t(:,1)];
         [t(:,3), t(:,1), (1:TrianglesTotal)', 2*ones(TrianglesTotal,1), t(:,2)]; 
         [t(:,1), t(:,2), (1:TrianglesTotal)', 3*ones(TrianglesTotal,1), t(:,3)]];
% Use a compression technique that sums the edge node indeces, shifts them 
% by 9 decimal places through multiplication with 1e9 and then adds the
% absolute value of their subtraction. The resulting integer is appended to
% the 'edges' as a sixth column.
% Compression formula: 
%     EC1 = (edge_node1 + edge_node2)*1e9 + abs(edge_node1 - edge_node2)
edges = [edges,(edges(:,1) + edges(:,2))*1e9 + abs(edges(:,1) - edges(:,2))];
% Use similar compression approach to keep the information of the starting
% and ending edge node indeces in a single integer number. Shift edge_node1
% index by 9 decimal places through multiplication with 1e9 and then add
% the edge_node2 index. Append the resulting integer to the 'edges' block
% as a seventh column.
% Compression formula:
%    EC2 = edge_node1*1e9 + edge_node2
edges = [edges,edges(:,1)*1e9+edges(:,2)];
% Repeat the same compression process to store the integer for the reverse
% case. Namely, shift the edge_node2 index by 9 decimal places and add the
% edge_node1 index. Append the resulting integer to the 'edges' block as an
% eighth column.
% Compression formula:
%    EC2 = edge_node2*1e9 + edge_node1
edges = [edges,edges(:,2)*1e9+edges(:,1)];
% Perform a logical indexing preparation to obtain the places where
% duplicate edges are found in the 'edges' block. Logical_Idx is created
% through dimension broadcasting the 7th and 8th columns and checking for 
% elemnt-wise equality.
Logical_Idx = edges(:,8)==edges(:,7)';
% Sum the Logical_Idx across its columns and append the result to the
% 'edges' block as a ninth column. (Col_09: common_edge_logical_idx)
edges = [edges,sum(Logical_Idx,2)];
% Col_10: other triangle sharing the edge
edges = [edges,sum(Logical_Idx.*edges(:,3)',2)];
% Col_11: other triangle local edge number
edges = [edges,sum(Logical_Idx.*edges(:,4)',2)];
% Col_12: other free node
edges = [edges,sum(Logical_Idx.*edges(:,5)',2)];
% isolate common edges
c_edges = edges(logical(edges(:,9)),:);
% keep unique common edges
[~,idx] = unique(c_edges(:,6));
% rwg = [Tp,Tm,Ep,Em,n1,n2,Vp,Vm,lbl,idx]
rwg = c_edges(idx,[3,10,4,11,1,2,5,12]);
% keep rwg index
rwg = [rwg,t(rwg(:,1),4).*t(rwg(:,2),4),(1:length(rwg(:,1)))'];
end
