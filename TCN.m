function [ tcn ] = TCN( rwg, mesh_size )
%TCN Connectivity cell array from the triangles to the RWG edges.
%--------------------------------------------------------------------------
% Usage:
%    tcn = TCN( rwg, mesh_size );
%--------------------------------------------------------------------------
% Input:
%    rwg - ( rwgTotal x 10 ) - RWG-to-Triangles Connectivity
%--------------------------------------------------------------------------
% Output:
%    tcn - { mesh_size x 3 }
% The connectivity cell array tcn consists of two columns of cells.
%    In tcn{:,1} we store the number of instances of the triangle used.
%    In tcn{:,2} we store the RWG index where the triangle is found.
%    In the case where the triangle is used as the T- then a minus sign is
%    assigned to the concerned RWG indices.
%    In tcn{:,3} we store the local edge index associated with the basis
%    function for which the triangle is acting as the support.
%--------------------------------------------------------------------------
%  Copyright: Klearchos A. Samaras, kleasama@gmail.com, December 2022
%--------------------------------------------------------------------------
tcn = {mesh_size, 3};
for i = 1: mesh_size
  locplus  =  rwg((rwg(:,1) == i),10);
  locminus = -rwg((rwg(:,2) == i),10);
  tcn{i,1} = length(locplus) + length(locminus);
  tcn{i,2} = [locplus; locminus];
  tcn{i,3} = [rwg(locplus,3); rwg(-locminus,4)];
end
