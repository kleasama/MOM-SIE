function [ GP ] = GPC( T )
%GPC - Geometric pre-calculations over planar triangles. Create a structure
% that acts as the container for a collection of geometric precalculated
% quantities for a planar triangle.
%--------------------------------------------------------------------------
% Usage:
%    GP = GPC( T );
%--------------------------------------------------------------------------
% Input:
%    T  - Triangle object.
%--------------------------------------------------------------------------
% Output:
%    GP - (MeshSize x 29) - [A,L1,L2,L3,n,s1,s2,s3,m1,m2,m3,c,SId]
%
% where the following geometric data are block-stored as follows:
%    A  - (MeshSize x 1) -> GP(:, 1)    - Triangle Area
%    L1 - (MeshSize x 1) -> GP(:, 2)    - Length of Edge 1
%    L2 - (MeshSize x 1) -> GP(:, 3)    - Length of Edge 2
%    L3 - (MeshSize x 1) -> GP(:, 4)    - Length of Edge 3
%    n  - (MeshSize x 3) -> GP(:, 5: 7) - Face Normal Unit Vector (Counterclockwise)
%    s1 - (MeshSize x 3) -> GP(:, 8:10) - Tangential Unit Vector of Edge 1
%    s2 - (MeshSize x 3) -> GP(:,11:13) - Tangential Unit Vector of Edge 2
%    s3 - (MeshSize x 3) -> GP(:,14:16) - Tangential Unit Vector of Edge 3
%    m1 - (MeshSize x 3) -> GP(:,17:19) - Outward Normal Unit Vector of Edge 1
%    m2 - (MeshSize x 3) -> GP(:,20:22) - Outward Normal Unit Vector of Edge 2
%    m3 - (MeshSize x 3) -> GP(:,23:25) - Outward Normal Unit Vector of Edge 3
%    c  - (MeshSize x 3) -> GP(:,26:28) - Triangle Centroid
%    SId- (MeshSize x 1) -> GP(:,29)    - Self Term Integration I = ∫∫{∫∫(1/R)ds'}ds
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------

  % Block all the computed data in a single array.
  GP  = [ T.getArea, T.getEdgeLength(1), T.getEdgeLength(2), T.getEdgeLength(3),...
          T.getFaceNormal, T.getEdgeUnit(1), T.getEdgeUnit(2), T.getEdgeUnit(3),...
          T.getEdgeNormal(1), T.getEdgeNormal(2), T.getEdgeNormal(3), T.getCentroid,...
          T.getSelfTermId ];
end
