% Node coordinate list.
p = [ 0.00000, 0.21651,-0.12500;...
     -0.18750, 0.10825,-0.12500;...
      0.00000, 0.21651, 0.12500;...
     -0.18750, 0.10825, 0.12500;...
      0.00000, 0.00000,-0.25000];
% Triangle-to-node connectivity.
t = [1, 2, 3, 2, 1;...
     3, 2, 4, 2, 1;...
     1, 5, 2, 2, 1];
% Build RWG connectivity.
%      T+   T-   le+  le-  gn1  gn2  gfv+ gfv- BC   basis_idx
rwg = [1,   2,   1,   3,   2,   3,   1,   4,   4,   1;...
       1,   3,   3,   2,   2,   1,   3,   5,   4,   2];
tcn = TCN(rwg, length(t(:,1)));
mesh_size = length(t(:,1));
% Wavenumber
k = 2*pi*23.86e6 / 299792458;