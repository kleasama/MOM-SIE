function [gradG] = dGF(k,R)
%%dGF Evaluation of the scalar part of the gradient of the free space
% Helmholtz Green's function that is employed as the integration kernel
% in the field integral equations for time-harmonic electromagnetic waves
% that propagate in a homogeneous, unbounded space.
%--------------------------------------------------------------------------
% Input:
%   k    - The wavenumber of operation.
%   R    - The distance between the field and source point.
%--------------------------------------------------------------------------
% Output:
%   dG   - The scalar part of the Green's function gradient evaluation for
%          the distance R and wavenumber k.
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------
  gradG = -(1 + 1i .* k .* R) .* GF(k, R) ./ R.^2;
end
