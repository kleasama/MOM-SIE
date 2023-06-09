function [ G ] = GF(k, R)
%%GF Evaluation of the free space Helmholtz Green's function that is
% employed as the integration kernel in the field integral equations
% for time-harmonic electromagnetic waves that propagate in a homogeneous,
% unbounded space.
%--------------------------------------------------------------------------
% Input:
%   k    - The wavenumber of operation.
%   R    - The distance between the field and source point.
%--------------------------------------------------------------------------
% Output:
%   G - The Green's function evaluation for the distance R and wavenumber k.
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------
  oneover4pi = 0.25 ./ pi;
  G = oneover4pi .* exp( -1i .* k .* R ) ./ R;
end

