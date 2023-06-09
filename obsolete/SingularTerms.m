function [ Im, K1, K2v1, K2v2, K2v3, K3, K4v1, K4v2, K4v3 ] = SingularTerms( T, r, N )
%SingularTerms Analytic evaluation of integrations containing the singuar
%terms of the Helmholtz free space Green's function Taylor series expansion
%around R = 0 using recursive formulas.
%
% INPUT
% T - Triangle Object.
% r - Field point vector with x,y,z components along its columns.
% N - Number of extracted singular terms.
%
% OUTPUT
% Im   - Contour integral ∫  Rⁿ dl'.
% K1   - Surface integral ∫∫ Rⁿ ds'.
% K2v1 - Surface integral ∫∫ Rⁿ(r' - v₁) ds'.
% K2v2 - Surface integral ∫∫ Rⁿ(r' - v₂) ds'.
% K2v3 - Surface integral ∫∫ Rⁿ(r' - v₃) ds'.
% K3   - Surface integral ∫∫ ∇Rⁿ ds'.
% K4v1 - Surface integral ∫∫ ∇Rⁿ(r' - v₁) ds'.
% K4v2 - Surface integral ∫∫ ∇Rⁿ(r' - v₂) ds'.
% K4v3 - Surface integral ∫∫ ∇Rⁿ(r' - v₃) ds'.
% where n = -1,1,3,...,(2N-1)
%--------------------------------------------------------------------------
% REFERENCES:
% [1] - P. Ylä-Oijala, M. Taskinen, "Calculation of CFIE Impedance Matrix
%       Elements With RWG and nxRWG Functions", IEEE Transactions on
%       Antennas and Propagation, Vol. 51, No. 8, August 2003.
%
%--------------------------------------------------------------------------
% Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
%--------------------------------------------------------------------------

  %% TRIANGLE PARAMETERS
  m1  = T.getEdgeNormal(1);
  m2  = T.getEdgeNormal(2);
  m3  = T.getEdgeNormal(3);
  n   = T. getFaceNormal;
  %% INTERACTION WITH FIELD POINT
  w0  = T.getW0(r);
  rpj = T.getProjection(r);
  s1m = T.getSm(r, 1);
  s1p = T.getSp(r, 1);
  s2m = T.getSm(r, 2);
  s2p = T.getSp(r, 2);
  s3m = T.getSm(r, 3);
  s3p = T.getSp(r, 3);
  t01 = T.getT0(r, 1);
  t02 = T.getT0(r, 2);
  t03 = T.getT0(r, 3);
  R01 = T.getR0(r, 1);
  R02 = T.getR0(r, 2);
  R03 = T.getR0(r, 3);
  R1p = T.getRp(r, 1);
  R2p = T.getRp(r, 2);
  R3p = T.getRp(r, 3);
  R1m = T.getRm(r, 1);
  R2m = T.getRm(r, 2);
  R3m = T.getRm(r, 3);
  %% SINGULAR INTEGRALS - INITIAL TERMS
  I1m1 = log( (R1p + s1p)./(R1m + s1m) );
  I2m1 = log( (R2p + s2p)./(R2m + s2m) );
  I3m1 = log( (R3p + s3p)./(R3m + s3m) );
  b1 = atan2( (t01.*s1p),(R01.^2 + abs(w0).*R1p) ) - atan2( (t01.*s1m),(R01.^2 + abs(w0).*R1m) );
  b2 = atan2( (t02.*s2p),(R02.^2 + abs(w0).*R2p) ) - atan2( (t02.*s2m),(R02.^2 + abs(w0).*R2m) );
  b3 = atan2( (t03.*s3p),(R03.^2 + abs(w0).*R3p) ) - atan2( (t03.*s3m),(R03.^2 + abs(w0).*R3m) );
  w0K1m3 = sign(w0).*(b1 + b2 + b3);
  w0K1m3(w0==0) = 0;
  %% SINGULAR INTEGRALS - RECURSIVE FORMULAS
  % Memory allocation
  I1   = zeros(T.getTrianglesTotal,1,N+1);
  I2   = zeros(T.getTrianglesTotal,1,N+1);
  I3   = zeros(T.getTrianglesTotal,1,N+1);
  K1   = zeros(T.getTrianglesTotal,1,N);
  K2v1 = zeros(T.getTrianglesTotal,3,N);
  K2v2 = zeros(T.getTrianglesTotal,3,N);
  K2v3 = zeros(T.getTrianglesTotal,3,N);
  K3   = zeros(T.getTrianglesTotal,3,N);

  % Contour Integrals
  I1(:,:,1) = I1m1;
  I2(:,:,1) = I2m1;
  I3(:,:,1) = I3m1;
  for i = 2:N+1
    ni = 2*(i-1) - 1;
    I1(:,1,i) = 1./(ni + 1) .* ( s1p.*R1p.^ni - s1m.*R1m.^ni + ni.*R01.^2.*I1(:,1,i-1) );
    I2(:,1,i) = 1./(ni + 1) .* ( s2p.*R2p.^ni - s2m.*R2m.^ni + ni.*R02.^2.*I2(:,1,i-1) );
    I3(:,1,i) = 1./(ni + 1) .* ( s3p.*R3p.^ni - s3m.*R3m.^ni + ni.*R03.^2.*I3(:,1,i-1) );
  end
  Im = m1.*I1 + m2.*I2 + m3.*I3;

  % Surface Integrals
  K1(:,:,1)     = t01.*I1m1 + t02.*I2m1 + t03.*I3m1;
  K1(w0~=0,:,1) = K1(w0~=0,:,1) - w0(w0~=0).*w0K1m3(w0~=0);
  K2v1(:,:,1)   = Im(:,:,2) + (rpj - T.getP(1)).*K1(:,:,1);
  K2v2(:,:,1)   = Im(:,:,2) + (rpj - T.getP(2)).*K1(:,:,1);
  K2v3(:,:,1)   = Im(:,:,2) + (rpj - T.getP(3)).*K1(:,:,1);
  K3(:,:,1)     = -w0K1m3.*n - Im(:,:,1);
  for i = 2:N
    ni = 2*(i-1) - 1;
    K1(:,:,i)     = 1./(ni + 2) .* ( t01.*I1(:,:,i) + t02.*I2(:,:,i) + t03.*I3(:,:,i) );
    K1(w0~=0,:,i) = K1(w0~=0,:,i) + 1./(ni + 2) .* ( ni.*(w0(w0~=0)).^2.*K1(w0~=0,:,i-1) );
    K2v1(:,:,i)   = 1./(ni + 2) .* Im(:,:,i+1) + T.getRho(rpj, 1).*K1(:,:,i);
    K2v2(:,:,i)   = 1./(ni + 2) .* Im(:,:,i+1) + T.getRho(rpj, 2).*K1(:,:,i);
    K2v3(:,:,i)   = 1./(ni + 2) .* Im(:,:,i+1) + T.getRho(rpj, 3).*K1(:,:,i);
    K3(:,:,i)     = ni.*w0.*n.*K1(:,:,i-1) - Im(:,:,i);
  end
  K4v1 = -cross(repmat((r - T.getP(1)),[1,1,N]),K3,2);
  K4v2 = -cross(repmat((r - T.getP(2)),[1,1,N]),K3,2);
  K4v3 = -cross(repmat((r - T.getP(3)),[1,1,N]),K3,2);
end
