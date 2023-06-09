function [ Cn ] = sphereCMs(ur, tandu, sigmau, er, tande, sigmae, radius, isPEC, fL, fU, nmax)
% sphereCMs computes the characteristic eigenvalues Cn for a multilayered sphere
% having a number of homogeneous dielectric layers, with or without a perfectly
% electrically conducting (PEC) core. It also plots the characteristic modes
% and stores the real and imaginary parts of Cn in the TXT comma delimited files
% "ev_re.txt" and "ev_im.txt".
%
% Inputs:
%   ur: Vector containing the relative permeability for each dielectric region
%        from outermost to innermost (the unbounded/free space region should
%        be the first entry).
%   tandu: Vector containing the magnetic loss tangent.
%   sigmau: Vector containing the magnetic conductivity.
%   er: Vector containing the relative permittivity for each dielectric region
%        from outermost to innermost (the unbounded/free space region should
%        be the first entry).
%   tande: Vector containing the dielectric loss tangent.
%   sigmae: Vector containing the electric conductivity.
%   radius: Vector containing the radius of each spherical dielectric interface,
%        from outermost to innermost.
%   isPEC: A flag for when the innermost sphere (core) is conducting.
%   fL: Starting frequency sample (Hz).
%   fU: Ending frequency sample (Hz).
%   nmax: Maximum mode for computing Bessel functions.
%
% Outputs:
%   Cn: Array of the characteristic eigenvalues for the TE and TM modes of the
%        multilayered spherical structure.
%
% Usage:
%-------------------------------------------------------------------------------
% >> EXAMPLE #1 - Lossy dielectric penetrable sphere with er = 12, tande = 0.02
%         and radius = 1m.
%-------------------------------------------------------------------------------
% ur = [1, 1];
% tandu = [0, 0];
% sigmau = [0, 0];
% er = [1, 12];
% tande = [0, 0.02];
% sigmae = [0, 0];
% radius = 1;
% isPEC = 0;
% fL = 23.86e6;
% fU = 71.57e6;
% nmax = 3;
%
% [Cn] = sphereCMs(ur, tandu, sigmau, er, tande, sigmae, radius, isPEC, fL, fU, nmax);
%
%-------------------------------------------------------------------------------
% >> EXAMPLE #2 - Lossy conducting spherical shell with surface resistance
%         rs = 30, thickness ts = 1e-3m and radius = 1m.
%-------------------------------------------------------------------------------
% rs = 30;
% ts = 1e-3;
% ur = [1, 1];
% tandu = [0, 0];
% sigmau = [0, 0];
% er = [1, 1];
% tande = [0, 0];
% sigmae = [0, 1/(ts*rs)];
% radius = [1 + 0.5*ts, 1 - 0.5*ts];
% isPEC = 0;
% fL = 40.00e6;
% fU = 150.00e6;
% nmax = 5;
%
% [Cn] = sphereCMs(ur, tandu, sigmau, er, tande, sigmae, radius, isPEC, fL, fU, nmax);
%
%-------------------------------------------------------------------------------
% >> EXAMPLE #3 - Multilayered spherical structure consisting of four lossy
%         dielectric and magnetic compartments.
%         ■ Compartment 1: er = 7, ur = 4, tande = 0.01, tandu = 0.01, r = 1.00m
%         ■ Compartment 2: er = 6, ur = 3, tande = 0.02, tandu = 0.02, r = 0.75m
%         ■ Compartment 3: er = 3, ur = 9, tande = 0.03, tandu = 0.05, r = 0.50m
%         ■ Compartment 4: er = 5, ur = 2, tande = 0.04, tandu = 0.03, r = 0.25m
%-------------------------------------------------------------------------------
% ur = [1, 4, 3, 9, 2];
% tandu = [0.00, 0.01, 0.02, 0.05, 0.03];
% sigmau = [0, 0, 0, 0, 0];
% er = [1, 7, 6, 3, 5];
% tande = [0.00, 0.01, 0.02, 0.03, 0.04];
% sigmae = [0, 0, 0, 0, 0];
% radius = [1.00, 0.75, 0.50, 0.25];
% isPEC = 0;
% fL = 23.86e6;
% fU = 71.57e6;
% nmax = 3;
%
% [Cn] = sphereCMs(ur, tandu, sigmau, er, tande, sigmae, radius, isPEC, fL, fU, nmax);
%
%-------------------------------------------------------------------------------
% >> EXAMPLE #4 - Lossy conducting spherical shell with surface resistance
%         rs = 30, thickness ts = 1e-3m and radius = 1.0m, coated with a thin
%         lossy dielectric and magnetic layer of thickness tl = 0.07m, ur = 2,
%         tandu = 0.13, er = 4, tande = 0.21.
%-------------------------------------------------------------------------------
% tl = 0.07;
% rs = 30;
% ts = 1e-3;
% ur = [1, 2, 1];
% tandu = [0.00, 0.13, 0.00];
% sigmau = [0, 0, 0];
% er = [1, 4, 1];
% tande = [0.00, 0.21, 0.00];
% sigmae = [0, 0, 1/(ts*rs)];
% radius = [1 + 0.5*ts + tl, 1 + 0.5*ts, 1 - 0.5*ts];
% isPEC = 0;
% fL = 40.00e6;
% fU = 150.00e6;
% nmax = 5;
%
% [Cn] = sphereCMs(ur, tandu, sigmau, er, tande, sigmae, radius, isPEC, fL, fU, nmax);
%

  % Create a lineraly spaced discrete frequency sample vector.
  frequency = linspace(fL,fU,501);
  % Compute the size parameter based on the maximu_rm radius.
  x = 2*pi*frequency/2.99792458e8*max(radius);
  % Compute the Mie scattering coefficients for the TE and TM modes. Call the
  % function mieScatCoeffs for every frequency sample.
  sigmauOVERtwopimu_r      = sigmau ./ ( 2 * pi * ur * 4 * pi * 1e-7 );
  sigmaeOVERtwopiepsilon_r = sigmae ./ ( 2 * pi * er * 8.8541878128e-12 );
  for i=1:length(frequency)
      f = frequency(i);
      mu_r      = ur .* ( 1 - 1i * (tandu + sigmauOVERtwopimu_r/f) );
      epsilon_r = er .* ( 1 - 1i * (tande + sigmaeOVERtwopiepsilon_r/f) );
      [An(i,:), Bn(i,:)] = mieScatCoeffs(mu_r, epsilon_r, radius, isPEC, f, nmax);
  end
  % Create a container array to store all the resulting eigenvalues starting
  % from the TE modes and then appending the TM modes.
  Cn = [An, Bn].';
  Cn = (1 - Cn)./(1i*Cn);

  % Output the real and imaginary parts of the eigenvalues in the delimited TXT
  % files "ev_re.txt" and "ev_im.txt".
  strTE = '';
  strTM = '';
  for i = 1:nmax
    strTE = [strTE,',TE',num2str(i)];
    strTM = [strTM,',TM',num2str(i)];
  end
  str = ['freq',strTE,strTM];
  dlmwrite('ev_re.txt', str, '');
  dlmwrite('ev_re.txt', [frequency',real(Cn).'],"-append");
  dlmwrite('ev_im.txt', str, '');
  dlmwrite('ev_im.txt', [frequency',imag(Cn).'],"-append");

  % Create a figure that will contain a 3x2 array of subplots. Plot the
  % eigenvalues and the derived quantities such as modal significance and
  % characteristic angle.
  figure('Name','Eigenvalue Spectrum');

  % Real part of the eigenvalues in linear scale.
  subplot(3,2,1);plot(x, real(Cn), 'LineWidth', 1.3);
  grid on;xlim([x(1),x(end)]);ylim([-15,15]);
  xlabel('Size parameter, x = k_0\alpha','FontSize',11);
  ylabel('real(\lambda)','FontSize',11);

  % IMaginary part of the eigenvalues in linear scale.
  subplot(3,2,2);plot(x, imag(Cn), 'LineWidth', 1.3);
  grid on;xlim([x(1),x(end)]);ylim([-15,15]);
  xlabel('Size parameter, x = k_0\alpha','FontSize',11);
  ylabel('imag(\lambda)','FontSize',11);

  % Real part of the eigenvalues in log scale.
  subplot(3,2,3);plot(x, abs(real(Cn)), 'LineWidth', 1.3);
  set(gca,'YScale','log');grid on;xlim([x(1),x(end)]);ylim([1e-2,1e4]);
  xlabel('Size parameter, x = k_0\alpha','FontSize',11);
  ylabel('|real(\lambda)|','FontSize',11);

  % Imaginary part of the eigenvalues in log scale.
  subplot(3,2,4);plot(x, abs(imag(Cn)), 'LineWidth', 1.3);
  set(gca,'YScale','log');grid on;xlim([x(1),x(end)]);ylim([1e-2,1e4]);
  xlabel('Size parameter, x = k_0\alpha','FontSize',11);
  ylabel('|imag(\lambda)|','FontSize',11);

  % Modal significance.
  subplot(3,2,5);plot(x, 1./abs(1 + 1i*Cn), 'LineWidth', 1.3);
  grid on;xlim([x(1),x(end)]);ylim([0,1]);
  xlabel('Size parameter, x = k_0\alpha','FontSize',11);
  ylabel('MS_n','FontSize',11);

  % Characteristic angle.
  subplot(3,2,6);plot(x, 180 - rad2deg(atan(real(Cn))), 'LineWidth', 1.3);
  grid on;xlim([x(1),x(end)]);ylim([90,270]);
  xlabel('Size parameter, x = k_0\alpha','FontSize',11);
  ylabel('\alpha_n','FontSize',11);
end

function [An, Bn] = mieScatCoeffs(mu_r, epsilon_r, radius, isPEC, frequency, nmax)
% Compute the sphere coefficients An and Bn for a sphere having a number
% of homogeneous dielectric layers, with or without a perfectly
% electrically conducting (PEC) core.
% Follows the treatment in Chapter 3 of
%
% Ruck, et. al. "Radar Cross Section Handbook", Plenum Press, 1970.
%
% The incident electric field is in the -z direction (theta = 0) and is
% theta-polarized. The time-harmonic convention exp(jwt) is assumed, and
% the Green's function is of the form exp(-jkr)/r.
%
% Inputs:
%   mu_r: Relative complex permeability for each dielectric region (the
%        unbounded/free space region should be the first entry)
%   epsilon_r: Relative complex permittivity for each dielectric region (the
%        unbounded/free space region should be the first entry)
%   radius: Vector containing the radius of each spherical dielectric interface,
%        from outermost to innermost
%   isPEC: A flag for when the innermost sphere (core) is conducting
%   frequency: Operating frequency (Hz)
%   nmax: Maximum mode for computing Bessel functions
% Outputs:
%   An: Array of Mie solution constants
%   Bn: Array of Mie solution constants

  % speed of light
  c = 299792458.0;

  % radian frequency
  w = 2.0*pi*frequency;

  % wavenumber
  k0 = w/c;

  % total number of dielectric interfaces
  numInterfaces = length(mu_r);

  % free space impedance
  eta0 = 376.7303134617707;

  % impedance of each dielectric region
  eta = eta0 * sqrt(mu_r./epsilon_r);

  % scale factor
  m = sqrt(mu_r.*epsilon_r);

  % mode numbers
  mode = 1:nmax;
  Z = zeros(numInterfaces, length(mode));
  Y = zeros(numInterfaces, length(mode));
  numIFC = numInterfaces - 1;
  for L = numIFC:-1:1
    if L == numIFC
      if isPEC == 0
        x = k0 * m(L + 1) * radius(L);
        [J JP] = SphericalBessel_JJP(mode, x);

        % Ruck, et. al. (3.4-23)
        P = (J ./ JP);

        % Ruck, et. al. (3.4-26), at innermost interface
        Z(L,:) = eta(L + 1) * P;

        % Ruck, et. al. (3.4-26), at innermost interface
        Y(L,:) = P / eta(L + 1);
      end
    else
      x1 = k0 * m(L + 1) * radius(L);
      x2 = k0 * m(L + 1) * radius(L + 1);
      [J1 JP1] = SphericalBessel_JJP(mode, x1);
      [H1 HP1] = SphericalHankel_HHP(mode, 2, x1);
      [J2 JP2] = SphericalBessel_JJP(mode, x2);
      [H2 HP2] = SphericalHankel_HHP(mode, 2, x2);

      % Ruck, et. al. (3.4-24)
      U = JP2 .* HP1 ./ (JP1 .* HP2);

      % Ruck, et. al. (3.4-25)
      V = (J2 .* H1) ./ (J1 .* H2);

      % Ruck, et. al. (3.4-23)
      P1 = J1 ./ JP1;

      % Ruck, et. al. (3.4-23)
      P2 = J2 ./ JP2;

      % Ruck, et. al. (3.4-23)
      Q2 = H2 ./ HP2;
      if (L == numIFC - 1) && (isPEC == 1)
        % Ruck, et. al. (3.4-27) for PEC boundary condition
        Z(L,:) = eta(L + 1) * P1 .* (1.0 - V) ./ (1.0 - U .*  P2 ./ Q2);

        % Ruck, et. al. (3.4-28) for PEC boundary condition
        Y(L,:) = (P1 / eta(L + 1)) .* (1.0 - V .*  Q2 ./ P2) ./ (1.0 - U);
      else
        % Ruck, et. al. (3.4-27)
        Z(L,:) = eta(L + 1) * P1 .* (1.0 - V .*((1.0 - Z(L+1,:)...
             ./(eta(L+1)*P2))./(1.0 - Z(L+1,:)./(eta(L+1)*Q2))))...
             ./(1.0 - U .* ((1.0 - eta(L+1)*P2./Z(L+1,:))...
             ./(1.0 - eta(L+1)*Q2./Z(L+1,:))));

        % Ruck, et. al. (3.4-28)
        Y(L,:) = (P1 / eta(L + 1)) .* (1.0 - V .*((1.0 - eta(L+1)* Y(L+1,:)...
             ./ P2)./(1.0 - eta(L+1)*Y(L+1,:)./Q2)))...
             ./(1.0 - U .* ((1.0 - P2./(eta(L+1)*Y(L+1,:)))...
             ./(1.0 - Q2./(eta(L + 1)*Y(L+1,:)))));
      end
    end
  end
  % Ruck, et. al. (3.4-29)
  Zn = i*Z(1,:)/eta0;

  % Ruck, et. al. (3.4-29)
  Yn = i*eta0*Y(1,:);
  x = k0*radius(1);
  [J JP] = SphericalBessel_JJP(mode, x);
  [H HP] = SphericalHankel_HHP(mode, 2, x);

  % Ruck, et. al. (3.4-1)
  %An = -((i).^(mode)) .* (2*mode + 1) ./ (mode.*(mode + 1))...
  %     .* (J + i*Zn .* JP) ./ (H + i*Zn .* HP);
  An = (J + i*Zn .* JP) ./ (H + i*Zn .* HP);

  % Ruck, et. al. (3.4-2) - there is an error in Ruck which is fixed here
  %Bn = ((i).^(mode+1)) .* (2*mode + 1) ./ (mode.*(mode + 1))...
  %     .* (J + i*Yn .* JP) ./ (H + i*Yn .* HP);
  Bn = (J + i*Yn .* JP) ./ (H + i*Yn .* HP);

  function [J JP] = SphericalBessel_JJP(mode, x)
    % compute spherical bessel functions and their derivatives
    s = sqrt(0.5*pi/x);
    [J] = besselj(mode + 1/2, x); J = J*s;
    [J2] = besselj(mode - 1/2, x); J2 = J2*s;
    JP = (x * J2 - mode .* J);
    J = x*J;
  end

  function [H HP] = SphericalHankel_HHP(mode, arg, x)
    % compute spherical hankel functions and their derivatives
    s = sqrt(0.5*pi/x);
    [H] = besselh(mode + 1/2, arg, x); H = H*s;
    [H2] = besselh(mode - 1/2, arg, x); H2 = H2*s;
    HP = (x * H2 - mode .* H);
    H = x*H;
  end
end
