classdef Kernel
  %KERNEL Integration kernel for the Helmholtz equations.
  % Complete set of computations for Green's function kernels
  % and its gradients in time-harmonic electromagnetic waves.
  %
  %      rf       o
  %  (field point) \
  %                 \
  %                  \ R = |D| = |rf - rs|
  %                   \
  %                k   \
  %                     o rs
  %                  (source point)
  %
  %
  %--------------------------------------------------------------------------
  % Copyright: Klearchos A. Samaras, kleasama@gmail.com, June 2023
  %--------------------------------------------------------------------------
  properties
    % The wavenumber of operation [rad/sec].
    k
    % The field point rf [m].
    rf
    % The source point rs [m].
    rs
    % Number of singular terms to be extracted.
    nl
    % The distance vector D = rf - rs [m].
    D
    % The distance R = |D| = |rf - rs| [m].
    R
    % The free space Green's function interaction.
    G
    % The extracted part of the free space Green's function.
    Gext
    % The scalar part of the free space Green's function gradient.
    dG
    % The extracted scalar part of the free space Green's function gradient.
    dGext
  end
  
  methods
    function [ obj ] = Kernel( wavenumber, field_point, source_point, num_extr_terms )
    %%Constructor of Kernel object.
      obj.k  = wavenumber;
      obj.rf = field_point;
      obj.rs = source_point;
      obj.nl = num_extr_terms;
      obj.D = obj.calcD;
      obj.R = obj.calcR;
      obj.G = obj.GF;
      obj.Gext = obj.GFext;
      obj.dG = obj.dGF;
      obj.dGext = obj.dGFext;
    end

    function [ obj ] = processKernel(obj)
    %%processKernel Wrapper of the basic computations for the kernel
    % object.
      obj.D = obj.calcD;
      obj.R = obj.calcR;
      obj.G = obj.GF;
      obj.Gext = obj.GFext;
      obj.dG = obj.dGF;
      obj.dGext = obj.dGFext;
    end
    
    function [ wavenumber ] = getWavenumber(obj)
    %%getWavenumber Get the wavenumber of operation k that is associated
    % with the kernel object.
      wavenumber = obj.k;
    end
    
    function [ obj ] = setWavenumber(wavenumber)
    %%setWavenumber Modify the wavenumber of operation k that is associated
    % with the kernel object.
      obj.k = wavenumber;
    end
    
    function [ field_point ] = getFieldPoint(obj)
    %%getFieldPoint Get the field point rf of the kernel object.
      field_point = obj.rf;
    end
    
    function [ obj ] = setFieldPoint(field_point)
    %%setFieldPoint Modify the field point rf of the kernel object.
      obj.rf = field_point;
    end
    
    function [ source_point ] = getSourcePoint(obj)
    %%getSourcePoint Get the source point rs of the kernel object.
      source_point = obj.rs;
    end
    
    function [ obj ] = setSourcePoint(source_point)
    %%setSourcePoint Modify the source point rs of the kernel object.
      obj.rs = source_point;
    end
    
    function [ num_extr ] = getNl(obj)
    %%getNl Get the number of the ectracted singular terms.
      num_extr = obj.nl;
    end
    
    function [ obj ] = setNl(num_extr)
    %%setNl Modify the number of the ectracted singular terms.
      obj.nl = num_extr;
    end
    
    function [ D ] = calcD(obj)
    %%calcD Evaluation of the distance vector between the field point rf
    % and the source point rs. D = rf - rs.
    %----------------------------------------------------------------------
    % Output:
    %   D - Distance vector rf - rs.
    %----------------------------------------------------------------------
      D = permute(obj.rf, [3, 2, 1]) - obj.rs;
    end
    
    function [d_vec] = getD(obj)
    %%getD Get the distance vector D of the kernel object.
      d_vec = obj.D;
    end
    
    function [obj] = setD(obj, d_vec)
    %%setD Modify the distance vector D of the kernel object.
      obj.D = d_vec;
    end
    
    function [ R ] = calcR(obj)
    %%calcR Evaluation of the distance between the field point rf
    % and the source point rs. R = |D| = |rf - rs|.
    %----------------------------------------------------------------------
    % Output:
    %   R - The distance between the field point rf and the source point.
    %----------------------------------------------------------------------
      R = sqrt(sum(obj.getD .* obj.getD, 2));
    end
    
    function [r_dist] = getR(obj)
    %%getR Get the distance R of the kernel object.
      r_dist = obj.R;
    end
    
    function [obj] = setR(obj, r_dist)
    %%setR Modify the distance R of the kernel object.
      obj.R = r_dist;
    end
    
    function [ G ] = GF(obj)
    %%GF Evaluation of the free space Helmholtz Green's function that is
    % employed as the integration kernel in the field integral equations
    % for time-harmonic electromagnetic waves that propagate in a
    % homogeneous, unbounded space.
    %----------------------------------------------------------------------
    % Output:
    %   G - The Green's function kernel interaction.
    %----------------------------------------------------------------------
      oneover4pi = 0.25 ./ pi;
      G = oneover4pi .* exp( -1i .* obj.k .* obj.getR ) ./ obj.getR;
    end
    
    function [ Gext ] = GFext(obj)
    %%GFext Extracted portion of the free space Helmholtz Green's function
    % intended for singularity extraction techniques.
    %----------------------------------------------------------------------
    % Output:
    %   Gext - The overall extracted part of the Green's function kernel.
    %----------------------------------------------------------------------
      Gext = 0;
      for l = 1 : obj.nl
        n = 2 * (l - 1);
        Gext = Gext + obj.GFsingCoeff(n) .* obj.getR.^(n - 1);
      end
    end
    
    function [ coeff ] = GFsingCoeff(obj, n)
    %%GFsingCoeff Coefficient of the extracted term of the free space
    % Helmholtz Green's function intended for singularity extraction
    % techniques.
    %----------------------------------------------------------------------
    % Input:
    %   n     - The exponent of the extrated term using MacLaurin series.
    %----------------------------------------------------------------------
    % Output:
    %   coeff - Coefficient of the singular term raised to the n-th power.
    %----------------------------------------------------------------------
      oneover4pi = 0.25 ./ pi;
      mjk = -1i .* obj.k;
      coeff = oneover4pi .* mjk.^(n) ./ factorial(n);
    end
    
    function [dG] = dGF(obj)
    %%dGF Evaluation of the scalar part of the gradient of the free space
    % Helmholtz Green's function that is employed as the integration kernel
    % in the field integral equations for time-harmonic electromagnetic
    % waves that propagate in a homogeneous, unbounded space.
    %----------------------------------------------------------------------
    % Output:
    %  dG - The scalar part of the Green's function kernel gradient.
    %----------------------------------------------------------------------
      dG = -(1 + 1i .* obj.k .* obj.getR) ...
         .* obj.G ./ obj.getR.^2;
    end
    
    function [ dGext ] = dGFext(obj)
    %%dGFext Extracted portion of the scalar part of the free space
    % Helmholtz Green's function gradient intended for singularity
    % extraction techniques.
    %----------------------------------------------------------------------
    % Output:
    %   dGext - The overall extracted part of the Green's function kernel.
    %----------------------------------------------------------------------
      dGext = 0;
      for l = 1 : obj.nl
        n = 2 * (l - 1);
        dGext = dGext + obj.dGFsingCoeff(n) .* obj.getR.^(n - 3);
      end
    end
    
    function [ coeff ] = dGFsingCoeff(obj, n)
    %%dGFsingCoeff Coefficient of the extracted term of the free space
    % Helmholtz Green's function gradient intended for singularity
    % extraction techniques.
    %----------------------------------------------------------------------
    % Input:
    %   n     - The exponent of the extrated term using MacLaurin series.
    %----------------------------------------------------------------------
    % Output:
    %   coeff - The coefficient of the singular term raised to the n-th
    %           power.
    %----------------------------------------------------------------------
      coeff = (n - 1) .* obj.GFsingCoeff(n);
    end
  end
end