clear;clc;close all;addpath('..\');
% Load a model.
airbody;
% load dipole
% mesh_size = length(t(:,1));
% rwg = RWG(t);
% tcn = TCN(rwg, mesh_size);
% k = 50;
% Input parameters.
% Number of field integration points.
Nf = 7;
% Number of source integration points.
Ns = 7;
% Singularity extraction options.
Nl = 2;
force_doextract = 0;
never_doextract = 0;
% Visualisation options.
visualise_mesh = 0;
show_progress = 1;
% Debug options.
enable_debug = 1;
debug_mg = 1;
debug_ng = 2;
debug_mb = 1;
debug_nb = 2;
visualise_debug_tris = 0;
% Plot mesh.
if visualise_mesh == 1; viewer(p, t); end
% Start timer.
tic;
% Set the thresholds for the singularity extraction.
doextract_thrsh_L = (0.1 * pi / real(k))^2;
doextract_thrsh_K = (0.4 * pi / real(k))^2;
% Memory allocation.
L  = zeros(rwg(end, end)) + 1i * zeros(rwg(end, end));
Ln = zeros(rwg(end, end)) + 1i * zeros(rwg(end, end));
K  = zeros(rwg(end, end)) + 1i * zeros(rwg(end, end));
Kn = zeros(rwg(end, end)) + 1i * zeros(rwg(end, end));
disp('>> Looping over triangles for the efficient matrix fill');
% Field triangle loop.
for mg = 1 : mesh_size
  % Create the field triangle object.
  [ Tf ] = Triangle(p(t(mg,1),:), p(t(mg,2),:), p(t(mg,3),:));
  % Create the field integration points.
  [ rf, wf ] = GQ(Tf.getP(1), Tf.getP(2), Tf.getP(3), Nf);
  % Pull the number of basis functions for the mg-th triangle.
  mg_basis_num = tcn{mg,1};
  % Source triangle loop.
  for ng = 1 : mesh_size
    % Create the source triangle object.
    [ Ts ] = Triangle(p(t(ng,1),:), p(t(ng,2),:), p(t(ng,3),:));
    % Centroid distance between field and source triangles.
    tfts_dist = sum((Tf.getCentroid - Ts.getCentroid).^2, 2);
    % Control the singularity extraction if needed.
    if (force_doextract == 1)
      doextract_L = 1;
      doextract_K = 1;
    elseif (never_doextract == 1)
      doextract_L = 0;
      doextract_K = 0;
    elseif (tfts_dist < doextract_thrsh_K)
      doextract_K = 1;
      if (tfts_dist < doextract_thrsh_L)
        doextract_L = 1;
      else
        doextract_L = 0;
      end
    else
      doextract_L = 0;
      doextract_K = 0;
    end
    % Create the source integration points.
    [ rs, ws ] = GQ(Ts.getP(1), Ts.getP(2), Ts.getP(3), Ns);
    % Pull the number of basis functions for the ng-th triangle.
    ng_basis_num = tcn{ng,1};
    if (enable_debug == 1)
      if and(mg == debug_mg, ng == debug_ng)
        disp('-------------------------------------------------------------------------------');
        disp('>> Interacting triangles:');
        disp(['   T_mg = ', num2str(mg), ' comprised of global nodes:']);
        disp(['   n1 = ', num2str(t(mg,1)), ':   ',num2str(Tf.getP(1),[' X = %5.5f,', ' Y = %5.5f,', ' Z = %5.5f'])]);
        disp(['   n2 = ', num2str(t(mg,2)), ':   ',num2str(Tf.getP(2),[' X = %5.5f,', ' Y = %5.5f,', ' Z = %5.5f'])]);
        disp(['   n3 = ', num2str(t(mg,3)), ':   ',num2str(Tf.getP(3),[' X = %5.5f,', ' Y = %5.5f,', ' Z = %5.5f'])]);
        disp(['   A total of ',  num2str(mg_basis_num),' basis functions is associated with T_mg = ', num2str(mg)]);
        disp(' ');
        disp(['   T_ng = ', num2str(ng), ' comprised of global nodes:']);
        disp(['   n1 = ', num2str(t(ng,1)), ':   ',num2str(Ts.getP(1),[' X = %5.5f,', ' Y = %5.5f,', ' Z = %5.5f'])]);
        disp(['   n2 = ', num2str(t(ng,2)), ':   ',num2str(Ts.getP(2),[' X = %5.5f,', ' Y = %5.5f,', ' Z = %5.5f'])]);
        disp(['   n3 = ', num2str(t(ng,3)), ':   ',num2str(Ts.getP(3),[' X = %5.5f,', ' Y = %5.5f,', ' Z = %5.5f'])]);
        disp(['   A total of ',  num2str(ng_basis_num),' basis functions is associated with T_ng = ', num2str(mg)]);
        disp('-------------------------------------------------------------------------------');
        if and(visualise_mesh == 1, visualise_debug_tris == 1)
          hold on;
          % Plot source and field integration points.
          plot3(rs(:,1), rs(:,2), rs(:,3), '*r');
          plot3(rf(:,1), rf(:,2), rf(:,3), '*b');
          hold off;
        end
      elseif or(debug_mg > mesh_size, debug_mg < 0)
        error('Incorrect values for debug_mg');
      elseif or(debug_ng > mesh_size, debug_ng < 0)
        error('Incorrect values for debug_ng');
      end
    end
    % Kernel object creation for the calculation of the Green's function
    % interactions and its gradient.
    KER = Kernel(k, rf, rs, Nl);
    % Loop over the associated local edges to the mg-th triangle.
    for local_mj = 1 : mg_basis_num
      % Associated field basis function.
      mg_basis_m   = tcn{mg,2}(local_mj);
      % Associated local edge with the field basis.
      mj = tcn{mg,3}(local_mj);
      % Loop over the associated local edges to the ng-th triangle.
      for local_nj = 1 : ng_basis_num
        % Associated source basis function.
        ng_basis_n = tcn{ng,2}(local_nj);
        % Associated local edge with the source basis.
        nj = tcn{ng,3}(local_nj);
        % Projection integral used in the testing procedure.
        Id  = Idmn(KER, Tf, Ts, mg, ng, mj, nj, wf, ws, Nf, Ns, doextract_L);
        Idn = 0;
        Is  = Ismn(KER, Tf, Ts, mg, ng, mj, nj, wf, ws, Nf, Ns, doextract_L);
        Isn = Isnmn(KER, Tf, Ts, mg, ng, mj, nj, wf, ws, Nf, Ns, doextract_L);
        Ic  = Icmn(KER, Tf, Ts, mg, ng, mj, nj, wf, ws, Nf, Ns, doextract_K);
        Icn = Icnmn(KER, Tf, Ts, mg, ng, mj, nj, wf, ws, Nf, Ns, doextract_K);
        Ii  = Iimn(Tf, mg, ng, mj, nj, +1);
        Iin = Iinmn(Tf, mg, ng, mj, nj, +1);
        % Debug output.
        if (enable_debug == 1)
          if and(mg == debug_mg, ng == debug_ng)
            disp(['   Basis-index = ', num2str(tcn{mg,2}(local_mj)), ...
              ' and local edge index = ', num2str(mj)]);
            disp(['   Basis-index = ', num2str(tcn{ng,2}(local_nj)), ...
              ' and local edge index = ', num2str(nj)]);
            disp(['   Scalar potential interaction:   Id  = ', ...
              num2str([real(Id), imag(Id)], ['%4.12e  ','%4.12ei'])]);
            disp(['   Vector potential interaction:   Is  = ', ...
              num2str([real(Is), imag(Is)], ['%4.12e  ','%4.12ei'])]);
            disp(['   nxVector potential interaction: Isn = ', ...
              num2str([real(Isn), imag(Isn)], ['%4.12e  ','%4.12ei'])]);
            disp(['   Cross potential interaction:    Ic  = ', ...
              num2str([real(Ic), imag(Ic)], ['%4.12e  ','%4.12ei'])]);
            disp(['   nxCross potential interaction:  Icn = ', ...
              num2str([real(Icn), imag(Icn)], ['%4.12e  ','%4.12ei'])]);
            disp(['   Principal Value interaction:    Ii  = ', ...
              num2str([real(Ii), imag(Ii)], ['%4.12e  ','%4.12ei'])]);
            disp(['   nxPrincipal Value interaction:  Iin = ', ...
              num2str([real(Iin), imag(Iin)], ['%4.12e  ','%4.12ei'])]);
            disp('-------------------------------------------------------------------------------');
          end
        end
        % L-operator assembly.
        L(abs(tcn{mg,2}(local_mj)),abs(tcn{ng,2}(local_nj))) = ...
             L(abs(tcn{mg,2}(local_mj)),abs(tcn{ng,2}(local_nj))) ...
             + sign(tcn{mg,2}(local_mj)) .* sign(tcn{ng,2}(local_nj))...
             .* 1i .* k .* (Is - 1 ./ k.^2 .* Id);
        % Ln-operator assembly.
        Ln(abs(tcn{mg,2}(local_mj)),abs(tcn{ng,2}(local_nj))) = ...
             Ln(abs(tcn{mg,2}(local_mj)),abs(tcn{ng,2}(local_nj))) ...
             + sign(tcn{mg,2}(local_mj)) .* sign(tcn{ng,2}(local_nj))...
             .* 1i .* k .* (-Isn + 1 ./ k.^2 .* Idn);
        % K-operator assembly.
        K(abs(tcn{mg,2}(local_mj)),abs(tcn{ng,2}(local_nj))) = ...
             K(abs(tcn{mg,2}(local_mj)),abs(tcn{ng,2}(local_nj))) ...
             + sign(tcn{mg,2}(local_mj)) .* sign(tcn{ng,2}(local_nj))...
             .* (0.5 * Ii + Ic);
        % Kn-operator assembly.
        Kn(abs(tcn{mg,2}(local_mj)),abs(tcn{ng,2}(local_nj))) = ...
             Kn(abs(tcn{mg,2}(local_mj)),abs(tcn{ng,2}(local_nj))) ...
             + sign(tcn{mg,2}(local_mj)) .* sign(tcn{ng,2}(local_nj))...
             .* (-0.5 * Iin - Icn);   
      end
    end
  end
  if (show_progress == 1)
    disp(['   Completion = ', num2str(mg / mesh_size * 100, ["%3.2f%%"])]);
  end
end
% Debug output.
if (enable_debug == 1)
  if (debug_mb ~= 0 && debug_nb ~= 0)
    disp('-------------------------------------------------------------------------------');
    disp(['>> Matrix entry (', ...
      num2str(debug_mb),',', num2str(debug_nb),') for L, Ln, K and Kn operators']);
    disp('-------------------------------------------------------------------------------');
    disp(['   L(', num2str(debug_mb), ',', num2str(debug_nb),') = ', ...
      num2str([real(L(debug_mb, debug_nb)), imag(L(debug_mb, debug_nb))], ...
      ['%4.12e  ','%4.12ei'])]);
    disp(['   Ln(', num2str(debug_mb), ',', num2str(debug_nb),') = ', ...
      num2str([real(Ln(debug_mb, debug_nb)), imag(Ln(debug_mb, debug_nb))], ...
      ['%4.12e  ','%4.12ei'])]);
    disp(['   K(', num2str(debug_mb), ',', num2str(debug_nb),') = ', ...
      num2str([real(K(debug_mb, debug_nb)), imag(K(debug_mb, debug_nb))], ...
      ['%4.12e  ','%4.12ei'])]);
    disp(['   Kn(', num2str(debug_mb), ',', num2str(debug_nb),') = ', ...
      num2str([real(Kn(debug_mb, debug_nb)), imag(Kn(debug_mb, debug_nb))], ...
      ['%4.12e  ','%4.12ei'])]);
    disp('-------------------------------------------------------------------------------');
  end
end
disp('>> Finished');
% Count elapsed time.
toc;