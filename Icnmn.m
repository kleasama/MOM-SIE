function [ Icn_entry ] = Icnmn(KER, Tf, Ts, mg, ng, mj, nj, wf, ws, Nf, Ns, doextract)
  % Proper treatment of the self terms.
  if mg == ng
    Icn_entry = 0;
  else
    % Allocate some memory.
    Ic_s = zeros(Nf, 3);
    for nf = 1 : Nf
      if (doextract == 1)
        % Source integral with singularity extraction.
        Ic_s(nf,:) = SubdomainIntegral(ws, ...
                                       Ns, ...
                                       Ts.getArea, ...
                                       cross(KER.D(:,:,nf), Ts.getF(KER.rs, nj), 2), ...
                                       KER.dG(:,:,nf), ...
                                       doextract, ...
                                       KER.dGext(:,:,nf), ...
                                       SingularSumK4(KER, Ts, nj, nf));
      else
        % Source integral without singularity extraction.
        Ic_s(nf,:) = SubdomainIntegral(ws, ...
                                       Ns, ...
                                       Ts.getArea, ...
                                       cross(KER.D(:,:,nf), Ts.getF(KER.rs, nj), 2), ...
                                       KER.dG(:,:,nf), ...
                                       doextract, ...
                                       0, ...
                                       0);
      end
    end
    % Projection integral used in the testing procedure using the singularity
    % extraction technique.
    Icn_entry = ProjectionIntegral(wf, ...
                                   Nf, ...
                                   Tf.getArea, ...
                                   Tf.getNxF(KER.rf, mj), ...
                                   Ic_s);
  end
end
