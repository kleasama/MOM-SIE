function [ Iin_entry ] = Iinmn(Tf, mg, ng, mj, nj, s)
  Iin_entry = 0;
  % Only compute the principal value for the self terms.
  if mg == ng
    Iin_entry = s .* Tf.getSelfTermIn(mj, nj);
  end
end
