function [ Ii_entry ] = Iimn(Tf, mg, ng, mj, nj, s)
  Ii_entry = 0;
  % Only compute the principal value for the self terms.
  if mg == ng
    Ii_entry = s .* Tf.getSelfTermIi(mj, nj);
  end
end
