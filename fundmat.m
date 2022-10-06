if isOctave
  pkg load symbolic;
end

syms t real;

## P = [2 -1 -1; 1 0 -1; 3 -1 -2];
## P = eye(9);
P = [2 -1; 1 2];

N = length(P);
P = sym(P)
b = sym('b', [N 1]);

% Находим собственные числа и их алгебраическую кратность
eig_vec = eig(P);
eigenvalues = unique(eig_vec).';
syms count;

A = sym('A', size(P));
e = sym('e', size(b));

ignored_eigenvalues = [];

index = 1;
for i = 1:numel(eigenvalues)
  eigenvalue = eigenvalues(i);
  if length(find(ignored_eigenvalues==eigenvalue)) ~= 0
    fprintf(" -> ignored:\n")
    eigenvalue
    continue
  end

  count = length(find(eig_vec==eigenvalue));

  eigenvalue
  count

  if imag(eigenvalue) == 0
    for row = 1:count
      e(index) = t^(row - 1) * exp(eigenvalue * t);
      for col = 1:N
        A(index, col) = getCoeff(row, col, eigenvalue);
      end

      index += 1;
    end
  else
    syms z;
    lhs = exp(z * t);

    for row = 1:count
      for col = 1:N
        rhs(col) = getCoeff(row, col, eigenvalue);
      end

      lhs_sub = subs(lhs, z, eigenvalue);
      e(index) = real(lhs_sub);
      A(index, 1:N) = real(rhs);
      index += 1;

      e(index) = imag(lhs_sub);
      A(index, 1:N) = imag(rhs);
      index += 1;

      lhs = diff(lhs, z);
    end

    fprintf("ignored_eigenvalues updated:\n");
    ignored_eigenvalues = [ignored_eigenvalues; conj(eigenvalue)]
  end

end

b = simplify(A \ e)

Y = zeros(size(P));
Y = sym(Y);

for k = 0:N-1
  Y += b(k+1) * P^k;
end

simplify(Y)
