if isOctave
  pkg load symbolic;
end

SHOULD_TEST = 1;
DEBUG_PRINTING = 0;

if DEBUG_PRINTING
  debugPrint = @(varargin) printf(varargin{:});
else
  debugPrint = @(varargin) 1;
end

syms t real;

## P = [2 -1 -1; 1 0 -1; 3 -1 -2];
## P = eye(9);
P = [2 -1; 1 2];

N = length(P);
P = sym(P);
b = sym('b', [N 1]);

debugPrint(" -> P:\n%s\n", pretty(P));

% Находим собственные числа и их алгебраическую кратность
eig_vec = eig(P);
eigenvalues = unique(eig_vec).';
debugPrint("\n -> eigenvalues:\n%s\n", pretty(eigenvalues));

syms count;

A = sym('A', size(P));
e = sym('e', size(b));

ignored_eigenvalues = [];

index = 1;
for i = 1:numel(eigenvalues)
  eigenvalue = eigenvalues(i);
  debugPrint("\n -> processing eigenvalue: %s", pretty(eigenvalue));
  if length(find(ignored_eigenvalues==eigenvalue)) ~= 0
    debugPrint("   -> ignoring\n");
    continue
  end

  count = length(find(eig_vec==eigenvalue));
  debugPrint("   -> algebraic multiplicity: %d\n", count);

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
    rhs = sym('rhs', [1 N]);

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

    ignored_eigenvalues = [ignored_eigenvalues; conj(eigenvalue)];
    debugPrint("   -> ignored_eigenvalues updated:\n%s\n", pretty(ignored_eigenvalues));
  end

end

debugPrint("\n -> calculating b...");
b = A \ e;
debugPrint("\n -> b:\n%s\n", pretty(b));
debugPrint(" -> simplifying b...");
b = simplify(b);
printf("\nb:\n%s\n", pretty(b));

Y = zeros(size(P));
Y = sym(Y);

debugPrint("\n -> calculating Y...");
for k = 0:N-1
  Y += b(k+1) * P^k;
end
debugPrint("\n -> Y:\n%s\n", pretty(Y));

debugPrint(" -> simplifying Y...");
Y = simplify(Y);
printf("\n\nY:\n%s\n", pretty(Y));

if SHOULD_TEST
  printf("\n--- TEST ---\n");
  printf(" 1. Y' = PY\n")
  printf("    Y' = PY: %s\n", mat2str(areIdentical(diff(Y), P * Y)))
  printf(" 1. Y(0) = E\n")
  printf("    Y(0) = E: %s\n", mat2str(areIdentical(subs(Y, 0), eye(N))))
end
