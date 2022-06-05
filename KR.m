function U = KR(A,B)
% usage: U = KR(A,B)
% Khatri-Rao product. A and B have to have the same number of columns

[a,c] = size(A);
[b,c1] = size(B);

if (c ~= c1), error('KR: A and B do not have same number of columns'); end

U = zeros(a*b,c);

for l=1:c
  U(:,l) = kron(A(:,l),B(:,l));
end

