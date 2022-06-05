function Cost = cost(A,MM)

[n,n,K] = size(MM);

% show modeling error.  First compute optimal Lambda
Av = KR(conj(A),A);

Y = [];
for k=1:K,
   Y = [Y vec(squeeze(MM(:,:,k)))];
end

L = pinv(Av)*Y;		% each col is optimal Lambda_k

Err = [];
for k=1:K,
   Yk = squeeze(MM(:,:,k));
   err = Yk - A* diag(L(:,k)) * A';
   Err = [Err; err];
end
Cost = norm(Err,'fro');

% same as norm(Y - Av*L,'fro')
