function [A,LL,Cost] = acdc(MM)
% function A = acdc(MM)
% MM = M(:,:,k) sequence of square matrices to be jointly diagonalized
% LL = L(:,:,k) sequence of diagonal factors
% M_k = A * L_k * A'
%
% see ICA-2000 Arie Yeredor

[n,n,K] = size(MM);
maxloop = 150;

% init
M1 = squeeze(MM(:,:,1));
M2 = squeeze(MM(:,:,2));
[A,junk] = eig(M1*pinv(M2));

LL = zeros(n,n,K);
for k = 1:K,
   LL(:,:,k) = diag(diag( inv(A) * squeeze(MM(:,:,k)) * inv(A') ));
end


Cost = cost(A,MM);

for loop=1:maxloop,
    % AC
    for l = 1:n,
       P = zeros(n,n);
       for k = 1:K
	   S = zeros(n,n);
	   for j = 1:n,
	       if (j ~= l),
		   S = S + LL(j,j,k) * A(:,j)*A(:,j)' ; 
	       end
	   end
	   P = P + LL(l,l,k)*( squeeze(MM(:,:,k)) -  S) ;
       end
       [V,Lambda] = eig(P);
       [maxi,ii] = max(diag(Lambda));
       if (maxi < 0),
	   A(:,l) = zeros(n,1);
       else
	   A(:,l) = V(:,ii) * sqrt(maxi) / sqrt(sum( (squeeze(LL(l,l,:)) ).^2) );
       end
    end

    % DC

    G = pinv( conj(A'*A) .* (A'*A) );

    for k = 1:K,
	LL(:,:,k) = diag( G* diag(A' * squeeze(MM(:,:,k)) * A));
    end

    % collect cost (model error)
    Cost = [Cost cost(A,MM)];
end
