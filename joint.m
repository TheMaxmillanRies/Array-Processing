function [theta, f] = joint(X, d, m)
    N = size(X, 2);
    
    Z = zeros(m*size(X, 1), N-m);
    
    for j = 1:N-m
        counter = 0;
        for i = 1:size(X, 1):m*size(X, 1)
            Z(i:i+size(X, 1)-1,j) = X(:,j+counter);
            counter = counter + 1;
        end
    end
   
    [U,S,V] = svd(Z,"econ");
    
    U = U(:,1:d);
    
    deltaX = [eye(size(X,1)) zeros(size(X,1),1)];
    deltaY = [zeros(size(X,1),1) eye(size(X,1))];
    
    tempX = transpose(kron(eye(m), deltaX));
    tempY = transpose(kron(eye(m), deltaY));
    
    Ux = tempX * U;
    Uy = tempY * U;
    
    pUx = pinv(Ux);
    pUxUy_theta = pUx * Uy;
    
    %[Vectors,Values] = eig(pUxUy);
    %theta = zeros(size(Values,1),1);
    %for i = 1:size(Values,1)
    %    theta(i) = 180/pi*asin(angle(Values(i,i))/pi);
    %end
    
    deltaX = transpose(kron([eye(m) zeros(m,1)], eye(size(X,1))));
    deltaY = transpose(kron([zeros(m,1) eye(m)], eye(size(X,1))));
    
    Ux = deltaX * U;
    Uy = deltaY * U;
    
    pUx = pinv(Ux);
    pUxUy_f = pUx * Uy;

    MM = zeros(2,2,2);
    MM(:,:,1) = pUxUy_theta;
    MM(:,:,2) = pUxUy_f;
    [A,LL,Cost] = acdc(MM);

%     [Vectors,Values] = eig(pUxUy);
%     f = zeros(size(Values,1),1);
%     for i = 1:size(Values,1)
%         f(i) = angle(Values(i,i)) / (2*pi);
%     end
end