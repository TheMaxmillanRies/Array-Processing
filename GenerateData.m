[X, A, S] = gendata(5, 20, 0.5, [-20, 30], [0.1, 0.3], 100000000);

Singular = svd(X);

%plot(Singular,'*r')

theta = esprit_direction(X, 2);

function [X, A, S] = gendata(M, N, Delta, theta, f, SNR)
    % Create empty matrix MxN -> ReceiverAntenna x SamplesMeasured
    X = zeros(M, N);

    % S = source vector
    % Create and initialize vector of sources (irrespective of receiver
    % antenna)
    num_sources = size(f, 2);
    S = zeros(num_sources, 1);
    for i = 1:num_sources
        S(i) = exp(1i*2*pi*f(i));
    end

    % NOTE: x(t) = A*s(t) + n(t)
    % A = attenuation caused by angle and antenna distance (delta)
    % k = M?
    A = zeros(M, num_sources);
    for i = 0:M-1
        for j = 1:num_sources
            A(i+1,j) = exp(1i*2*pi*i*Delta*sind(theta(j)));
        end
    end

    % Add Noise - AWGN
    for i = 1:N
        S_prime = S.^i;
        temp = A * S_prime;
        X(:,i) = awgn(temp, SNR);
    end
end

function theta = esprit_direction(X, d)
    X_top = X(1:size(X,1) - 1, :);
    X_bottom = X(2:size(X,1), :);
    
    Z = [X_top; X_bottom];
        
    [U,S,V] = svd(Z,"econ");
    
    %Instead of this for loop to cut down small values, I think d needs to
    %be used since we assume the number of sources
    for i = size(S,1):-1:1
        if S(i,i) < 0.00005
            U(:,i) = [];
            S(i,:) = [];
            S(:,i) = [];
        end
    end
    
    Ux = U(1:size(U,1)/2, :);
    Uy = U(size(U,1)/2+1:size(U,1), :);
    
    pUx = pinv(Ux);
    pUxUy = pUx * Uy;

    [Vectors,Values] = eig(pUxUy);
    
    theta = zeros(size(Values,1),1);
    
    for i = 1:size(Values,1)
        theta(i) = 180/pi*asin(acos(real(Values(i,i)))/pi);
    end
end