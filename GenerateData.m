clc; clear all
[X, A, S] = gendata(5, 20, 0.5, [-20, 30], [0.1, 0.3], 100000000);

Singular = svd(X);

%plot(Singular,'*r')

theta = esprit(X, 2);

f = espritfreq(X, 2);

%theta, f = joint(X, 2, 10);

%% Comparision
d = 2; M = 3; N = 20; theta = [-20, 30]; f = [0.1, 0.12];
SNR = 0:4:20; test_num = 1000;
theta1 = zeros(d, size(SNR,2),test_num);
f2 = zeros(d, size(SNR,2),test_num);
theta3 = zeros(d, size(SNR,2),test_num);
f3 = zeros(d, size(SNR,2),test_num);

for i=1:size(SNR,2)
    for j=1:test_num
        snr = SNR(i);
        [X, A, S] = gendata(M, N, 0.5, theta, f, snr);

        theta1(:,i,j) = esprit(X, 2);
        f2(:,i,j) = espritfreq(X, 2);
%         [theta3(:,i,j), f3(:,i,j)] = joint(X, 2, 10);
    end

    % Mean
    theta1_mean = mean(theta1,3);
    f2_mean = mean(f2,3);
%     theta3_mean = mean(theta3,3);
%     f3_mean = mean(f3,3);

    % Variance
    theta1_std = std(theta1,0,3);
    f2_std = std(f2,0,3);
%     theta3_std = std(theta3,0,3);
%     f3_std = std(f3,0,3);
end

%% Functions

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

function theta = esprit(X, d)
    X_top = X(1:size(X,1) - 1, :);
    X_bottom = X(2:size(X,1), :);
    
    Z = [X_top; X_bottom];
        
    [U,S,V] = svd(Z,"econ");

    S = S(1:d,1:d);
    U = U(:,1:d);
    
    %Instead of this for loop to cut down small values, I think d needs to
    %be used since we assume the number of sources
    %for i = size(S,1):-1:1
    %    if S(i,i) < 0.00005
    %        U(:,i) = [];
    %        S(i,:) = [];
    %        S(:,i) = [];
    %    end
    %end
    
    Ux = U(1:size(U,1)/2, :);
    Uy = U(size(U,1)/2+1:size(U,1), :);
    
    pUx = pinv(Ux);
    pUxUy = pUx * Uy;

    [Vectors,Values] = eig(pUxUy);
    
    theta = zeros(size(Values,1),1);
    
    for i = 1:size(Values,1)
        theta(i) = 180/pi*asin(angle(Values(i,i))/pi);
    end
end

function f = espritfreq(X, d)
    x_t = X(1,:);

    N = size(x_t, 2);
    m = N/2;

    Z = zeros(m, N-m);

    for j = 1:N-m
        Z(:,j) = x_t(1+(j-1):m+(j-1));
    end

    [U,S,V] = svd(Z,"econ");

    S = S(1:d,1:d);
    U = U(:,1:d);

    Ux = U(1:size(U,1) - 1, :);
    Uy = U(2:size(U,1), :);
    
    pUx = pinv(Ux);
    pUxUy = pUx * Uy;

    [Vectors,Values] = eig(pUxUy);

    f = zeros(size(Values,1),1);

    for i = 1:size(Values,1)
        f(i) = angle(Values(i,i)) / (2*pi);
    end
end

% Joint works but, but I calculate freq and angle separately, normally stop
% at the pseudoinverse multiplication and do the joint solver.
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
    
    S = S(1:d,1:d);
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
    x = 0;
%     [Vectors,Values] = eig(pUxUy);
%     f = zeros(size(Values,1),1);
%     for i = 1:size(Values,1)
%         f(i) = angle(Values(i,i)) / (2*pi);
%     end
end