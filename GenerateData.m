close all;

[X, A, S] = gendata(5, 20, 0.5, [20, 30], [0.20, 0.3], 20);

%plot(Singular,'*r')

%theta = esprit(X, 2);

%f = espritfreq(X, 2);

%[theta, f] = joint(X, 2, 10);

%[U,S,V] = svd(X,"econ");
%plot(diag(S), 'o');





esp_theta = zeros(2, 6);
esp_f = zeros(2, 6);
jt_theta = zeros(2, 6);
jt_f = zeros(2, 6);

for snr = 0:4:20
    esprit_theta = zeros(2, 1000);
    esprit_f = zeros(2, 1000);
    joint_theta = zeros(2, 1000);
    joint_f = zeros(2,1000);

    for i = 1:1000
        [X, A, S] = gendata(3, 20, 0.5, [-20, 30], [0.1, 0.12], snr);
        
        theta = esprit(X, 2);
        theta = sort(theta);
        esprit_theta(1,i) = theta(1);
        esprit_theta(2,i) = theta(2);

        f = espritfreq(X, 2);
        f = sort(f);
        esprit_f(1,i) = f(1);
        esprit_f(2,i) = f(2);

        [theta, f] = joint(X, 2, 10);
        theta = sort(theta);
        f = sort(f);
        joint_theta(1,i) = theta(1);
        joint_theta(2,i) = theta(2);
        
        joint_f(1,i) = f(1);
        joint_f(2,i) = f(2);

    end
    

    esp_theta(:,snr/4+1) = mean(esprit_theta, 2); % TODO: fix sign on the values
    esp_f(:,snr/4+1) = mean(esprit_f, 2);
    jt_theta(:,snr/4+1) = mean(joint_theta, 2);
    jt_f(:,snr/4+1) = mean(joint_f, 2);
       
end

% f1 = figure('Name','ESPRIT THETA');
% hold on
% plot((0:4:20), esp_theta(1,:),'r.');
% plot((0:4:20), esp_theta(2,:),'ro');
% % Add joint
% hold off
% 
% f2 = figure('Name','ESPRIT F');
% hold on
% plot((0:4:20), esp_f(1,:),'g.');
% plot((0:4:20), esp_f(2,:),'go');
% % Add joint
% hold off
% 
% f3 = figure('Name','JOINT THETA');
% hold on
% plot((0:4:20), jt_theta(1,:),'g.');
% plot((0:4:20), jt_theta(2,:),'go');
% % Add joint
% hold off
% 
% f4 = figure('Name','JOINT F');
% hold on
% plot((0:4:20), jt_f(1,:),'g.');
% plot((0:4:20), jt_f(2,:),'go');
% % Add joint
% hold off

[X, A, S] = gendata(3, 20, 0.5, [-20, 30], [0.1, 0.12], 10);

theta = esprit(X, 2);
f = espritfreq(X, 2);

[S_theta, w_H_theta] = zero_forcing_theta(X, 0.5, theta);
[S_f, w_H_f] = zero_forcing_freq(X, f);

plot_spatial_response_theta(w_H_theta, 0.5); %TODO: Add magnitude

plot_spatial_response_f(w_H_f, 0.5); %TODO: Add magnitude

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
    
%     [Vectors,Values] = eig(pUxUy_theta);
%     theta = zeros(size(Values,1),1);
%     for i = 1:size(Values,1)
%         theta(i) = -180/pi*asin(angle(Values(i,i))/pi);
%     end
    
    deltaX = transpose(kron([eye(m) zeros(m,1)], eye(size(X,1))));
    deltaY = transpose(kron([zeros(m,1) eye(m)], eye(size(X,1))));
    
    Ux = deltaX * U;
    Uy = deltaY * U;
    
    pUx = pinv(Ux);
    pUxUy_f = pUx * Uy;

    % Solving Joint Diagonalization
    M = [pUxUy_theta pUxUy_f];
    [V,D] =  joint_diag(M,1e-8);
    
    D1 = D(:,1:d);
    D2 = D(:,d+1:2*d);

    theta_tmp=diag(D1);
    theta = -asin(angle(theta_tmp)./(pi))*180/pi;

    phi = diag(D2);
    f=-angle(phi)/(2*pi);


    [theta, index] = sort(theta);

    f = f(index);

%     [Vectors,Values] = eig(pUxUy_f);
%     f = zeros(size(Values,1),1);
%     for i = 1:size(Values,1)
%         f(i) = -angle(Values(i,i)) / (2*pi);
%     end
end

function [S_estimate, w_H] = zero_forcing_theta(X, Delta, theta)
    M = size(X, 1);
    num_sources = size(theta,1);

    A = zeros(M, num_sources);
    for i = 0:M-1
        for j = 1:num_sources
            A(i+1,j) = exp(1i*2*pi*i*Delta*sind(theta(j)));
        end
    end

    tempa = A';
    tempb = A' * A;
    tempc = eye(num_sources)/tempb;

    w_H = tempc * tempa;

    S_estimate = w_H * X;
end

function [S_estimate, w_H] = zero_forcing_freq(X, f)
    num_sources = size(f, 1);
    num_samples = size(X, 2);

    S = zeros(num_sources, num_samples);
    for i = 1:num_sources
        for j = 1:num_samples
            S(i,j) = exp(1i*2*pi*f(i)*j);
        end
    end

    tempa = S';
    tempb = S*tempa;
    tempc = eye(num_sources) / tempb;
    A = (X * tempa) * tempc;
    
    tempa = A';
    tempb = A' * A;
    tempc = eye(num_sources)/tempb;

    w_H = tempc * tempa;
    S_estimate = w_H * X;
end

function spatial_responses = plot_spatial_response_theta(w_H, Delta)
    M = size(w_H, 2);
    num_sources = size(w_H, 1);

    y = zeros(2,180);

    for angle = -90:90
        A = zeros(M, num_sources);
        for i = 0:M-1
            for j = 1:num_sources
                A(i+1,j) = exp(1i*2*pi*i*Delta*sind(angle));
            end
        end
        temp = w_H * A;
        y(1,angle+91) = abs(temp(1,1));
        y(2,angle+91) = abs(temp(2,1));
    end

    x_axis = (-90:90);
    f3 = figure;
    plot(x_axis, y(1,:));

    f4 = figure;
    plot(x_axis, y(2,:));
end

function spatial_responses = plot_spatial_response_f(w_H, Delta)
    M = size(w_H, 2);
    num_sources = size(w_H, 1);

    y = zeros(2,180);

    for angle = -90:90
        A = zeros(M, num_sources);
        for i = 0:M-1
            for j = 1:num_sources
                A(i+1,j) = exp(1i*2*pi*i*Delta*sind(angle));
            end
        end
        temp = w_H * A;
        y(1,angle+91) = abs(temp(1,1));
        y(2,angle+91) = abs(temp(2,1));
    end

    x_axis = (-90:90);
    f7 = figure;
    plot(x_axis, y(1,:));

    f8 = figure;
    plot(x_axis, y(2,:));
end





