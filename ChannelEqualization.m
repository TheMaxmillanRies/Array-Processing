
temp = 1/sqrt(2) + 1i*1/sqrt(2);
s = zeros(500, 1);
for i = 1:500
    a = rand(1);
    b = rand(1);
    
    if a >= 0.5 && b >= 0.5
        s(i) = 1/sqrt(2) + 1i*1/sqrt(2);
    elseif a >= 0.5 && b < 0.5
        s(i) = 1/sqrt(2) - 1i*1/sqrt(2);
    elseif a < 0.5 && b >= 0.5
        s(i) = -1/sqrt(2) + 1i*1/sqrt(2);
    else
        s(i) = -1/sqrt(2) - 1i*1/sqrt(2);
    end
end

%s(1) = 0 + 1i * 0;

[X, H] = gendata_conv(s, 8, 500, 0);
%X2 = gendata_conv2(s, 4, 500, 0);

% Zero Forcing Receiver
wH_ZF = pinv(H);
wH_ZF = wH_ZF(1,:);
ZF_S = wH_ZF * X;

% Wiener Receiver
E = (1/size(X,1))*(X * X');
Rx = sum(E, 'all');
Rxs = H;
w_W = Rxs / Rx;
wH_W = w_W';
wH_W = wH_W(1,:);
W_S = wH_W * X;

% Plotting
f1 = figure();
plot(ZF_S, '.');

f2 = figure();
plot(W_S, '.');

f3 = figure();
plot(s, '.')


function [X, H] = gendata_conv(s, P, N, sigma)

    X = zeros(2*P, N-1);
    for n = 1:N-1
        H = zeros(P, 2);
        
        for i = 1:P
            H(i,1) = getH(0 + (i-1)/P);
            H(i+P,2) = getH(0 + (i-1)/P);
        end

        X(:,n) = H * s(n:n+1);
        
        for i = 1:2*P
            X(i,n) = X(i,n) + (sigma.*rand(1, 1) + 1i*sigma.*rand(1,1));
        end
    end
end

function x = gendata_conv2(s, P, N, sigma)
    x = zeros(2*P, N);

    for i = 1:2*P
        for j = 1:N
            x(i,j) = sigma.*rand(1, 1) + 1i*sigma.*rand(1,1);
        end
    end
    w = x;

    for i = 1:2*P
        for j = 1:N
            h = zeros(N,1);
            for k  = 1:N
                h(k) = getH((j-1) - (k-1) + (i-1)/P);
            end

            conv_res = conv(s, h, 'valid');
            x(i,j) = x(i,j) + conv_res;
        end
    end
end

function h = getH(t)
    if t >= 0 && t < 0.25
        h = 1;
    elseif t >= 0.25 && t < 0.5
        h = -1;
    elseif t >= 0.5 && t < 0.75
        h = 1;
    elseif t >= 0.75 && t < 1
        h = -1;
    else
        h = 0;
    end
end




