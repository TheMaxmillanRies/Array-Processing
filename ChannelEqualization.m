
temp = 1/sqrt(2) + 1i*1/sqrt(2);
s = zeros(500, 1);
for i = 1:500
    s(i) = temp;
end

%s(1) = 0 + 1i * 0;

[X, H] = gendata_conv(s, 4, 500, 0.5);
%X2 = gendata_conv2(s, 4, 500, 0);

% Zero Forcing Receiver
wH_ZF = pinv(H);
wH_ZF = wH_ZF(1,:);
ZF_S = wH_ZF * X;

% Wiener Receiver
E = (1/size(X,1))*(X * X');
Rx = sum(E, 'all');
Rxs = H;
w_W = H / Rx;
wH_W = w_W';
W_S = wH_ZF * X;


function [X, H_c] = gendata_conv(s, P, N, sigma)

    H_c = zeros(2*P*N, N);
    for n = 0:N-1
        H = zeros(2*P,1);
        
        for i = 1:P
            H(i) = getH(0 + (i-1)/P);
        end
        H(P+1:end) = H(1:P);
        
        H_c(2*P*n+1:2*P + n*2*P, n+1) = H;
        
        %X(:,n+1) = H * s(n+1);
        
    end
    
    X = H_c * s;
    for i = 1:2*P*N
        X(i) = X(i) + (sigma.*rand(1, 1) + 1i*sigma.*rand(1,1));
    end
    %X = reshape(X, [N, 2*P]);
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




