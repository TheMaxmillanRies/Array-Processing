
temp = 1/sqrt(2) + 1i*1/sqrt(2);
s = [temp, temp, temp, temp];

x = gendata_conv(s, 1, 4, 0.5);

function x = gendata_conv(s, P, N, sigma)
    x = zeros(1, N);
    
    for i = 0:N-1
        for j = 1:size(s,2)
            x(i+1) = x(i+1) + getH(i-j) * s(j) + normrnd(0,sigma) + 1i*normrnd(0,sigma);
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

k + n/p
n = 0 -> p-1