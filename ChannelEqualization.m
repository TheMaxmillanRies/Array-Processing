
temp = 1/sqrt(2) + 1i*1/sqrt(2);
s = zeros(500, 1);
for i = 1:500
    s(i) = temp;
end

x = gendata_conv(s, 4, 500, 0);

function x = gendata_conv(s, P, N, sigma)
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
                h(j) = getH((j-1) - (k-1) + (i-1)/P);
            end

            conv_res = conv(s, h, 'valid');
            x(i,j) = x(i,j) + conv_res;
        end
    end

    del = 0;
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




