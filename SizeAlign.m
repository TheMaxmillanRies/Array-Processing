function [y] = SizeAlign(y, len)

% len is target length
% y is the signal needs to be clapped or padded

if len > length(y)
    % Zero-Padding
    y = padarray(y.', [0, len - size(y,1)], 0, 'post');
else
    % Clip
    y = y(1:len)';
end

