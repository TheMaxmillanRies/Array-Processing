kernel = [1,2,1,0,0,0,0,0,0];
input = [1,1,1,1,1,0,0,0,0]; %M + L - 1

conv_out = conv(input, kernel);

fftkernel = fft(kernel);
fftinput = fft(input);

fftconv_out = ifft(fftinput .* fftkernel);

%x = [1 2 3 4 0 0]; y = [-3 5 -4 0 0 0];
x = input; y = kernel;
con_xy1 = conv(x,y);
con_xy2 = ifft(fft(x).*fft(y));

%         fft_seg = zeros(4, window_length + 399);
% 
%         fft_seg(1,:) = [segment(1,:) zeros(1,399)];
%         fft_seg(2,:) = [segment(2,:) zeros(1,399)];
%         fft_seg(3,:) = [segment(3,:) zeros(1,399)];
%         fft_seg(4,:) = [segment(4,:) zeros(1,399)];
% 
%         fft_seg(1,:) = fft(fft_seg(1,:));
%         fft_seg(2,:) = fft(fft_seg(2,:));
%         fft_seg(3,:) = fft(fft_seg(3,:));
%         fft_seg(4,:) = fft(fft_seg(4,:));
%     
%         add processing here
%         a = zeros(4,1423); % a_1(k,l)
%         for j = 1:4
%             a(j,:) = fft([impulse_responses.h_target(j,:) zeros(1,1023)]);
%         end
% 
%         tempa = a';
%         tempb = a' * a;
%         tempc = eye(size(tempb,1)) / tempb;
%         tempd = tempc * a';
%         tempd = tempd.';
%         source = tempd .* fft_seg;
% 
%         fft_seg(1,:) = ifft(source(1,:));
%         fft_seg(2,:) = ifft(source(2,:));
%         fft_seg(3,:) = ifft(source(3,:));
%         fft_seg(4,:) = ifft(source(4,:));