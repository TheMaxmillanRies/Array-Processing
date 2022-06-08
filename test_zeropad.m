clc; clear all
% Load signal files
[clean_1,fs] = audioread("Datasets\clean_speech.wav");

% Load impuse responses
impulse_responses = load("impulse_responses.mat");

x = conv(clean_1, impulse_responses.h_target(1,:), 'same');
pad_length = size(clean_1,1);
h_zeropad = padarray(impulse_responses.h_target(1,:),[0 pad_length-400],0,'post');
x_pad = conv(clean_1, h_zeropad, 'full')';
x_new = x_pad(1:pad_length);

% t1=(0:length(x)-1)/fs;
% figure;plot(t1,x);xlabel('t/s');ylabel('Amplitude');
% t2=(0:length(x_pad)-1)/fs;
% figure;plot(t1,x_new);xlabel('t/s');ylabel('Amplitude');

nfft=128;
win = hamming(nfft,'periodic');
[S,f,t] = stft(clean_1, fs, 'Window', win,'OverlapLength',0.5*nfft);
X = stft(x, fs, 'Window', win,'OverlapLength',0.5*nfft);
H = stft(h_zeropad', fs, 'Window', win,'OverlapLength',0.5*nfft);

% for i = 1:size(microphone_data,2)
%     A(:,:,i) = X(:,:,i) ./ S;
% end
% A(isnan(A)) = 0;
% 
% figure;
% imagesc(t, f, 20*log10((abs(A(1024:end,:,1)))));xlabel('Samples'); ylabel('Freqency');
% colorbar;