% Load signal files
[clean_1,Fs_clean_1] = audioread("clean_speech.wav");
[clean_2,Fs_clean_2] = audioread("clean_speech_2.wav");
[babble,Fs_babble] = audioread("babble_noise.wav");
[artif_nonstat,Fs_artif_nonstat] = audioread("aritificial_nonstat_noise.wav");
[speech_shaped,Fs_speech_shaped] = audioread("Speech_shaped_noise.wav");

% Load impuse responses
impulse_responses = load("impulse_responses.mat");

% Find largest file
padded_size = max([size(clean_1,1), ...
                   size(clean_2,1), ...
                   size(babble,1), ...
                   size(artif_nonstat,1), ...
                   size(speech_shaped,1)]);

% Pad all files to match the largest
clean_1 = padarray(clean_1.', [0, padded_size - size(clean_1,1)], 0, 'post');
clean_2 = padarray(clean_2.', [0, padded_size - size(clean_2,1)], 0, 'post');
babble = padarray(babble.', [0, padded_size - size(babble,1)], 0, 'post');
artif_nonstat = padarray(artif_nonstat.', [0, padded_size - size(artif_nonstat,1)], 0, 'post');
speech_shaped = padarray(speech_shaped.', [0, padded_size - size(speech_shaped,1)], 0, 'post');

%% 

% Microphone data array
microphone_data = zeros(4, padded_size-399);

% Add clean_1 data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(clean_1, impulse_responses.h_target(i,:), 'valid');
end

% Add clean_2 data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(clean_2, impulse_responses.h_inter1(i,:), 'valid');
end

% Add babble data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(babble, impulse_responses.h_inter2(i,:), 'valid');
end

% Add artif_nonstat data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(artif_nonstat, impulse_responses.h_inter3(i,:), 'valid');
end

% Add artif_nonstat data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(speech_shaped, impulse_responses.h_inter4(i,:), 'valid');
end

%% 

% Overlapp-add method
window_length = 400;
hop_size = window_length / 2;
window = hamming(window_length).';
new_data = [];

% FFT of IR
a = zeros(4,800);
for j = 1:4
    a(j,:) = fft(impulse_responses.h_target(j,:), 800);
end

norm = a(1,:);
% Make relative IR
for j = 1:4
    a(j,:) = a(j,:) / norm;
end

% Create windows and process in STFT space
test1 = zeros(4, size(microphone_data,2));
test2 = zeros(4, size(microphone_data,2));

new_data = zeros(4, size(microphone_data,2)); % x
for i = 1:hop_size:size(microphone_data,2)
        segment = [];
        % Final window is incomplete due to size % window_length != 0 but
        % doesn't matter since it's mostly 0's due to padding being used
        if i + window_length > size(microphone_data,2)
%             segment = zeros(4, size(microphone_data,2) - i);
%             segment(1,1:window_length) = microphone_data(1,i:size(microphone_data,2) - 1);
%             segment(2,1:window_length) = microphone_data(2,i:size(microphone_data,2) - 1);
%             segment(3,1:window_length) = microphone_data(3,i:size(microphone_data,2) - 1);
%             segment(4,1:window_length) = microphone_data(4,i:size(microphone_data,2) - 1);
% 
%             window = window(1:size(segment,2));
              continue;
        else
            segment = zeros(4, window_length);
            segment(1,1:window_length) = microphone_data(1,i:i+window_length - 1);
            segment(2,1:window_length) = microphone_data(2,i:i+window_length - 1);
            segment(3,1:window_length) = microphone_data(3,i:i+window_length - 1);
            segment(4,1:window_length) = microphone_data(4,i:i+window_length - 1);
        end

        segment(1,:) = segment(1,:) .* window;
        segment(2,:) = segment(2,:) .* window;
        segment(3,:) = segment(3,:) .* window;
        segment(4,:) = segment(4,:) .* window;

%         fft_seg = zeros(4, 800);
%         for j = 1:4
%             fft_seg(j,:) = fft(segment(j,:), 800);
%         end 
% 
%         a_p = pinv(a);
%         fft_s_estimate = a_p.' .* fft_seg;
%         ifft_s_estimate = ifft(fft_s_estimate);
%         ifft_s_estimate = ifft_s_estimate(:,1:400);
% 
%         w = (1/4)*a.';
%         fft_s_estimate2 = w.' .* fft_seg;
%         ifft_s_estimate2 = ifft(fft_s_estimate2);
%         ifft_s_estimate2 = ifft_s_estimate2(:,1:400);
%         
%         test1(:, i:i+size(ifft_s_estimate,2) - 1) = test1(:, i:i+size(segment,2) - 1) + ifft_s_estimate;
%         test2(:, i:i+size(ifft_s_estimate2,2) - 1) = test2(:, i:i+size(segment,2) - 1) + ifft_s_estimate2;


        new_data(1, i:i+size(segment,2) - 1) = new_data(1, i:i+size(segment,2) - 1) + segment(1,:);
        new_data(2, i:i+size(segment,2) - 1) = new_data(2, i:i+size(segment,2) - 1) + segment(2,:);
        new_data(3, i:i+size(segment,2) - 1) = new_data(3, i:i+size(segment,2) - 1) + segment(3,:);
        new_data(4, i:i+size(segment,2) - 1) = new_data(4, i:i+size(segment,2) - 1) + segment(4,:);

        disp(i);
end

f1 = figure;
plot(test1(1, 1:6000));

f2 = figure;
plot(test2(1, 1:6000));

