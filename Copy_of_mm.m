clc; clear all
% Load signal files
[clean_1,fs] = audioread("Datasets\clean_speech.wav");
[clean_2,~] = audioread("Datasets\clean_speech_2.wav");
[babble,~] = audioread("Datasets\babble_noise.wav");
[artif_nonstat,~] = audioread("Datasets\aritificial_nonstat_noise.wav");
[speech_shaped,~] = audioread("Datasets\Speech_shaped_noise.wav");

% Load impuse responses
impulse_responses = load("impulse_responses.mat");

speech_size = size(clean_1,1);

% Clip all files to match the target speech file
clean_1 = clean_1';
clean_2 = padarray(clean_2.', [0, speech_size - size(clean_2,1)], 0, 'post');
babble = babble(1:speech_size,1)';
artif_nonstat = artif_nonstat';
speech_shaped = speech_shaped(1:speech_size,1)';

%% Microphone data array
M = 4;  % number of microphones
target = conv(clean_1, impulse_responses.h_target(1,:), 'same');
microphone_data = zeros(M, speech_size);

% Add noise
for i = 1:M
    microphone_data(i,:) = microphone_data(i,:) + conv(clean_2, impulse_responses.h_inter1(i,:), 'same');
    microphone_data(i,:) = microphone_data(i,:) + conv(babble, impulse_responses.h_inter2(i,:), 'same');
    microphone_data(i,:) = microphone_data(i,:) + conv(artif_nonstat, impulse_responses.h_inter3(i,:), 'same');
    microphone_data(i,:) = microphone_data(i,:) + conv(speech_shaped, impulse_responses.h_inter4(i,:), 'same');
end
noise_orig = microphone_data;

% Add clean_1 data
for i = 1:M
    microphone_data(i,:) = microphone_data(i,:) + conv(clean_1, impulse_responses.h_target(i,:), 'same');
end

%% Overlapp-add method
r = 1;
nfft = 512;
window_length = 256;
hop_size = window_length / 2;
window = hamming(window_length)';

das_data = zeros(1, size(microphone_data,2)); % x
mvdr_data = zeros(1, size(microphone_data,2)); % x
mcw_data = zeros(1, size(microphone_data,2)); % x

% STFT
l=1;
window = repmat(window, M, 1);
for i = 1:hop_size:size(microphone_data,2)
%     segment = [];
    if i + window_length*2 > size(microphone_data,2)
        continue;
    else
        segment = zeros(M, window_length);
        seg_noise = zeros(M, window_length);
        segment(:,1:window_length) = microphone_data(:,i:i+window_length - 1);
        seg_noise(:,1:window_length) = noise_orig(:,i:i+window_length - 1);
    end

    segment = segment .* window;
    seg_noise = seg_noise .* window;

    fft_seg = fft(segment',nfft)';
    fft_N = fft(seg_noise',nfft)';

    source_das = zeros(1, window_length);
    source_mvdr = zeros(1, window_length);
    source_mcw = zeros(1, window_length);
    a_est = zeros(M,nfft);

    for k = 1:size(fft_seg,2)
        % Estimate ATF and Rs
        Rx = fft_seg(:,k) * fft_seg(:,k)' ./ M +  eps*eye(M);
        Rn = fft_N(:,k) * fft_N(:,k)' ./ M +  eps*eye(M);
        [V,D] = eig(Rn\Rx);
        [d,ind] = sort(diag(D),'descend');
        Ds = D(ind,ind);
        Vs = V(:,ind);
        Q = Vs^(-1)';
        a_est(:,k) = Q(:,1);
        % Normalize a_est
%         a_est(:,k) = a_est(:,k) ./ a_est(1,k);
        sigma_s = d(1) - 1;

        % Delay & Sum Beamformer
        temp = (a_est(:,k)' * a_est(:,k));
        w_das = a_est(:,k) ./ temp;
        source_das(k) = w_das' * fft_seg(:,k);

        % MVDR Beamfomer
        temp1 = Rx \ a_est(:,k);
        %         temp2 = a_est(:,k)' * Rx^(-1) * a_est(:,k);
        temp2 = a_est(:,k)' * temp1;

        w_mvdr = temp1 ./ temp2;
        source_mvdr(k) = w_mvdr' * fft_seg(:,k);

        % Multi-Channel Wiener filter
        % = Single-Channel Wiener filter + MVDR
        w_scw = sigma_s^2 / (sigma_s^2 + temp2^(-1));
        w_mcw = w_scw * w_mvdr;
        %         w_mcw = sigma_s * temp1;
        source_mcw(k) = w_mcw' * fft_seg(:,k);

    end

    ifft_source_das = (ifft(source_das));
    ifft_source_mvdr = (ifft(source_mvdr));
    ifft_source_mcw = (ifft(source_mcw));

    das_data(i:i+size(ifft_source_das,2) - 1) = das_data(i:i+size(ifft_source_das,2) - 1) + ifft_source_das;
    mvdr_data(i:i+size(ifft_source_mvdr,2) - 1) = mvdr_data(i:i+size(ifft_source_mvdr,2) - 1) + ifft_source_mvdr;
    mcw_data(i:i+size(ifft_source_mcw,2) - 1) = mcw_data(i:i+size(ifft_source_mcw,2) - 1) + ifft_source_mcw;

    %     disp(i);
end

% sound(das_new_data,fs)
das_data = real(das_data);
mvdr_data = real(mvdr_data);
mcw_data = real(mcw_data);
%% Evaluation
stoi_before = stoi(target,microphone_data(1,:),fs);
stoi_mvdr = stoi(target,mvdr_data,fs);
stoi_mcw = stoi(target,mcw_data,fs);

%% Plot
len=500000;
t=(0:len-1)/fs;
figure(1)
subplot(511);
plot(t,target(1:len));title('clean speech');xlabel('t/s');ylabel('Amplitude');
subplot(512);
plot(t,microphone_data(1,1:len));title('noisy speech (mic 1)');xlabel('t/s');ylabel('Amplitude');
subplot(513);
plot(t,das_data(1:len));title('Delay & Sum enhanced speech');xlabel('t/s');ylabel('Amplitude');
subplot(514);
plot(t,mvdr_data(1:len));title('MVDR enhanced speech');xlabel('t/s');ylabel('Amplitude');
subplot(515);
plot(t,mcw_data(1:len));title('Multi-Channel Wiener filter enhanced speech');xlabel('t/s');ylabel('Amplitude');

