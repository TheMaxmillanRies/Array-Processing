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

%%
target = conv(clean_1, impulse_responses.h_target(1,:), 'same');
% Microphone data array
microphone_data = zeros(4, speech_size);

% Add clean_2 data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(clean_2, impulse_responses.h_inter1(i,:), 'same');
end

% Add babble data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(babble, impulse_responses.h_inter2(i,:), 'same');
end

% Add artif_nonstat data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(artif_nonstat, impulse_responses.h_inter3(i,:), 'same');
end

% Add artif_nonstat data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(speech_shaped, impulse_responses.h_inter4(i,:), 'same');
end

noise_orig = microphone_data;

% Add clean_1 data
for i = 1:size(microphone_data,1)
    microphone_data(i,:) = microphone_data(i,:) + conv(clean_1, impulse_responses.h_target(i,:), 'same');
end

%% Beamforming
% Overlapp-add method
r = 1;
M = size(microphone_data,1);
nfft = 800;
window_length = 400;
hop_size = window_length / 2;
window = hamming(window_length).';

% FFT of IR
a = zeros(M,nfft);
for j = 1:M
    a(j,:) = fft(impulse_responses.h_target(j,:), 800);
end

norm = a(1,:);
% Make relative IR
for j = 1:M
    a(j,:) = a(j,:) ./ norm;
end

das_data = zeros(1, size(microphone_data,2)); % x
mvdr_data = zeros(1, size(microphone_data,2)); % x
mcw_data = zeros(1, size(microphone_data,2)); % x
Rn = zeros(M,M,nfft);

% Estimate Rn
noise_only_len = 0.5 * fs;
for i = 1:hop_size:noise_only_len+1
        segment = [];
        if i + window_length*2 > size(microphone_data,2)
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

        fft_seg = zeros(M, nfft);
        for j = 1:M
            fft_seg(j,:) = fft(segment(j,:), nfft);
        end 

        for k = 1:size(fft_seg,2)
            Rn(:,:,k) = Rn(:,:,k) + fft_seg(:,k) * fft_seg(:,k)' ./ M;
        end
end
Rn = Rn / (noise_only_len/hop_size);

for i = 1:hop_size:size(microphone_data,2)
    segment = [];
    if i + window_length*2 > size(microphone_data,2)
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

    fft_seg = zeros(M, nfft);
    for j = 1:M
        fft_seg(j,:) = fft(segment(j,:), nfft);
    end

    % Estimate Rn
%     if i <= 8001 % first 0.5s as noise-only segment
%         Rn = zeros(M,M,nfft);
%         for k = 1:size(fft_seg,2)
%             Rn(:,:,k) = fft_seg(:,k) * fft_seg(:,k)' ./ M + eps*eye(M);
%         end
%     end

    source_das = zeros(1, window_length);
    source_mvdr = zeros(1, window_length);
    source_mcw = zeros(1, window_length);
    a_est = zeros(M,nfft);

    for k = 1:size(fft_seg,2)
        % Estimate ATF and Rs
        Rx = fft_seg(:,k) * fft_seg(:,k)' ./ M +  eps*eye(M);
        %         [V,D] = eig(Rx,Rn(:,:,k));
        [V,D] = eig(Rn(:,:,k)\Rx);
        %         if find(isinf(D)) ~= 0
        %             Rn(:,:,k) = Rn(:,:,k) + 1.00e-10;
        %             [V,D] = eig(Rx,Rn(:,:,k));
        %         end
        [d,ind] = sort(diag(D),'descend');
        Ds = D(ind,ind);
        Vs = V(:,ind);
        Q = Vs^(-1)';
        a_est(:,k) = Q(:,1);
        D1 = Ds(1:r,1:r);
        sigma_s = D1 - 1;

        % Delay & Sum Beamformer
        temp = (a_est(:,k)' * a_est(:,k));
        w_das = a_est(:,k) ./ temp;
        source_das(k) = w_das' * fft_seg(:,k);

        % MVDR Beamfomer
        temp1 = Rx \ a_est(:,k);
%         temp2 = a_est(:,k)' * Rx^(-1) * a_est(:,k);
        temp2 = a_est(:,k)' * temp1;

        %             temp1 = Rn(:,:,k) \ aa(:,k);
        %             temp2 = aa(:,k)' * Rn(:,:,k)^(-1) * aa(:,k);

        w_mvdr = temp1 ./ temp2;
        source_mvdr(k) = w_mvdr' * fft_seg(:,k);

        % Multi-Channel Wiener filter
        % = Single-Channel Wiener filter + MVDR
        w_scw = sigma_s^2 / (sigma_s^2 + temp2^(-1));
        w_mcw = w_scw * w_mvdr;
%         w_mcw = sigma_s * temp1;
        source_mcw(k) = w_mcw' * fft_seg(:,k);

    end

    ifft_source_das = real(ifft(source_das));
    ifft_source_mvdr = real(ifft(source_mvdr));
    ifft_source_mcw = real(ifft(source_mcw));

    das_data(i:i+size(ifft_source_das,2) - 1) = das_data(i:i+size(ifft_source_das,2) - 1) + ifft_source_das;
    mvdr_data(i:i+size(ifft_source_mvdr,2) - 1) = mvdr_data(i:i+size(ifft_source_mvdr,2) - 1) + ifft_source_mvdr;
    mcw_data(i:i+size(ifft_source_mcw,2) - 1) = mcw_data(i:i+size(ifft_source_mcw,2) - 1) + ifft_source_mcw;

    %     disp(i);
end

% sound(das_new_data,fs)

%% Evaluation


%% Plot
t=(0:100000-1)/fs;
figure(1)
subplot(511);
plot(t,target(1:100000));ylim([-0.04,0.04]);title('clean speech');xlabel('t/s');ylabel('Amplitude');
subplot(512);
plot(t,microphone_data(1,1:100000));ylim([-0.04,0.04]);title('noisy speech (mic 1)');xlabel('t/s');ylabel('Amplitude');
subplot(513);
plot(t,das_data(1:100000));ylim([-0.04,0.04]);title('Delay & Sum enhanced speech');xlabel('t/s');ylabel('Amplitude');
subplot(514);
plot(t,mvdr_data(1:100000));title('MVDR enhanced speech');xlabel('t/s');ylabel('Amplitude');
subplot(515);
plot(t,mcw_data(1:100000));title('Multi-Channel Wiener filter enhanced speech');xlabel('t/s');ylabel('Amplitude');
