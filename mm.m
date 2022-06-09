clc; clear all
% Load signal files
[clean_1,fs] = audioread("Datasets\clean_speech.wav");
[clean_2,~] = audioread("Datasets\clean_speech_2.wav");
[babble,~] = audioread("Datasets\babble_noise.wav");
[artif_nonstat,~] = audioread("Datasets\aritificial_nonstat_noise.wav");
[speech_shaped,~] = audioread("Datasets\Speech_shaped_noise.wav");

% Load impuse responses
impulse_responses = load("impulse_responses.mat");

% Find largest file
% padded_size = max([size(clean_1,1), ...
%                    size(clean_2,1), ...
%                    size(babble,1), ...
%                    size(artif_nonstat,1), ...
%                    size(speech_shaped,1)]);
speech_size = size(clean_1,1);

% Clip all files to match the target speech file
clean_2 = padarray(clean_2.', [0, speech_size - size(clean_2,1)], 0, 'post')';
babble = babble(1:speech_size,1);
speech_shaped = speech_shaped(1:speech_size,1);

%% 

% Microphone data array
microphone_data = zeros(speech_size,4);

% Add clean_1 data
for i = 1:size(microphone_data,2)
    microphone_data(:,i) = microphone_data(:,i) + conv(clean_1, impulse_responses.h_target(i,:), 'same');
end

% Add clean_2 data
for i = 1:size(microphone_data,2)
    microphone_data(:,i) = microphone_data(:,i) + conv(clean_2, impulse_responses.h_inter1(i,:), 'same');
end

% Add babble data
for i = 1:size(microphone_data,2)
    microphone_data(:,i) = microphone_data(:,i) + conv(babble, impulse_responses.h_inter2(i,:), 'same');
end

% Add artif_nonstat data
for i = 1:size(microphone_data,2)
    microphone_data(:,i) = microphone_data(:,i) + conv(artif_nonstat, impulse_responses.h_inter3(i,:), 'same');
end

% Add artif_nonstat data
for i = 1:size(microphone_data,2)
    microphone_data(:,i) = microphone_data(:,i) + conv(speech_shaped, impulse_responses.h_inter4(i,:), 'same');
end

%% Beamforming

% Delay & Sum beamformer
nfft=256;
win = hamming(nfft,'periodic');
S = stft(clean_1, fs, 'Window', win,'OverlapLength',0.5*nfft);
X = stft(microphone_data, fs, 'Window', win,'OverlapLength',0.5*nfft);

%% Estimate ATF
noise_only = microphone_data(1:8000,:); % only one reference microphone or all microphones?
N = stft(noise_only, fs, 'Window', win,'OverlapLength',0.5*nfft);
t_repeat = ceil(size(X,2)/size(N,2));
N_rm = repmat(N,1,t_repeat,1);

for i = 1:size(X,1)
    for j = 1:size(X,2)
        Rn = squeeze(N_rm(i,j,:)) * squeeze(N_rm(i,j,:))';
%         [Vectors,Values] = eig(Rn);
        pre_whiten = Rn^(-1/2);
        I = pre_whiten * Rn * pre_whiten;
        X_pw = pre_whiten * squeeze(X(i,j,:));
        Rx_pw = X_pw * X_pw';
        [V,D] = eig(Rx_pw);
        [d,ind] = sort(diag(D),'descend');
        Ds = D(ind,ind);
        Vs = V(:,ind);
        a(i,j,:) = Vs(:,1);
    end
end

% pre_whiten = sqrtm(Rn);
noise_whiten = pre_whiten * noise_only;
Rn_w = noise_whiten * noise_whiten';

% for i = 1:size(X,1)
%     for j = 1:size(X,2)
%         Rx = squeeze(X(i,j,:)) * squeeze(X(i,j,:))';
%         [U,V,Q] = eig(Rx);
%         if j == 2000
%             print(j);
%         end
%        
%     end
% end