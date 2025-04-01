clear
close all
clc
Te = 10;
level = 1;

%% Chargement du fichier contenant le signal reçu

load 'signal_recu.mat';
signal_recu = signal_recu(20000:200000) -  mean(signal_recu(20000:200000));

% if level > 1
%   signal_recu = signal_recu(1:5:end);
% end

% Paramètres pour Welch
Nfft = 512; 
window = hamming(1024);
noverlap = length(window) / 10;

% Periodogramme de Welch du signal_recu
[Px, F] = pwelch(signal_recu, window, noverlap, Nfft, 1, 'centered');

Px_db = 10*log10(Px);

% sigma = 2; % Écart-type du filtre gaussien
% windowSize = 10000;
% gaussFilter = fspecial('gaussian', [1, windowSize], sigma);
% Px_db_2 = conv(Px_db, gaussFilter, 'same');

% Affichage du periodogramme
figure;
plot(F, Px_db);
xlabel('Fréquence (Hz)');
ylabel('Densité spectrale de puissance (dB/Hz)');
title('Periodogramme de Welch du signal reçu');
grid on;


% Estimation du rolloff

% Détection de la bande à -13 dB
threshold = - 13; % Seuil à -13 dB
indices = find(Px_db >= threshold); % Indices des fréquences au-dessus du seuil

f_min = F(indices(1));
f_max = F(indices(end));
B_Ts = f_max - f_min;

Ts = 1/B_Ts
% Détection de la bande à -40 dB
threshold_2 = - 55; % Seuil à -40 dB
indices_2 = find(Px_db >= threshold_2); % Indices des fréquences au-dessus du seuil
f_min_2 = F(indices_2(1));
f_max_2 = F(indices_2(end));
B_w = f_max_2 - f_min_2

rolloff =  B_w*Ts-1



scatterplot(signal_recu(1:200))

