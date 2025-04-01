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
[Px, F] = pwelch(signal_recu, window, noverlap, Nfft, 1);

Px_db = 10*log10(fftshift(Px));

% sigma = 2; % Écart-type du filtre gaussien
% windowSize = 10000;
% gaussFilter = fspecial('gaussian', [1, windowSize], sigma);
% Px_db_2 = conv(Px_db, gaussFilter, 'same');

% Affichage du periodogramme
figure;
hold on
plot(F, Px_db_2);
plot(F, Px_db);
xlabel('Fréquence (Hz)');
ylabel('Densité spectrale de puissance (dB/Hz)');
title('Periodogramme de Welch du signal reçu');
grid on;
hold off

% Estimation du Temps Symbole
Ts = 10^(-11/10);

% Estimation du rolloff

% Détection de la bande à -40 dB
threshold = - 40; % Seuil à -40 dB
indices = find(Px_db >= threshold); % Indices des fréquences au-dessus du seuil

f_min = F(indices(1));
f_max = F(indices(end));
B_db = f_max - f_min


rolloff =  B_db*Ts-1;

h = rcosdesign(rolloff, Nfft, 1);
H = fftshift(fft(h, Nfft));
plot(H);



