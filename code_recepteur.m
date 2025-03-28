clear
close all
clc
Te = 10;
level = 1;

%% Chargement du fichier contenant le signal reçu

load 'signal_recu.mat';
if level > 1
  signal_recu = signal_recu(1:5:end);
end

% Paramètres pour Welch
Nfft = 1024; % Taille de la FFT
window = hamming(256); % Fenêtre de Hamming
noverlap = length(window) / 2; % Chevauchement à 50%

% Calcul du periodogramme de Welch
[Px, F] = pwelch(signal_recu, window, noverlap, Nfft, 1/Te); % Te est le temps d'échantillonnage

% Affichage du periodogramme
figure;
plot(F, 10*log10(Px)); % Échelle logarithmique
xlabel('Fréquence (Hz)');
ylabel('Densité spectrale de puissance (dB/Hz)');
title('Periodogramme de Welch du signal reçu');
grid on;

%% Votre récepteur 
% En entrée : signal_recu, signal équivalent à rl(kTe) avec Te le temps
% d'échantillonnage

% hatB doit être une matrice de log2(M) lignes et Ns
% générer la matrice hat(B)
% Définir M (ordre de modulation)
M = 4; % Exemple pour QPSK
log2M = log2(M); % Nombre de bits par symbole

% Définir Ns (nombre de symboles)
Ns = floor(length(signal_recu) / log2M); % Prend la partie entière

% Générer une séquence binaire aléatoire pour hatB (exemple)
hatB = randi([0 1], log2M, Ns);

% calculé grace à la fonction de2bi(foo,2) foo étant ici une représentation entière des étiquettes
%% Décodage de source
hatMatBitImg = reshape(hatB(:),[],8);
matImg = bi2de(hatMatBitImg);
T = 1 % Changer ici la taille de l'image
Img = reshape(matImg,T,T);

%% Affichage
figure
imagesc(Img)
colormap gray
