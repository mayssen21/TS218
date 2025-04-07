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


%% === Synchronisation par boucle de Costas ===
% Paramètres de la boucle
alpha = 0.01;    % Coefficient proportionnel
beta = 0.0001;   % Coefficient intégral

theta_hat = 0;         % Phase estimée initiale
v_n = 0;               % Sortie du filtre (correction)
N = length(signal_recu);
signal_corrige = zeros(1, N);  % Signal corrigé en phase

for n = 2:N
    % Correction du signal avec estimation de phase
    signal_corrige(n) = signal_recu(n) * exp(-1j * theta_hat);
    
    % Détecteur de phase :
    e_n = imag(signal_corrige(n)^2);  % Erreur de phase
    
    % Filtre de boucle : intégration avec alpha et beta
    v_n = v_n + beta * e_n;
    theta_hat = theta_hat + alpha * e_n + v_n;  % mise à jour de la phase
end


% Periodogramme de Welch du signal_recu
[Px, F] = pwelch(signal_corrige, window, noverlap, Nfft, 1);

Px_db = 10*log10(fftshift(Px));

% sigma = 2; % Écart-type du filtre gaussien
% windowSize = 10000;
% gaussFilter = fspecial('gaussian', [1, windowSize], sigma);
% Px_db_2 = conv(Px_db, gaussFilter, 'same');

% Décision 
constellation = [exp(1i*pi/4) exp(3*1i*pi/4) exp(5*1i*pi/4) exp(7*1i*pi/4)];
bits_toestimate = [0 0; 0 1; 1 1; 1 0];
bits_recus = [];
for k = 1:N_symb
    [~, idx] = min(abs(symbols(k) - constellation));
    bits_recus = [bits_recus bits_toestimate(idx,:)];
end


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
threshold = - 13; 
indices = find(Px_db >= threshold); % Indices des fréquences au-dessus du seuil

f_min = F(indices(1));
f_max = F(indices(end));
B_db = f_max - f_min


rolloff =  B_db*Ts-1;

h = rcosdesign(rolloff, Nfft, 1);
H = fftshift(fft(h, Nfft));
plot(H);





%% Votre récepteur 
% En entrée : signal_recu, signal équivalent à rl(kTe) avec Te le temps
% d'échantillonnage

% hatB doit être une matrice de log2(M) lignes et Ns
% calculé grace à la fonction de2bi(foo,2) foo étant ici une représentation entière des étiquettes
Ns = length(bits_recus)/2;
hatB = reshape(bits_recus, 2, Ns);  

%% Décodage de source
hatMatBitImg = reshape(hatB(:),[],8);
matImg = bi2de(hatMatBitImg);
T = 1 % Changer ici la taille de l'image
Img = reshape(matImg,T,T);

%% Affichage
figure
imagesc(Img)
colormap gray