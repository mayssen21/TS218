clear
close all
clc
level = 1;

%Fse = 1;                    % Facteur de suréchantillonnage

%Fe = Ds * Fse;              % Fréquence d’échantillonnage

%% Chargement du fichier contenant le signal reçu

load 'signal_recu.mat'; %655407

%signal_recu = signal_recu(20000:200000) -  mean(signal_recu(20000:200000));

signal_recu = signal_recu(41:655407);

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




% Estimation du rolloff

% Détection de la bande à -13 dB
threshold = - 13; % Seuil à -13 dB
indices = find(Px_db >= threshold); % Indices des fréquences au-dessus du seuil

f_min = F(indices(1));
f_max = F(indices(end));
B_Ts = f_max - f_min;

Fse = round(1/B_Ts)
% Détection de la bande à -40 dB
threshold_2 = - 40; % Seuil à -40 dB
indices_2 = find(Px_db >= threshold_2); % Indices des fréquences au-dessus du seuil
f_min_2 = F(indices_2(1));
f_max_2 = F(indices_2(end));
B_w = f_max_2 - f_min_2

rolloff =  B_w*Fse-1

Fse = 10;
%round(5/0.35)
h = rcosdesign(0.35, 128, Fse);

H = 10*log10(fftshift(fft(h, Nfft)))-15;

% Recepteur

% Filtrage adapté
s_l = conv(signal_recu, conj(flip(h)));

symbols = downsample(s_l, Fse, 8); 

% Programmation de la boucle à verrouillage de phase (Youpi)
symbolsp4 = symbols.^4;
N_symb = length(symbols);
phi_k = 2; 
delta_k = 2;
alpha = 1;
beta = 1;
phi = zeros(1, N_symb);
delta = zeros(1, N_symb);

for k = 2:N_symb
    delta(k) = delta(k-1) + beta * imag(symbolsp4(k) * exp(-1j * phi(k-1)));
    phi(k) = phi(k-1) + delta(k) +alpha * imag(symbolsp4(k) * exp(-1j * phi(k-1)));
end
M=4;
symbols = symbols .* exp(-1j * (phi.'/M + 3*pi/4));



% Décision 
constellation = [exp(1i*pi/4) exp(3*1i*pi/4) exp(5*1i*pi/4) exp(7*1i*pi/4)];
bits_toestimate = [0 0; 0 1; 1 1; 1 0];
bits_recus = [];
iddx = [];

for k = 1:N_symb
    [~, idx] = min(abs(symbols(k) - constellation));
    bits_recus = [bits_recus bits_toestimate(idx,:)'];
    iddx = [iddx idx];
end

iddx


% Affichage du periodogramme
figure;
hold on
plot(F, Px_db);
plot(F, H);
xlabel('Fréquence (Hz)');
ylabel('Densité spectrale de puissance (dB/Hz)');
title('Periodogramme de Welch du signal reçu');
grid on;
hold off

% hatB doit être une matrice de log2(M) lignes et Ns
% calculé grace à la fonction de2bi(foo,2) foo étant ici une représentation entière des étiquettes
%% Décodage de source

hatB = bits_recus(:,1:65536);

hatMatBitImg = reshape(hatB(:),[],8);
matImg = bi2de(hatMatBitImg);
T = 128; 
Img = reshape(matImg,T,T);

%% Affichage
figure;
imagesc(Img);
colormap gray;
title('Image reconstruite');

%% Affichage
figure
imagesc(Img)
colormap gray