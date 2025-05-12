%% CODI PROJECT PART 1
%---------------------------------------------
% Data: 25/04/2025
% Membres: Antonio Luque, Pau Cornudella i Alex Aleñà
% Assignatura: AERODINÀMICA, MECÀNICA DE VOL I ORBITAL
% Grup: 6
% ---------------------------------------------
clc; clear; close all;


fprintf('Introduiu el número del fitxer per llegir: \n');
disp('1. NACA_22112_N_16.txt');
disp('2. NACA_22112_N_32.txt');
disp('3. NACA_22112_N_64.txt');
disp('4. NACA_22112_N_128.txt');
disp('5. NACA_22112_N_256.txt');
disp('6. NACA_22112_N_512.txt');

file_name = input('Introduïu el número (1, 2, 3, 4, 5, 6): ');

switch file_name
    case 1
        file_name = 'NACA_22112_N_16.txt';
    case 2
        file_name = 'NACA_22112_N_32.txt';
    case 3
        file_name = 'NACA_22112_N_64.txt';
    case 4
        file_name = 'NACA_22112_N_128.txt'; 
    case 5
        file_name = 'NACA_22112_N_256.txt';
    case 6
        file_name = 'NACA_22112_N_512.txt'; 
    otherwise
        disp('Elecció no vàlida. S assigna el fitxer per defecte: NACA_0012_N_16.txt');
        file_name = 'NACA_22112_N_16.txt';
end

% Cridar la funció per llegir les coordenades X i Z
[X, Z] = llegir_punts_perfils(file_name);

% Rotar el perfil
alpha = input('Introduïu l angle datac (en graus): ');
[X, Z] = rotar_perfil(X, Z, alpha);

% Discretitzar el perfil
[X_c, Z_c, N_panells,nx,nz,l] = metode_panells(X, Z);



% Calcular la sustentació (C_L) i la circulació (gamma)
[Cl, Cm14, Cm0, Cp] = calcular_coeficients(gamma_vec, l, Q_inf_m, alpha_rad, N, X, X_c);
