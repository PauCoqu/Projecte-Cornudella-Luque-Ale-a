%% CODI PROJECT PART 1
%---------------------------------------------
% Data: 25/04/2025
% Membres: Antonio Luque, Pau Cornudella i Alex Aleñà
% Assignatura: AERODINÀMICA, MECÀNICA DE VOL I ORBITAL
% Grup: 6
% ---------------------------------------------
clc; clear; close all;

% Nom del fitxer a llegir:
file_name = 'NACA_22112_N_16.txt';

% Rotar el perfil
alpha = 8; 
alpha_rad = deg2rad(alpha);

% Cridar la funció per llegir les coordenades X i Z
[X, Z] = llegir_punts_perfil(file_name); %Funció que llegeix els punts del perfil i el grafica



% [X, Z] = rotar_perfil(X, Z, alpha);
% 
% % Discretitzar el perfil
% [X_c, Z_c, N_panells,nx,nz,l] = discretitzar_perfil(X, Z);
% 
% % Calcular la sustentació (C_L) i la circulació (gamma)
% [gamma, CL] = calcular_gamma_CL(X, Z, X_c, Z_c, nx, nz, N_panells, alpha, l);

