%% CODI PROJECTE PART 1
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
        disp('Elecció no vàlida. S''assigna el fitxer per defecte: NACA_22112_N_16.txt');
        file_name = 'NACA_22112_N_16.txt';
end

% Llegir les coordenades X i Z
[X, Z] = llegir_punts_perfils(file_name);
alpha = input('Introduïu l''angle datac (en graus): ');
alpha_rad = deg2rad(alpha);

% Discretitzar el perfil. Calcul gamma.
[X_c, Z_c, l, Tangent, Normal, gamma, Q_inf, Q_inf_modul, N] = metode_panells(X, Z, alpha_rad);

% Calcular la sustentació (C_l) i Cm_1/4
[Cl, Cm14, Cm0, Cp] = calcular_CL_CM14(gamma, l, Q_inf_modul, X, Z, X_c, Z_c);

% Rotar el perfil i gràficar Cp
rotar_perfil(X, Z,alpha_rad, Normal,Cp)


%2. Mach crític per a alpha=[0, 2, 4]º%%
gamma_aire=1.4;
%[Cp_0, Cp_kt, Cp_star, Cp_crit, M_crit] = M_critic(Cp, N, gamma_aire);

%3. C_l per a alpha=2 i per diferents M_inf%%
rho_aire=1.225;
T_inf=288.15;
R=287;

%[Mach_vector, Q_inf_3, Cl_3, Cl_3_corregit] = apartat3(M_crit, gamma_aire, R, T_inf, rho_aire, gamma, l);

