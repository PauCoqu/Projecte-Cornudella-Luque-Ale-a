%% CODI PROJECTE PART 2
%---------------------------------------------
% Data: 20/05/2025
% Membres: Antonio Luque, Pau Cornudella i Alex Aleñà
% Assignatura: AERODINÀMICA, MECÀNICA DE VOL I ORBITAL
% Grup: 6
% ---------------------------------------------
clc; clear; close all;

%Variables ala
Nw = 512; % Number of spanwise segments for the wing
b = 24;
c_r = 1.8;
c_t = 1.2;
lambda = c_t/c_r; %(pag. 24)
c_mitja = (2/3)*c_r*(1+lambda+lambda^2)/(1+lambda);%(pag. 24)
S = c_mitja*b;

%Variables canard
Nc = 512; % Number of spanwise segments for the canard
b_h = 6; 
c_rh = 1;
c_th = 0.6;
l_h = 7.5;
lambda_h = c_th/c_rh; %(pag. 24)
c_mitjah = (2/3)*c_rh*(1+lambda_h+lambda_h^2)/(1+lambda_h);%(pag. 24)
S_h = c_mitjah*b_h;

%Variables VTP
Q_inf = 1;
Cd_VTP = 0.0062;
i_w = deg2rad(0);
i_h = deg2rad(4); 
S_v = 2.1;
rho = 1.225;  %EEEEEEEPPPP
Re = (rho*Q_inf*c_mitjah)/(1.81e-5);


%Rectes Cl vs alpha de cada perfil
[Cl_alpha_22112,Cl_alpha_0012, Cl_0_0012, Cl_0_22112] = parametres_perfils ();

%Geometria avió
%[A] = Geometria_avio ();





%FEM PRIMER UN CAS AMB TWIST=0º I LLAVORS HO FIQUEM DINS UN FOR I GRAFIQUEM




