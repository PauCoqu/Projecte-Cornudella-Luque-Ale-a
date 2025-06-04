%% CODI PROJECTE PART 2
%---------------------------------------------
% Data: 20/05/2025
% Membres: Antonio Luque, Pau Cornudella i Alex Aleñà
% Assignatura: AERODINÀMICA, MECÀNICA DE VOL I ORBITAL
% Grup: 6
% ---------------------------------------------
clc; clear; close all;

N = 512;
%Variables ala
b = 24;
c_r = 1.8;
c_t = 1.2;
lambda = c_t/c_r; %(pag. 24)
c_mitja = (2/3)*c_r*(1+lambda+lambda^2)/(1+lambda);%(pag. 24)
S = c_mitja*b;

%Variables canard
b_h = 6; 
c_rh = 1;
c_th = 0.6;
l_h = 7.5;
lambda_h = c_th/c_rh; %(pag. 24)
c_mitjah = (2/3)*c_rh*(1+lambda_h+lambda_h^2)/(1+lambda_h);%(pag. 24)
S_h = c_mitjah*b_h;

%Variables VTP
Cd_VTP = 0.0062;
S_v = 2.1;

Q_inf = 1;
i_w = deg2rad(0);
i_h = deg2rad(4); 
rho = 1.225;  %EEEEEEEPPPP
Re = (rho*Q_inf*c_mitjah)/(1.81e-5);


%Rectes Cl vs alpha de cada perfil
[Cl_alpha_22112,Cl_alpha_0012, Cl_0_0012, Cl_0_22112] = parametres_perfils ();

%Discretitzem l'ala i el canard, i la coorda en cada punt
[Coords_ala, Coords_centre_ala, c_ala, Coords_canard, Coords_centre_canard, c_canard] = geometria_avio (N,b,c_r,c_t,b_h,c_rh,c_th,l_h);

%1)Definir l'angle de twist adequat+distribució sustentació
%Definim un angle de twist=0º i calculem la distribució de lift.
%Un cop ho tiguem ho grafiquem també per altres angles de twist i el que
%millor distribució de lift doni serà el que triem.

twist_tip = [-10,-8,-6,-4,-2,0,2,4,6,8,10]; %angle de twist a la punta de l'ala
twist_tip = deg2rad(twist_tip);
alpha_ala = deg2rad(4); %enunciat

for i  = 1: length(twist_tip)
[twist_centre_panell] = calcul_twist(twist_tip(i), N);
[gama_centre_panell] = calcul_gama(c_ala, alpha_ala, Cl_0_22112, Cl_alpha_22112, N, Coords_centre_ala, Coords_ala, i_w,Q_inf,twist_centre_panell);
end



