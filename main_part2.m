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
b_a = 24;
c_r = 1.8;
c_t = 1.2;
lambda = c_t/c_r; %(pag. 24)
c_mitja = (2/3)*c_r*(1+lambda+lambda^2)/(1+lambda);%(pag. 24)
S = c_mitja*b_a;
Cd_visc_coeffs_ala = [0.0080, - 0.0013, 0.0063]; %enunciat (ala)

%Variables canard
b_h = 6; 
c_rh = 1;
c_th = 0.6;
l_h = 7.5;
lambda_h = c_th/c_rh; %(pag. 24)
c_mitjah = (2/3)*c_rh*(1+lambda_h+lambda_h^2)/(1+lambda_h);%(pag. 24)
S_h = c_mitjah*b_h;
Cd_visc_coeffs_can = [0.0052 , 0.0071 ]; %enunciat (canard)

%Variables VTP
Cd_VTP = 0.0062;
S_v = 2.1;

%Altres variables
Q_inf = 1;
i_w = deg2rad(0);
i_h = deg2rad(3); 
rho = 1.225; 
Re = (rho*Q_inf*c_mitjah)/(1.81e-5);

%%

%Rectes Cl vs alpha de cada perfil
[Cl_alpha_22112,Cl_alpha_0012, Cl_0_0012, Cl_0_22112, Cl_0_0012_flap12, Cl_alpha_0012_flap12, Cm_0_0012_flap12, Cm_0_22112, Cm_0_0012] = parametres_perfils ();

%Discretitzem l'ala i el canard, i la coorda en cada punt
[Coords_ala, Coords_centre_ala, c_ala, Coords_canard, Coords_centre_canard, c_canard] = geometria_avio (N,b_a,c_r,c_t,b_h,c_rh,c_th,l_h);

%1)Definir l'angle de twist adequat+distribució sustentació
%Definim un angle de twist=0º i calculem la distribució de lift.
%Un cop ho tiguem ho grafiquem també per altres angles de twist i el que
%millor distribució de lift doni serà el que triem.

%twist_tip = [-1,-0.99,-0.98,-0.97,-0.96,-0.95,-0.94,-0.93,-0.92,-0.91,-0.9,-0.89,-0.88,-0.87,-0.86]; %angle de twist a la punta de l'ala
twist_tip = [-10,-8,-6,-4,-2,0,2,4,6,8,10];
twist_tip = deg2rad(twist_tip);
alpha_ala = deg2rad(4); %enunciat


% Inicialització de les matrius per emmagatzemar els resultats
CL_ala_ap1 = zeros(1, length(twist_tip));
CL_pan_ap1 = zeros(N, length(twist_tip));
alpha_ind_ap1 = zeros(N, length(twist_tip));
Cd_visc_pan_ap1 = zeros(N, length(twist_tip));
Cd_ind_ap1 = zeros(N, length(twist_tip));
Cd_tot = zeros(N, length(twist_tip));
CD_ap1 = zeros(1, length(twist_tip)); 
Eff_ap1 = zeros(1, length(twist_tip));
Lift_ap1 = zeros(1, length(twist_tip));

u_r = [-cos(alpha_ala); 0; sin(alpha_ala)]; %vector unitari de la velocitat incident
Cd_visc_coeffs = Cd_visc_coeffs_ala;
% Bucle per recórrer cada valor de twist_tip
for i = 1:length(twist_tip)
    % Càlculs per a cada iteració
    [twist_centre_panell] = calcul_twist(twist_tip(i), N);
    [gamma_centre_panell] = calcul_gama(c_ala, alpha_ala, Cl_0_22112, Cl_alpha_22112, N, Coords_centre_ala, Coords_ala, i_w, Q_inf, twist_centre_panell, u_r);
    [CL_ala_ap1(i), CL_pan_ap1(:, i), alpha_ind_ap1(:, i), Cd_visc_pan_ap1(:, i), Cd_ind_ap1(:, i),Cd_tot(:, i), CD_ap1(i), Eff_ap1(i), Lift_ap1(i)] = calcul_coef(N, gamma_centre_panell, alpha_ala, Cl_alpha_22112, Cl_0_22112, i_w, twist_centre_panell, Coords_ala, rho, Q_inf, S, c_ala, Cd_visc_coeffs);
end

[eff_max, filaMax] = max(Eff_ap1);
fprintf('El valor de màxima eficiència de twisting és %g i es troba a un valor de %d\n graus de twist', eff_max, twist_tip(filaMax));


%% GRAFICS

%1) CL al llarg de l'envergadura per diferents angles de twist
figure;
spanwise_pos = Coords_centre_ala(:, 2);
hold on;
for i = 1:length(twist_tip)
    plot(spanwise_pos, CL_pan_ap1(:, i), 'LineWidth', 1.5, 'DisplayName', sprintf('$\\theta = %.2f^\\circ$', rad2deg(twist_tip(i))));
end

title('(C_{L}) front envergadura');
xlabel('Posició en l''envergadura (m)');
ylabel('C_{L}');
legend('show', 'Interpreter', 'latex');
grid on;
xlim([min(spanwise_pos), max(spanwise_pos)]);
ylim([min(CL_pan_ap1(:)), max(CL_pan_ap1(:))]);
hold off;


%2) CD al llarg de l'envergadura per diferents angles de twist
figure;
hold on;
for i = 1:length(twist_tip)
    plot(spanwise_pos, Cd_tot(:, i), 'LineWidth', 1.5, 'DisplayName', sprintf('$\\theta = %.2f^\\circ$', rad2deg(twist_tip(i))));
end

title('C_{D} front envergadura');
xlabel('Posició en l''envergadura (m)');
ylabel('C_{D}');
legend('show', 'Interpreter', 'latex');
grid on;
xlim([min(spanwise_pos), max(spanwise_pos)]);
ylim([min(Cd_tot(:)), max(Cd_tot(:))]);
hold off;

%3) Eficiència al llarg de l'envergadura per diferents angles de twist
theta_deg = rad2deg(twist_tip);    % 1×11
figure;
plot(theta_deg, Eff_ap1, '-o', 'LineWidth',1.5)
text(theta_deg, Eff_ap1, compose('%.1f', Eff_ap1), ...
     'VerticalAlignment','bottom','HorizontalAlignment','right')
grid on
xlabel('$\theta\;(^\circ)$','Interpreter','latex')
ylabel('$C_L/C_D$','Interpreter','latex')
title('Eficiència vs. angle de twist','Interpreter','latex')

%2) CD_viscós al llarg de l'envergadura per diferents angles de twist
figure;
hold on;
for i = 1:length(twist_tip)
    plot(spanwise_pos, Cd_visc_pan_ap1(:, i), 'LineWidth', 1.5, 'DisplayName', sprintf('$\\theta = %.2f^\\circ$', rad2deg(twist_tip(i))));
end

title('C_{D,visc} front envergadura');
xlabel('Posició en l''envergadura (m)');
ylabel('C_{D,visc}');
legend('show', 'Interpreter', 'latex');
grid on;
xlim([min(spanwise_pos), max(spanwise_pos)]);
ylim([min(Cd_visc_pan_ap1(:)), max(Cd_visc_pan_ap1(:))]);
hold off;


%2) CD_induit al llarg de l'envergadura per diferents angles de twist
figure;
hold on;
for i = 1:length(twist_tip)
    plot(spanwise_pos, Cd_ind_ap1(:, i), 'LineWidth', 1.5, 'DisplayName', sprintf('$\\theta = %.2f^\\circ$', rad2deg(twist_tip(i))));
end

title('C_{D,ind} front envergadura');
xlabel('Posició en l''envergadura (m)');
ylabel('C_{D,ind}');
legend('show', 'Interpreter', 'latex');
grid on;
xlim([min(spanwise_pos), max(spanwise_pos)]);
ylim([min(Cd_ind_ap1(:)), max(Cd_ind_ap1(:))]);
hold off;


%% 2.2

%Definim un twist de màxim rendiment

Twist_ala = deg2rad(-0.89);
Twist_flap = 0;
[twist_centre_panell_actualitzat] = calcul_twist(Twist_ala, N);
[twist_centre_flap] = calcul_twist(Twist_flap, N);
alpha_ala=deg2rad(4);

[gamma, A, b] = gamma_ala_can(alpha_ala, i_w, i_h, Q_inf, N, Coords_ala, Coords_centre_ala, c_ala, twist_centre_panell_actualitzat, ...
    Cl_alpha_22112, Cl_0_22112, Coords_canard, Coords_centre_canard, c_canard, Cl_alpha_0012, Cl_0_0012, u_r);

gamma_ala=gamma(1:N,1);
gamma_can=gamma(N+1:end,1);

[CL_ala, CL_panala, alpha_ind_ala, Cd_visc_panala, Cd_ind_ala, Cd_tot_ala, CD_ala, Eff_ala, Lift_ala] = ...
    calcul_coef(N, gamma_ala, alpha_ala, Cl_alpha_22112, Cl_0_22112, i_w, twist_centre_panell_actualitzat, Coords_ala, rho, Q_inf, S, c_ala, Cd_visc_coeffs_ala);

[CL_can, Cl_pancan, alpha_ind_can, Cd_visc_pancan, Cd_ind_can, Cd_tot_can, CD_can, Eff_can, Lift_can] = ...
    calcul_coef(N, gamma_can, alpha_ala, Cl_alpha_0012, Cl_0_0012, i_h, twist_centre_flap, Coords_canard, rho, Q_inf, S_h, c_canard, Cd_visc_coeffs_can);

CM = calculCentreMassa(Coords_ala, Coords_canard, N, Lift_ala, Lift_can, Cm_0_0012, Cm_0_22112, c_mitja, c_mitjah);


%% 3. Aerodynamic polar curve per angles d'atac de 0 a 6.

alpha_ala =  0:0.5:6;
alpha_ala = deg2rad(alpha_ala);

CL_ap3=zeros(length(alpha_ala),1);
CD_ap3=zeros(length(alpha_ala),1);

for i=1:length(alpha_ala)
u_r = [-cos(alpha_ala(i)); 0; sin(alpha_ala(i))]; %vector unitari de la velocitat incident
[gamma_ap3, A_ap3, b_ap3] = gamma_ala_can(alpha_ala(i), i_w, i_h, Q_inf, N, Coords_ala, Coords_centre_ala, c_ala, twist_centre_panell_actualitzat, ...
    Cl_alpha_22112, Cl_0_22112, Coords_canard, Coords_centre_canard, c_canard, Cl_alpha_0012, Cl_0_0012, u_r);
gamma_ala_ap3=gamma_ap3(1:N,1);
gamma_can_ap3=gamma_ap3(N+1:end,1);


[CL_ala_ap3, CL_panala, alpha_ind_ala, Cd_visc_panala, Cd_ind_ala, Cd_tot_ala, CD_ala_ap3, Eff_ala, Lift_ala] = ...
    calcul_coef(N, gamma_ala_ap3, alpha_ala(i), Cl_alpha_22112, Cl_0_22112, i_w, twist_centre_panell_actualitzat, Coords_ala, rho, Q_inf, S, c_ala, Cd_visc_coeffs_ala);


[CL_can_ap3, Cl_pancan, alpha_ind_can, Cd_visc_pancan, Cd_ind_can, Cd_tot_can, CD_can_ap3, Eff_can, Lift_can] = ...
    calcul_coef(N, gamma_can_ap3, alpha_ala(i), Cl_alpha_0012, Cl_0_0012, i_h, twist_centre_flap, Coords_canard, rho, Q_inf, S_h, c_canard, Cd_visc_coeffs_can);

CL_ap3(i,1) = CL_ala_ap3 + CL_can_ap3;
CD_ap3(i,1) = CD_ala_ap3 + CD_can_ap3 + Cd_VTP;
end


%Corba polar aerodinàmica
figure;
hold on;
plot(CL_ap3, CD_ap3, 'LineWidth', 1.5);
title('Corba polar aerodinàmica');
xlabel('C_{L}');
ylabel('C_{D}');
grid on;
box on
xlim([min(CL_ap3), max(CL_ap3)]);
ylim([min(CD_ap3), max(CD_ap3)]);
hold off;






%% 4. Calculs per alpha = 4.

alpha_ala =  4;
alpha_ala = deg2rad(alpha_ala);

u_r = [-cos(alpha_ala); 0; sin(alpha_ala)]; %vector unitari de la velocitat incident
[gamma_ap4, A_ap4, b_ap4] = gamma_ala_can(alpha_ala, i_w, i_h, Q_inf, N, Coords_ala, Coords_centre_ala, c_ala, twist_centre_panell_actualitzat, ...
    Cl_alpha_22112, Cl_0_22112, Coords_canard, Coords_centre_canard, c_canard, Cl_alpha_0012_flap12, Cl_0_0012_flap12, u_r);
gamma_ala_ap4=gamma_ap4(1:N,1);
gamma_can_ap4=gamma_ap4(N+1:end,1);

[CL_ala_ap4, ~, ~, ~, ~, ~, CD_ala_ap4, ~, ~] = ...
    calcul_coef(N, gamma_ala_ap4, alpha_ala, Cl_alpha_22112, Cl_0_22112, i_w, twist_centre_panell_actualitzat, Coords_ala, rho, Q_inf, S, c_ala, Cd_visc_coeffs_ala);

[CL_can_ap4, ~, ~, ~, ~, ~, CD_can_ap4, ~, ~] = ...
    calcul_coef(N, gamma_can_ap4, alpha_ala, Cl_alpha_0012_flap12, Cl_0_0012_flap12, i_h, twist_centre_flap, Coords_canard, rho, Q_inf, S_h, c_canard, Cd_visc_coeffs_can);

CL_ap4 = CL_ala_ap4 + CL_can_ap4;
CD_ap4 = CD_ala_ap4 + CD_can_ap4 + Cd_VTP;


Cm_CM = Cm_0_0012_flap12 + Cm_0_22112 + CL_can_ap4*(CM-Coords_canard(N,1)) + CL_ala_ap4*(CM-Coords_ala(N,1));
