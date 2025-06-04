function [Cl_alpha_22112,Cl_alpha_0012, Cl_0_0012, Cl_0_22112] = parametres_perfils ()

%Emmagatzemem els coeficients aerodinamics Cl i de moment Cm 1/4 de l'apartat 1 per cada angle d'atac
alpha = [0,1,2,3,4,5,6,7,8,9,10];
alpha_rad = deg2rad(alpha);

%Parametres NACA 0012 (Canard) --> Perfil simetric (Cl_0 = 0;)
%Cl_0012 = [0, -0.0285, 0.1829, 0.4492, 0.5431, 0.6194, 0.6919, 0.7691]; %Airfool tools (de 0-8º)
Cl_0012 = [0, 0.1198, 0.2395, 0.3591, 0.4787, 0.5981, 0.7173, 0.8363, 0.9550, 1.0735, 1.1916]; %Nostres valors (de 0-10º)
Cm_1I4_0012 = [0, -0.0017, -0.0035, -0.0051, -0.0066, -0.0081, -0.0093, -0.0104, -0.0113, -0.0119, -0.0123]; %Nostres valors (de 0-10º)

%Parametres NACA 22112 (Ala)  
%Cl_22112 = [-0.0329, 0.3280, 0.4460, 0.5233, 0.6058, 0.6908, 0.7685, 0.8397, 0.8556]; %Airfool tools (de 0-8º)
Cl_22112 = [0.0896, 0.2074 ,0.3252, 0.4429 , 0.5604, 0.6778 , 0.7950, 0.9119 , 1.0285, 1.1449 , 1.2608]; %Nostres valors (de 0-10º)
Cm_1I4_22112 = [9.5e-4, -5.3067e-4, -0.0020, -0.0035, -0.0049, -0.0062, -0.0074, -0.0084, -0.0093, -0.01, -0.0104]; %Nostres valors (de 0-10º)


%Fem un ployfit; Cl = Cl_0 + Cl_alpha*alpha
pCl_0012 = polyfit(alpha_rad, Cl_0012, 1); %p(1) = Cl_alpha ; p(2) = Cl_0;
pCl_22112 = polyfit(alpha_rad,Cl_22112,1); %p(1) = Cl_alpha ; p(2) = Cl_0;

%Amb polyval podem saber qualsevol valor del polyfit Cl_05 = ployval(p,0.5);
Cl_alpha_22112 = pCl_22112(1);
Cl_alpha_0012 = pCl_0012(1);
Cl_0_22112 = pCl_22112(2);
Cl_0_0012 = pCl_0012(2);

Cm_mig_22112 = mean(Cm_1I4_22112);
Cm_mig_0012 = mean(Cm_1I4_0012);

% figure
% plot(alpha_rad, Cl_22112, '-o', 'MarkerSize', 6, 'DisplayName', 'NACA 22112'); 
% hold on;
% plot(alpha_rad, Cl_0012, '-s', 'MarkerSize', 6, 'DisplayName', 'NACA 0012'); 
% xlabel('\alpha [rad]');
% ylabel('C_l');
% title('C_l vs \alpha');
% grid on;
% legend('Location', 'best');
% hold off;

% figure
% plot(alpha_rad, Cm_1I4_22112, '-o', 'MarkerSize', 6, 'DisplayName', 'NACA 22112'); 
% hold on;
% plot(alpha_rad, Cm_1I4_0012, '-s', 'MarkerSize', 6, 'DisplayName', 'NACA 0012'); 
% yline(Cm_mig_22112, '--', 'Color', [0, 0.5, 0], 'LineWidth', 2, 'DisplayName', 'C_{m_{mig}} 22112'); 
% yline(Cm_mig_0012, ':', 'Color', [0.7, 0, 0], 'LineWidth', 2, 'DisplayName', 'C_{m_{mig}} 0012');
% xlabel('\alpha [rad]');
% ylabel('C_{m_{1/4}}');
% title('C_{m_{1/4}} vs \alpha');
% grid on;
% legend('Location', 'best');
% hold off;

end

