function [twist_centre_panell] = calcul_twist(twist_tip, N);
%Cl_w_iso_t = zeros(size(theta));
%E_t_w = zeros(size(theta));
%Cl_w_iso_dis = zeros(Nw, length(theta));
%alpha_ind_w = zeros(Nw, length(theta));
%Cd_visc_w = zeros(Nw, length(theta));
%Cd_visc_t_w = zeros(1, length(theta));
%Cd_induced_w = zeros(Nw, length(theta));

  
%El twist varia linealment des de la punta de l'ala al centre i del centre
twist_semiala_1 = linspace(twist_tip, 0, N/2 + 1); 
twist_semiala_2 = linspace(0, twist_tip, N/2 + 1);
twist_ala = [twist_semiala_1, twist_semiala_2(2:end)]; % Així no repetim el 0

for i = 1:N
       twist_centre_panell(i) = (twist_ala(i) + twist_ala(i + 1)) / 2; %Twist al centre del panell
end

twist_centre_panell = twist_centre_panell';

end