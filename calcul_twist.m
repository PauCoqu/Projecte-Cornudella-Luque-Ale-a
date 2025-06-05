function [twist_centre_panell] = calcul_twist(twist_tip, N)

%El twist varia linealment des de la punta de l'ala al centre i del centre
twist_semiala_1 = linspace(twist_tip, 0, N/2 + 1); 
twist_semiala_2 = linspace(0, twist_tip, N/2 + 1);
twist_ala = [twist_semiala_1, twist_semiala_2(2:end)]; % Així no repetim el 0
twist_centre_panell = zeros(N,1);
for i = 1:N
       twist_centre_panell(i) = (twist_ala(i) + twist_ala(i + 1)) / 2; %Twist al centre del panell
end
twist_centre_panell = twist_centre_panell'; %transposem en columna
end