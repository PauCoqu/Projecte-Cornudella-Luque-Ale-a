function [gama_centre_panell] = calcul_gama(alpha_ala, N,Coords_ala)
 
a_ij = zeros(N, N); %pag 19
b_i = zeros(N, 1); %pag 19

%Referencia_k = Vector perpendicular a Q_inf. Referencia_k = (-a,0,b) 
Referencia_k = [-sin(alpha_ala + inc); 0; cos(alpha_ala + inc)]; %(pag 16)

for i = 1:N

b_i(i) = 0.5 * c(i) * Q_inf * (Cl0 + Cl_alpha_ala * (alpha_ala + t(i) + inc));

    for j = 1:N
        if j == i
             % Velocidad inducida por el propio vórtice (self-induced)
                V_ii = self_vortex(Xc(:, i), X(:, j), X(:, j + 1), Ur);
                a_ww(i, j) = -0.5 * Cl_alpha_ala * c(i) * dot(V_ii, k) + 1;
        else
                % Velocidad inducida por otros vórtices (interacción)
                V_ij = horseshoe_vortex(Xc(:, i), X(:, j), X(:, j + 1), Ur);
                a_ww(i, j) = -0.5 * Cl_alpha_ala * c(i) * dot(V_ij, k);
        end
    end
end

% Resolución del sistema lineal A * gamma = b
gamma_isol = a_ww \ b_ww;


%gama_centre_panell = 1;

end