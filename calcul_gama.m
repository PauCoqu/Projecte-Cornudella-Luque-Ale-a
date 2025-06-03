function [gama_centre_panell] = calcul_gama(c_ala, alpha_ala, Cl_0_ala, Cl_alpha_ala, N, i_ala)
 
%c_ala(i) = corda ala o canard (c_ala, c_canard....)
%alpha_ala = angle ala o canard (alpha_ala, .....)
%Cl_0_ala = Cl_0_0012, Cl_0_22112 ....
%Cl_alpha_ala = Cl_alpha_0012, Cl_alpha_22112 ....
%i_ala = incidencia ala o canard (i_w, i_h ....)
%twist(i) = twist_centre_panell....

a_ij = zeros(N, N); %pag 19
b_i = zeros(N, 1); %pag 19

%Referencia_k = Vector perpendicular a Q_inf. Referencia_k = (-a,0,b) 
Referencia_k = [-sin(alpha_ala + i_ala); 0; cos(alpha_ala + i_ala)]; %(pag 16)

for i = 1:N

b_i(i) = (1/2)*c_ala(i)*Q_inf*(Cl_0_ala + Cl_alpha_ala*(alpha_ala + i_ala + twist(i))); %(pag 19)

    for j = 1:N
        if j == i
                V_ii = self_vortex(Xc(:, i), X(:, j), X(:, j + 1), Ur); %Velocitat autoinduida
                a_ij(i, j) = (-1/2)*Cl_alpha_ala*c_ala(i)*dot(V_ii, k)+1; %(pag 19) Vel
        else
                V_ij = horseshoe_vortex(Xc(:, i), X(:, j), X(:, j + 1), Ur); %Velocitat induida
                a_ij(i, j) = (-1/2)*Cl_alpha_ala * c(i) * dot(V_ij, k); %(pag 19) 
        end
    end
end

% A*gamma = b
gama_centre_panell = a_ij \ b_ij;

end