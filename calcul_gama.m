function [gama_centre_panell] = calcul_gama(c_ala, alpha_ala, Cl_0_ala, Cl_alpha_ala, N, Coords_centre_ala, Coords_ala, i_ala,Q_inf,twist)
 
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
b_i(i) = (1/2)*c_ala(i)*Q_inf*(Cl_0_ala + Cl_alpha_ala*(alpha_ala + i_ala + twist(i))); % sempre constant (pag 19)

    for j = 1:N
        if j == i %cas vortex autoinduit  
            %Per cada panell calculem la velocitat del vortex autoinduit (V_ij = VinfA + V_AB - V_infB)
            r1 = (Coords_centre_ala(i,:) - Coords_ala(j,:))'; %(pag 9)
            r2 = (Coords_centre_ala(i,:) - Coords_ala(j+1,:))'; %(pag 9)
            u_r = [-cos(alpha_ala); 0; sin(alpha_ala)]; %vector unitari de la velocitat incident
            u_r1 = r1/norm(r1); %vector unitari u_r1
            u_r2 = r2/norm(r2); %vector unitari u_r1
            V_infA = 1/(4*pi)*(1-dot(u_r,u_r1))/(norm(cross(u_r,r1))).^2 *cross(u_r,r1); %(pag 10)
            V_infB = 1/(4*pi)*(1-dot(u_r,u_r2))/(norm(cross(u_r,r2))).^2 *cross(u_r,r2); %(pag 10)
            %V_AB = (1/(4*pi))*(norm(r1)+norm(r2))/(norm(r1)*norm(r2)*(norm(r1)*norm(r2) + dot(r1,r2)))*cross(r1,r2); %(pag 9) hauria de donar 0 en aquest cas
            V_ij = V_infA - V_infB; %(pag 11)    
            a_ij(i, j) = (-1/2)*Cl_alpha_ala*c_ala(i)*dot(V_ij, Referencia_k) + 1; %(pag 19)
        else
            % Per cada panell calculem la velocitat del vortex (V_ij = VinfA + V_AB - V_infB)
            r1 = (Coords_centre_ala(i,:) - Coords_ala(j,:))'; %(pag 9)
            r2 = (Coords_centre_ala(i,:) - Coords_ala(j+1,:))'; %(pag 9)
            u_r = [-cos(alpha_ala); 0; sin(alpha_ala)]; %vector unitari de la velocitat incident
            u_r1 = r1/norm(r1); %vector unitari u_r1
            u_r2 = r2/norm(r2); %vector unitari u_r1
            V_infA = 1/(4*pi)*(1-dot(u_r,u_r1))/(norm(cross(u_r,r1)))^2 *cross(u_r,r1); %(pag 10)
            V_infB = 1/(4*pi)*(1-dot(u_r,u_r2))/(norm(cross(u_r,r2)))^2 *cross(u_r,r2); %(pag 10)
            V_AB = (1/(4*pi))*(norm(r1)+norm(r2))/(norm(r1)*norm(r2)*(norm(r1)*norm(r2) + dot(r1,r2)))*cross(r1,r2); %(pag 9) hauria de donar 0 en aquest cas
            V_ij = V_infA + V_AB - V_infB; %(pag 11)    
            a_ij(i, j) = (-1/2)*Cl_alpha_ala*c_ala(i)*dot(V_ij, Referencia_k) + 1; %(pag 19) 
        end
    end
end

% A*gamma = b
gama_centre_panell = a_ij \ b_i;

end