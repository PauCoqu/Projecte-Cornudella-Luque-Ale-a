function [gamma, A, b] = gamma_ala_can(alpha, incid_a, incid_can, Q_inf, N, Coords_a, Coords_c_a, c_a, twist_a, ...
    Cl_alpha_22112, Cl_0_22112, Coords_can, Coords_c_can, c_can, Cl_alpha_0012, Cl_0_0012, u_r)

Na=N;
Nc=N;
k=[0; 0; 1];
A= zeros(N,N);
b=zeros(N,1);

% A11, b11 (ala sobre ala)
for i=1:Na
    b_11(i)=0.5*c_a(i)*Q_inf*(Cl_0_22112+Cl_alpha_22112*(alpha+twist_a(i)+incid_a));
    for j=1:Na
        if j==i
            V_ii = auto_velocinduida(Coords_c_a(i,:), Coords_a(j,:), Coords_a(j+1,:), u_r);
            A_11(i,j)=-0.5*Cl_alpha_22112*c_a(i)*dot(V_ii,k)+1;

            
        else
            V_ij = velocinduida(Coords_c_a(i,:), Coords_a(j,:), Coords_a(j+1,:), u_r);
            A_11(i,j)=-0.5*Cl_alpha_22112*c_a(i)*dot(V_ij,k);
        end
    end
end


% A12, b12 (canard sobre ala)
for i=1:Na
    for j=1:Nc
        V_ij_comp= velocinduida(Coords_c_a(i,:), Coords_can(j,:), Coords_can(j+1,:), u_r);
        A_12(i,j)=-0.5*Cl_alpha_22112*c_a(i)*dot(V_ij_comp,k);
    end
end

% A21, b21 (ala sobre canard)
for i=1:Nc
    for j=1:Na
        V_ij_comp= velocinduida(Coords_c_can(i,:), Coords_a(j,:), Coords_a(j+1,:), u_r);
        A_21(i,j)=-0.5*Cl_alpha_0012*c_a(i)*dot(V_ij_comp,k);
    end
end


% A22, b22 (can sobre can)
for i=1:Nc
    b_22(i)=0.5*c_a(i)*Q_inf*(Cl_0_0012+Cl_alpha_0012*(alpha+twist_a(i)+incid_can));
    for j=1:Nc
        if j==i
            V_ii = auto_velocinduida(Coords_c_can(i,:), Coords_can(j,:), Coords_can(j+1,:), u_r);
            A_22(i,j)=-0.5*Cl_alpha_0012*c_can(i)*dot(V_ii,k)+1;

            
        else
            V_ij = velocinduida(Coords_c_can, Coords_can(j,:), Coords_can(j+1,:), u_r);
            A_22(i,j)=-0.5*Cl_alpha_0012*c_can(i)*dot(V_ij,k);
        end
    end
end



%NO CAL FER LA INCIDÈNCIA DEL VTP SOBRE L'ALA I EL FLAP PERQUÈ NO ES
%GENEREN ENTRE ELLS, (SÓN PERPENDICULARS)


A = [A11, A12; A21, A22];
b = [b11;b22];

gamma=A\b;