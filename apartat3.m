function [Mach_vector, Q_inf_3, Cl_3, Cl_3_corregit] = apartat3(M_crit, gamma_aire, R, T_inf, rho_aire, gamma, l)

%Vector llista dels machs que ens diu l'enunciat. Calcular a=sqrt(gamma_air*R*T)

Mach_vector = M_crit - [0.16, 0.12, 0.08, 0.04, 0];
a=sqrt(gamma_aire*R*T_inf);

%Calculem vel. corrent lliure de l'apartat 3 Q_inf_3
Q_inf_3=Mach_vector*a;

%Calcular Lift i Cl per a cada Mach. 

L=rho_aire*Q_inf_3*sum(gamma .* l); %sum(gamma .* l) calcula gamma_1*l_1+gamma_2*l_2+...+gamma_N*l_N

Cl_3=L./(0.5*rho_aire*Q_inf_3.^2*1); %M1_L1 p.46

%Aplicar les correccions per règims compressibles:(M1_L3 a partir de la p.10)

beta=sqrt(1-Mach_vector.^2);

Cl_3_corregit=(1./beta) .*Cl_3;

end