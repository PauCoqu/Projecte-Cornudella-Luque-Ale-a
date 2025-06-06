function [CL, Cl_pan, alpha_ind, Cd_visc_pan, Cd_ind, CD, Eff, Lift] = .... 
    calcul_coef(N, gamma_centre_panell, alpha, Cl_alpha, Cl_0, incid, twist, Coords, rho, Q_inf, S, c, component)

circ_pan = zeros(N,1);
Cl_pan = zeros(N,1);
alpha_ind = zeros(N,1);
D_ind = zeros(N,1);
Cd_ind = zeros(N,1);
Cd_visc_pan = zeros(N,1);
Cd_tot = zeros(N,1);

for i=1:N
    delta_y = abs(Coords(i+1,2)-Coords(i,2));
    circ_pan(i) = gamma_centre_panell(i)*delta_y;

    %Calculem l'angle d'atac induit
    Cl_pan(i) = (2*gamma_centre_panell(i))/(c(i)*Q_inf); %pag 21
    alpha_ind(i) = (Cl_pan(i) - Cl_0) / Cl_alpha - (alpha + incid) - twist(i);


    %Càlcul de l'arrossegament induit
    D_ind(i) = -rho * norm(Q_inf) * circ_pan(i) * alpha_ind(i); %(pag. 22)
    Cd_ind(i) = 2 * D_ind(i) / (rho * norm(Q_inf)^2 * S); %(pag. 22)
     

    %Càlcul arrossegament viscós
    if component=="ala"
        Cd_visc=0.0080*Cl_pan(i)^2 - 0.0013*Cl_pan(i) + 0.0063;
    else
        Cd_visc=0.0052*Cl_pan(i) + 0.0071;
    end

    Cd_visc_pan(i) = (Cd_visc * c(i) * delta_y)/S;

    %Cd total
    Cd_tot(i) = Cd_visc_pan(i) + Cd_ind(i);
    
end
circ_tot = sum(circ_pan);
CD = sum(Cd_tot);
Lift = circ_tot * rho * Q_inf; %pag 20
CL = 2*Lift / (rho*Q_inf^2*S);
Eff = CL/CD;
