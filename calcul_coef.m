function [CL, CL_pan, alpha_ind, CD_visc_pan, CD_ind,CD_tot, CD, Eff, Lift] = calcul_coef(N, gamma_centre, alpha, Cl_alpha, Cl_0, incid, twist, Coords, rho, Q_inf, S, c, Cd_visc_coeffs)

circ_pan = zeros(N,1);
CL_pan = zeros(N,1);
alpha_ind = zeros(N,1);
D_ind = zeros(N,1);
CD_ind = zeros(N,1);
CD_visc_pan = zeros(N,1);
CD_tot = zeros(N,1);

for i=1:N
    Delta_y = abs(Coords(i+1,2)-Coords(i,2));
    circ_pan(i) = gamma_centre(i)*Delta_y;

    %Calcul del Drag induit
    CL_pan(i) = (2*gamma_centre(i))/(c(i)*norm(Q_inf)); %(pag 21)
    alpha_ind(i) = (CL_pan(i) - Cl_0)/(Cl_alpha) - (alpha + incid) - twist(i); %(pag 21)
    D_ind(i) = -rho * norm(Q_inf) * circ_pan(i) * alpha_ind(i); %(pag. 22)
    CD_ind(i) = 2 * D_ind(i) / (rho * norm(Q_inf)^2 * S); %(pag. 22)
     

    %Calcul arrossegament viscós
    
    Cd_visc = polyval(Cd_visc_coeffs, CL_pan(i)); % Avalua polinomi: a*CL^2 + b*CL + c
    CD_visc_pan(i) = (Cd_visc*c(i)*Delta_y)/S;
    CD_tot(i) = CD_visc_pan(i) + CD_ind(i); %CD total
       
    
end

circ_tot = sum(circ_pan);
CD = sum(CD_tot);
Lift = circ_tot*rho*norm(Q_inf); %(pag 20)
CL = 2*Lift / (rho*norm(Q_inf)^2*S);%(pag 20)
Eff = CL/CD;
