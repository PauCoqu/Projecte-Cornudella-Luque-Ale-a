function [Cl, Cm14, Cm0, Cp] = calcular_CL_CM14(gamma, l, Q_inf_modul, X, Z, X_c, Z_c)

    N = length(gamma);
    Cp = zeros(N,1);
    Vi_ext = zeros(N,1);
    F = 0;          
    Cm0_sum = 0;    
    c = 1; % Longitud de la corda (assumida constant)

    for i = 1:N
        Vi_ext(i) = abs(gamma(i));  % Intensitat de la velocitat induïda (pag. 49)
        Cp(i) = 1 - (gamma(i) / Q_inf_modul)^2;          
        F = F + gamma(i) * l(i); %sumatori de contribucions dels panells sobre el Lift (pag. 46)
        
        deltaX = X(i+1) - X(i);
        deltaZ = Z(i+1) - Z(i);
        Cm0_sum = Cm0_sum + Cp(i) * ((X_c(i) * deltaX + Z_c(i) * deltaZ) / c^2);
    end

    Cl = 2 * F / (Q_inf_modul * c);
    Cm0 = Cm0_sum;
    Cm14 = Cm0 + 0.25 * Cl;

    %Corregim el Cp al perfil N/4 (mitjana dels veins)
    Cp(N/4) = (Cp((N/4)+1) + Cp((N/4)-1))/2;
