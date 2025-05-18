function [Cl, Cm14, Cm0, Cp] = calcular_CL_CM14(gamma, l, Q_inf_modul, X, Z, X_c, Z_c, Normal)

    N = length(gamma);
    Cp = zeros(N,1);
    Vi_ext = zeros(N,1);
    F = 0;          
    Cm0_sum = 0;    
    c = 1; % Longitud de la corda (assumida constant)

    for i = 1:N
        Vi_ext(i) = abs(gamma(i));  % Intensitat de la velocitat induïda (pag. 49)
        Cp(i) = 1 - (gamma(i) / Q_inf_modul)^2;          
        F = F + gamma(i) * l(i); %sumatori de contribucions dels panells sobre el Lift
        
        deltaX = X(i+1) - X(i);
        deltaZ = Z(i+1) - Z(i);
        Cm0_sum = Cm0_sum + Cp(i) * ((X_c(i) * deltaX + Z_c(i) * deltaZ) / c^2);
    end

    Cl = 2 * F / (Q_inf_modul * c);
    Cm0 = Cm0_sum;
    Cm14 = Cm0 + 0.25 * Cl;


    %Corregim el Cp al perfil N/4
    Cp(N/4) = (Cp((N/4)+1) + Cp((N/4)-1))/2;


%Plot perfil+ coeficient pressions:
figure;
hold on;

% Dibuixar el perfil
plot(X, Z, 'k-', 'DisplayName', 'Perfil');

% Escala dels vectors
scale = 0.5;

for i = 1:length(Cp)
    % Escalem la normal amb el valor absolut de Cp
    dx = scale * abs(Cp(i)) * Normal(i,1);
    dz = scale * abs(Cp(i)) * Normal(i,2);

    % Assignem color segons el signe de Cp (però el vector sempre surt cap a fora)
    if Cp(i) < 0
        color = 'b';
    else
        color = 'r';
    end

    % Dibuix del vector
    quiver(X_c(i), Z_c(i), dx, dz, 0, color, 'LineWidth', 1);
end

axis equal;
xlabel('X');
ylabel('Z');
title('Vectors del coeficient de pressió C_p');
grid on;
hold off;

end

