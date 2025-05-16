function [Cl, Cm14, Cm0, Cp] = calcular_CL_CM14(gamma, l, Q_inf_modul, X, Z, X_c, Z_c, Normal)

    N = length(gamma);

    Cl = 0;
    Cm0 = 0;
    Cp = zeros(N,1);
    Vi_ext = zeros(N,1);

    F = 0;          % Força de sustentació (proporcional a la circulació total)
    Cm0_sum = 0;    % Suma per calcular el moment respecte l'origen

    c = 1; % Longitud de la corda (assumida constant)

    for i = 1:N
        partf = gamma(i) * l(i);                       % Contribució de cada panell a la circulació
        Vi_ext(i) = abs(gamma(i));                     % Intensitat de la velocitat induïda
        Cp(i) = 1 - (gamma(i) / Q_inf_modul)^2;  %Hi havia Q_inf_modul          
        F = F + partf;
        
        deltaX = X(i+1) - X(i);
        deltaZ = Z(i+1) - Z(i);
        Cm0_sum = Cm0_sum + Cp(i) * ((X_c(i) * deltaX + Z_c(i) * deltaZ) / c^2);
    end

    Cl = 2 * F / (Q_inf_modul * c);
    Cm0 = Cm0_sum;
    Cm14 = Cm0 + 0.25 * Cl;


%%Plot perfil+ coeficient pressions:
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
title('Vectors del coeficient de pressió C_p (cap a fora)');
grid on;
hold off;

end
