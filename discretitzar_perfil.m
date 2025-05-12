function [X_c, Z_c, N_panells,nx,nz,l] = discretitzar_perfil(X, Z)
    % Discretitza el perfil en panells i calcula els punts de control i els punts centrals 
   
    N_panells = length(X) - 1; %Numero de panells     
    nx = zeros(N_panells, 1);
    nz = zeros(N_panells, 1);
    l = zeros(N_panells, 1);
    
    % Calcular longituds i normals dels panells
    for j = 1:N_panells
        dx = X(j+1) - X(j);
        dz = Z(j+1) - Z(j);
        l(j) = sqrt(dx^2 + dz^2);
        nx(j) = dz / l(j);   %ALTANTO
        nz(j) = -dx / l(j);  %ALTANTO
    end
    
    % Punts de control per panell (el punt mitjà de cada panell)
    X_c = (X(1:end-1) + X(2:end)) / 2;
    Z_c = (Z(1:end-1) + Z(2:end)) / 2;

% Visualitzar el perfil i les normals
figure;
plot(X, Z, '-o', 'DisplayName', 'Perfil');
hold on;

% Punts de control
plot(X_c, Z_c, 'k.', 'MarkerSize', 10, 'DisplayName', 'Punts de control');

% Dibuixar normals cap a fora
scale = 0.2;  % Escala de les normals
for i = 1:length(X_c)
    x_start = X_c(i);
    z_start = Z_c(i);

    % Invertim la direcció de la normal
    x_end = x_start - scale * nx(i);
    z_end = z_start - scale * nz(i);

    plot([x_start, x_end], [z_start, z_end], 'r-', 'LineWidth', 1.5);
end

% Ajustos visuals
title('Perfil amb normals cap a fora');
xlabel('X'); ylabel('Z');
axis equal; grid on;

end
