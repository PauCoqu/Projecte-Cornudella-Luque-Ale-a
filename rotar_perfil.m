function rotar_perfil(X, Z,alpha_rad, Normal,Cp)

    % Matriu de rotació
    R = [cos(alpha_rad), sin(alpha_rad);
         -sin(alpha_rad), cos(alpha_rad)];
     
    % Coordenades rotades
    coords = R * [X; Z];
    
    % Retornar les noves coordenades
    X = coords(1, :);
    Z = coords(2, :);

    % Punts de control per panell
    X_c = (X(1:end-1) + X(2:end)) / 2;
    Z_c = (Z(1:end-1) + Z(2:end)) / 2;

    %Grafic
    figure;
    hold on;

    % Dibuixar el perfil
    plot(X, Z, 'k-', 'DisplayName', 'Perfil');

    % Escala dels vectors
    scale = 0.5;

    % Colors
    blau_suau = [0.4 0.6 1];
    vermell_suau = [1 0.4 0.4];

    for i = 1:length(Cp)
    n = Normal(i,:) / norm(Normal(i,:));  % normalitzem la normal
    dx = scale * abs(Cp(i)) * n(1);
    dz = scale * abs(Cp(i)) * n(2); 

    % Assignem color segons signe de Cp
    if Cp(i) < 0
        color = blau_suau;
    else
        color = vermell_suau;
    end

    % Dibuix del vector
    quiver(X_c(i), Z_c(i), dx, dz, 0, 'Color', color, 'LineWidth', 0.5);
    end

    axis equal;
    xlabel('X');
    ylabel('Z');
    title('Vectors del coeficient de pressió C_p');
    grid on;
    hold off;


end
