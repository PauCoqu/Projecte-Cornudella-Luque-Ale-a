function [X, Z, alpha_rad] = rotar_perfil(X, Z, alpha, alpha_rad)

    % Convertim l'angle a radians
    alpha_rad = deg2rad(alpha);
    
    % Matriu de rotació
    R = [cos(alpha_rad), sin(alpha_rad);
         -sin(alpha_rad), cos(alpha_rad)];
     
    % Coordenades rotades
    coords = R * [X; Z];
    
    % Retornar les noves coordenades
    X = coords(1, :);
    Z = coords(2, :);

    % Visualitzar el perfil
    figure;
    plot(X, Z, '-o');
    title(['Perfil NACA 0012 amb angle d''atac \alpha = ', num2str(alpha), '°']);
    xlabel('X');
    ylabel('Z');
    grid on;
    axis equal;
end
