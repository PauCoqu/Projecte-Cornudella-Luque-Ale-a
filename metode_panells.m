function [X_c, Z_c, l, Tangent, Normal, gamma, Q_inf, Q_inf_modul] = metode_panells(X, Z, alpha_rad);

    N = length(X) - 1; %Numero de panells     
    l = zeros(N, 1);
    Tangent = zeros(N, 2);
    Normal = zeros(N, 2);   
    calphaj = zeros (1,N);
    salphaj = zeros (1,N);
    
    % Calcular longituds i normals dels panells
    for j = 1:N
        dx = X(j+1) - X(j);
        dz = Z(j+1) - Z(j);
        l(j) = sqrt(dx^2 + dz^2);
        calphaj(j) = dx / l(j);         % cos(α) del panell
        salphaj(j) = -dz / l(j);        % -sin(α) del panell

        % Vector tangent i normal (sentit matemàtic estàndard)
        Tangent(j, :) = [calphaj(j), -salphaj(j)];
        Normal(j, :) = [salphaj(j),  calphaj(j)];
    
    end
    
    % Punts de control per panell
    X_c = (X(1:end-1) + X(2:end)) / 2;
    Z_c = (Z(1:end-1) + Z(2:end)) / 2;

    % Cálcul de gamma

    U_inf=1;
    Q_inf = [U_inf*cos(alpha_rad), U_inf*sin(alpha_rad)];
    Q_inf_modul=norm(Q_inf);

    % Inicialització de matrius
    A = zeros(N, N);
    b = zeros(N, 1);

    % Construir matriu A i vector b
    for i = 1:N
        b(i) = -dot(Q_inf, Tangent(i, :)); %producte escalar (pag. 40)
        for j = 1:N
            if j == i
                A(i,j) = -0.5; %autoinducció (pag. 43)
            else
                % Coordenades del punt de control i respecte al panell j (pag. 41)
                X_diff = X_c(i) - X(j);
                Z_diff = Z_c(i) - Z(j);
                X_cpanj = X_diff * calphaj(j) - Z_diff * salphaj(j);
                Z_cpanj = X_diff * salphaj(j) + Z_diff * calphaj(j);

                r1 = sqrt(X_cpanj^2 + Z_cpanj^2);
                r2 = sqrt((X_cpanj - l(j))^2 + Z_cpanj^2);

                theta1 = atan2(Z_cpanj, X_cpanj);
                theta2 = atan2(Z_cpanj, X_cpanj - l(j));

                % Velocitat induida (pag 41)
                w = (1/(4*pi)) * log(r2^2 / r1^2); 
                u = (theta2 - theta1) / (2*pi); 
                
                % Velocitat coordenades globals (pag 43)
                u_ij = u * calphaj(j) + w * salphaj(j);
                w_ij = -u * salphaj(j) + w * calphaj(j);
                A(i,j) = u_ij * calphaj(1,i) - w_ij *salphaj(1,i); %producte escalar

            end
        end
    end

    % Condició de Kutta al panell de sortida (aproximadament 1/2 del perfil)
    kutta_idx = round(N/2); %
    A(kutta_idx, :) = 0;
    A(kutta_idx, 1) = 1;
    A(kutta_idx, N) = 1;
    b(kutta_idx) = 0;

    % Resolució del sistema lineal
    gamma = A \ b;

end
