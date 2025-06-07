function [A, b, li, X_c, Z_c, Normali] = metode_panells_ap4(Xi, Zi, Xj, Zj, Ni, Nj, Q_inf)

    lj = zeros(Nj, 1);  
    calphaj = zeros (Nj,1);
    salphaj = zeros (Nj,1);
    Tangentj = zeros(Nj, 2);
    Normalj = zeros(Nj, 2); 

    li = zeros(Ni, 1);
    calphai = zeros (Ni,1);
    salphai = zeros (Ni,1);
    Tangenti = zeros(Ni, 2);
    Normali = zeros(Ni, 2); 
    
    % Calcular longituds, tangencials i normals dels panells j
    for j = 1:Nj
        dx = Xj(j+1) - Xj(j);
        dz = Zj(j+1) - Zj(j);
        lj(j) = sqrt(dx^2 + dz^2);
        calphaj(j) = dx / lj(j);         % cos(α) del panell
        salphaj(j) = -dz / lj(j);        % -sin(α) del panell
        Tangentj(j, :) = [calphaj(j), -salphaj(j)];
        Normalj(j, :) = [salphaj(j),  calphaj(j)];
    end

        % Calcular longituds, tangencials i normals dels panells i
    for i = 1:Ni
        dx = Xi(i+1) - Xi(i);
        dz = Zi(i+1) - Zi(i);
        li(i) = sqrt(dx^2 + dz^2);
        calphai(i) = dx / li(i);         % cos(α) del panell
        salphai(i) = -dz / li(i);        % -sin(α) del panell
        Tangenti(i, :) = [calphai(i), -salphai(i)];
        Normali(i, :) = [salphai(i),  calphai(i)];
    end

    % Punts de control per panell
    X_c = (Xi(1:end-1) + Xi(2:end)) / 2;
    Z_c = (Zi(1:end-1) + Zi(2:end)) / 2;

    
    % Càlcul de gamma

    % Inicialització de matrius
    A = zeros(Ni, Nj);
    b = zeros(Ni, 1);

    % Construir matriu A i vector b
    for i = 1:Ni
        b(i) = -dot(Q_inf, Tangenti(i, :)); %producte escalar (pag. 40)
        for j = 1:Nj
            if abs(Xi(1,1) - Xj(1,1)) < eps && j == i % Condició específica pels casos A11 i A22
                A(i,j) = -0.5; %autoinducció (pag. 43)
            else
                % Coordenades del punt de control i respecte al panell j (pag. 41)
                X_diff = X_c(i) - Xj(j);
                Z_diff = Z_c(i) - Zj(j);
                X_cpanj = X_diff * calphaj(j) - Z_diff * salphaj(j);
                Z_cpanj = X_diff * salphaj(j) + Z_diff * calphaj(j);
                
                r1 = sqrt(X_cpanj^2 + Z_cpanj^2);
                r2 = sqrt((X_cpanj - lj(j))^2 + Z_cpanj^2);
                
                theta1 = atan2(Z_cpanj, X_cpanj);
                theta2 = atan2(Z_cpanj, X_cpanj - lj(j));

                % Velocitat induida (pag. 41)
                w = (1/(4*pi)) * log(r2^2 / r1^2); 
                u = (theta2 - theta1) / (2*pi); 
                
                % Velocitat coordenades globals (pag. 43)
                u_ij = u * calphaj(j) + w * salphaj(j);
                w_ij = -u * salphaj(j) + w * calphaj(j);
                A(i,j) = u_ij * calphai(i,1) - w_ij *salphai(i,1); %producte escalar

            end
        end
    end


    % Condició de Kutta (aproximadament 1/2 de la corda inferior)
    kutta_idx = round(Ni/4); % 

    if abs(Xi(1,1) - Xj(1,1)) < eps % Condició específica pels casos A11 i A22
        A(kutta_idx, :) = 0;
        A(kutta_idx, 1) = 1;
        A(kutta_idx, Nj) = 1;
        b(kutta_idx) = 0;
    else
        % Condició específica pels casos A12 i A21
        A(kutta_idx, :) = 0;
        b(kutta_idx) = 0;
    end
   
end