function [X, Z] = llegir_punts_perfils(file_name)
    % Llegeix les coordenades X i Z del perfil NACA des d'un fitxer .txt
    % file_name: nom del fitxer que conté les dades de les coordenades
    
    % Obrim el fitxer
    fileID = fopen(file_name, 'r');
    if fileID == -1
        error('No s''ha pogut obrir el fitxer');
    end
    
    % Llegir les coordenades
    data = fscanf(fileID, '%f %f %f\n', [3, Inf]);  % Llegeix tres columnes de dades
    fclose(fileID);
    
    % Extraiem les coordenades X i Z (segona i tercera columna)
    X = data(2, :);  % Coordenades X
    Z = data(3, :);  % Coordenades Z

    % Visualitzar el perfil
    %figure;
    %plot(X, Z, '-o');
    %title('Perfil NACA 0012 sense rotar');
    %xlabel('X');
    %ylabel('Z');
    %grid on;
    %axis equal;
end
