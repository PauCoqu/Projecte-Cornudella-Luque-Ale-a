function [X, Z] = llegir_punts_perfils(file_name)
        
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

end
