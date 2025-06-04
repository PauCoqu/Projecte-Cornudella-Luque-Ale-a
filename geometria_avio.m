function [Coords_ala, Coords_centre_panell, c_ala, Coords_canard, Coords_centre_canard, c_canard] = geometria_avio (N,b,c_r,c_t,b_h,c_rh,c_th,l_h)

% Discretitzem l'ala i el canard en N panells per tal de poder aplicar la teoria de lifting line de Prandtl
Ala_coords_x = zeros(N+1, 1);                 
Ala_coords_y = linspace(-b/2, b/2, N+1)';
Ala_coords_z = zeros(N+1, 1);
Coords_ala = [Ala_coords_x, Ala_coords_y, Ala_coords_z];

% Punts de control
Coords_centre_panell = (Coords_ala(1:end-1, :) + Coords_ala(2:end, :)) / 2;

% Formula d'evolució de la corda considerant ala trapezoidal (varia en l'eix Y)
%corda = c_r - coord * pendent (y=0 ; c=c_r | y=b/2 ; c=c_h)
c_ala = c_r - (abs(Coords_centre_panell(:,2))*(c_r-c_t))/(b/2);

%Ara discretitzem el canard (CORREGIR!)
Canard_coords_x = zeros(N+1, 1);                 
Canard_coords_y = linspace(-b_h/2, b_h/2, N+1)';
Canard_coords_z = zeros(N+1, 1);
Coords_canard = [Canard_coords_x, Canard_coords_y, Canard_coords_z];

% Punts de control
Coords_centre_canard = (Coords_canard(1:end-1, :) + Coords_canard(2:end, :)) / 2;

% Formula evolució corda considerant ala trapezoidal (varia en l'eix Y)
%corda = c_r - coord * pendent (y=0 ; c=c_r | y=b/2 ; c=c_h)
c_canard = c_r - (abs(Coords_centre_canard(:,2))*(c_rh-c_th))/(b_h/2);

end