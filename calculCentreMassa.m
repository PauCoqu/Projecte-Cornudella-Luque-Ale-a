function [CentreMassa] = calculCentreMassa(Coords_ala, Coords_canard, N, Lift_ala, Lift_can)

%Calcul CM; tant l'ala com el canard no tenen fletxa, per la qual cosa
%considerem el CL concentrat en un punt al llarg del fuselatge

centreL_ala = Coords_ala(N,1);
centreL_can=Coords_canard(N,1);

CentreMassa = (Lift_ala*centreL_ala+Lift_can*centreL_can)/(Lift_ala+Lift_can);
