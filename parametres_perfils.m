function [Cl_alpha_0012, Cl_alpha_22112, Cl_0_0012, Cl_0_22112] = parametres_perfils ()

% parametres Cl_alpha, Cl_0 i Cm_1/4 de cada peerfil (dee moment de Google per+o s'hauria de treure del nostre treball apartat 1)
alpha = linspace(0,10,2);
alpha_rad = deg2rad(alpha_deg);


%Parametres NACA 0012 (Canard) --> Perfil simetric (Cl_0 = 0;)
Cl_0012 = [0,0,0,0,0]; %introduir Cl
Cm_1I4_0012 = [0,0,0,0,0];
Cl_0_0012 = Cl_0012(1,1);

%Parametres NACA 22112 (Ala)  
Cl_22112 = [0,0,0,0,0];
Cm_1I4_22112 = [0,0,0,0,0];
Cl_0_22112 = Cl_0012(1,1);


%Extreiem la pendent Cl_alpha

Cl_alpha_0012 = %....
Cl_alpha_22112 = %....

end

