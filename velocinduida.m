function [V_ij] = velocinduida(Coords_centre, Coords_1, Coords_2, u_r)

% Per cada panell calculem la velocitat del vortex (V_ij = VinfA + V_AB - V_infB)
r1 = (Coords_centre - Coords_1)'; %(pag 9)
r2 = (Coords_centre - Coords_2)'; %(pag 9)

u_r1 = r1/norm(r1); %vector unitari u_r1
u_r2 = r2/norm(r2); %vector unitari u_r1

V_infA = 1/(4*pi)*(1-dot(u_r,u_r1))/(norm(cross(u_r,r1))).^2 *cross(u_r,r1); %(pag 10)
V_infB = 1/(4*pi)*(1-dot(u_r,u_r2))/(norm(cross(u_r,r2))).^2 *cross(u_r,r2); %(pag 10)
V_AB = (1/(4*pi))*(norm(r1)+norm(r2))/(norm(r1)*norm(r2)*(norm(r1)*norm(r2) + dot(r1,r2)))*cross(r1,r2);

V_ij = V_infA + V_AB - V_infB; %(pag 11)    



