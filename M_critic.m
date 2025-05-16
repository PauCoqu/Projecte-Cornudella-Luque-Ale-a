function [Cp_0, Cp_kt, Cp_star, Cp_crit, M_crit] = M_critic(Cp, N, gamma_aire)

M_inf = linspace(0.1, 1, N); %Discretització del Mach per realitzar el plot
Cp_0 = min(Cp); %Cp mínim

%Càlcul del Cp_kt i Cp_crit i veure on es creuen per determinar M_crit
Cp_kt = Cp_0 ./ (sqrt(1 - M_inf.^2) + (Cp_0 / 2) .* ((M_inf.^2) ./ (1 + sqrt(1 - M_inf.^2))));
% sqrt_term = sqrt(1 - M_inf.^2);
% Cp_kt = Cp_0 ./ (sqrt_term + (Cp_0 / 2) .* (M_inf.^2) ./ sqrt_term .* (1 + (gamma_aire - 1)/2 .* M_inf.^2));
% 


    Cp_star = (2 ./ (gamma_aire * M_inf.^2)) .* ( ((2 + (gamma_aire - 1) * M_inf.^2) ./ (1 + gamma_aire)).^(gamma_aire / (gamma_aire - 1)) - 1 );   

%Buscar quan es creuen les gràfiques:
%--------------
    delta = Cp_kt - Cp_star;

    % Buscar el primer canvi de signe
    idx = find(delta(1:end-1) .* delta(2:end) < 0, 1, 'first');

    if isempty(idx)
        M_crit = NaN;
        Cp_crit = NaN;
        disp('No hi ha intersecció entre les corbes.');
    else
        % Interpolació per trobar la intersecció precisa
        M1 = M_inf(idx);
        M2 = M_inf(idx+1);
        d1 = delta(idx);
        d2 = delta(idx+1);

        M_crit = interp1([d1, d2], [M1, M2], 0); % M_inf per a delta = 0
        Cp_crit = interp1(M_inf, Cp_kt, M_crit); % Cp corresponent
    end
%--------------

figure;
    plot(M_inf, Cp_kt, 'b--', 'LineWidth', 1); hold on;
    plot(M_inf, Cp_star, 'r-', 'LineWidth', 1);
    xlim([0.2 0.8]);
    xlabel('M_{\infty}');
    ylabel('C_p');
    title('Comparació entre C_p Karman-Tsien i C_p^*');
    legend('C_p (Karman-Tsien)', 'C_p^* (Crític)', 'Location', 'best');
    grid on;
    set(gca, 'YDir','reverse') % Per seguir convenció aerodinàmica
end


