function P_rad = Prad_Assign(U, dth, dph, th, ph)

        P_rad = 0;
  
        P_rad_i =  U .* sin(th) .*dth .* dph;
        
        for j = 1:size(U, 2)
            for k = 1:size(U, 1)
                P_rad_f(k) = sum(P_rad_i(:, k));
            end
            P_rad = P_rad + sum(P_rad_f(:)); 
        end
end