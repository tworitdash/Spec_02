function P_rad = Prad(P_rad, U, dth, dph, th, ph)
  
        P_rad_i =  U .* sin(th) .*dth .* dph;
        for j = 1:size(U,3)
            for k = 1:180
                P_rad_f(k, j) = sum(P_rad_i(:, k, j));
            end
            P_rad(j) = sum(P_rad_f(:,j)); 
        end
end