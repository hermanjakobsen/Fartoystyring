function [chi_d] = guidance(x1, y1, x2, y2, ship_x, ship_y)
    % path-tangential angle with respect to the North axis
    pi_p = atan2(y2-y1, x2-x1); 
    
    % cross-track error expressed in NED 
    y_e = -(ship_x-x1) * sin(pi_p) + (ship_y-y1) * cos(pi_p);   % crosstrackWpt.m
    
    % proportional gain (eq. 12.79)
    Kp = 1 / 1; 
   
    % desired course angle
    chi_d = pi_p - atan(Kp*y_e); 
end