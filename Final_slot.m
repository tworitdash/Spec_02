 s = slot('Length', 0.014, 'Width', 0.001, 'GroundPlaneLength', 0.1, 'GroundPlaneWidth', 0.1);
 show(s);
 
impedance(s,10.7e9, 0, 1:1:360);

% s = dipole('Length', 0.014, 'Width', 0.002);
% show(s);
% 
% impedance(s, 5e9:0.01e9:15e9);