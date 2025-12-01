function rot_mat = Rx_roll(gamma)
    rot_mat = [1 0 0; 0 cos(gamma) -sin(gamma); 0 sin(gamma) cos(gamma)];
end

