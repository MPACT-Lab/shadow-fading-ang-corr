function rot_mat = Ry_pitch(beta)
    rot_mat = [cos(beta) 0 sin(beta); 0 1 0; -sin(beta) 0 cos(beta)];
end