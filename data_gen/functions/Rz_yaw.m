function rot_mat = Rz_yaw(alpha)
    rot_mat = [cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
end