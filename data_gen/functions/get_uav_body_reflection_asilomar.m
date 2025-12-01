function pitch_towards_src = get_uav_body_reflection_asilomar(azim,elev,all_ori,all_pitch,all_roll,dis2d,hs)
    % input: azimuth and elevation angle of the UAV with respect to
    % Transmitter (Tx), the UAV's Euler orientation angles, its distance from
    % Tx, UAV altitude
    % output: UAV tilt angle towards soruce (\delta) = difference between
    % elevation angles of the Tx from the viewpoint of the UAV, before and
    % after incorporating 3D orientation of the UAV (roll, pitch, yaw)

    azim = azim - pi .* (azim>0) + pi.* (azim<=0); % -pi to pi
    elev = - elev; 
    N = length(azim);
    dis3d = dis2d.^2 + hs.^2;
    pitch_towards_src = zeros(size(azim)); % pitch_towards_src>0 means rcv elevation angle will increase (UAV is leaning up to the src)

    for ii = 1:N
        % calculate the inverse rotation matrix
        alpha = all_ori(ii)*pi/180;% yaw
        gamma = all_roll(ii)*pi/180; % roll
        beta = all_pitch(ii)*pi/180; % pitch
        % 3D rotation matrix (book reference: J. J. Craig, Introduction to Robotics: Mechanics and Control, 4th ed. Pearson Education International, 2018.
        R3 = Rx_roll(pi)*Ry_pitch(-gamma)*Rx_roll(beta)*Rz_yaw(-alpha);

        [input_x,input_y,input_z] = sph2cart(azim(ii),elev(ii),dis3d(ii)); 
        output_cartesian = R3 * [input_x; input_y; input_z];

        [rcv_azimuth,rcv_elevation,r_r] = cart2sph(output_cartesian(1),output_cartesian(2),output_cartesian(3));

        % calculate tilt towards source
        pitch_towards_src(ii) = rcv_elevation + elev(ii); % at zero tilt angle rcv_elevation would be negative of elev(ii)
        % so this is effectively calculating differece between the
        % elevation angles
    end


end
