function meas_correlation_all = meas_correlation_spatial(lat_meas,lon_meas, h_meas, corr_coef, var_sha)
    % input: 3D coordinates of the measurement samples, and correlation
    % profiling parameters
    % output: measured correlation among measured samples as a matrix
    num_meas = length(lat_meas);
    R_earth = 6378.137 * 10^3 ; % [m]

    meas_correlation_all = ones(num_meas, num_meas);

    for i=1:num_meas
        for j=1:num_meas
            dist_ver_ij = 0; %abs(h_meas(i) - h_meas(j));
            dist_hor_ij = R_earth * acos(sin(lat_meas(i) * pi/180) * sin(lat_meas(j) * pi/180)...
                + cos(lat_meas(i) * pi/180) * cos(lat_meas(j) * pi/180) * cos(...
                (lon_meas(i) - lon_meas(j))  * pi/180 )) ;
            if (i==j)
                dist_hor_ij = 0;
            end
            meas_correlation_all(i,j) = get_3D_corr_profile(dist_hor_ij,dist_ver_ij,corr_coef) * (var_sha); 
        end
    end

    function r=get_3D_corr_profile(dis_hor, dis_ver, corr_coef)
        % correlation_f
        r = corr_coef(1)*exp(-(corr_coef(2)*dis_hor)) + (1-corr_coef(1))*exp(-(corr_coef(3)*dis_hor)); %using data points hor 250 ver 50
    end
end