function meas_correlation_all = meas_correlation_spatial_angular(lat_meas,lon_meas, h_meas, corr_coef, var_sha,...
    elev_meas, tilt_meas, tilt_sensitivity_exp, elev_sensitivity_exp, elev_sep, tilt_src_sep, k_val_tilt, k_val_elev)
    
    % input: 3D coordinates of the measurement samples, and correlation
    % profiling parameters (both spatial and angular)
    % output: measured correlation among measured samples as a matrix

    num_meas = length(lat_meas);
    R_earth = 6378.137 * 10^3 ; % [m]

    meas_correlation_all = ones(num_meas, num_meas);

    for i=1:num_meas
        % assuming elev_tar(i) is between elev_sep(1) and elev_sep(end)
        % (0 to 90 degree) or (-30 to 30) for delta theta
        [~, temp] = find(elev_meas(i) > elev_sep);
        if ~numel(temp)
            elev_indx1 = 1;
        else
            elev_indx1 = min(max(temp), length(elev_sep)-1);
        end
        
        [~, temp] = find(tilt_meas(i) > tilt_src_sep);
        if ~numel(temp)
            tilt_indx1 = 1;
        else
            tilt_indx1 = min(max(temp), length(tilt_src_sep)-1);
        end
        increase_expo_tilt = tilt_sensitivity_exp(elev_indx1, tilt_indx1, 1);
        decrease_expo_tilt = tilt_sensitivity_exp(elev_indx1, tilt_indx1, 2);

        increase_expo_elev = elev_sensitivity_exp(tilt_indx1, elev_indx1, 1);
        decrease_expo_elev = elev_sensitivity_exp(tilt_indx1, elev_indx1, 2);

        for j=1:num_meas
            if j<i
                meas_correlation_all(i,j) = meas_correlation_all(j,i);
            end
            dist_ver_ij = 0; %abs(h_meas(i) - h_meas(j));
            dist_hor_ij = R_earth * acos(sin(lat_meas(i) * pi/180) * sin(lat_meas(j) * pi/180)...
                + cos(lat_meas(i) * pi/180) * cos(lat_meas(j) * pi/180) * cos(...
                (lon_meas(i) - lon_meas(j))  * pi/180 )) ;
            meas_correlation_all(i,j) = get_3D_corr_profile(dist_hor_ij,dist_ver_ij,corr_coef) * (var_sha);
            
            delta_tilt = abs(tilt_meas(i) - tilt_meas(j));
            delta_elev = abs(elev_meas(i) - elev_meas(j));
            
            increase_tilt = tilt_meas(j) > tilt_meas(i);
            increase_elev = elev_meas(j) > elev_meas(i);
            expo_tilt = increase_tilt*increase_expo_tilt + (1-increase_tilt)*decrease_expo_tilt;
            expo_elev = increase_elev*increase_expo_elev + (1-increase_elev)*decrease_expo_elev;

            tilt_corr = exp(-delta_tilt*expo_tilt*k_val_tilt);
            elev_corr = exp(-delta_elev*expo_elev*k_val_elev);
            meas_correlation_all(i,j) = meas_correlation_all(i,j) * tilt_corr * elev_corr; %asilomar_mat(indx1, indx2);
        end
    end

    function r=get_3D_corr_profile(dis_hor, dis_ver, corr_coef)
        % correlation_f
        r = corr_coef(1)*exp(-(corr_coef(2)*dis_hor)) + (1-corr_coef(1))*exp(-(corr_coef(3)*dis_hor)); %using data points hor 250 ver 50
    end
end