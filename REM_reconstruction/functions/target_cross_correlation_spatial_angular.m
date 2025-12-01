function target_xcorr_all = target_cross_correlation_spatial_angular(lat_meas,lon_meas, h_meas, lat_tar, lon_tar, h_tar, corr_coef, var_sha, elev_meas,...
    elev_tar, tilt_meas, tilt_tar, tilt_sensitivity_exp, elev_sensitivity_exp, elev_sep, tilt_src_sep, k_val_tilt, k_val_elev)
    
    % input: 3D coordinates of the measurement samples and target samples, and correlation
    % profiling parameters (both spatial and angular)
    % output: measured correlation among measured and target samples as a matrix...
    
    num_tar = length(lat_tar);
    num_meas = length(lat_meas);
    R_earth = 6378.137 * 10^3 ; % [m]

    target_xcorr_all = ones(num_tar, num_meas);

    for i=1:num_meas
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
    end

    for i=1:num_tar
        % assuming elev_tar(i) is between elev_sep(1) and elev_sep(end)
        % (0 to 90 degree) or (-30 to 30) for delta theta
        [~, temp] = find(elev_tar(i) > elev_sep);
        if ~numel(temp)
            elev_indx1 = 1;
        else
            elev_indx1 = min(max(temp), length(elev_sep)-1);
        end

        [~, temp] = find(tilt_tar(i) > tilt_src_sep);
        if ~numel(temp)
            tilt_indx1 = 1;
        else
            tilt_indx1 = min(max(temp), length(tilt_src_sep)-1);
        end

        increase_expo_tilt = tilt_sensitivity_exp(elev_indx1, tilt_indx1, 1);
        decrease_expo_tilt = tilt_sensitivity_exp(elev_indx1, tilt_indx1, 2);

        increase_expo_elev = elev_sensitivity_exp(tilt_indx1, elev_indx1, 1);
        decrease_expo_elev = elev_sensitivity_exp(tilt_indx1, elev_indx1, 2);

        indx1 = (elev_indx1-1)*(length(tilt_src_sep)-1) + tilt_indx1;
        for j=1:num_meas
            dist_hor_ij = R_earth * acos(sin(lat_tar(i) * pi/180) * sin(lat_meas(j) * pi/180)...
                + cos(lat_tar(i) * pi/180) * cos(lat_meas(j) * pi/180) * cos(...
                (lon_tar(i) - lon_meas(j))  * pi/180 )) ;
            target_xcorr_all(i,j) = get_2D_corr_profile(dist_hor_ij,corr_coef) * (var_sha); 

            delta_tilt = abs(tilt_tar(i) - tilt_meas(j));
            delta_elev = abs(elev_tar(i) - elev_meas(j));
            
            increase_tilt = tilt_meas(j) > tilt_tar(i);
            increase_elev = elev_meas(j) > elev_tar(i);
            expo_tilt = increase_tilt*increase_expo_tilt + (1-increase_tilt)*decrease_expo_tilt;
            expo_elev = increase_elev*increase_expo_elev + (1-increase_elev)*decrease_expo_elev;

            tilt_corr = exp(-delta_tilt*expo_tilt*k_val_tilt);
            elev_corr = exp(-delta_elev*expo_elev*k_val_elev);
            target_xcorr_all(i,j) = target_xcorr_all(i,j) * tilt_corr * elev_corr; %asilomar_mat(indx1, indx2);
        end
    end

    function r=get_2D_corr_profile(dis_hor, corr_coef)
        % correlation_f
        r = corr_coef(1)*exp(-(corr_coef(2)*dis_hor)) + (1-corr_coef(1))*exp(-(corr_coef(3)*dis_hor)); %using data points hor 250 ver 50
    end
end