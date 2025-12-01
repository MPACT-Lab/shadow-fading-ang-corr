function [RSRP_PL_free, RSRP_PL_two] = baseline_flipx_anechoic_test_asilomar(power_all, ori_all, roll_all, pitch_all,...
    azim, elev, dist_3D, dist_2D_new, h_UAV_new, h_BS, ugv_angle, G_t_dB, G_r_dB, theta2, phi2, pl_expo, P_Tx_dBm, lamda)
    
    % input: PL calculation inputs: distance, azimuth and elevation angles,
    % antenna gains
    % output: predicted RSRP with PL model: 1) FSPL, 2) Two-ray path loss

    N_sam = length(power_all); 
    G_r_max = 0 ; % [dBi]
    phi_os = 0; % 0 ~ 119 % 88

    % Transmit power
    unknown_offset = 0;
    P_Tx_2_dBm = P_Tx_dBm + unknown_offset ;
    P_Tx_2 = 10.^(P_Tx_2_dBm/10) ;
 
     
    % Elevation angle index
    [~,idx_theta] = min(abs(elev - theta2'*pi/180)') ;
    % Azimuth angle (of UAV w.r.t BS) index
    theta_tx = azim - ugv_angle;
    theta_tx = theta_tx + (theta_tx<-pi).*2*pi - (theta_tx>pi).*2*pi;
    [~,idx_phi] = min(abs(theta_tx - circshift(phi2, phi_os)*pi/180)') ;


    % Azimuth angle (of BS w.r.t UAV using UAV yaw angle) index
    azim_wrt_uav_ori = azim + ori_all*pi/180;
    azim_wrt_uav_ori = azim_wrt_uav_ori - (azim_wrt_uav_ori>=0).*pi + (azim_wrt_uav_ori<0).*pi;
    [~,idx_phi2] = min(abs(-azim_wrt_uav_ori  - circshift(phi2, phi_os)*pi/180)') ;

    % antenna gains at Tx and Rx antenna
    for ii=1:N_sam
        G_t(ii,1)= G_t_dB(idx_theta(ii), idx_phi(ii)) ; % [dBi]
    end
    for ii=1:N_sam
        G_r(ii,1)= G_r_dB(idx_theta(ii), idx_phi2(ii)) + G_r_max ; % [dBi]
    end
    
    G_t_l = 10.^(G_t/10) ;
    G_r_l = 10.^(G_r/10) ;
    
    RSRP_PL_free_lin = P_Tx_2 * (lamda/(4*pi))^2 * abs( sqrt(G_t_l) .* sqrt(G_r_l) ./ ( dist_3D ) ).^pl_expo ;
    RSRP_PL_free = 10*log10(RSRP_PL_free_lin) ;
    
    % Two-ray ground reflection pathloss model
    rho = 0.005 ;
    eps = 15 ;
    eps_0 = eps - j*60*rho*lamda ;

    th_r = atan( ( h_BS + h_UAV_new ) ./ (dist_2D_new) ) ; % from 0 to pi/2
    th_r_2 = 2*pi - th_r ;
    th_r_22 = - th_r ;
    
    r_1 = h_BS ./ sin(th_r) ;
    r_2 = h_UAV_new ./ sin(th_r) ;
    r_3 = dist_3D ;
    
    phi_del = 2*pi*( r_1 + r_2 - r_3)/lamda ;
    
    z_V = sqrt(eps_0 - cos(th_r).^2) / eps_0 ;
    R_V = (sin(th_r) - z_V) ./ (sin(th_r) + z_V) ;
    
    % Azimuth and elvation angle of the reflected ray
    [~,idx_theta_t_r] = min(abs(th_r_22 - theta2'*pi/180)') ; % transmitter ray downwards
    [~,idx_theta_r_r] = min(abs(th_r - theta2'*pi/180)') ; % receiver ray recives upward ray
    
    for ii=1:N_sam
        G_t_r_dB(ii,1)= G_t_dB(idx_theta_t_r(ii), idx_phi(ii)) ; % [dBi]
    end
    for ii=1:N_sam
        G_r_r_dB(ii,1)= G_r_dB(idx_theta_r_r(ii), idx_phi2(ii)) + G_r_max ; % [dBi]
    end
    
    G_t_r = 10.^(G_t_r_dB/10) ;
    G_r_r = 10.^(G_r_r_dB/10) ;
    
    RSRP_PL_two_lin = P_Tx_2 * (lamda/(4*pi))^2 * abs( sqrt(G_t_l) .* sqrt(G_r_l) ./ ( r_3 ) + R_V .* sqrt(G_t_r) .* sqrt(G_r_r) ./ ( r_1 + r_2) .* exp( -1i * phi_del ) ).^pl_expo ;
    RSRP_PL_two = 10*log10(RSRP_PL_two_lin);
end
