clc;
clear;
close all;

addpath('functions\')

% constants related to earth and geocoding
R_earth = 6378.137 * 10^3 ; % [m]

% Dataset main source: https://aerpaw.org/dataset/aerpaw-find-a-rover-afar-challenge-in-december-2023/
% this includes experiments from five teams, denoted by 288, 300, 301, 309,
% and 328
% For each team three seprate experiments were conducted for three
% locations of Transmitter, Tx (kept on a UGV)
% for 309, experiment for Tx location 1 was not conducted
% so, there were 14 experiments in total


exps = [328];
locations = [3];

% antenna pattern (measured in anechoic chamber)
load("./antenna_patterns/rad_pat.mat")
G_t_dB = SA3300;
G_r_dB = SA3300;

% The unknown azimuth orientation of Rx and Tx antennas were calculated finding opmital
% angle that minimized PL calculation error
unknown_angle = -50;
unknown_angle2 = 25;

% constants related to PL calculation
pl_expo = 2;
h_BS = 1.5;
P_Tx_dBm = 100 ; % [dBm]
c = 3.0e8;
f0 = 3.32e9; % what is the actual f0 here
lamda = c/f0;

% use one exp and one location at a time, to avoid unexpected errors,
% run 14 times for 14 experiments, it is easy
for iii=1:length(exps)
    for jjj=1:length(locations)
        exp_no = exps(iii);
        location = locations(jjj);
        if exp_no==309 && location==1
            continue
        end
        
        % UGV (Tx carrier) azimuth angle taken into account
        if location < 2 % 1
            ugv_angle = 60.7648 + unknown_angle2; % 60.7648 + unknown_angle2
        elseif location < 3 % 2
            ugv_angle = 60.7648 + unknown_angle2;
        else % 3
            ugv_angle = -70.4993 + unknown_angle2; % -70.4993 + unknown_angle2
        end
        
        % data source: https://aerpaw.org/dataset/aerpaw-find-a-rover-afar-challenge-in-december-2023/
        matdata = load("raw_data/results_useful_ori_exp_"+num2str(exp_no)+"_loc_"+num2str(location)+".mat");
        
        % rename variables
        lat_all = matdata.mY;
        lon_all = matdata.mX;
        h_all = matdata.mZ;
        power_all = matdata.mpower;
        vz_all = matdata.mVz;
        q_all = matdata.mquality;
        speed_all = matdata.mSpeed;
        ori_all = matdata.mYaw+unknown_angle; % 0 to 360 (UAV angle)
        roll_all = matdata.mRoll;
        pitch_all = matdata.mPitch;
        
        % UAV orientation angles (in degree)
        ori_all = ori_all - 360 .* (ori_all>=180) + 360 .* (ori_all<-180);
        roll_all = roll_all - 360 .* (roll_all>=180) + 360 .* (roll_all<-180);
        pitch_all = pitch_all - 360 .* (pitch_all>=180) + 360 .* (pitch_all<-180);
        
        height_avg = mode(round(h_all));
        
        % filter measurements during take off and lading
        valid_height = ((h_all>height_avg-5) & (h_all<height_avg+5));
        lat_all = lat_all(valid_height);
        lon_all = lon_all(valid_height);
        h_all = h_all(valid_height);
        power_all = power_all(valid_height);
        vz_all = vz_all(valid_height);
        q_all = q_all(valid_height);
        speed_all = speed_all(valid_height);
        ori_all = ori_all(valid_height);
        roll_all = roll_all(valid_height);
        pitch_all = pitch_all(valid_height);
        
        N_sam = length(power_all);
        
        % get Tx coordinate (based on UGV location number 1/2/3)
        if location ==1
        origin_y=35.72806709;
        origin_x=-78.69730398;
        end
        if location==2
        origin_y=35.72911779;
        origin_x=-78.69918128;
        end
        if location==3
        origin_y=35.72985129;
        origin_x=-78.69711002;
        end
        lat_BS = origin_y;
        lon_BS = origin_x;
        
        % calculate the UAV's elevation and azimuth angle, and 2D and 3D distance of UAV from Tx, 
        h_UAV_new = h_all;
        dist_2D_new = R_earth .* acos(sin(lat_BS * pi/180) .* sin(lat_all * pi/180) + cos(lat_BS * pi/180) .* cos(lat_all * pi/180) .* cos( (lon_all - lon_BS)  * pi/180 )) ;
        dist_3D = sqrt( (dist_2D_new).^2 + (h_UAV_new - h_BS).^2 ) ;
        elev = atan( ( h_UAV_new - h_BS )./dist_2D_new) ;
        for ii=1:length(lon_all)
            azim(ii,1) = great_circle_azimuth(lat_BS,lon_BS,lat_all(ii),lon_all(ii));
        end
        
        % calculate \delta (tilt angle towards source) from 3D orientation
        % angle and LoS direction
        pitch_towards_src = get_uav_body_reflection_asilomar(azim,elev,ori_all,pitch_all,roll_all,dist_2D_new,h_all);

        % PL calculation
        [RSRP_PL_free, RSRP_PL_two] = baseline_flipx_anechoic_test_asilomar(power_all, ori_all, roll_all, pitch_all,...
    azim, elev, dist_3D, dist_2D_new, h_UAV_new, h_BS, ugv_angle, G_t_dB, G_r_dB, theta2, phi2, pl_expo, P_Tx_dBm, lamda);
        % save the 3D coordinates of UAV, actual RSRP, calculated PL, and
        % tilt angles
        save("processed_data/results_useful_with_tworay_tilt_src_exp_"+num2str(exp_no)+"_loc_"+num2str(location)+".mat", 'power_all', 'ori_all', 'roll_all', 'pitch_all', 'azim', 'elev', 'dist_3D', 'dist_2D_new', 'RSRP_PL_free', 'RSRP_PL_two', 'lat_all', 'lon_all', 'h_all', 'speed_all',...
            'pitch_towards_src');
    end
end