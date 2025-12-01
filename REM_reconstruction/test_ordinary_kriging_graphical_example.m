clc;
clear;
close all;

addpath('functions\')

% constant used in geocoding
R_earth = 6378.137 * 10^3 ; % [m]

%% initialize and load correlation profiles
% variants are for testing different correlation models:
% corr_dis2d = coventional spatial correlation model (bi exponential
% function of horizontal distance)
% corr_elev_tilt = proposed correlation model with spatial and angular
% correlation
% corr_elev and corr_tilt = ablation studies for proposed correlation
% model, using only elevation angle or using tilt angle
variants_all = ["corr_dis2d", "corr_elev", "corr_tilt", "corr_elev_tilt"];

% using spatial correlation profile from experiment 309
% Using UAV data from a single trajectory, i.e, trajectory of experiment
% 309 so that the effects of 2D distance are less hampered by changed
% trajectory
load("../corr_profile/learned_profile/correlation_profile_asilomar_afar_309_loc1_2_3.mat")
% parameters of horizontal biexponential correlation fitting
corr_coef2 = [correlation_f.a, correlation_f.b, correlation_f.c];
% using angular correlation profile from experiment 309 and 328
load("../corr_profile/learned_profile/correlation_profile_asilomar_afar_309_328_loc1_2_3.mat")
% correlation profile obtained: tilt_sensitivity_exp,
% tilt_only_sensitivity_exp (ablation), elev_only_sensitivity_exp (ablation)

%% load test data from Kriging
exp_nos = [301]; % test experiment
exp_txt = num2str(exp_nos(1));
exp_no = exp_nos(1);
location = 2;

load("../data_gen/processed_data/results_useful_with_tworay_tilt_src_exp_"+num2str(exp_no)+"_loc_"+num2str(location)+".mat")
sha_all = power_all - RSRP_PL_two;
angle_all = pitch_towards_src*180/pi;
elev_all = elev; 
% equalize distribution (so that Kriging performance not hampered by too
% many nearby points as the known data)
dis_travelled = cumsum(speed_all);
equalized_dist = linspace(0, dis_travelled(end), length(dis_travelled));
lat_all = interp1(dis_travelled, lat_all, equalized_dist)';
lon_all = interp1(dis_travelled, lon_all, equalized_dist)';
sha_all = interp1(dis_travelled, sha_all, equalized_dist)';
angle_all = interp1(dis_travelled, angle_all, equalized_dist)';
elev_all = interp1(dis_travelled, elev_all, equalized_dist)';


%% Kriging interpolation setting

all_heihgts = [27];
added_term = var_sha*2;
num_repeat = 1;

num_known_samples = [150]; % per height
len_num_known_samples = 1;
radius_all = [200];
len_radius_all = 1;

radius_colors = ["r", "g", "b", "c", "m", "y", "k"];
total_points = length(lat_all); % 43k
train_points = num_known_samples(1); % 350
eligible_test_points = total_points - train_points;
test_point_selection = 1:40:eligible_test_points;
num_test_points = length(test_point_selection); % per height

rng(32322); % 623 32323 32322

%% Kriging interpolation

error_all_2d = ones(num_repeat, length(all_heihgts), length(num_known_samples), length(radius_all), num_test_points) * NaN;
sha_tar_2d = ones(num_repeat, length(all_heihgts), length(num_known_samples), length(radius_all), num_test_points) * NaN;
est_tar_2d = ones(num_repeat, length(all_heihgts), length(num_known_samples), length(radius_all), num_test_points) * NaN;
num_meas_used_2d = zeros(num_repeat, length(all_heihgts), length(num_known_samples), length(radius_all), num_test_points);
ang_tar_2d = ones(num_repeat, length(all_heihgts), length(num_known_samples), length(radius_all), num_test_points) * NaN;
elev_tar_2d = ones(num_repeat, length(all_heihgts), length(num_known_samples), length(radius_all), num_test_points) * NaN;

variance_all_2d = ones(num_repeat, length(all_heihgts), length(num_known_samples), length(radius_all), num_test_points) * NaN;

error_all_2d_ablation1 =  error_all_2d;
est_tar_2d_ablation1 =  est_tar_2d;
variance_all_2d_ablation1 =  variance_all_2d;

error_all_2d_ablation2 =  error_all_2d;
est_tar_2d_ablation2 =  est_tar_2d;
variance_all_2d_ablation2 =  variance_all_2d;

error_all_2d_conventional =  error_all_2d;
est_tar_2d_conventional =  est_tar_2d;
variance_all_2d_conventional =  variance_all_2d;


lat_calc = zeros(num_repeat, length(all_heihgts), num_test_points);
lon_calc = zeros(num_repeat, length(all_heihgts), num_test_points);
angle_calc = zeros(num_repeat, length(all_heihgts), num_test_points);
elev_calc = zeros(num_repeat, length(all_heihgts), num_test_points);
h_calc = zeros(num_repeat, length(all_heihgts), num_test_points);


plot_fig = 0;

for repeat_id = 1:num_repeat
    fprintf('iter %d|\n',repeat_id);
    for height_id = 1:length(all_heihgts) % all data in same height for a experiment in afar
        % target point set
        meas_ids_temp = randi([1 total_points],1,ceil(train_points*1.1));
        meas_ids_temp = unique(meas_ids_temp);
        meas_ids_temp_idx = randperm(length(meas_ids_temp));
        meas_ids = meas_ids_temp(meas_ids_temp_idx(1:train_points));

        lat_meas = lat_all(meas_ids) ;
        lon_meas = lon_all(meas_ids) ;
        h_meas = h_all(meas_ids) ;
        elev_meas = elev_all(meas_ids) ;
        angle_meas = angle_all(meas_ids) ;
        sha_meas = sha_all(meas_ids) ;
        
        tar_ids_probable_set = setdiff(1:total_points, meas_ids);
        target_ids = tar_ids_probable_set(test_point_selection);

        lat_tar = lat_all(target_ids) ;
        lon_tar = lon_all(target_ids) ;
        h_tar = h_all(target_ids) ;
        elev_tar = elev_all(target_ids) ;
        angle_tar = angle_all(target_ids) ;

        sha_tar = sha_all(target_ids) ;
        
        lat_calc(repeat_id, :, :) = lat_tar;
        lon_calc(repeat_id, :, :) = lon_tar;
        elev_calc(repeat_id, :, :) = elev_tar;
        angle_calc(repeat_id, :, :) = angle_tar;
        h_calc(repeat_id, :, :) = h_tar;
    
        
    
        for meas_i = 1:len_num_known_samples % length(num_known_samples)
            % meas point set
            num_meas = num_known_samples(meas_i);
            

            %dist_ver_tar = abs(h_tar - h_meas'); %(num_tar, num_meas)
            dist_hor_tar = R_earth .* acos(sin(lat_tar * pi/180) .* sin(lat_meas' * pi/180) + cos(lat_tar * pi/180) .* cos(lat_meas' * pi/180) .* cos( (lon_tar - lon_meas')  * pi/180 )) ;
            %dist_target = sqrt(dist_ver_tar.^2 + dist_hor_tar.^2);

            save_data = 1;
            for variant = variants_all
                if variant == "corr_elev_tilt"
                     k_val_tilt = 1;
                     k_val_elev = 1;
                     meas_correlation_all = meas_correlation_spatial_angular(lat_meas,lon_meas, h_meas, corr_coef2, var_sha, elev_meas,...
    angle_meas, tilt_sensitivity_exp, elev_sensitivity_exp, elev_sep, tilt_ang_sep, k_val_tilt, k_val_elev); %(num_meas, num_meas)
                    target_cross_correlation_all = target_cross_correlation_spatial_angular(lat_meas,lon_meas, h_meas, lat_tar, lon_tar, h_tar, corr_coef2, var_sha, elev_meas,...
    elev_tar, angle_meas, angle_tar, tilt_sensitivity_exp, elev_sensitivity_exp, elev_sep, tilt_ang_sep, k_val_tilt, k_val_elev); %(num_tar, num_meas)
                    sha_tar_kriging = sha_tar;
                    sha_meas_kriging = sha_meas;
                elseif variant == "corr_elev"
                     k_val_tilt = 0;
                     k_val_elev = 1;
                     % pass the elevation sensitivity irrespective of tilt
                     %that is elev_only_sensitivity_exp with size (1, len(elev_sep)-1)
                     % we need to make it (len(tilt_sep)-1,
                     % len(elev_sep)-1) to use within
                     % meas_correlation_spatial_angular function
                     elev_only_sensitivity_exp_compatible = ones(size(elev_sensitivity_exp)) .* elev_only_sensitivity_exp;
                     meas_correlation_all = meas_correlation_spatial_angular(lat_meas,lon_meas, h_meas, corr_coef2, var_sha, elev_meas,...
    angle_meas, tilt_sensitivity_exp, elev_only_sensitivity_exp_compatible, elev_sep, tilt_ang_sep, k_val_tilt, k_val_elev); %(num_meas, num_meas)
                    target_cross_correlation_all = target_cross_correlation_spatial_angular(lat_meas,lon_meas, h_meas, lat_tar, lon_tar, h_tar, corr_coef2, var_sha, elev_meas,...
    elev_tar, angle_meas, angle_tar, tilt_sensitivity_exp, elev_only_sensitivity_exp_compatible, elev_sep, tilt_ang_sep, k_val_tilt, k_val_elev); %(num_tar, num_meas)
                    sha_tar_kriging = sha_tar;
                    sha_meas_kriging = sha_meas;
                elseif variant == "corr_tilt"
                     k_val_tilt = 1;
                     k_val_elev = 0;
                     % pass the tilt sensitivity irrespective of elevation
                     % that is tilt_only_sensitivity_exp with size (1, len(tilt_sep)-1)
                     % we need to make it (len(elev_sep)-1, len(tilt_sep)-1
                     % ) to use within
                     % meas_correlation_spatial_angular function
                     tilt_only_sensitivity_exp_compatible = ones(size(tilt_sensitivity_exp)) .* tilt_only_sensitivity_exp;
                     meas_correlation_all = meas_correlation_spatial_angular(lat_meas,lon_meas, h_meas, corr_coef2, var_sha, elev_meas,...
    angle_meas, tilt_only_sensitivity_exp_compatible, elev_sensitivity_exp, elev_sep, tilt_ang_sep, k_val_tilt, k_val_elev); %(num_meas, num_meas)
                    target_cross_correlation_all = target_cross_correlation_spatial_angular(lat_meas,lon_meas, h_meas, lat_tar, lon_tar, h_tar, corr_coef2, var_sha, elev_meas,...
    elev_tar, angle_meas, angle_tar, tilt_only_sensitivity_exp_compatible, elev_sensitivity_exp, elev_sep, tilt_ang_sep, k_val_tilt, k_val_elev); %(num_tar, num_meas)
                    sha_tar_kriging = sha_tar;
                    sha_meas_kriging = sha_meas;
                else % corr_dis2d
                    meas_correlation_all = meas_correlation_spatial(lat_meas,lon_meas, h_meas, corr_coef2, var_sha); %(num_meas, num_meas)
                    target_cross_correlation_all = target_cross_correlation_spatial(lat_meas,lon_meas, h_meas, lat_tar, lon_tar, h_tar, corr_coef2, var_sha); %(num_tar, num_meas)
                    sha_tar_kriging = sha_tar;
                    sha_meas_kriging = sha_meas;
                end
       
             
                for radi_i = 1:len_radius_all % length(radius_all) parfor
                    index_valid_radi_2d = (dist_hor_tar < radius_all(radi_i));% & (B>0); % & (dist_ver_tar<10) ;
        

                    index_valid_radi = index_valid_radi_2d;
                    
                   [est_tar_kriging, est_variance, num_points] = my_oridinary_kriging(meas_correlation_all,...
        target_cross_correlation_all, index_valid_radi, sha_meas_kriging, var_sha);


                   if (save_data)
                       sha_tar_2d(repeat_id, :, meas_i, radi_i, :) = sha_tar;
                       ang_tar_2d(repeat_id, :, meas_i, radi_i, :) = angle_tar;
                       elev_tar_2d(repeat_id, :, meas_i, radi_i, :) = elev_tar;
                       num_meas_used_2d(repeat_id, :, meas_i, radi_i, :) = num_points;
                   end

                    if variant == "corr_elev_tilt"
                       est_tar = est_tar_kriging;
                       error_all_2d(repeat_id, :, meas_i, radi_i, :) = est_tar - sha_tar;
                       est_tar_2d(repeat_id, :, meas_i, radi_i, :) = est_tar;
                       variance_all_2d(repeat_id, :, meas_i, radi_i, :) = est_variance;
                    elseif variant == "corr_elev"
                       est_tar = est_tar_kriging;
                       error_all_2d_ablation1(repeat_id, :, meas_i, radi_i, :) = est_tar - sha_tar;
                       est_tar_2d_ablation1(repeat_id, :, meas_i, radi_i, :) = est_tar;
                       variance_all_2d_ablation1(repeat_id, :, meas_i, radi_i, :) = est_variance;
                    elseif variant == "corr_tilt"
                       est_tar = est_tar_kriging;
                       error_all_2d_ablation2(repeat_id, :, meas_i, radi_i, :) = est_tar - sha_tar;
                       est_tar_2d_ablation2(repeat_id, :, meas_i, radi_i, :) = est_tar;
                       variance_all_2d_ablation2(repeat_id, :, meas_i, radi_i, :) = est_variance;
                    else % corr_dis2d
                       est_tar = est_tar_kriging;
                       error_all_2d_conventional(repeat_id, :, meas_i, radi_i, :) = est_tar - sha_tar;
                       est_tar_2d_conventional(repeat_id, :, meas_i, radi_i, :) = est_tar;
                       variance_all_2d_conventional(repeat_id, :, meas_i, radi_i, :) = est_variance;
                    end
                end
            end
        end
    end
end

%% Calculate RMSE of predictions

true_value = sha_tar;
pred_conventional = est_tar_2d_conventional(:);
pred_proposed = est_tar_2d(:);
pred_ablation1_elev = est_tar_2d_ablation1(:);
pred_ablation2_tilt = est_tar_2d_ablation2(:);

% RMSE is done for 100 samples
err1 = pred_conventional-true_value;
rmse_set = [];
for i=1:100:length(err1)
    rmse_set = [rmse_set, sqrt(mean(err1.^2, 'omitnan'))];
end
rmse_conventional = median(rmse_set);

err1 = pred_proposed-true_value;
rmse_set = [];
for i=1:100:length(err1)
    rmse_set = [rmse_set, sqrt(mean(err1.^2, 'omitnan'))];
end
rmse_proposed = median(rmse_set);

err1 = pred_ablation1_elev-true_value;
rmse_set = [];
for i=1:100:length(err1)
    rmse_set = [rmse_set, sqrt(mean(err1.^2, 'omitnan'))];
end
rmse_ablation1 = median(rmse_set);

err1 = pred_ablation2_tilt-true_value;
rmse_set = [];
for i=1:100:length(err1)
    rmse_set = [rmse_set, sqrt(mean(err1.^2, 'omitnan'))];
end
rmse_ablation2 = median(rmse_set);

%% SF prediction maps illustration

close all;
offset = [-78.698 35.727] ;
scaler = 111139 ;
%figure(repeat_id)
%FigH = figure('Position', get(0, 'Screensize'), 'visible','off');
mX_range = [-78.6998, -78.6955];
mY_range = [35.7270, 35.7295];
new_offset = [-78.6998,35.7270];
color_range = [-40, 40];
font_size = 16;
lat_eNBs = [35.72911779];
lon_eNBs = [-78.69918128];

lon_plot = lon_tar;
lat_plot = lat_tar;

% all true values
rsrp_plot = true_value;
FigH = figure(1);
hold on
scatter((lon_plot-new_offset(1))*scaler,(lat_plot-new_offset(2))*scaler,[],rsrp_plot,'filled')
scatter((lon_eNBs(1)-new_offset(1))*scaler,(lat_eNBs(1)-new_offset(2))*scaler,'r^',"filled")
colormap jet
ylim([0, (mY_range(2)-new_offset(2))*scaler])
xlim([0, (mX_range(2)-new_offset(1))*scaler])
clim(color_range);
h = colorbar;
title(h, "SF [dB]")
xlabel('X [m]')
ylabel('Y [m]')
fontsize(FigH,font_size,"points")
text(50,220,sprintf('Tx'), 'FontSize', 16);
grid on
box on

% Tuning samples (150 values used for Kriging)
rsrp_plot = sha_meas;
FigH = figure(10);
hold on
scatter((lon_meas-new_offset(1))*scaler,(lat_meas-new_offset(2))*scaler,[],rsrp_plot,'filled')
scatter((lon_eNBs(1)-new_offset(1))*scaler,(lat_eNBs(1)-new_offset(2))*scaler,'r^',"filled")
colormap jet
ylim([0, (mY_range(2)-new_offset(2))*scaler])
xlim([0, (mX_range(2)-new_offset(1))*scaler])
clim(color_range);
h = colorbar;
title(h, "SF [dB]")
xlabel('X [m]')
ylabel('Y [m]')
fontsize(FigH,font_size,"points")
text(50,220,sprintf('Tx'), 'FontSize', 16);
grid on
box on

% R(dis_2d)
rsrp_plot = pred_conventional;
FigH = figure(2);
hold on
scatter((lon_plot-new_offset(1))*scaler,(lat_plot-new_offset(2))*scaler,[],rsrp_plot,'filled')
scatter((lon_eNBs(1)-new_offset(1))*scaler,(lat_eNBs(1)-new_offset(2))*scaler,'r^',"filled")
colormap jet
ylim([0, (mY_range(2)-new_offset(2))*scaler])
xlim([0, (mX_range(2)-new_offset(1))*scaler])
clim(color_range);
h = colorbar;
title(h, "SF [dB]")
xlabel('X [m]')
ylabel('Y [m]')
fontsize(FigH,font_size,"points")
%text(160,40,sprintf('Median of RMSE = %.1f dB', rmse_conventional), 'FontSize', 16);
text(50,220,sprintf('Tx'), 'FontSize', 16);
grid on
box on

% R(dis_2d, theta, delta)
rsrp_plot = pred_proposed;
FigH = figure(3);
hold on
scatter((lon_plot-new_offset(1))*scaler,(lat_plot-new_offset(2))*scaler,[],rsrp_plot,'filled')
scatter((lon_eNBs(1)-new_offset(1))*scaler,(lat_eNBs(1)-new_offset(2))*scaler,'r^',"filled")
colormap jet
ylim([0, (mY_range(2)-new_offset(2))*scaler])
xlim([0, (mX_range(2)-new_offset(1))*scaler])
clim(color_range);
h = colorbar;
title(h, "SF [dB]")
xlabel('X [m]')
ylabel('Y [m]')
fontsize(FigH,font_size,"points")
%text(160,40,sprintf('Median of RMSE = %.1f dB', rmse_proposed), 'FontSize', 16);
text(50,220,sprintf('Tx'), 'FontSize', 16);
grid on
box on

% R(dis_2d, theta)
rsrp_plot = pred_ablation1_elev;
FigH = figure(5);
hold on
scatter((lon_plot-new_offset(1))*scaler,(lat_plot-new_offset(2))*scaler,[],rsrp_plot,'filled')
scatter((lon_eNBs(1)-new_offset(1))*scaler,(lat_eNBs(1)-new_offset(2))*scaler,'r^',"filled")
colormap jet
ylim([0, (mY_range(2)-new_offset(2))*scaler])
xlim([0, (mX_range(2)-new_offset(1))*scaler])
clim(color_range);
h = colorbar;
title(h, "SF [dB]")
xlabel('X [m]')
ylabel('Y [m]')
fontsize(FigH,font_size,"points")
%text(160,40,sprintf('Median of RMSE = %.1f dB', rmse_ablation1), 'FontSize', 16);
text(50,220,sprintf('Tx'), 'FontSize', 16);
grid on
box on

% R(dis_2d, delta)
rsrp_plot = pred_ablation2_tilt;
FigH = figure(6);
hold on
scatter((lon_plot-new_offset(1))*scaler,(lat_plot-new_offset(2))*scaler,[],rsrp_plot,'filled')
scatter((lon_eNBs(1)-new_offset(1))*scaler,(lat_eNBs(1)-new_offset(2))*scaler,'r^',"filled")
colormap jet
ylim([0, (mY_range(2)-new_offset(2))*scaler])
xlim([0, (mX_range(2)-new_offset(1))*scaler])
clim(color_range);
h = colorbar;
title(h, "SF [dB]")
xlabel('X [m]')
ylabel('Y [m]')
fontsize(FigH,font_size,"points")
%text(160,40,sprintf('Median of RMSE = %.1f dB', rmse_ablation2), 'FontSize', 16);
text(50,220,sprintf('Tx'), 'FontSize', 16);
grid on
box on