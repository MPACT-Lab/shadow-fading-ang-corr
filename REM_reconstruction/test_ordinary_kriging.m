clc;
clear;
close all;

addpath('functions\')

% constant used in geocoding
R_earth = 6378.137 * 10^3 ; % [m]

out_folder = "ok_output";

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

%% load test data for Kriging
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

all_heights = [27];
added_term = var_sha*2;
num_repeat = 10;
num_known_samples = [50, 150, 250, 350, 450];
len_num_known_samples = 5;
radius_all = [ 70, 200];
len_radius_all = 2;

radius_colors = ["r", "g", "b", "c", "m", "y", "k"];
num_test_points = 100;

%% Kriging interpolation

error_all_2d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;
sha_tar_2d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;
est_tar_2d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;
num_meas_used_2d = zeros(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points);
ang_tar_2d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;
elev_tar_2d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;

variance_all_2d = ones(num_repeat, length(all_heights), length(num_known_samples), length(radius_all), num_test_points) * NaN;

error_all_2d_ablation1 =  error_all_2d;
est_tar_2d_ablation1 =  est_tar_2d;
variance_all_2d_ablation1 =  variance_all_2d;

error_all_2d_ablation2 =  error_all_2d;
est_tar_2d_ablation2 =  est_tar_2d;
variance_all_2d_ablation2 =  variance_all_2d;

error_all_2d_conventional =  error_all_2d;
est_tar_2d_conventional =  est_tar_2d;
variance_all_2d_conventional =  variance_all_2d;

lat_calc = zeros(num_repeat, length(all_heights), num_test_points);
lon_calc = zeros(num_repeat, length(all_heights), num_test_points);
angle_calc = zeros(num_repeat, length(all_heights), num_test_points);
elev_calc = zeros(num_repeat, length(all_heights), num_test_points);
h_calc = zeros(num_repeat, length(all_heights), num_test_points);


for repeat_id = 1:num_repeat
    fprintf('iter %d|\n',repeat_id);
    for height_id = 1:length(all_heights) % all data in same height for a experiment in afar
        % target point set
        tot_points = length(sha_all);
        target_ids = randi([1 tot_points],1,num_test_points);
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
    
        meas_ids_probable_set = setdiff(1:tot_points, target_ids);
    
        for meas_i = 1:len_num_known_samples % length(num_known_samples)
            % meas point set
            num_meas = num_known_samples(meas_i);
            meas_ids_temp = unique(randi([1, tot_points-num_test_points],1,ceil(num_meas*2)));
            randomize_idx = randperm(length(meas_ids_temp));
            meas_ids = meas_ids_probable_set(meas_ids_temp(randomize_idx(1:num_meas)));
            
            lat_meas = lat_all(meas_ids) ;
            lon_meas = lon_all(meas_ids) ;
            h_meas = h_all(meas_ids) ;
            elev_meas = elev_all(meas_ids) ;
            angle_meas = angle_all(meas_ids) ;
            sha_meas = sha_all(meas_ids) ;

            dist_hor_tar = R_earth .* acos(sin(lat_tar * pi/180) .* sin(lat_meas' * pi/180) + cos(lat_tar * pi/180) .* cos(lat_meas' * pi/180) .* cos( (lon_tar - lon_meas')  * pi/180 )) ;

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
                     %that is tilt_only_sensitivity_exp with size (1, len(tilt_sep)-1)
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
                       est_tar = est_tar_kriging; % + mu_test;
                       error_all_2d_ablation1(repeat_id, :, meas_i, radi_i, :) = est_tar - sha_tar;
                       est_tar_2d_ablation1(repeat_id, :, meas_i, radi_i, :) = est_tar;
                       variance_all_2d_ablation1(repeat_id, :, meas_i, radi_i, :) = est_variance;
                    elseif variant == "corr_tilt"
                       est_tar = est_tar_kriging; % + mu_test;
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
cd(out_folder)
eval("save res_kriging_ordinary_afar_"+num2str(exp_no)+"_"+num2str(location)+"_asilomar.mat error_all_2d error_all_2d_ablation1 error_all_2d_ablation2 error_all_2d_conventional num_meas_used_2d lat_calc lon_calc h_calc elev_calc angle_calc radius_all num_known_samples variance_all_2d variance_all_2d_ablation1 variance_all_2d_ablation2 variance_all_2d_conventional all_heights est_tar_2d est_tar_2d_ablation1 est_tar_2d_ablation2 est_tar_2d_conventional sha_tar_2d ang_tar_2d elev_tar_2d")
cd("../")

%% Load and Plot RMSE results

% load ordinary kriging resuts
load(sprintf("ok_output/res_kriging_ordinary_afar_301_2_asilomar.mat"))

%% CDF plots

radius_colors = ["r", "g", "b", "c", "m", "y", "k" ];
method_colors_1 = "k";
method_colors_2 = "m";
method_colors_3 = [1,0.5,0];
method_colors_4 = "b";
radius_plot_type = ["-^", "-s", "-v", "-o", "-<", "-v", "->"];
markertypes = ['^','s','v','o','<'];

font_size = 16;

data1_legend = 'Conventional, \boldmath$R\left(d_{\mathrm{2D}}\right)$';
data2_legend = 'Proposed, \boldmath$R\left(d_{\mathrm{2D}},\theta,\delta\right)$'; % change it
data3_legend = 'Ablation 1, \boldmath$R\left(d_{\mathrm{2D}},\theta\right)$';
data4_legend = 'Ablation 2, \boldmath$R\left(d_{\mathrm{2D}},\delta\right)$';

all_legends = [data1_legend+"", data2_legend+"", data3_legend+"", data4_legend+""];

FigH = figure(28);
hold on
set(FigH, 'Position', ceil([100,100,800,650]), 'visible','on');
h = cdfplot(sqrt(mean(error_all_2d_conventional.^2,5)));
h.LineWidth = 2; h.LineStyle = '--'; h.Color = method_colors_1; set(h, 'HandleVisibility', 'off');
x_h = get(h, 'XData'); y_h = get(h, 'YData');
h = scatter(x_h(10:4000:end),y_h(10:4000:end),50, method_colors_1, 'Marker', markertypes(1), 'LineWidth', 2);
set(h, 'HandleVisibility', 'off');
% dummy plot for legend
plot([25,26],[0.5,0.75],"--"+markertypes(1),'LineWidth',2,'Color',method_colors_1)

h = cdfplot(sqrt(mean(error_all_2d.^2,5)));
h.LineWidth = 2; h.LineStyle = '-'; h.Color = method_colors_2; set(h, 'HandleVisibility', 'off');
x_h = get(h, 'XData'); y_h = get(h, 'YData');
h = scatter(x_h(10:4000:end),y_h(10:4000:end),50, method_colors_2, 'Marker', markertypes(2), 'LineWidth', 2);
set(h, 'HandleVisibility', 'off');
plot([25,26],[0.5,0.75],"-"+markertypes(2),'LineWidth',2,'Color',method_colors_2)

h = cdfplot(sqrt(mean(error_all_2d_ablation1.^2,5)));
h.LineWidth = 2; h.LineStyle = ':'; h.Color = method_colors_3; set(h, 'HandleVisibility', 'off');
x_h = get(h, 'XData'); y_h = get(h, 'YData');
h = scatter(x_h(10:4000:end),y_h(10:4000:end),50, method_colors_3, 'Marker', markertypes(3), 'LineWidth', 2);
set(h, 'HandleVisibility', 'off');
plot([25,26],[0.5,0.75],":"+markertypes(3),'LineWidth',2,'Color',method_colors_3)

h = cdfplot(sqrt(mean(error_all_2d_ablation2.^2,5)));
h.LineWidth = 2; h.LineStyle = '-.'; h.Color = method_colors_4; set(h, 'HandleVisibility', 'off');
x_h = get(h, 'XData'); y_h = get(h, 'YData');
h = scatter(x_h(10:4000:end),y_h(10:4000:end),50, method_colors_4, 'Marker', markertypes(4), 'LineWidth', 2);
set(h, 'HandleVisibility', 'off');
plot([25,26],[0.5,0.75],"-."+markertypes(4),'LineWidth',2,'Color',method_colors_4)
%cdfplot(sqrt(mean(error_all_2d.^2,5)))
%legend(["R(d)", "R(pro,k=1)", "R(pro(theta))", "R(pro(delta))", "R(pro,k=2)"])
legend(all_legends, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 16)
xlim([8,16])
box on
title('')
xlabel('RMSE [dB]')
ylabel('CDF')
fontsize(FigH,font_size,"points")

%% Median plots

num_known_samples1 = num_known_samples;
data1_all_heights =  error_all_2d_conventional; %error_all_2d error_all_3d num_meas_used_2d num_meas_used_3d
if exist('sha_tar_2d', 'var')
    data1_sha_tar_all_heights =  sha_tar_2d;
    data1_est_tar_all_heights =  est_tar_2d;
end

data1_legend = 'Conventional, \boldmath$R\left(d_{\mathrm{2D}}\right)$';
keep_radius = [2]; % use only 200 m radius instead of 70 m

data2_all_heights = error_all_2d; % num_meas_used_2d error_all_2d
if exist('sha_tar_2d', 'var')
    data2_sha_tar_all_heights =  sha_tar_2d;
    data2_est_tar_all_heights =  est_tar_2d;
end
data2_legend = 'Proposed, \boldmath$R\left(d_{\mathrm{2D}},\theta,\delta\right)$'; % change it

y_label = 'Median of RMSE (dB)'; % change it
y_label2 = 'Valid estimation ratio'; %

for height_i = 1:length(all_heights)
    current_heihgt = all_heights(height_i);
    data1 = squeeze(data1_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,radius,targets)
    data2 = squeeze(data2_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,num_fixed_points,targets)

    if exist('data1_sha_tar_all_heights', 'var')
        data1_est_tar_2d = squeeze(data1_est_tar_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,radius,targets)
        data1_sha_tar_2d = squeeze(data1_sha_tar_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,num_fixed_points,targets)

        data2_est_tar_2d = squeeze(data2_est_tar_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,radius,targets)
        data2_sha_tar_2d = squeeze(data2_sha_tar_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,num_fixed_points,targets)
    end

    [rmse_meas_radius_data1, rmse_meas_radius_data2, valid_est1, valid_est2] = calc_rmse_median(data1,data2);
    FigH = figure('Position', ceil([100,100,800,650]), 'visible','on');
    hold on
    all_legends = [];
    for radi_i = keep_radius
        all_legends = [all_legends, data1_legend + ""]; % + ", R = " + num2str(radius_all(radi_i))+" m"
        plot(num_known_samples1, rmse_meas_radius_data1(:,radi_i),method_colors_1+'-'+radius_plot_type(1),'LineWidth',2)
    end
    for radi_i = keep_radius
        plot(num_known_samples1, rmse_meas_radius_data2(:,radi_i),method_colors_2+radius_plot_type(2),'LineWidth',2)
        all_legends = [all_legends, data2_legend + "" ]; % + ", R = " + num2str(radius_all(radi_i))+" m"
    end
    xlabel('Number of tuning samples, M')
    ylabel(y_label)
    grid on
    font_size = 16;
    fontsize(FigH,font_size,"points")
    box on
    %total 47243 samples , 100 * [50:100:450] / 47243
    xticks(num_known_samples);                     % positions
    xticklabels({'50 (0.1%)','150 (0.3%)','250 (0.5%)','350 (0.7%)','450 (1%)'}); % custom text
    legend(all_legends, 'Location', 'best', 'Interpreter', 'latex', 'FontSize',18);
end


data1_all_heights =  error_all_2d_ablation1; %error_all_2d error_all_3d num_meas_used_2d num_meas_used_3d
if exist('sha_tar_2d', 'var')
    data1_sha_tar_all_heights =  sha_tar_2d;
    data1_est_tar_all_heights =  est_tar_2d;
end
data1_legend = 'Ablation 1, \boldmath$R\left(d_{\mathrm{2D}},\theta\right)$';
data2_all_heights = error_all_2d_ablation2; % num_meas_used_2d error_all_2d
if exist('sha_tar_2d', 'var')
    data2_sha_tar_all_heights =  sha_tar_2d;
    data2_est_tar_all_heights =  est_tar_2d;
end
data2_legend = 'Ablation 2, \boldmath$R\left(d_{\mathrm{2D}},\delta\right)$';
for height_i = 1:length(all_heights)
    current_heihgt = all_heights(height_i);
    data1 = squeeze(data1_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,radius,targets)
    data2 = squeeze(data2_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,num_fixed_points,targets)

    if exist('data1_sha_tar_all_heights', 'var')
        data1_est_tar_2d = squeeze(data1_est_tar_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,radius,targets)
        data1_sha_tar_2d = squeeze(data1_sha_tar_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,num_fixed_points,targets)

        data2_est_tar_2d = squeeze(data2_est_tar_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,radius,targets)
        data2_sha_tar_2d = squeeze(data2_sha_tar_all_heights(:,height_i,:,:,:)); % (iter,height,num_meas,num_fixed_points,targets)
    end
    

    [rmse_meas_radius_data1, rmse_meas_radius_data2, valid_est1, valid_est2] = calc_rmse_median(data1,data2);
    for radi_i = keep_radius
        all_legends = [all_legends, data1_legend + ""]; % + ", R = " + num2str(radius_all(radi_i))+ " m"
        plot(num_known_samples1, rmse_meas_radius_data1(:,radi_i),radius_plot_type(3),'LineWidth',2,'Color',method_colors_3)
    end
    for radi_i = keep_radius
        plot(num_known_samples1, rmse_meas_radius_data2(:,radi_i),radius_plot_type(4),'LineWidth',2,'Color',method_colors_4)
        all_legends = [all_legends, data2_legend + ""]; % + ", R = " + num2str(radius_all(radi_i))+ " m"
    end
    fontsize(FigH,font_size,"points")
    legend(all_legends, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 20);
end