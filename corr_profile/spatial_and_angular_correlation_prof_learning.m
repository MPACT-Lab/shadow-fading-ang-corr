clc;
clear;
close all;

% Constant used in geocoding
R_earth = 6378.137 * 10^3 ; % [m]

% Define custom colors
colors = [
    0.20 0.60 0.80;  % blue
    0.90 0.40 0.20;  % orange
    0.30 0.80 0.30;  % green
    0.70 0.30 0.80;  % purple
    0.10 0.45 0.75   % deep blue,
    0.75 0.20 0.10;  % red
    0.10 0.70 0.70;  % teal
    0.95 0.80 0.25;  % yellow-gold
    0.85 0.33 0.10;  % reddish-brown
    0.50 0.50 0.50;  % gray
    0.60 0.30 0.10;  % brown
    0.10 0.45 0.85   % deep blue,
    0.00 1.00 0.00;  % deep green
];

%% load data
% variables for storing UAV x,y,z coordinates, tilt angle \delta,
% elevation angle \theta
lat_all_collection = [];
lon_all_collection = [];
h_all_collection = [];
sha_all_collection = [];
arm_all_collection = [];
elev_all_collection = [];

data_str = "results_useful_with_tworay_tilt_src_exp_";
experiments = [309]; % should use experiments of same altitudes only
% 288 -> 36 m (UAV altitude)
% 300 -> 48 m (UAV altitude)
% 301 -> 19 m (UAV altitude)
% 309, 328 -> 28 m (UAV altitude)
exp_txt = "309_loc1_2_3"; % change accordingly
locations = [1,2,3];

% load UAV x,y,z, 3D orientation and corresponding shadow fading (SF)
for i = 1:length(experiments)
    for j = 1:length(locations)
        exp_no = experiments(i);
        location = locations(j);
        if exp_no==309 && location==1
            continue
        end
        load("../data_gen/processed_data/"+data_str +num2str(exp_no)+"_loc_"+num2str(location)+".mat")
        lat_all_collection = [lat_all_collection; lat_all];
        lon_all_collection = [lon_all_collection; lon_all];
        h_all_collection = [h_all_collection; h_all];
        sha_all_collection = [sha_all_collection; power_all - RSRP_PL_two]; % sha_all was power
        arm_all_collection = [arm_all_collection; pitch_towards_src*180/pi];
        elev_all_collection = [elev_all_collection; elev];
    end
end
lat_all = lat_all_collection;
lon_all = lon_all_collection;
h_all = h_all_collection;
sha_all = sha_all_collection;
arm_all = arm_all_collection;
elev_all = elev_all_collection;

%% histogram / cdf plot of shadow fading (unconditional and conditional cases)

% tilt angle group separation: i.e., group 1: delta < -7 degree
tilt_ang_sep = [-7, -3.0, 3.0, 7.0];
% elevation angle group separation: i.e., group 1: 0 < theta < 8 degree
elev_sep = [0.0, 8.0, 16.0, 24.0, 90] *pi /180;

font_size = 20;

% unconditional case (for all elevation and tilt angles)
FigH = figure(21);
set(gcf, 'Position', [200 200 565 480]);   % square: 500×500 pixels
h = cdfplot(sha_all);
xlabel('Shadowing [dB]', 'FontSize', font_size); % 'Interpreter', 'latex'
ylabel('CDF', 'FontSize', font_size); % 'Interpreter', 'latex'
xlim([-50, 50])
set(gca, 'FontSize', font_size)
h.LineWidth = 3;        % Set the line width
title('')
grid on;
box on;
hold off;

% conditional on elevation angle (histogram plot)
FigH = figure(23); clf; hold on;
legend_entries = cell(1, length(elev_sep)-1);
for ii = 1:length(elev_sep)-1
    % Extract data in current separation range
    sha_all_sep = sha_all(elev_all >= elev_sep(ii) & elev_all < elev_sep(ii+1));
    % Plot histogram for current range
    h = histogram(sha_all_sep,'Normalization', 'probability', ...   % optional: normalize
        'FaceColor', colors(ii,:), ...
        'BinWidth', 1.0, ...
        'FaceAlpha', 0.7);
    legend_entries{ii} = sprintf('%.1f',elev_sep(ii)*180/pi)+"$^\circ\,\leq \theta <$"+sprintf('%.1f',elev_sep(ii+1)*180/pi)+"$^\circ$";%sprintf('%.1f',elev_sep(ii)*180/pi) + "$\leq elev <$" + sprintf('%.1f$', elev_sep(ii+1)*180/pi);
end
xlabel('Shadowing [dB]', 'Interpreter', 'latex', 'FontSize', font_size); 
ylabel('Probability', 'Interpreter', 'latex', 'FontSize', font_size);
legend(legend_entries, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
xlim([-50, 50])
grid on;
box on;
hold off;

% conditional on elevation angle (cdf plot)
FigH = figure(231); clf; hold on;
set(gcf, 'Position', [200 200 565 480]);   % square: 500×500 pixels
legend_entries = cell(1, length(elev_sep)-1);
line_styles = ["-", "--",":","-."];
for ii = 1:length(elev_sep)-1
    % Extract data in current separation range
    sha_all_sep = sha_all(elev_all >= elev_sep(ii) & elev_all < elev_sep(ii+1));
    h = cdfplot(sha_all_sep);
    h.LineWidth = 3;
    h.LineStyle = line_styles(ii);
    h.Color = colors(ii,:);
    legend_entries{ii} = sprintf('%d',elev_sep(ii)*180/pi)+"$^\circ\,\leq \theta <$"+sprintf('%d',elev_sep(ii+1)*180/pi)+"$^\circ$";%sprintf('%.1f',elev_sep(ii)*180/pi) + "$\leq elev <$" + sprintf('%.1f$', elev_sep(ii+1)*180/pi);
end
xlabel('Shadowing [dB]', 'FontSize', font_size); %, 'Interpreter', 'latex'
ylabel('CDF', 'FontSize', font_size); % 'Interpreter', 'latex'
xlim([-50, 50])
set(gca, 'FontSize', font_size)
legend(legend_entries, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 22, 'FontWeight', 'bold');
grid on;
title('')
box on;
hold off;

% conditional on both elevation angle and tilt angle (cdf plot)
line_styles = ["-", "--",":","-."];
for ii = 1:length(elev_sep)-1
    FigH = figure(100+ii); clf; hold on;
    set(gcf, 'Position', [200 200 565 480]);   % square: 500×500 pixels
    legend_entries = cell(1, length(tilt_ang_sep)-1);
    for jj = 1:length(tilt_ang_sep)-1
        % Extract data in current separation range
        sha_all_sep = sha_all(arm_all >= tilt_ang_sep(jj) & arm_all < tilt_ang_sep(jj+1)...
            & elev_all >= elev_sep(ii) & elev_all < elev_sep(ii+1));
        % Plot cdf for current range
        h = cdfplot(sha_all_sep);
        h.Color = colors(jj,:);
        h.LineWidth = 3;
        h.LineStyle = line_styles(jj);
        legend_entries{jj} = sprintf('%d',tilt_ang_sep(jj))+"$^\circ\leq \delta <$"+sprintf('%d',tilt_ang_sep(jj+1))+"$^\circ$";
    end
    xlabel('Shadowing [dB]', 'FontSize', font_size);  % 'Interpreter', 'latex'
    ylabel('CDF', 'FontSize', font_size); % 'Interpreter', 'latex'
    set(gca, 'FontSize', font_size)
    legend(legend_entries, 'Location', 'best', 'FontSize', 22, 'Interpreter', 'latex', 'FontWeight', 'bold'); % 
    xlim([-50, 50])
    title('')
    grid on;
    box on;
    hold off;
end

%% Spatial correlation of SF (Modeling SF as a biexponential function of 2D distance)

% to reduce computation, get ~8k samples randomly from the 109k total samples
rng(100)
N_seed = 8000;
s_idx = randi(length(sha_all),round(N_seed*1.5),1) ;
s_idx = shuffle(unique(s_idx)) ; % caution unique sorts the numbers
N_seed = length(s_idx);

% horizontal distance grid: from 0 m to 250 m
d_int = 2 ;
d_max = 250 ;
d_grid = [0:d_int:d_max]' ;

% mean and variance
mu_sha = mean(sha_all(s_idx)) ;
var_sha = var(sha_all(s_idx)) ;

% pair-wise 2D distance and correlation calculation
dist_hor = zeros(N_seed, N_seed);
corr_all = zeros(N_seed, N_seed);
for ii=1:N_seed
    dist_hor(:,ii) = real(R_earth .* acos(sin(lat_all(s_idx(ii)) * pi/180) .* sin(lat_all(s_idx) * pi/180) + cos(lat_all(s_idx(ii)) * pi/180) .* cos(lat_all(s_idx) * pi/180) .* cos( (lon_all(s_idx) - lon_all(s_idx(ii)))  * pi/180 ))) ;
    corr_all(:,ii) = ( sha_all(s_idx(ii)) - mu_sha ) .* ( sha_all(s_idx) - mu_sha ) / var_sha ;
end
dist_hor2 = dist_hor(:);
corr_all2 = corr_all(:);

% correlation for different 2D distances
corr_grid = zeros(1, length(d_grid)-1);
N_grid = zeros(1, length(d_grid)-1);
hor_dist_grid = zeros(1, length(d_grid)-1);
for ii = 1:length(d_grid)-1
    corr_tmp1 = corr_all2(dist_hor2>=d_grid(ii) & dist_hor2<d_grid(ii+1));
    N_grid(1,ii) = length(corr_tmp1) ;
    corr_grid(1,ii) = mean(corr_tmp1) ;
    hor_dist_grid(1,ii) = (d_grid(ii) + d_grid(ii+1))/2;
end
hor_dist_list = hor_dist_grid(:);
corr_list = corr_grid(:);

% filter nan values
valid_idx = ~isnan(corr_list);
hor_dist_list = hor_dist_list(valid_idx);
corr_list = corr_list(valid_idx); 

% fit with a biexponential model
ft = fittype( '(a*exp(-(b*x)) + (1-a)*exp(-(c*x)))', 'independent', {'x'}, 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', [1,4,4], 'Lower',[0,0,0]);
correlation_f = fit([hor_dist_list],corr_list,ft,opts);
% Evaluate the fitted model at the input data points
z_fitted = correlation_f(hor_dist_list);
% Calculate residuals
residuals = corr_list - z_fitted;
% Calculate summary statistics of residuals
mean_residual = mean(residuals);
std_residual = std(residuals);
sum_squared_residuals = sum(residuals.^2);
disp(['Mean Residual: ', num2str(mean_residual)]);
disp(['Standard Deviation of Residuals: ', num2str(std_residual)]);
disp(['Sum of Squared Residuals: ', num2str(sum_squared_residuals)]);

% plot the fitted curve across scatter correlation values
figure(3)
hold on
scatter(hor_dist_list, corr_list, 5, corr_list) % Use scatter instead of scatter3 for 2D plot
colormap(jet);
colorbar;
caxislim = [0 1];
caxis(caxislim)
xlabel('Horizontal distance [m]')
ylabel('Correlation')
hcb = colorbar;
colorTitleHandle = get(hcb, 'Title');
titleString = 'Correlation';
set(colorTitleHandle, 'String', titleString);
set(gcf, 'color', 'w');
grid on;

% Define the fitted function
fitted_f = @(x) (correlation_f.a * exp(-(correlation_f.b * x)) + (1 - correlation_f.a) * exp(-(correlation_f.c * x)));

% Plot the fitted curve using fplot (for 2D)
fplot(fitted_f, [min(hor_dist_list), max(hor_dist_list)], 'k-', 'LineWidth', 2);

%% Angular correlation of SF

% tilt angle group separation: i.e., group 1: , 30 < delta < -7 degree
tilt_ang_sep = [-30.0 -7.0 -3.0, 3.0, 7.0, 30.0];
tilt_ang_sep_avg = [-10.0 -5.0 0.0, 5.0, 10.0]; % length reduces by 1

% elevation angle group separation: i.e., group 1: 0 < theta < 10 degree
elev_sep = [0.0, 10.0, 30.0, 50.0, 90] *pi /180;
elev_sep_avg = [5.0, 20.0, 40, 70] *pi / 180;


% calculate pairwise correlation among joint groups (elevation, tilt),
% i.e., between (\theta_1, \tilt_1) and (\theta_2, \tilt_2) 
total_sep = (length(tilt_ang_sep)-1) * (length(elev_sep)-1);
corr_elev_tilt = zeros(total_sep, total_sep);
corr_tilt_elev = zeros(total_sep, total_sep);
N_elev_tilt = zeros(total_sep, total_sep);
N_tilt_elev = zeros(total_sep, total_sep);
for ii = 1:length(elev_sep)-1
    for jj = 1:length(tilt_ang_sep)-1
        ii_jj = (ii-1)*(length(tilt_ang_sep)-1)+jj;
        ii_jj_p_e = (jj-1)*(length(elev_sep)-1)+ii;
        % Extract data in current separation range
        idx_sep = find(arm_all >= tilt_ang_sep(jj) & arm_all < tilt_ang_sep(jj+1)...
            & elev_all >= elev_sep(ii) & elev_all < elev_sep(ii+1));
        if length(idx_sep) > 4000
            s_idx = randi(length(idx_sep),round(4000*1.2),1) ;
            s_idx = shuffle(unique(s_idx)) ; % caution unique sorts the numbers
            idx_sep = idx_sep(s_idx);
        end
        sha_all_sep = sha_all(idx_sep);
        lat_sep = lat_all(idx_sep);
        lon_sep = lon_all(idx_sep);

        for ii2 = 1:length(elev_sep)-1
            for jj2 = 1:length(tilt_ang_sep)-1
                ii_jj2 = (ii2-1)*(length(tilt_ang_sep)-1)+jj2;
                ii_jj2_p_e = (jj2-1)*(length(elev_sep)-1)+ii2;
                if ii_jj2 < ii_jj
                    N_elev_tilt(ii_jj, ii_jj2) = N_elev_tilt(ii_jj2, ii_jj);
                    corr_elev_tilt(ii_jj, ii_jj2) = corr_elev_tilt(ii_jj2, ii_jj);
                    continue
                end
                % Extract data in current separation range
                idx_sep2 = find(arm_all >= tilt_ang_sep(jj2) & arm_all < tilt_ang_sep(jj2+1)...
                    & elev_all >= elev_sep(ii2) & elev_all < elev_sep(ii2+1));
                if length(idx_sep2) > 4000
                    s_idx = randi(length(idx_sep2),round(4000*1.2),1) ;
                    s_idx = shuffle(unique(s_idx)) ; % caution unique sorts the numbers
                    idx_sep2 = idx_sep2(s_idx);
                end
                sha_all_sep2 = sha_all(idx_sep2);
                lat_sep2 = lat_all(idx_sep2);
                lon_sep2 = lon_all(idx_sep2);

                sha_all_sep = sort(sha_all_sep);
                sha_all_sep2 = sort(sha_all_sep2);
                min_length = min(length(sha_all_sep), length(sha_all_sep2));
                if (min_length<10)
                    corr_elev_tilt(ii_jj, ii_jj2) = NaN;
                    N_elev_tilt(ii_jj, ii_jj2) = min_length;

                    corr_tilt_elev(ii_jj_p_e, ii_jj2_p_e) = NaN;
                    N_tilt_elev(ii_jj_p_e, ii_jj2_p_e) = min_length;
                    corr_tilt_elev(ii_jj2_p_e, ii_jj_p_e) = NaN;
                    N_tilt_elev(ii_jj2_p_e, ii_jj_p_e) = min_length;
                    continue
                end
                % equalize number of samples in sha_all_sep and sha_all_sep2
                % keeping the distribution same
                if length(sha_all_sep) ~= length(sha_all_sep2)
                    if length(sha_all_sep)> length(sha_all_sep2)
                        k_factor = length(sha_all_sep)/ length(sha_all_sep2);
                        sha_all_sep2_resampled = interp1((1:length(sha_all_sep2))*k_factor, sha_all_sep2, 1:length(sha_all_sep));
                        sha_all_sep2 = sha_all_sep2_resampled';
                    else
                        k_factor = length(sha_all_sep2)/ length(sha_all_sep);
                        sha_all_sep_resampled = interp1((1:length(sha_all_sep))*k_factor, sha_all_sep, 1:length(sha_all_sep2));
                        sha_all_sep = sha_all_sep_resampled';
                    end
                end
                valid_indx = (~isnan(sha_all_sep)) & (~isnan(sha_all_sep2));
                sha_all_sep = sha_all_sep(valid_indx);
                sha_all_sep2 = sha_all_sep2(valid_indx);
 
                var1 = mean((sha_all_sep - mu_sha).^2);
                var2 = mean((sha_all_sep2 - mu_sha).^2);
                raw_correlation = (sha_all_sep - mu_sha) .* (sha_all_sep2 - mu_sha) ./ sqrt(var1*var2);

                corr_elev_tilt(ii_jj, ii_jj2) = mean(raw_correlation(:),'omitnan');
                N_elev_tilt(ii_jj, ii_jj2) = min_length;

                corr_tilt_elev(ii_jj_p_e, ii_jj2_p_e) = corr_elev_tilt(ii_jj, ii_jj2);
                N_tilt_elev(ii_jj_p_e, ii_jj2_p_e) = min_length;
                corr_tilt_elev(ii_jj2_p_e, ii_jj_p_e) = corr_elev_tilt(ii_jj, ii_jj2);
                N_tilt_elev(ii_jj2_p_e, ii_jj_p_e) = min_length;
            end
        end
    end
end

% show it as 4x4 (for elevation) or 5x5 matrix (for tilt)
% for 5 different tilt angles will plot 5 plots, each having 4x4
% elevation-wise correlations
font_size = 20;
for i = 1:length(tilt_ang_sep)-1
    figH = figure('Position', [500 500 400 300]);
    submatrix = corr_tilt_elev((length(elev_sep)-1)*(i-1)+(1:length(elev_sep)-1), (length(elev_sep)-1)*(i-1)+(1:length(elev_sep)-1));
    imagesc(submatrix)
    set(gca,'Color','white');          % background color
    alphaData = ~isnan(submatrix);            % NaNs become transparent
    set(gca().Children, 'AlphaData', alphaData);
    submatrix;
    axis xy
    colormap(jet);
    colorbar;
    caxis([0 1])
    hcb=colorbar;
    hcb.Ticks = [0 1];   % set tick positions
    xticks(1:4)
    xticklabels({'5^\circ','20^\circ','40^\circ','70^\circ'})
    yticks(1:4)
    yticklabels({'5^\circ','20^\circ','40^\circ','70^\circ'})
    set(gca, 'FontSize',18);
    xlabel('Elevation angle, $\theta_1$', 'Interpreter', 'latex', 'FontSize', font_size);
    ylabel('Elevation angle, $\theta_2$', 'Interpreter', 'latex', 'FontSize', font_size);

end

% for 4 different elevation angles will plot 4 plots, each having 5x5
% tilt-wise correlations
font_size = 18;
for i = 1:length(elev_sep)-1
    figH = figure('Position', [500 500 400 300]);
    submatrix = corr_elev_tilt((length(tilt_ang_sep)-1)*(i-1)+(1:length(tilt_ang_sep)-1), (length(tilt_ang_sep)-1)*(i-1)+(1:length(tilt_ang_sep)-1));
    imagesc(submatrix)
    set(gca,'Color','white');          % background color
    alphaData = ~isnan(submatrix);            % NaNs become transparent
    set(gca().Children, 'AlphaData', alphaData);
    submatrix;
    axis xy
    colormap(jet);
    caxis([0.8, 1])
    hcb=colorbar;
    hcb.Ticks = [0.8 1];   % set tick positions
    xticks(1:5)
    xticklabels({'-10^\circ','-5^\circ','0^\circ','5^\circ','10^\circ'})
    yticks(1:5)
    yticklabels({'-10^\circ','-5^\circ','0^\circ','5^\circ','10^\circ'})
    %set(gca, 'TickLabelInterpreter','latex');
    set(gca, 'FontSize',16);
    xlabel('Tilt angle, $\delta_1$', 'Interpreter', 'latex', 'FontSize', font_size);
    ylabel('Tilt angle, $\delta_2$', 'Interpreter', 'latex', 'FontSize', font_size);

end

%% Fitting the SF Angular correlation with exponential models

% we will process tilt difference for identical elevation angles (theta_i =
% theta_j), and we fit an exponential fit across tilt difference for each
% of the elevations (5), for each tilt angles (4)
% we also consider 2 directions separately (increasing tilt and decreasing tilt) we have
% another dimension of 2, so we will fit total 4 x 5 x 2 curves with
% exponentials

% similarly we will also have correlation across elevation deifference, 
% then we will fit have another 5 x 4 x 2 curves with exponentials

% Ablation study: if we merge all the elevation angles and fit across tilt
% angle differences for each tilt angles (5) , we will fit 5 x 2 curves

% Ablation study: if we merge all the tilt angles and fit across elevation
% angle differences for each elevation angles (4) , we will fit 4 x 2 curves

% exponential fitting of (delta_1, delta_2, theta)
tilt_sensitivity_exp = zeros(length(elev_sep)-1, length(tilt_ang_sep)-1, 2);
for ii = 1:length(elev_sep)-1
    for jj = 1:length(tilt_ang_sep)-1
        ii_jj = (ii-1)*(length(tilt_ang_sep)-1)+jj;
        ii_jj_p_e = (jj-1)*(length(elev_sep)-1)+ii;  % not used in this subsection
        % ii2 = ii
        left_set = (ii-1)*(length(tilt_ang_sep)-1) + (1:jj);
        right_set = (ii-1)*(length(tilt_ang_sep)-1) + (jj:(length(tilt_ang_sep)-1));
        left_set_data = corr_elev_tilt(ii_jj,left_set);
        right_set_data = corr_elev_tilt(ii_jj,right_set);

        left_set_delta = abs(tilt_ang_sep_avg(1:jj) - tilt_ang_sep_avg(jj));
        right_set_delta = abs(tilt_ang_sep_avg(jj:end) - tilt_ang_sep_avg(jj));
        
        valid_index = ~isnan(left_set_data);
        left_set_data = left_set_data(valid_index);
        left_set_delta = left_set_delta(valid_index);
        valid_index = ~isnan(right_set_data);
        right_set_data = right_set_data(valid_index);
        right_set_delta = right_set_delta(valid_index);

        if length(left_set_data) > 1
            ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 1, 'Lower',0);
            correlation_temp = fit([left_set_delta'],left_set_data',ft,opts);
            tilt_sensitivity_exp(ii, jj, 1) = correlation_temp.a;
        end

        if length(right_set_data) > 1
            ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 1, 'Lower',0);
            correlation_temp = fit([right_set_delta'],right_set_data',ft,opts);
            tilt_sensitivity_exp(ii, jj, 2) = correlation_temp.a;
        end
    end
end

% figure
% imagesc(tilt_sensitivity_exp(:,:,1))
% colormap(jet);
% colorbar;
% 
% figure
% imagesc(tilt_sensitivity_exp(:,:,2))
% colormap(jet);
% colorbar;

% exponential fitting of (theta_1, theta_2, delta)
elev_sensitivity_exp = zeros(length(tilt_ang_sep)-1, length(elev_sep)-1, 2);
for ii = 1:length(elev_sep)-1
    for jj = 1:length(tilt_ang_sep)-1
        ii_jj = (ii-1)*(length(tilt_ang_sep)-1)+jj; % not used in this subsection
        ii_jj_p_e = (jj-1)*(length(elev_sep)-1)+ii;
        % ii2 = ii
        left_set = (jj-1)*(length(elev_sep)-1) + (1:ii);
        right_set = (jj-1)*(length(elev_sep)-1) + (ii:(length(elev_sep)-1));
        left_set_data = corr_tilt_elev(ii_jj_p_e,left_set);
        right_set_data = corr_tilt_elev(ii_jj_p_e,right_set);

        left_set_delta = abs(elev_sep_avg(1:ii) - elev_sep_avg(ii));
        right_set_delta = abs(elev_sep_avg(ii:end) - elev_sep_avg(ii));
        
        valid_index = ~isnan(left_set_data);
        left_set_data = left_set_data(valid_index);
        left_set_delta = left_set_delta(valid_index);
        valid_index = ~isnan(right_set_data);
        right_set_data = right_set_data(valid_index);
        right_set_delta = right_set_delta(valid_index);

        if length(left_set_data) > 1
            ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 6, 'Lower',0);
            correlation_temp = fit([left_set_delta'],left_set_data',ft,opts);
            elev_sensitivity_exp(jj, ii, 1) = correlation_temp.a;
        end

        if length(right_set_data) > 1
            ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 6, 'Lower',0);
            correlation_temp = fit([right_set_delta'],right_set_data',ft,opts);
            elev_sensitivity_exp(jj, ii, 2) = correlation_temp.a;
        end
    end
end

% figure
% imagesc(elev_sensitivity_exp(:,:,1))
% colormap(jet);
% colorbar;
% 
% figure
% imagesc(elev_sensitivity_exp(:,:,2))
% colormap(jet);
% colorbar;

%% exponential fitting example for (delta_1, delta_2, theta)

font_size = 14;
markertypes = ['^','s','v','o','<'];

fig = figure(320);
set(gcf, 'Position', [200 200 565 480]);   % square: 500×500 pixels
hold on
% for 3 (3rd index of \delta, -3< >3)
for_what  = 3;
legend_entries = cell(1, length(elev_sep)-1);
for i = 1:length(elev_sep)-1
%scatter((1:length(tilt_ang_sep)-1),corr_elev_tilt((i-1)*(length(tilt_ang_sep)-1) + for_what, (i-1)*(length(tilt_ang_sep)-1) + (1: (length(tilt_ang_sep)-1))) , 80, colors(i,:), 'Marker', markertypes(i),'LineWidth',2 )
plot(tilt_ang_sep_avg, corr_elev_tilt((i-1)*(length(tilt_ang_sep)-1) + for_what, (i-1)*(length(tilt_ang_sep)-1) + (1: (length(tilt_ang_sep)-1))) , "--"+markertypes(i), 'Color', colors(i,:), 'LineWidth', 2 )
left_set_data = corr_elev_tilt((i-1)*(length(tilt_ang_sep)-1) + for_what, (i-1)*(length(tilt_ang_sep)-1) + (1: for_what));
right_set_data = corr_elev_tilt((i-1)*(length(tilt_ang_sep)-1) + for_what, (i-1)*(length(tilt_ang_sep)-1) + (for_what: length(tilt_ang_sep)-1));

left_set_delta = abs(tilt_ang_sep_avg(1:for_what) - tilt_ang_sep(for_what));
right_set_delta = abs(tilt_ang_sep_avg(for_what:end) - tilt_ang_sep(for_what));

valid_index = ~isnan(left_set_data);
left_set_data = left_set_data(valid_index);
left_set_delta = left_set_delta(valid_index);
valid_index = ~isnan(right_set_data);
right_set_data = right_set_data(valid_index);
right_set_delta = right_set_delta(valid_index);

if length(left_set_data) > 1
    ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 6, 'Lower',0);
    correlation_temp = fit([left_set_delta'],left_set_data',ft,opts);
    plot_x_data = -10:0.01:tilt_ang_sep_avg(for_what);
    plot_x_data_scatter = -10+1.25*(i-1):5:tilt_ang_sep_avg(for_what);
    plot_y_data = exp(-correlation_temp.a*abs(plot_x_data-tilt_ang_sep_avg(for_what)));
    plot_y_data_scatter = exp(-correlation_temp.a*abs(plot_x_data_scatter-tilt_ang_sep_avg(for_what)));
    h = plot(plot_x_data, plot_y_data, 'LineWidth', 2, 'Color', colors(i,:));
    set(h, 'HandleVisibility', 'off');
    h = scatter(plot_x_data_scatter, plot_y_data_scatter, 50, colors(i,:), 'Marker', markertypes(i), 'LineWidth', 2);
    set(h, 'HandleVisibility', 'off');
end

if length(right_set_data) > 1
    ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 6, 'Lower',0);
    correlation_temp = fit([right_set_delta'],right_set_data',ft,opts);
    plot_x_data = tilt_ang_sep_avg(for_what):0.01:10;
    plot_x_data_scatter = 1.25*(i-1)+tilt_ang_sep_avg(for_what):5:10;
    plot_y_data = exp(-correlation_temp.a*abs(plot_x_data-tilt_ang_sep_avg(for_what)));
    plot_y_data_scatter = exp(-correlation_temp.a*abs(plot_x_data_scatter-tilt_ang_sep_avg(for_what)));
    h = plot(plot_x_data, plot_y_data, 'LineWidth', 2, 'Color', colors(i,:));
    set(h, 'HandleVisibility', 'off');
    h = scatter(plot_x_data_scatter, plot_y_data_scatter, 50, colors(i,:), 'Marker', markertypes(i), 'LineWidth', 2);
    set(h, 'HandleVisibility', 'off');
end
 legend_entries{i} = sprintf('%.1f',elev_sep(i)*180/pi)+"$^\circ\,\leq \theta <$"+sprintf('%.1f',elev_sep(i+1)*180/pi)+"$^\circ$";%sprintf('%.1f',elev_sep(ii)*180/pi) + "$\leq elev <$" + sprintf('%.1f$', elev_sep(ii+1)*180/pi);
 end
%xlabel('Elevation, $\theta$ [degree]', 'Interpreter', 'latex'); 
xlabel('Tilt angle, $\delta_2$ [degree]', 'Interpreter', 'latex', 'FontSize', font_size); 
xticks(-10:5:10)
xticklabels({'-10^\circ','-5^\circ','0^\circ','5^\circ','10^\circ'})
ylabel('Correlation, $R_{\mathrm{tlt}}\left(\delta_1,\delta_2,\theta\right)$', 'Interpreter', 'latex', 'FontSize', font_size);
legend(legend_entries, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
grid on
box on
ylim([0.9, 1.0])
set(gca, 'FontSize', font_size)

%% exponential fitting example for (theta_1, theta_2, delta)
fig = figure(321);
set(gcf, 'Position', [200 200 565 480]);   % square: 500×500 pixels
hold on
% for 3 (3rd index of \theta, 30< >50)
for_what  = 2;
legend_entries = cell(1, length(tilt_ang_sep)-1);
for i = 1:length(tilt_ang_sep)-1
%scatter((1:length(tilt_ang_sep)-1),corr_elev_tilt((i-1)*(length(tilt_ang_sep)-1) + for_what, (i-1)*(length(tilt_ang_sep)-1) + (1: (length(tilt_ang_sep)-1))) , 80, colors(i,:), 'Marker', markertypes(i),'LineWidth',2 )
plot(elev_sep_avg, corr_tilt_elev((i-1)*(length(elev_sep)-1) + for_what, (i-1)*(length(elev_sep)-1) + (1: (length(elev_sep)-1))) , "--"+markertypes(i), 'Color', colors(i,:), 'LineWidth', 2 )
left_set_data = corr_tilt_elev((i-1)*(length(elev_sep)-1) + for_what, (i-1)*(length(elev_sep)-1) + (1: for_what));
right_set_data = corr_tilt_elev((i-1)*(length(elev_sep)-1) + for_what, (i-1)*(length(elev_sep)-1) + (for_what: length(elev_sep)-1));

left_set_delta = abs(elev_sep_avg(1:for_what) - elev_sep(for_what));
right_set_delta = abs(elev_sep_avg(for_what:end) - elev_sep(for_what));

valid_index = ~isnan(left_set_data);
left_set_data = left_set_data(valid_index);
left_set_delta = left_set_delta(valid_index);
valid_index = ~isnan(right_set_data);
right_set_data = right_set_data(valid_index);
right_set_delta = right_set_delta(valid_index);

if length(left_set_data) > 1
    ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 6, 'Lower',0);
    correlation_temp = fit([left_set_delta'],left_set_data',ft,opts);
    plot_x_data = 0:0.01*pi/180:elev_sep_avg(for_what);
    plot_x_data_scatter = 5*pi/180*(i-1):20*pi/180:elev_sep_avg(for_what);
    plot_y_data = exp(-correlation_temp.a*abs(plot_x_data-elev_sep_avg(for_what)));
    plot_y_data_scatter = exp(-correlation_temp.a*abs(plot_x_data_scatter-elev_sep_avg(for_what)));
    h = plot(plot_x_data, plot_y_data, 'LineWidth', 2, 'Color', colors(i,:));
    set(h, 'HandleVisibility', 'off');
    h = scatter(plot_x_data_scatter, plot_y_data_scatter, 50, colors(i,:), 'Marker', markertypes(i), 'LineWidth', 2);
    set(h, 'HandleVisibility', 'off');
end

if length(right_set_data) > 1
    ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
    opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 6, 'Lower',0);
    correlation_temp = fit([right_set_delta'],right_set_data',ft,opts);
    plot_x_data = elev_sep_avg(for_what):0.01*pi/180:90*pi/180;
    plot_x_data_scatter = elev_sep_avg(for_what)+5*pi/180*(i-1):20*pi/180:90*pi/180;
    plot_y_data = exp(-correlation_temp.a*abs(plot_x_data-elev_sep_avg(for_what)));
    plot_y_data_scatter = exp(-correlation_temp.a*abs(plot_x_data_scatter-elev_sep_avg(for_what)));
    h = plot(plot_x_data, plot_y_data, 'LineWidth', 2, 'Color', colors(i,:));
    set(h, 'HandleVisibility', 'off');
    h = scatter(plot_x_data_scatter, plot_y_data_scatter, 50, colors(i,:), 'Marker', markertypes(i), 'LineWidth', 2);
    set(h, 'HandleVisibility', 'off');
end
if i==1
    legend_entries{i} = "$\delta <$"+sprintf('%.1f',tilt_ang_sep(i+1))+"$^\circ$";%sprintf('%.1f',elev_sep(ii)*180/pi) + "$\leq elev <$" + sprintf('%.1f$', elev_sep(ii+1)*180/pi);
elseif i==5
    legend_entries{i} = "$\delta >$" + sprintf('%.1f',tilt_ang_sep(i)) + "$^\circ$";
else
    legend_entries{i} = sprintf('%.1f',tilt_ang_sep(i))+"$^\circ\,\leq \delta <$"+sprintf('%.1f',tilt_ang_sep(i+1))+"$^\circ$";%sprintf('%.1f',elev_sep(ii)*180/pi) + "$\leq elev <$" + sprintf('%.1f$', elev_sep(ii+1)*180/pi);
    
end
 end
xlabel('Elevation angle, $\theta_2$ [degree]', 'Interpreter', 'latex', 'FontSize', font_size); 
xticks([0, 5, 20, 40, 70, 90]*pi/180)
xticklabels({'0^\circ','5^\circ','20^\circ','40^\circ','70^\circ','90^\circ'})
ylabel('Correlation, $R_{\mathrm{elv}}\left(\theta_1,\theta_2,\delta\right)$', 'Interpreter', 'latex', 'FontSize', font_size);
legend(legend_entries, 'Location', 'best', 'Interpreter', 'latex', 'FontSize', 14, 'FontWeight', 'bold');
grid on
box on
%ylim([0.9, 1.0])
set(gca, 'FontSize', font_size)

%% ablation study (using only elevation)

total_sep = length(elev_sep)-1;
corr_elev = zeros(total_sep, total_sep);
N_elev = zeros(total_sep, total_sep);

for ii = 1:length(elev_sep)-1
        ii_jj = ii;
        % Extract data in current separation range
        idx_sep = find(elev_all >= elev_sep(ii) & elev_all < elev_sep(ii+1));
        if length(idx_sep) > 4000
            s_idx = randi(length(idx_sep),round(4000*1.2),1) ;
            s_idx = shuffle(unique(s_idx)) ; % caution unique sorts the numbers
            idx_sep = idx_sep(s_idx);
        end
        sha_all_sep = sha_all(idx_sep);
        lat_sep = lat_all(idx_sep);
        lon_sep = lon_all(idx_sep);

        for ii2 = 1:length(elev_sep)-1
                ii_jj2 = ii2;
                if ii_jj2 < ii_jj
                    N_elev(ii_jj, ii_jj2) = N_elev(ii_jj2, ii_jj);
                    corr_elev(ii_jj, ii_jj2) = corr_elev(ii_jj2, ii_jj);
                    continue
                end
                % Extract data in current separation range
                idx_sep2 = find(elev_all >= elev_sep(ii2) & elev_all < elev_sep(ii2+1));
                if length(idx_sep2) > 4000
                    s_idx = randi(length(idx_sep2),round(4000*1.2),1) ;
                    s_idx = shuffle(unique(s_idx)) ; % caution unique sorts the numbers
                    idx_sep2 = idx_sep2(s_idx);
                end
                sha_all_sep2 = sha_all(idx_sep2);
                lat_sep2 = lat_all(idx_sep2);
                lon_sep2 = lon_all(idx_sep2);
                
                sha_all_sep = sort(sha_all_sep);
                sha_all_sep2 = sort(sha_all_sep2);
                min_length = min(length(sha_all_sep), length(sha_all_sep2));
                % equalize number of samples in sha_all_sep and sha_all_sep2
                % keeping the distribution same
                if length(sha_all_sep) ~= length(sha_all_sep2)
                    if length(sha_all_sep)> length(sha_all_sep2)
                        k_factor = length(sha_all_sep)/ length(sha_all_sep2);
                        sha_all_sep2_resampled = interp1((1:length(sha_all_sep2))*k_factor, sha_all_sep2, 1:length(sha_all_sep));
                        sha_all_sep2 = sha_all_sep2_resampled';
                    else
                        k_factor = length(sha_all_sep2)/ length(sha_all_sep);
                        sha_all_sep_resampled = interp1((1:length(sha_all_sep))*k_factor, sha_all_sep, 1:length(sha_all_sep2));
                        sha_all_sep = sha_all_sep_resampled';
                    end
                end
                valid_indx = (~isnan(sha_all_sep)) & (~isnan(sha_all_sep2));
                sha_all_sep = sha_all_sep(valid_indx);
                sha_all_sep2 = sha_all_sep2(valid_indx);
                var1 = mean((sha_all_sep - mu_sha).^2);
                var2 = mean((sha_all_sep2 - mu_sha).^2);
                raw_correlation = (sha_all_sep - mu_sha) .* (sha_all_sep2 - mu_sha) ./ sqrt(var1*var2);

                corr_elev(ii_jj, ii_jj2) = mean(raw_correlation(:),'omitnan');
                N_elev(ii_jj, ii_jj2) = min_length;


        end
end

fig = figure(33);
imagesc(corr_elev)
axis xy
colormap(jet);
colorbar;
caxis([0 max(corr_elev(:))])
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = 'Correlation';
set(colorTitleHandle ,'String',titleString);
xlabel('Pitch [degree]'); ylabel('Pitch [degree]'); zlabel('Correlation');

% fitting exponential
elev_only_sensitivity_exp = zeros(1, length(elev_sep)-1, 2);
for ii = 1:length(elev_sep)-1
        left_set = (1:ii);
        right_set = (ii:(length(elev_sep)-1));
        left_set_data = corr_elev(ii,left_set);
        right_set_data = corr_elev(ii,right_set);

        left_set_delta = abs(elev_sep_avg(1:ii) - elev_sep_avg(ii));
        right_set_delta = abs(elev_sep_avg(ii:end) - elev_sep_avg(ii));
        
        valid_index = ~isnan(left_set_data);
        left_set_data = left_set_data(valid_index);
        left_set_delta = left_set_delta(valid_index);
        valid_index = ~isnan(right_set_data);
        right_set_data = right_set_data(valid_index);
        right_set_delta = right_set_delta(valid_index);

        if length(left_set_data) > 1
            ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 6, 'Lower',0);
            correlation_temp = fit([left_set_delta'],left_set_data',ft,opts);
            elev_only_sensitivity_exp(1, ii, 1) = correlation_temp.a;
        end

        if length(right_set_data) > 1
            ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
            opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 6, 'Lower',0);
            correlation_temp = fit([right_set_delta'],right_set_data',ft,opts);
            elev_only_sensitivity_exp(1, ii, 2) = correlation_temp.a;
        end
end

% figure
% imagesc(elev_only_sensitivity_exp(:,:,1))
% colormap(jet);
% colorbar;
% 
% figure
% imagesc(elev_only_sensitivity_exp(:,:,2))
% colormap(jet);
% colorbar;

%% Ablation study 2 (using only tilt)

total_sep = length(tilt_ang_sep)-1;
corr_tilt = zeros(total_sep, total_sep);
N_tilt = zeros(total_sep, total_sep);

for jj = 1:length(tilt_ang_sep)-1
    ii_jj = jj;
    % Extract data in current separation range
    idx_sep = find(arm_all >= tilt_ang_sep(jj) & arm_all < tilt_ang_sep(jj+1));
    if length(idx_sep) > 4000
        s_idx = randi(length(idx_sep),round(4000*1.2),1) ;
        s_idx = shuffle(unique(s_idx)) ; % caution unique sorts the numbers
        idx_sep = idx_sep(s_idx);
    end
    sha_all_sep = sha_all(idx_sep);
    lat_sep = lat_all(idx_sep);
    lon_sep = lon_all(idx_sep);

    for jj2 = 1:length(tilt_ang_sep)-1
        ii_jj2 = jj2;
        if ii_jj2 < ii_jj
            N_tilt(ii_jj, ii_jj2) = N_tilt(ii_jj2, ii_jj);
            corr_tilt(ii_jj, ii_jj2) = corr_tilt(ii_jj2, ii_jj);
            continue
        end
        % Extract data in current separation range
        idx_sep2 = find(arm_all >= tilt_ang_sep(jj2) & arm_all < tilt_ang_sep(jj2+1));
        if length(idx_sep2) > 4000
            s_idx = randi(length(idx_sep2),round(4000*1.2),1) ;
            s_idx = shuffle(unique(s_idx)) ; % caution unique sorts the numbers
            idx_sep2 = idx_sep2(s_idx);
        end
        sha_all_sep2 = sha_all(idx_sep2);
        lat_sep2 = lat_all(idx_sep2);
        lon_sep2 = lon_all(idx_sep2);

        sha_all_sep = sort(sha_all_sep);
        sha_all_sep2 = sort(sha_all_sep2);
        min_length = min(length(sha_all_sep), length(sha_all_sep2));
        % equalize number of samples in sha_all_sep and sha_all_sep2
        % keeping the distribution same
      
        if length(sha_all_sep) ~= length(sha_all_sep2)
            if length(sha_all_sep)> length(sha_all_sep2)
                k_factor = length(sha_all_sep)/ length(sha_all_sep2);
                sha_all_sep2_resampled = interp1((1:length(sha_all_sep2))*k_factor, sha_all_sep2, 1:length(sha_all_sep));
                sha_all_sep2 = sha_all_sep2_resampled';
            else
                k_factor = length(sha_all_sep2)/ length(sha_all_sep);
                sha_all_sep_resampled = interp1((1:length(sha_all_sep))*k_factor, sha_all_sep, 1:length(sha_all_sep2));
                sha_all_sep = sha_all_sep_resampled';
            end
        end
        valid_indx = (~isnan(sha_all_sep)) & (~isnan(sha_all_sep2));
        sha_all_sep = sha_all_sep(valid_indx);
        sha_all_sep2 = sha_all_sep2(valid_indx);
        var1 = mean((sha_all_sep - mu_sha).^2);
        var2 = mean((sha_all_sep2 - mu_sha).^2);
        raw_correlation = (sha_all_sep - mu_sha) .* (sha_all_sep2 - mu_sha) ./ sqrt(var1*var2);

        corr_tilt(ii_jj, ii_jj2) = mean(raw_correlation(:), 'omitnan');%trimmean(raw_correlation(:),10);
        N_tilt(ii_jj, ii_jj2) = min_length;
    end
end

fig = figure(34);
imagesc(corr_tilt)
axis xy
colormap(jet);
colorbar;
caxis([0.8 max(corr_tilt(:))])
hcb=colorbar;
colorTitleHandle = get(hcb,'Title');
titleString = 'Correlation';
set(colorTitleHandle ,'String',titleString);
%set(gcf,'color','w')
xlabel('Pitch [degree]'); ylabel('Pitch [degree]'); zlabel('Correlation');

% exponential fit
tilt_only_sensitivity_exp = zeros(1, length(tilt_ang_sep)-1, 2);
for jj = 1:length(tilt_ang_sep)-1
    left_set = (1:jj);
    right_set = (jj:(length(tilt_ang_sep)-1));
    left_set_data = corr_tilt(jj,left_set);
    right_set_data = corr_tilt(jj,right_set);

    left_set_delta = abs(tilt_ang_sep_avg(1:jj) - tilt_ang_sep_avg(jj));
    right_set_delta = abs(tilt_ang_sep_avg(jj:end) - tilt_ang_sep_avg(jj));
    
    valid_index = ~isnan(left_set_data);
    left_set_data = left_set_data(valid_index);
    left_set_delta = left_set_delta(valid_index);
    valid_index = ~isnan(right_set_data);
    right_set_data = right_set_data(valid_index);
    right_set_delta = right_set_delta(valid_index);

    if length(left_set_data) > 1
        ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 1, 'Lower',0);
        correlation_temp = fit([left_set_delta'],left_set_data',ft,opts);
        tilt_only_sensitivity_exp(1, jj, 1) = correlation_temp.a;
    end

    if length(right_set_data) > 1
        ft = fittype( 'exp(-(a*x))', 'independent', {'x'}, 'dependent', 'y' );
        opts = fitoptions( 'Method', 'NonlinearLeastSquares', 'Upper', 1, 'Lower',0);
        correlation_temp = fit([right_set_delta'],right_set_data',ft,opts);
        tilt_only_sensitivity_exp(1, jj, 2) = correlation_temp.a;
    end
end

% figure
% imagesc(tilt_only_sensitivity_exp(:,:,1))
% colormap(jet);
% colorbar;
% 
% figure
% imagesc(tilt_only_sensitivity_exp(:,:,2))
% colormap(jet);
% colorbar;

%% save correlation profile
filename_intro = "learned_profile/correlation_profile_asilomar_afar";
save(sprintf(filename_intro +"_"+exp_txt+".mat"),'hor_dist_list',...
    'corr_list', 'mu_sha', "var_sha", 'd_int', 'correlation_f', 'tilt_ang_sep' ,...
    'elev_sep', 'corr_elev_tilt', 'N_elev_tilt', 'corr_elev', 'N_elev', 'corr_tilt', 'N_tilt',...
    'tilt_only_sensitivity_exp', 'elev_only_sensitivity_exp', 'tilt_sensitivity_exp', 'elev_sensitivity_exp')

function v=shuffle(v)
v=v(randperm(length(v)));
end
