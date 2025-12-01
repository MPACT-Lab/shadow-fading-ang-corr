function [rmse_meas_radius_2d, rmse_meas_radius_3d, valid_est1, valid_est2] = calc_rmse_median(error_all_2d,error_all_3d)
    rmse_iter_meas_radius_2d = zeros(size(error_all_2d,1), size(error_all_2d,2), size(error_all_2d,3));
    rmse_iter_meas_radius_3d = zeros(size(error_all_3d,1), size(error_all_3d,2), size(error_all_3d,3));


    for repeat_i = 1:size(error_all_2d,1)
        for meas_i = 1:size(error_all_2d,2)
            for radi_i = 1:size(error_all_2d,3)
                % 2D
                all_values = squeeze(error_all_2d(repeat_i,meas_i,radi_i,:));
                valid_value = ~isnan(all_values);
                valid_values = all_values(valid_value(:));
                if (sum(valid_value(:))) < 1 % minimum 10 out of num_test_points/5 (hieghts)
                    rmse_iter_meas_radius_2d(repeat_i, meas_i, radi_i) = NaN;
                else
                    rmse_iter_meas_radius_2d(repeat_i, meas_i, radi_i) = sqrt(mean(valid_values.^2));
                end
            end
        end
    end
    disp('data1')
    for repeat_i = 1:size(error_all_3d,1)
        for meas_i = 1:size(error_all_3d,2)
            for radi_i = 1:size(error_all_3d,3)
                % 3D
                all_values = squeeze(error_all_3d(repeat_i,meas_i,radi_i,:)); 
                valid_value = ~isnan(all_values);
                valid_values = all_values(valid_value(:));
                if (sum(valid_value(:))) < 1 % minimum 10 out of num_test_points/5 (hieghts)
                    rmse_iter_meas_radius_3d(repeat_i, meas_i, radi_i) = NaN;
                    disp('got nan')
                    disp(sum(valid_value(:)))
                else
                    rmse_iter_meas_radius_3d(repeat_i, meas_i, radi_i) = sqrt(mean(valid_values.^2));
                end
            end
        end
    end

    rmse_meas_radius_2d = zeros(size(error_all_2d,2), size(error_all_2d,3));
    rmse_meas_radius_3d = zeros(size(error_all_3d,2), size(error_all_3d,3));
    for meas_i = 1:size(error_all_2d,2)
        for radi_i = 1:size(error_all_2d,3)
            % 2D
            rmse_samples = squeeze(rmse_iter_meas_radius_2d(:, meas_i, radi_i));%
            rmse_samples = rmse_samples(~isnan(rmse_samples));
            %disp(length(rmse_samples))
            if length(rmse_samples) < 1
                rmse_meas_radius_2d(meas_i,radi_i) = NaN;
            else
                rmse_meas_radius_2d(meas_i,radi_i) = median(rmse_samples);
            end
        end
    end

    for meas_i = 1:size(error_all_3d,2)
        for radi_i = 1:size(error_all_3d,3)
            % 3D
            rmse_samples = squeeze(rmse_iter_meas_radius_3d(:, meas_i, radi_i));%
            rmse_samples = rmse_samples(~isnan(rmse_samples));
            %sv1 = rmse_samples;
            %disp(length(rmse_samples))
            if length(rmse_samples) < 1
                rmse_meas_radius_3d(meas_i,radi_i) = NaN;
            else
                rmse_meas_radius_3d(meas_i,radi_i) = median(rmse_samples);
                if rmse_meas_radius_3d(meas_i,radi_i)>100
                    disp(mean(rmse_samples))
                    disp(max(rmse_samples))
                    %dfsg = ohgoi + 1;
                end
            end
        end
    end

    valid_est_meas_radius_2d = zeros(size(error_all_2d,2), size(error_all_2d,3));
    valid_est_meas_radius_3d = zeros(size(error_all_3d,2), size(error_all_3d,3));

    for meas_i = 1:size(error_all_2d,2)
        for radi_i = 1:size(error_all_2d,3)
            % 2D
            all_values = squeeze(error_all_2d(:,meas_i,radi_i,:));
            valid_value = ~isnan(all_values);
            valid_est_meas_radius_2d(meas_i,radi_i) = valid_est_meas_radius_2d(meas_i,radi_i) + sum(valid_value(:));
        end
    end
    valid_est1 = valid_est_meas_radius_2d/(size(error_all_2d,1)*size(error_all_2d,4));
    for meas_i = 1:size(error_all_3d,2)
        for radi_i = 1:size(error_all_3d,3)
            % 3D
            all_values = squeeze(error_all_3d(:,meas_i,radi_i,:));
            valid_value = ~isnan(all_values);
            valid_est_meas_radius_3d(meas_i,radi_i) = valid_est_meas_radius_3d(meas_i,radi_i) + sum(valid_value(:));
        end
    end
    valid_est2 = valid_est_meas_radius_3d/(size(error_all_3d,1)*size(error_all_3d,4));
end

