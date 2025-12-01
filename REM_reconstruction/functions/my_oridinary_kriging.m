function [est_tar, est_variance, num_points] = my_oridinary_kriging(meas_correlation_all,...
    target_cross_correlation_all, index_valid_radi, meas_sha, var_sha)

    % input
    % meas_correlation_all (num_meas, num_meas)
    % target_cross_correlation_all (num_tar, num_meas)
    % index_valid_radi (num_tar, num_meas)

    % output
    % Estimated value using Kriging
    % Expected variance of estimation
    % number of training samples used for the prediction


    num_test_points = size(target_cross_correlation_all,1);

    est_tar = ones(num_test_points,1) * NaN;
    est_variance = ones(num_test_points,1) * NaN;
    num_points = zeros(num_test_points,1);

    added_term = 2 * var_sha;

    for tar_idx = 1 : num_test_points

        index_activated = squeeze(index_valid_radi(tar_idx, :));
        num_meas_valid = sum(index_activated>0);

        if num_meas_valid > 0
            meas_sha_activated = meas_sha(index_activated); % column matrix
    
            left_mat = ones(num_meas_valid + 1);
            left_mat(num_meas_valid + 1, num_meas_valid + 1) = 0;
            % semivariogram
            left_mat(1:num_meas_valid, 1:num_meas_valid) = added_term + var_sha - meas_correlation_all(index_activated,index_activated);
    
            tar_semivariogram = var_sha - target_cross_correlation_all(tar_idx,index_activated)'; % column matrix
    
            right_mat = [tar_semivariogram + added_term; 1];
            lammda0 = left_mat \ right_mat ; %left_mat \ right_mat ;
            lammda = lammda0(1:end-1) ;
            %lagrange_multi = lammda0(end) ; % lagrange multiplier
        
            est_tar(tar_idx, 1) = real(lammda')*meas_sha_activated ;
%             disp('weights')
%             fprintf('%.2f ', lammda')
%             fprintf('\n')
            est_variance(tar_idx, 1) = real(lammda0')*[tar_semivariogram; 1] ;
            num_points(tar_idx, 1) = num_meas_valid;

            %disp(est_tar)

            % if abs(est_tar(tar_idx, 1)) > 100
            %     figure
            %     corrs = squeeze(target_cross_correlation_all(tar_idx,index_activated))/var_sha;
            %     [~,srt_idx] = sort(corrs, 'ascend');
            %     sorted_meass = meas_sha_activated(srt_idx);
            %     sorted_lamda = lammda(srt_idx);
            %     plot(corrs(srt_idx), sorted_meass, '-o')
            %     xlabel('correlation')
            %     ylabel('RSRP')
            %     figure
            %     plot(corrs(srt_idx), sorted_lamda, '-o')
            %     xlabel('correlation')
            %     ylabel('Weight')
            %     res_out = cumsum(sorted_meass.*sorted_lamda);
            %     %disp(sorted_meass.*sorted_lamda)
            %     figure
            %     plot(corrs(srt_idx), res_out, '-o')
            %     xlabel('correlation')
            %     ylabel('Cum Result')
            %     %dsdg = jgioerg + 1;
            % 
            %     break
            % end
        end
    end
end



