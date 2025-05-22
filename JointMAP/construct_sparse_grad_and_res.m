function [grad_res, res] = construct_sparse_grad_and_res(num_ori, R_hat_Gauss, g, sigma_acc, sigma_mag, D_hat_mat, m_hat, acc_bias_hat, mag_bias_hat, dip_angle_hat, B_array, deltaR_ij, dRij_db, gyro_bias_est_err, y_acc_preInt, y_mag_preInt)
    % Calculate dimensions
    rows = 9*num_ori - 3;
    cols = 3*num_ori + 9 + 9 + 1;
    
    % Initialize res vector
    res = zeros(rows, 1);
    
    % Pre-allocate arrays for row indices, column indices, and values
    % Estimate the number of non-zero elements based on the structure
    nnz_per_ori_regular = 6*3 + 3*3 + 3*9 + 3*3 + 3*1;  % Regular blocks for each orientation
    nnz_per_ori_extra = 3*3 + 3*3 + 3*3;                % Extra blocks for j < num_ori
    max_nnz = num_ori * nnz_per_ori_regular + (num_ori-1) * nnz_per_ori_extra;
    
    row_indices = zeros(max_nnz, 1);
    col_indices = zeros(max_nnz, 1);
    values = zeros(max_nnz, 1);
    
    % Counter for the current position in our sparse arrays
    count = 1;
    
    for j = 1:num_ori
        % Block 1: rows 9j-8:9j-3, columns 3j-2:3j
        block1 = -[1/sigma_acc * vect2skew(-R_hat_Gauss(:,:,j)'*-g)*-R_hat_Gauss(:,:,j)'; ...
                  1/sigma_mag * D_hat_mat*vect2skew(-R_hat_Gauss(:,:,j)'*m_hat)*-R_hat_Gauss(:,:,j)'];
        
        for r = 1:6
            for c = 1:3
                row_indices(count) = (9*j-8) + (r-1);
                col_indices(count) = (3*j-2) + (c-1);
                values(count) = block1(r, c);
                count = count + 1;
            end
        end
        
        % Block 2: rows 9j-8:9j-6, columns 3*num_ori+9+1:3*num_ori+9+3
        block2 = -1/sigma_acc * eye(3);
        
        for r = 1:3
            for c = 1:3
                row_indices(count) = (9*j-8) + (r-1);
                col_indices(count) = (3*num_ori+9+1) + (c-1);
                values(count) = block2(r, c);
                count = count + 1;
            end
        end
        
        % Block 3: rows 9j-5:9j-3, columns 3*num_ori+1:3*num_ori+9
        block3 = -1/sigma_mag * kron((R_hat_Gauss(:,:,j)'*m_hat)', eye(3));
        
        for r = 1:3
            for c = 1:9
                row_indices(count) = (9*j-5) + (r-1);
                col_indices(count) = (3*num_ori+1) + (c-1);
                values(count) = block3(r, c);
                count = count + 1;
            end
        end
        
        % Block 4: rows 9j-5:9j-3, columns 3*num_ori+9+4:3*num_ori+9+6
        block4 = -1/sigma_mag * eye(3);
        
        for r = 1:3
            for c = 1:3
                row_indices(count) = (9*j-5) + (r-1);
                col_indices(count) = (3*num_ori+9+4) + (c-1);
                values(count) = block4(r, c);
                count = count + 1;
            end
        end
        
        % Block 5: rows 9j-5:9j-3, last column
        block5 = -1/sigma_mag * D_hat_mat * R_hat_Gauss(:,:,j)' * [0; -sin(dip_angle_hat); -cos(dip_angle_hat)];
        
        for r = 1:3
            row_indices(count) = (9*j-5) + (r-1);
            col_indices(count) = cols;  % Last column
            values(count) = block5(r);
            count = count + 1;
        end
        
        % Calculate residual vector 'res' for the current orientation
        res(9*j-8:9*j-3) = [1/sigma_acc*(y_acc_preInt(:,j) - R_hat_Gauss(:,:,j)' * -g - acc_bias_hat);...
                           1/sigma_mag*(y_mag_preInt(:,j) - D_hat_mat * R_hat_Gauss(:,:,j)' * m_hat - mag_bias_hat)];
        
        % Blocks for j < num_ori
        if j < num_ori
            B = squeeze(B_array(:,:,j));
            temp = log_SO3((deltaR_ij(:,:,j) * SO3_exp(dRij_db(:,:,j) * gyro_bias_est_err))' * R_hat_Gauss(:,:,j)'*R_hat_Gauss(:,:,j+1));
            
            % Calculate residual for the rotation constraint
            res(9*j-2:9*j) = B\temp;
            
            % Block 6: rows 9j-2:9j, columns 3j-2:3j
            block6 = -B\(invRightSO3Jacobian(temp)* R_hat_Gauss(:,:,j+1)');
            
            for r = 1:3
                for c = 1:3
                    row_indices(count) = (9*j-2) + (r-1);
                    col_indices(count) = (3*j-2) + (c-1);
                    values(count) = block6(r, c);
                    count = count + 1;
                end
            end
            
            % Block 7: rows 9j-2:9j, columns 3j+1:3j+3
            block7 = B\(invRightSO3Jacobian(temp)* R_hat_Gauss(:,:,j+1)');
            
            for r = 1:3
                for c = 1:3
                    row_indices(count) = (9*j-2) + (r-1);
                    col_indices(count) = (3*j+1) + (c-1);
                    values(count) = block7(r, c);
                    count = count + 1;
                end
            end
            
            % Block 8: rows 9j-2:9j, columns 3*num_ori+9+7:3*num_ori+9+9
            block8 = B\(-invRightSO3Jacobian(temp)* SO3_exp(temp)' * rightSO3Jacobian(dRij_db(:,:,j) * gyro_bias_est_err) * dRij_db(:,:,j));
            
            for r = 1:3
                for c = 1:3
                    row_indices(count) = (9*j-2) + (r-1);
                    col_indices(count) = (3*num_ori+9+7) + (c-1);
                    values(count) = block8(r, c);
                    count = count + 1;
                end
            end
        end
    end
    
    % Trim arrays to actual size used
    row_indices = row_indices(1:count-1);
    col_indices = col_indices(1:count-1);
    values = values(1:count-1);
    
    % Create sparse matrix
    grad_res = sparse(row_indices, col_indices, values, rows, cols);
end
