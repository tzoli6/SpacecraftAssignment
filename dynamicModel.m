
clc; clear; close all;
% Transformations for kinematic model
function C = C_k(theta)
    % This function computes the rotation matrix C_1 based on the angles theta
    % theta - 3x1 vector of angles

    % Rotation matrix
    C = [ cos(theta(2)) sin(theta(1))*sin(theta(2)) cos(theta(1))*sin(theta(2));
            0 cos(theta(1))*cos(theta(2)) -sin(theta(1))*cos(theta(2));
            0 sin(theta(1)) cos(theta(1))];
end

function t = t_k(theta)

    t = [sin(theta(3)); cos(theta(2))*cos(theta(3)); sin(theta(2))*sin(theta(3))];
end

function C = C_d(theta)
    % This function computes the rotation matrix C_2 based on the angles theta
    % theta - 3x1 vector of angles

    % Rotation matrix
    C = [ 0 -cos(theta(1))*cos(theta(2)) sin(theta(1))*cos(theta(2));
         cos(theta(1))*cos(theta(2)) 0 sin(theta(2));
          -sin(theta(1))*cos(theta(2)) -sin(theta(2)) 0];
end

function t = t_d(theta)

    t = [-sin(theta(2)); sin(theta(1))*cos(theta(2)); cos(theta(1))*cos(theta(2))];
end

function S = skew(v)
    % This function computes the skew-symmetric matrix of a 3x1 vector
    S = [  0    -v(3)  v(2);
          v(3)   0    -v(1);
         -v(2)  v(1)   0   ];
end

function d_x = dynamic_model(t, x, u)
    % This function computes the dynamic model of the system
    % x - state vector
    % u - control input
    
    % State variables
    theta = x(1:3);
    omega = x(4:6);
    % Disturbance 

    J = diag([2500, 2300, 3000]);
    n = sqrt(398600/(6378 + 700)^3);
    T_d = [0.001; 0.001; 0.001]; 
     
    d_theta = 1/cos(theta(2)) * C_k(theta)*omega + n/cos(theta(2)).*t_k(theta);
    d_omega = inv(J) * (-skew(omega)*J*omega + 3*n^2.*C_d(theta)*J*t_d(theta) + u + T_d);

    % Return state derivative
    d_x = [d_theta; d_omega];

end

function d_x = dynamic_model_noise(t, x, u, w)
    % This function computes the dynamic model of the system
    % x - state vector
    % u - control input
    % w - noise vector
    
    % State variables
    theta = x(1:3) + w(1:3);
    omega = x(4:6) + w(4:6);

    J = diag([2500, 2300, 3000]);
    n = sqrt(398600/(6378 + 700)^3);

    d_theta = 1/cos(theta(2))*C_k(theta)*omega + n/cos(theta(2)).*t_k(theta);
    d_omega = inv(J) * (-skew(omega)*J*omega + 3*n^2.*C_d(theta)*J*t_d(theta) + u);

    % Return state derivative
    d_x = [d_theta; d_omega];

end

function z = measurement_model(t, x, v)
    % Sensor model with simple additive Gaussian noise
    % x - state vector
    % v - measurement noise vector

    % State variables
    theta = x(1:3, 1);
    omega = x(4:6, 1);

    theta_m = theta + v(1:3);

    omega_m = omega + v(4:6);
 

    z = [theta_m; omega_m];

end

function [A, B, K] = linearise_system(x0, u0)

    syms x [6, 1] real;
    syms u [3, 1] real;

    f = dynamic_model(0, x, u);

    A_sym = jacobian(f, x);
    B_sym = jacobian(f, u);

    A = double(subs(A_sym, [x; u], [x0; u0]));
    B = double(subs(B_sym, [x; u], [x0; u0]));

    model = ss(A, B, eye(6, 6), []);

    G = tf(model);
    G.InputName = {'T_x', 'T_y', 'T_z'};
    G.OutputName = {'\theta_1', '\theta_2', '\theta_3', '\omega_1', '\omega_2', '\omega_3'};

    K = lqr(model,1e10*eye(6), eye(3));

    aug_model = ss(A-B*K, B, eye(6, 6), []);

    % impulse(aug_model);
    
end

function [] = fixed_structure_design()

    [A, B, K_init] = linearise_system(0.01*ones(6,1), zeros(3,1));
    
    C = eye(6, 6);
    G = ss(A, B, C, 0);
    eig(A-B*K_init)
    G.InputName = {'T_x', 'T_y', 'T_z'};
    G.OutputName = {'theta_1', 'theta_2', 'theta_3', 'omega_1', 'omega_2', 'omega_3'};

    % PID controllers for each axis
    C_x = tunablePID('C_x', 'pi');
    C_y = tunablePID('C_y', 'pi');
    C_z = tunablePID('C_z', 'pi');

    C_x.InputName = 'theta_1';  C_x.OutputName = 'p_T_x';
    C_y.InputName = 'theta_2';  C_y.OutputName = 'p_T_y';
    C_z.InputName = 'theta_3';  C_z.OutputName = 'p_T_z';


    p_x = tunablePID('p_x', 'p');
    p_y = tunablePID('p_y', 'p');
    p_z = tunablePID('p_z', 'p');
    p_x.InputName = 'omega_1';  p_x.OutputName = 'p_pmega_1';
    p_y.InputName = 'omega_2';  p_y.OutputName = 'p_pmega_2';
    p_z.InputName = 'omega_3';  p_z.OutputName = 'p_pmega_3';

    sum_x = sumblk('T_x = - p_T_x - p_pmega_1');
    sum_y = sumblk('T_y = - p_T_y - p_pmega_2');
    sum_z = sumblk('T_z = - p_T_z - p_pmega_3');

    % Connect all blocks
    C0 = connect(C_x, C_y, C_z, p_x, p_y, p_z, sum_x, sum_y, sum_z, {'theta_1', 'theta_2', 'theta_3', 'omega_1', 'omega_2', 'omega_3'}, {'T_x', 'T_y', 'T_z'});

    % Open in Control System Tuner
    wc = [0.6, 1];
    [G,C,gam,info] = looptune(G, C0, wc)
    %cl = connect(G, C, {}, {'theta_1', 'theta_2', 'theta_3', 'omega_1', 'omega_2', 'omega_3'});
    showTunable(C);
    tf(C)

end

function [x, t] = simulate_system(x0, u, dt)
   % This function simulates the system using the Extended Kalman Filter (EKF)
    % x0 - initial state vector
    % u - control input
    % dt - time step

    % This function simulates the system for a given time step
    % x0 - initial state vector
    % u - control input
    % dt - time step
    
    measurement_noise = [0.1, 0.1, 0.1, 0, 0, 0]';
    b = [0; 0; 0; 0.2; -0.2; 0.15];

    x_real = [x0];
    v = (randn(6, 1) .* measurement_noise) + b;
    z_k1 = measurement_model(0, x0, v);
    x_measured =[z_k1];
    t_real = [0];   
    
    % Initialize controller variables
    e_i = zeros(3, 1);

    for i = 1:3599
        % Calculate control input
        % K_p = [222; 258; 166];
        % K_i = [0.00186; 2.76e-05; 0.000945];
        % K_d = [2.67e+03; 2.46e+03; 3.21e+03];

        K_p = [187; 194; 161]*100;  
        K_i = [0.00172; 6.11e-05; 0.00094];
        K_d = [2.67e+03; 2.46e+03; 3.21e+03];

        e_i = e_i + z_k1(1:3) * dt  % Integral of error
        u_kk =  -K_p .* z_k1(1:3) - K_d .* z_k1(4:6) - K_i .* e_i; 

        % Simulate the system for one time step
        [t, x] = ode45(@(t, x) dynamic_model(t, x, u_kk), ...
                        [t_real(end), t_real(end)+dt], ...
                         x_real(:, end));
        x = x(end, :)';
        x_real = [x_real x];
        t = t(end);
        t_real = [t_real t];

        % Simulate measurement with noise
        v = (randn(6, 1) .* measurement_noise) + b;
        z_k1 = measurement_model(t, x, v);
        x_measured = [x_measured z_k1]; %#ok<*AGROW>

    end
    figure;
    set(gcf, 'Position', [100, 100, 900, 1800]); % Taller, elongated in y
    state_names = {'$\theta_1$', '$\theta_2$', '$\theta_3$', '$\omega_1$', '$\omega_2$', '$\omega_3$'};
    colors = lines(3); % Use distinguishable colors
    for i = 1:6
        subplot(6,1,i); % 6 rows, 1 column: all subplots stacked vertically
        hold on;
        if i <= 3
            % Convert angles from radians to degrees
            h1 = plot(t_real, rad2deg(x_real(i, :)), '-', 'Color', colors(1,:), 'LineWidth', 2);
            h3 = plot(t_real(1:5:end), rad2deg(x_measured(i, 1:5:end)), '--o', ...
                'Color', colors(3,:), 'MarkerSize', 3, 'LineWidth', 0.7);
        else
            h1 = plot(t_real, x_real(i, :), '-', 'Color', colors(1,:), 'LineWidth', 2);
            h3 = plot(t_real(1:5:end), x_measured(i, 1:5:end), '--o', ...
                'Color', colors(3,:), 'MarkerSize', 3, 'LineWidth', 0.7);
        end
        hold off;
        grid on;
        xlabel('Time [s]', 'FontSize', 18, 'Interpreter', 'latex');
        if i <= 3
            ylabel([state_names{i} ' [deg]'], 'FontSize', 18, 'Interpreter', 'latex');
        else
            ylabel([state_names{i} ' [rad/s]'], 'FontSize', 18, 'Interpreter', 'latex');
        end
        set(gca, 'FontSize', 11, 'LineWidth', 1.1, 'Box', 'on');
        if i == 1
            lgd = legend([h1 h3], {'Real', 'Measurement'}, ...
                'Location', 'northoutside', 'Orientation', 'horizontal', ...
                'FontSize', 12, 'Interpreter', 'latex');
            legend boxoff
        end
    end

        set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto');
        % Remove whitespace by tightly fitting the figure to the subplots
        set(gca, 'LooseInset', get(gca, 'TightInset'));
        % Use exportgraphics for tight bounding box (R2020a+)
        exportgraphics(gcf, 'noise_simulation_results.png', 'Resolution', 300, 'BackgroundColor', 'white', 'ContentType', 'image');


        % --- Performance Metrics Calculation ---
        threshold = 0.1; % Example threshold for settling (deg for angles, rad/s for rates)
        fprintf('\nPerformance Metrics:\n');
        for i = 1:6
            if i <= 3
                % Angles: convert to degrees
                y = rad2deg(x_real(i, :));
                ref = 0; % Reference is zero
            else
                y = x_real(i, :);
                ref = 0;
            end
            % Peak response and time
            [peak_val, peak_idx] = max(abs(y - ref));
            peak_time = t_real(peak_idx);

            % Settling time: time when |y-ref| stays below threshold
            above_thresh = abs(y - ref) > threshold;
            if any(~above_thresh)
                last_above = find(above_thresh, 1, 'last');
                if isempty(last_above) || last_above == length(y)
                    settling_time = NaN;
                    steady_state_val = NaN;
                else
                    settling_time = t_real(last_above+1);
                    % Calculate steady state value as average after settling time
                    idx_ss = last_above+1:length(y);
                    steady_state_val = mean(y(idx_ss));
                end
            else
                settling_time = NaN;
                steady_state_val = NaN;
            end

            fprintf('%s: Peak=%.4f at t=%.2f s, Settling Time=%.2f s, Steady State=%.4f\n', ...
                state_names{i}, peak_val, peak_time, settling_time, steady_state_val);
        end

end

function [] = simulate_system_ekf(x0, u, dt)
    % This function simulates the system using the Extended Kalman Filter (EKF)
    % x0 - initial state vector
    % u - control input
    % dt - time step

    % This function simulates the system for a given time step
    % x0 - initial state vector
    % u - control input
    % dt - time step

    syms x1 x2 x3 x4 x5 x6 w1 w2 w3 w4 w5 w6 v1 v2 v3 v4 v5 v6 real
    x = [x1; x2; x3; x4; x5; x6];
    w = [w1; w2; w3; w4; w5; w6];
    v = [v1; v2; v3; v4; v5; v6];
    u = sym('u', [3 1], 'real');
    t = sym('t');

    dx = dynamic_model_noise(t, x, u, w); 
    z = measurement_model(t, x, v);  

    F_x_sym = jacobian(dx, x);           
    F_w_sym = jacobian(dx, w);        

    H_x_sym = jacobian(z, x);          
    H_v_sym = jacobian(z, v);            

    % Dynamics Jacobian
    F_x = matlabFunction(F_x_sym, 'Vars', {t, x, u, w});
    F_w = matlabFunction(F_w_sym, 'Vars', {t, x, u, w});

    % Measurement Jacobian
    H_x = matlabFunction(H_x_sym, 'Vars', {t, x, v});
    H_v = matlabFunction(H_v_sym, 'Vars', {t, x, v});

    epsilon = 1e-10;
    epsilon2 = 0.1^2;
    m = 1;
    n = 1;
    % Process noise covariance (of w)
    Q = diag([epsilon*m, epsilon*m, epsilon*m, epsilon*m, epsilon*m, epsilon*m]);

    % Measurement noise covariance (of v)
    R = diag([epsilon2*n, epsilon2*n, epsilon2*n, epsilon2*n, epsilon2*n, epsilon2*n]);
    measurement_noise = [0.1, 0.1, 0.1, 0, 0, 0]';
    b = [0; 0; 0; 0.2; -0.2; 0.15];

    P_kk =  R;  % Initial covariance estimate

    x_real = [x0];
    x_kk = x0*2;  % Initial estimate
    x_predicted = [x_kk]; %#ok<*NBRAK2>
    x_measured =[x_kk];
    t_real = [0];
    
    % Initialize controller variables
    e_i = zeros(3, 1);

    for i = 1:3599

        % K_p = [222; 258; 166]*10;
        % K_i = [0.00186; 2.76e-05; 0.000945];
        % K_d = [2.67e+03; 2.46e+03; 3.21e+03];

        K_p = [187; 194; 161]*10;  
        K_i = [0.00172; 6.11e-05; 0.00094]*100;
        K_d = [2.67e+03; 2.46e+03; 3.21e+03];

        e_i = e_i + x_kk(1:3) * dt;  % Integral of error
        %
        u_kk =  -K_p .* x_kk(1:3) - K_d .* x_kk(4:6) - K_i .* e_i;  % Control input

        % Simulate the system for one time step
        [t, x] = ode45(@(t, x) dynamic_model(t, x, u_kk), ...
                        [t_real(end), t_real(end)+dt], ...
                         x_real(:, end));
        x = x(end, :)';
        x_real = [x_real x];
        t = t(end);
        t_real = [t_real t];

        % Simulate measurement with noise
        v = (randn(6, 1) .* measurement_noise) + b;
        z_k1 = measurement_model(t, x, v);
        x_measured = [x_measured z_k1]; %#ok<*AGROW>

        % Get best estimate from EKF
        [x_kk, P_kk] = ekf(x_kk, u_kk, z_k1, dt, P_kk, F_x, F_w, H_x, H_v, Q, R);
        x_predicted = [x_predicted x_kk];

    end
    figure;
    set(gcf, 'Position', [100, 100, 900, 1800]); % Taller, elongated in y
    state_names = {'$\theta_1$', '$\theta_2$', '$\theta_3$', '$\omega_1$', '$\omega_2$', '$\omega_3$'};
    colors = lines(3); % Use distinguishable colors
    for i = 1:6
        subplot(6,1,i); % 6 rows, 1 column: all subplots stacked vertically
        hold on;
        if i <= 3
            % Convert angles from radians to degrees
            h1 = plot(t_real, rad2deg(x_real(i, :)), '-', 'Color', colors(1,:), 'LineWidth', 2);
            h2 = plot(t_real, rad2deg(x_predicted(i, :)), '-', 'Color', colors(2,:), 'LineWidth', 2);
            h3 = plot(t_real(1:5:end), rad2deg(x_measured(i, 1:5:end)), '--o', ...
                'Color', colors(3,:), 'MarkerSize', 3, 'LineWidth', 0.7);
        else
            h1 = plot(t_real, x_real(i, :), '-', 'Color', colors(1,:), 'LineWidth', 2);
            h2 = plot(t_real, x_predicted(i, :), '-', 'Color', colors(2,:), 'LineWidth', 2);
            h3 = plot(t_real(1:5:end), x_measured(i, 1:5:end), '--o', ...
                'Color', colors(3,:), 'MarkerSize', 3, 'LineWidth', 0.7);
        end
        hold off;
        grid on;
        xlabel('Time [s]', 'FontSize', 18, 'Interpreter', 'latex');
        if i <= 3
            ylabel([state_names{i} ' [deg]'], 'FontSize', 18, 'Interpreter', 'latex');
        else
            ylabel([state_names{i} ' [rad/s]'], 'FontSize', 18, 'Interpreter', 'latex');
        end
        set(gca, 'FontSize', 11, 'LineWidth', 1.1, 'Box', 'on');
        if i == 1
            lgd = legend([h1 h2 h3], {'Real', 'EKF Predicted', 'Measurement'}, ...
                'Location', 'northoutside', 'Orientation', 'horizontal', ...
                'FontSize', 12, 'Interpreter', 'latex');
            legend boxoff
        end
    end

        set(gcf, 'Color', 'w', 'PaperPositionMode', 'auto');
        % Remove whitespace by tightly fitting the figure to the subplots
        set(gca, 'LooseInset', get(gca, 'TightInset'));
        % Use exportgraphics for tight bounding box (R2020a+)
        exportgraphics(gcf, 'ekf_simulation_results.png', 'Resolution', 300, 'BackgroundColor', 'white', 'ContentType', 'image');

    
     % --- Performance Metrics Calculation ---
    threshold = 1.0; % Example threshold for settling (deg for angles, rad/s for rates)
    fprintf('\nPerformance Metrics:\n');
    for i = 1:6
        if i <= 3
            % Angles: convert to degrees
            y = rad2deg(x_real(i, :));
            ref = 0; % Reference is zero
        else
            y = x_real(i, :);
            ref = 0;
        end
        % Peak response and time
        [peak_val, peak_idx] = max(abs(y - ref));
        peak_time = t_real(peak_idx);

        % Settling time: time when |y-ref| stays below threshold
        above_thresh = abs(y - ref) > threshold;
        if any(~above_thresh)
            last_above = find(above_thresh, 1, 'last');
            if isempty(last_above) || last_above == length(y)
                settling_time = NaN;
                steady_state_val = NaN;
            else
                settling_time = t_real(last_above+1);
                % Calculate steady state value as average after settling time
                idx_ss = last_above+1:length(y);
                steady_state_val = mean(y(idx_ss));
            end
        else
            settling_time = NaN;
            steady_state_val = NaN;
        end

        fprintf('%s: Peak=%.4f at t=%.2f s, Settling Time=%.2f s, Steady State=%.4f\n', ...
            state_names{i}, peak_val, peak_time, settling_time, steady_state_val);
    end
end

function [x_kk, P_kk] = ekf(x_kk, u_k1, z_k1, dt, P_kk, F_x, F_w, H_x, H_v, Q, R)
    % This function implements the Extended Kalman Filter (EKF) for the system
    % x_m - measurement vector
    % u_m - control input
    % dt - time step

    % PREDICTION STEP
    [x_k1k, z_k1k, P_k1k] = ekf_predict(dt, x_kk, u_k1, P_kk, F_x, F_w, Q);

    
    % CORRECTION STEP
    [x_kk, P_kk] = ekf_correction(dt, x_k1k, z_k1k, P_k1k, H_x, H_v, z_k1, R);

end

function [x_k1k, z_k1k, P_k1k] = ekf_predict(dt, x_kk, u_k1, P_kk, F_x, F_w, Q)
    % Predict next state
    [~, x_ode] = ode45(@(t, x) dynamic_model(t, x, u_k1), [0, dt], x_kk);
    x_k1k = x_ode(end, :)';
    % Predict covariance estimate
    Fx_k = dt.*F_x(0, x_k1k, u_k1, zeros(6, 1)) + eye(6, 6);
    Fw_k = dt.*F_w(0, x_k1k, u_k1, zeros(6, 1));
    P_k1k = Fx_k * P_kk * Fx_k' + Fw_k * Q * Fw_k';
    % Predict measurement
    z_k1k = measurement_model(0, x_k1k, zeros(6, 1));

end

function [x_k1k1, P_k1k1] = ekf_correction(dt, x_k1k, z_k1k, P_k1k, H_x, H_v, z_k1, R)
    % This function performs the correction step of the Extended Kalman Filter (EKF)
    % x_k1k - predicted state vector
    % z_k1k - predicted measurement vector
    % P_k1k - predicted covariance matrix
    % z_k1 - actual measurement vector
    % R - measurement noise covariance matrix

    H_x_k1 = H_x(0, x_k1k, zeros(6, 1));
    H_w_k1 = H_v(0, x_k1k, zeros(6, 1));
    S_k1 = H_x_k1 * P_k1k * H_x_k1' +  H_w_k1 * R * H_w_k1'; 
    K_k1 = P_k1k * H_x_k1' * pinv(S_k1);
    x_k1k1= x_k1k + K_k1 * (z_k1 - z_k1k);
    tmp = K_k1*H_x_k1;
    I = eye(size(tmp));
    P_k1k1 = (I - tmp) * P_k1k * (I - tmp)' + K_k1 * R * K_k1';

end

angle_0 = 10/180*pi;  % Initial angle in radians
simulate_system([angle_0; angle_0; angle_0; 0.01; 0.01; 0.01], [0; 0; 0], 0.05);
% fixed_structure_design();
% linearise_system([0; 0; 0; 0; 0; 0], [0; 0; 0]);