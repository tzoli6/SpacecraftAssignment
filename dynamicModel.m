
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

    J = diag([2500, 2300, 3000]);
    n = sqrt(398600/(6378 + 700)^3);
     
    d_theta = 1/cos(theta(2))*C_k(theta)*omega + n/cos(theta(2)).*t_k(theta);
    d_omega = eye(size(J))/J * (-skew(omega)*J*omega + 3*n^2.*C_d(theta)*J*t_d(theta) - u);

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
    d_omega = eye(size(J))/J * (-skew(omega)*J*omega + 3*n^2.*C_d(theta)*J*t_d(theta) - u);

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

    % Gyro model, nonzero mean noise, additive
    theta_m = theta + v(1:3);

    % General attitude senor model, zero mean noise, additive
    omega_m = omega + v(4:6);
 

    z = [theta_m; omega_m];

end

function [A, B, K] = linearise_system(x0, u0)

    syms x [6, 1] real;
    syms u [3, 1] real;

    f = dynamic_model(0, x, u);

    A_sym = jacobian(f, x);
    B_sym = jacobian(f, u);

    A = double(subs(A_sym, [x; u], [x0; u0]))
    B = double(subs(B_sym, [x; u], [x0; u0]));

    model = ss(A, B, eye(6, 6), []);

    G = tf(model);
    G.InputName = {'T_x', 'T_y', 'T_z'};
    G.OutputName = {'\theta_1', '\theta_2', '\theta_3', '\omega_1', '\omega_2', '\omega_3'};

    K = place(A, B, [-1 -2 -3 -4 -5 -6]);

    aug_model = ss(A-B*K, B, eye(6, 6), []);

    impulse(aug_model);
    
end

function [x, t] = simulate_system(x0, u, dt)
    % This function simulates the system for a given time step
    % x0 - initial state vector
    % u - control input
    % dt - time step

    % Initialize state vector
    x = x0;

    e_i = 0;
    prev_error = zeros(6, 1);

    K = [-4.7010    0.0150   -1.3278   -2.2329    0.0091 -0.2894;
    0.1349   -0.4609    0.2190    0.0307   -0.6905  0.0492;
   -0.9835    0.0097   -6.0186   -0.2255    0.0059 -2.7198];



    function d_x = dynamic_model_wrapper(t, x)
        dt = 0.01;
        r = zeros(6, 1);
        Kp = zeros(3, 6);
        Kd = zeros(3, 6);
        Ki = zeros(3, 6);
        Ki(1, 1) = 1000;
        Ki(2, 2) = 1000;
        Ki(3, 3) = 1000;

        e = r - x;
        e_d = (e-prev_error)/dt;
        e_i = e_i + e*dt;

        u = Kp*e + Ki*e_i + Kd*e_d;

        prev_error = e;
        
        u = -100.*K*x;
            
        d_x = dynamic_model(x, u);
    end

    [x, t] = ode45(@dynamic_model_wrapper, [0, 30000*dt], x0);

    plot(x, t)

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

     K = [-4.7010    0.0150   -1.3278   -2.2329    0.0091 -0.2894;
    0.1349   -0.4609    0.2190    0.0307   -0.6905  0.0492;
-0.9835    0.0097   -6.0186   -0.2255    0.0059 -2.7198];

    syms x1 x2 x3 x4 x5 x6 w1 w2 w3 w4 w5 w6 v1 v2 v3 v4 v5 v6 real
    x = [x1; x2; x3; x4; x5; x6];
    w = [w1; w2; w3; w4; w5; w6];
    v = [v1; v2; v3; v4; v5; v6];
    u = sym('u', [3 1], 'real');
    t = sym('t');

    dx = dynamic_model_noise(t, x, u, w);  % Symbolic dynamics model
    z = measurement_model(t, x, v);  % Symbolic measurement model

    F_x_sym = jacobian(dx, x);            % Symbolic Jacobian wrt x
    F_w_sym = jacobian(dx, w);            % Symbolic Jacobian wrt w

    H_x_sym = jacobian(z, x);            % Symbolic Jacobian wrt x
    H_v_sym = jacobian(z, v);            % Symbolic Jacobian wrt w

    % Dynamics Jacobian
    F_x = matlabFunction(F_x_sym, 'Vars', {t, x, u, w});
    F_w = matlabFunction(F_w_sym, 'Vars', {t, x, u, w});

    % Measurement Jacobian
    H_x = matlabFunction(H_x_sym, 'Vars', {t, x, v});
    H_v = matlabFunction(H_v_sym, 'Vars', {t, x, v});

    epsilon = 1e-10;
    epsilon2 = 1e-1;
    m = 0;  % Number of states
    n = 1e-1;  % Number of control inputs
    % Process noise covariance (of w)
    Q = diag([epsilon*m, epsilon*m, epsilon*m, epsilon*m, epsilon*m, epsilon*m]);

    % Measurement noise covariance (of v)
    R = diag([epsilon2*n, epsilon2*n, epsilon2*n, epsilon2*n, epsilon2*n, epsilon2*n]);
    b = [0; 0; 0; 0.2; -0.2; 0.15];

    P_kk =  0.1*eye(6, 6);  % Initial covariance estimate

    % Initialize state vector
    x = x0;

    e_i = 0;
    prev_error = zeros(6, 1);

    K = [-4.7010    0.0150   -1.3278   -2.2329    0.0091 -0.2894;
    0.1349   -0.4609    0.2190    0.0307   -0.6905  0.0492;
-0.9835    0.0097   -6.0186   -0.2255    0.0059 -2.7198];

    x_real = [x0];
    x_kk = x0*2;  % Initial estimate
    x_predicted = [x_kk]; %#ok<*NBRAK2>
    x_measured =[x_kk];
    t_real = [0];
    

    for i = 1:1600
        % Calculate control input
        u_kk = -1000.*K*x_kk;

        % Simulate the system for one time step
        [t, x] = ode45(@(t, x) dynamic_model(t, x, u_kk), ...
                        [t_real(end), t_real(end)+dt], ...
                         x_real(:, end));
        x = x(end, :)';
        x_real = [x_real x];
        t = t(end);
        t_real = [t_real t];

        % Simulate measurement with noise
        v = (randn(6, 1) .* sqrt([epsilon2, epsilon2, epsilon2, epsilon2, epsilon2, epsilon2]')) + b;
        z_k1 = measurement_model(t, x, v);
        x_measured = [x_measured z_k1]; %#ok<*AGROW>

        % Get best estimate from EKF
        [x_kk, P_kk] = ekf(x_kk, u_kk, z_k1, dt, P_kk, F_x, F_w, H_x, H_v, Q, R);
        x_predicted = [x_predicted x_kk];

    end
    figure;
    state_names = {'\theta_1', '\theta_2', '\theta_3', '\omega_1', '\omega_2', '\omega_3'};
    for i = 1:6
        if i <= 3
            subplot(3,2, (i-1)*2+1);
        else
            subplot(3,2, (i-4)*2+2);
        end
        plot(t_real, x_real(i, :), ...
            t_real, x_predicted(i, :), ...
            t_real(1:5:end), x_measured(i, 1:5:end), '--');
        legend('Real', 'EKF Predicted', 'Measurement');
        grid on;
        xlabel('Time [s]');
        if i <= 3
            ylabel([state_names{i} ' [rad]']);
        else
            ylabel([state_names{i} ' [rad/s]']);
        end
        grid on;
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
    x_k1k = x_ode(end, :)'
    % Predict covariance estimate
    Fx_k = dt.*F_x(0, x_k1k, u_k1, zeros(6, 1)) + eye(6, 6);
    Fw_k = dt.*F_w(0, x_k1k, u_k1, zeros(6, 1)) + eye(6, 6);
    % Fw_k = F_w(0, x_k1k, u_k1, zeros(6, 1));
    P_k1k = Fx_k * P_kk * Fx_k' + Q;
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

    H_x_k1 = dt * H_x(0, x_k1k, zeros(6, 1)) + eye(6, 6);
    H_w_k1 = dt * H_v(0, x_k1k, zeros(6, 1)) + eye(6, 6);
    S_k1 = H_x_k1 * P_k1k * H_x_k1' +  H_w_k1 * R * H_w_k1'; 
    K_k1 = P_k1k * H_x_k1' * pinv(S_k1);
    x_k1k1= x_k1k + K_k1 * (z_k1 - z_k1k);
    %P_k1k1 = (eye(size(K_k1*H_x_k1)) - K_k1*H_x_k1)*P_k1k;
    tmp = K_k1*H_x_k1;
    I = eye(size(tmp));
    P_k1k1 = (I - tmp) * P_k1k * (I - tmp)' + K_k1 * R * K_k1';

end

angle_0 = 10/180*pi;  % Initial angle in radians
simulate_system_ekf([angle_0; angle_0; angle_0; 0.01; 0.01; 0.01], [0; 0; 0], 0.01);

% linearise_system([0; 0; 0; 0; 0; 0], [0; 0; 0]);