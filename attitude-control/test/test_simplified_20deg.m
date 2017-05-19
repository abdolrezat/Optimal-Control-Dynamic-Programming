function test_simplified_20deg()
%% default parameters
skip_calc = 1;
simplified = 0;
q0 = [0;0;0.0174524064372835;0.999847695156391]; %angle2quat(0,0,roll = deg2rad(2))
X0 = [0;0;0;...
    ...  %quaternions equal to quat(deg2rad(-15), deg2rad(20), deg2rad(-10))
   q0];%0.999660006156261;0.00841930262082080;0.0176013597667272;0.0172968080698774];

w_min = deg2rad(-0.7);
w_max = deg2rad(0.7);
n_mesh_w = 10;
yaw_min = -5; %angles
yaw_max = 5;
pitch_min = -4;
pitch_max = 4;
roll_min = -5.5;
roll_max = 5.5;
n_mesh_t = 10;

% u = linspace(-1, 0, 1000);
% U_vector = [u(1:end-1), linspace(0, 1, 1000) ];
U_vector = [-0.01 0 0.01];
% n_mesh_u = length(U_vector);
Qw1 = 6;
Qw2 = 6;
Qw3 = 6;
Qt1 = 6;
Qt2 = 6;
Qt3 = 6;
R1 = 0.1;
R2 = 0.1;
R3 = 0.1;

T_final = 45;
h = 0.01;
N_stage = T_final/h;

if(~isinteger( N_stage))
    N_stage = ceil(N_stage);
    T_final = h*N_stage;
    warning('T_final is not a factor of h (dt), increasing T_final to %.2f\n',T_final)
end


J1 = 2;
J2 = 2.5;
J3 = 3;


%% mesh generation
s_w1 = linspace(w_min, w_max, n_mesh_w);
s_w2 = linspace(w_min, w_max, n_mesh_w);
s_w3 = linspace(w_min, w_max, n_mesh_w);

s_t1 = linspace(deg2rad(roll_min), deg2rad(roll_max), n_mesh_t);
s_t2 = linspace(deg2rad(pitch_min), deg2rad(pitch_max), n_mesh_t);
s_t3 = linspace(deg2rad(yaw_min), deg2rad(yaw_max), n_mesh_t);


%% initialization
%states, controls
[X1,X4,U1] = ndgrid(s_w1,s_t1,U_vector);
[X2,X5,U2] = ndgrid(s_w2,s_t2,U_vector);
[X3,X6,U3] = ndgrid(s_w3,s_t3,U_vector);

%u_opt
size_Umat = [n_mesh_w, n_mesh_t, N_stage];
U1_Opt = zeros(size_Umat,'single');
U2_Opt = zeros(size_Umat,'single');
U3_Opt = zeros(size_Umat,'single');

%J matrices, current and optimal next stage

J_current_1 = J_current(X1,X4,U1,Qw1,Qt1,R1);
F1 = griddedInterpolant({s_w1,s_t1},zeros(n_mesh_w,n_mesh_t),'linear');

J_current_2 = J_current(X2,X5,U2,Qw2,Qt2,R2);
F2 = griddedInterpolant({s_w2,s_t2},zeros(n_mesh_w,n_mesh_t),'linear');

J_current_3 = J_current(X3,X6,U3,Qw3,Qt3,R3);
F3 = griddedInterpolant({s_w3,s_t3},zeros(n_mesh_w,n_mesh_t),'linear');
%
J11 = J_current_1(:,:,1);
X11 = X1(:,:,1);
X44 = X4(:,:,1);
%
[w1_next,t1_next] = next_stage_states(X1, X4, U1, J1, h);
[w2_next,t2_next] = next_stage_states(X2, X5, U2, J2, h);
[w3_next,t3_next] = next_stage_states(X3, X6, U3, J3, h);
% keyboard
%% loop
if (skip_calc == 0 )
whandle = waitbar(0,'Calculation in Progress...');
for k_s = N_stage-1:-1:1
    tic
%     C = F1(w1_next,t1_next);
    [F1.Values, U1_idx] = min( J_current_1 + F1(w1_next,t1_next), [], 3);
    [F2.Values, U2_idx] = min( J_current_2 + F2(w2_next,t2_next), [], 3);
    [F3.Values, U3_idx] = min( J_current_3 + F3(w3_next,t3_next), [], 3);
    
    U1_Opt(:,:,k_s) = U_vector(U1_idx);
    U2_Opt(:,:,k_s) = U_vector(U2_idx);
    U3_Opt(:,:,k_s) = U_vector(U3_idx);
    %     QQ = U1_idx;
    waitbar( 1 - k_s/N_stage, whandle);
    fprintf('step %d - %f seconds\n', k_s, toc)
end
close(whandle)
fprintf('stage calculation complete... cleaning up...\n')
% clear J_current_1 J_current_2 J_current_3 F1 F2 F3 X1 X2 X3 X4 X5 X6 U1 U2 U3 ...
%     w1_next w2_next w3_next t1_next t2_next t3_next
fprintf('...Done!\n')
else
    load('simplified-test-20deg.mat')
    fprintf('Calculation skipped, loaded the required data\n')
end

keyboard

%% Simulate results
if(simplified == 1)
    %% test on Simplified system dynamics -----------------------------------------------------
    X_Simplified = zeros(6, N_stage);
    [t3,t2,t1] = quat2angle([X0(7),X0(6),X0(5),X0(4)]);
    X_Simplified(:,1) = [0;0;0;t1;t2;t3];
    U = zeros(3,N_stage);
    X_ANGLES = zeros(3,N_stage);
    tic
    for k_stage=1:N_stage-1
        t1 = X_Simplified(4,k_stage);
        t2 = X_Simplified(5,k_stage);
        t3 = X_Simplified(6,k_stage);
        %         t3 = 2*asin(X(4,k_stage));
        %         t2 = 2*asin(X(5,k_stage));
        %         t1 = 2*asin(X(6,k_stage));
        %     t1,t2,t3
        FU1 = griddedInterpolant({s_w1,s_t1}, single(U1_Opt(:,:,k_stage)),'nearest');
        FU2 = griddedInterpolant({s_w2,s_t2}, single(U2_Opt(:,:,k_stage)),'nearest');
        FU3 = griddedInterpolant({s_w3,s_t3}, single(U3_Opt(:,:,k_stage)),'nearest');
        
        U(1,k_stage) = FU1(X_Simplified(1,k_stage), t1);
        U(2,k_stage) = FU2(X_Simplified(2,k_stage), t2);
        U(3,k_stage) = FU3(X_Simplified(3,k_stage), t3);
        
        [w1_new,t1_new] = next_stage_states(X_Simplified(1,k_stage), t1, U(1,k_stage), J1, h);
        [w2_new,t2_new] = next_stage_states(X_Simplified(2,k_stage), t2, U(2,k_stage), J2, h);
        [w3_new,t3_new] = next_stage_states(X_Simplified(3,k_stage), t3, U(3,k_stage), J3, h);
        X_Simplified(:,k_stage+1) = [w1_new;w2_new;w3_new;t1_new;t2_new;t3_new];
        
        X_ANGLES(:,k_stage) = [rad2deg(t1);rad2deg(t2);rad2deg(t3)];
    end
    toc
    
      %% Plot Simplified System
    time_v = linspace(0, T_final, N_stage); %plot time vector
    % plot controls
    figure
    plot(time_v, U(1,:),'--')
    hold on
    plot(time_v, U(2,:),'--')
    plot(time_v, U(3,:),'--')
    grid on
    title('Control Inputs')
    legend('u1','u2','u3')
    
    %plot w
    figure
    hold on
    grid on
    for n_state = [1 2 3]
        plot(time_v, X_Simplified(n_state, :).*180/pi)
    end
    legend('\omega_1','\omega_2','\omega_3')
    xlabel('time (s)')
    ylabel('\omega (deg/sec)')
    
    figure
    hold on
    for i=1:3
        plot(time_v, X_ANGLES(i,:))
    end
    title('rotation angles')
    legend('\theta_1','\theta_2','\theta_3')
    grid on
    xlabel('time (s)')
    
    
    %% test on REAL SYSTEM DYNAMICS ----------------------------------------------------------
else
    X = zeros(7, N_stage);
    X(:,1) = X0;
    U = zeros(3,N_stage);
    X_ANGLES = zeros(9,N_stage);
    tic
    for k_stage=1:N_stage-1
%         t3 = x_roll;
%         t2 = x_pitch;
%         t1 = x_yaw;
        x_roll = 2*asin(X(4,k_stage));
        x_pitch = 2*asin(X(5,k_stage));
        x_yaw = 2*asin(X(6,k_stage));
%         
        t1 = x_roll;
        t2 = x_pitch;
        t3 = x_yaw;
        
        %     t1,t2,t3
        FU1 = griddedInterpolant({s_w1,s_t1}, single(U1_Opt(:,:,1)),'nearest');
        FU2 = griddedInterpolant({s_w2,s_t2}, single(U2_Opt(:,:,1)),'nearest');
        FU3 = griddedInterpolant({s_w3,s_t3}, single(U3_Opt(:,:,1)),'nearest');
        
        U(1,k_stage) = FU1(X(1,k_stage), t1);
        U(2,k_stage) = FU2(X(2,k_stage), t2);
        U(3,k_stage) = FU3(X(3,k_stage), t3);
        %% test on Real system dynamics

           X(:,k_stage+1) = real_system_dynamics(X(:,k_stage)', U(:,k_stage)', J1, J2, J3, h);
        [t1,t2,t3] = quat2angle([X(7,k_stage),X(6,k_stage),X(5,k_stage),X(4,k_stage)]);
    
        X_ANGLES(:,k_stage) = [X(1,k_stage);X(2,k_stage);X(3,k_stage);... % w1,w2,w3
            rad2deg(t1);rad2deg(t2);rad2deg(t3);... % angles
            U(:,k_stage)]; % controls
    end
    toc
    
     %% Plot REAL SYSTEM
    time_v = linspace(0, T_final, N_stage); %plot time vector
    % plot controls
    figure
    plot(time_v, U(1,:),'--')
    hold on
    plot(time_v, U(2,:),'--')
    plot(time_v, U(3,:),'--')
    grid on
    title('Control Inputs')
    legend('u1','u2','u3')
    
    %plot states
    figure
    hold on
    grid on
    for n_state = [1 2 3]
        plot(time_v, X(n_state, :).*180/pi)
    end
    legend('\omega_1','\omega_2','\omega_3')
    xlabel('time (s)')
    ylabel('\omega (deg/sec)')
    
    figure
    hold on
    for i=4:6
        plot(time_v, X_ANGLES(i,:))
    end
    title('rotation angles')
    legend('\theta_1','\theta_2','\theta_3')
    grid on
    xlabel('time (s)')
    
    
    figure
    hold on
    grid on
    for n_state = [4 5 6 7]
        plot(time_v, X(n_state, :))
    end
    legend('q1','q2','q3','q4')
    xlabel('time (s)')
    
end %end if

keyboard


k_s = 1:N_stage-2;
keyboard
this.X1_mesh = X1(:,:,1);
this.X2_mesh = X4(:,:,1);
%
this.u_star = U1_Opt;
%
if length(k_s) == 1
    figure
    plot3( this.X1_mesh, this.X2_mesh, this.u_star(:,:,k_s) )
else
    figure
    p = mesh(this.X1_mesh, this.X2_mesh, this.u_star(:,:,k_s(1)) );
    colormap winter
    % not allowing axis limits to change automatically
    axis manual
    for i=2:length(k_s)
        k_temp = k_s(i);
        p.ZData =  this.u_star(:,:,k_s(i));
        title(['Stage ',num2str(k_temp)]);
        pause(0.001)
    end
end
keyboard

function J_ = J_current(w,theta,U,Qw,Qt,R) 
J_ = (Qw * w.^2 + Qt * theta.^2 + R*U.^2);

function [w_next,t_next] = next_stage_states(W, T, U, J, h)
w_next = RK4_w(W, U, J, h);
t_next = RK4_t(T, W, h);


function W2 = RK4_w(w, U, J, h)
%calculates next stage (k+1) states
% X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
%first order taylor expansion
%X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);

% Runge-Kutta - 4th order
% h = dt;
k1 = wdynamics(w , U, J);
k2 = wdynamics((w + k1*h/2), U, J);
k3 = wdynamics((w + k2*h/2), U, J);
k4 = wdynamics((w + k3*h), U, J);

W2 = w + h*(k1 + 2*k2 + 2*k3 + k4)/6;


function T2 = RK4_t(T1, W1, h)
%calculates next stage (k+1) states
% X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
%first order taylor expansion
%X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);

% Runge-Kutta - 4th order
% h = dt;
k1 = tdynamics(W1);
k2 = tdynamics((W1 + k1*h/2));
k3 = tdynamics((W1 + k2*h/2));
k4 = tdynamics((W1 + k3*h));

T2 = T1 + h*(k1 + 2*k2 + 2*k3 + k4)/6;


function w_dot = wdynamics(~, U, J)
w_dot = U/J;

function t_dot = tdynamics(w)
t_dot = w;


function X2 = real_system_dynamics(X1, U, J1, J2, J3, h)
%calculates next stage (k+1) states
% X2 = X1 + dt*a_d(X,u) where a_d is the spacecraft dynamics
%first order taylor expansion
%X2 = X1 + dt*spacecraft_dynamics(spacecraft, X1, U);

% Runge-Kutta - 4th order
% h = dt;
k1 = spacecraft_dynamics_list( X1 , U, J1, J2, J3);
k2 = spacecraft_dynamics_list((X1 + k1*h/2), U, J1, J2, J3);
k3 = spacecraft_dynamics_list((X1 + k2*h/2), U, J1, J2, J3);
k4 = spacecraft_dynamics_list( (X1 + k3*h), U, J1, J2, J3);

X2 = X1 + h*(k1 + 2*k2 + 2*k3 + k4)/6;
q_sqsum = sqrt(X2(4)^2 +  X2(5)^2 +  X2(6)^2 +  X2(7)^2);
X2(4:7) = X2(4:7)/q_sqsum;


function X_dot = spacecraft_dynamics_list(X, U, J1, J2, J3)
%returns the derivatives x_dot = f(X,u)
% FIXED q~ convention
X_dot = zeros(size(X));
x1 = X(:,1);
x2 = X(:,2);
x3 = X(:,3);
x4 = X(:,4);
x5 = X(:,5);
x6 = X(:,6);
x7 = X(:,7);
u1 = U(:,1);
u2 = U(:,2);
u3 = U(:,3);
%J is assumed to be diagonal, J12 = J23 = ... = 0
X_dot(:,1) = (J2-J3)/J1*x2.*x3 + u1/J1;
X_dot(:,2) = (J3-J1)/J2*x3.*x1 + u2/J2;
X_dot(:,3) = (J1-J2)/J3*x1.*x2 + u3/J3;
X_dot(:,4) = 0.5*(x3.*x5 -x2.*x6 +x1.*x7);
X_dot(:,5) = 0.5*(-x3.*x4 +x1.*x6 +x2.*x7);
X_dot(:,6) = 0.5*(x2.*x4 -x1.*x5 +x3.*x7);
X_dot(:,7) = 0.5*(-x1.*x4 -x2.*x5 -x3.*x6);
