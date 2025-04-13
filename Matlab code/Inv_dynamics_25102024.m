clc;
clear all;
close all;
%%  inverse dyanmics calculation (to find joints torque using corresponding position, velocity, and accelaration)
% for 2- links manupulator L1=L2=L, theta_1 and theta_2
% steps stated below
% mass m1, m2
% length l
% g = 9.81 m/sec^2
% find initial and final postion of theta_1 and theta_2
% time of motion 
% trajectory (history)
% relation ship betbeen position, velocity, accelaration for each joints
% compute the polynomial theta_1 and theta_2  coefficients [a0,a1, a2, a3]
% at every instant of time (based on smapling time), find (theta, thetadot, thetadoubledot) for each joints anagles theta_1 and theta_2
% Put in dynamics equation
% tau = D(theta)*thetadoubledot+H(thetadot,theta)*thetadot+C(theta)
% plot tau1, tau2, joints torques , D,H, C, different compenents of tau for 3 different time instants.
% 




%% links mass
m1=1;   % unit kg
m2=1;  % unit kg


%% links length
l=1;  % unit length
%%gravational 
g=9.81; % m/sec^2


%% Declaring initail and final position of joint -1
theta0_1=0;
thetaf_1=40; % in deg
thetaf_1=thetaf_1*(pi/ 180);  % rad

%% Declaring initail and final position of joint -2
theta0_2=0;
thetaf_2=40; % in deg
thetaf_2=thetaf_2*(pi/180);   % rad

%%  storing the time for motion in a vector
t_f=[5,1,0.1];
%% declaring empty variables for for storing torque components due to inertia, coriolis, 
%and gravity terms and total torque
D=[];
H=[];
C=[];
tau=[];
time=[];
for i=1:length(t_f)
    % calculating the coefficients of the cubic polynomial for trajectory
    % genration for joint-1
    a0_1=theta0_1;
    a1_1=0;
    a2_1=(3/(t_f(i)^2))*(thetaf_1-theta0_1);
    a3_1=(-2/(t_f(i)^3))*(thetaf_1-theta0_1);
    % calculating the coefficients of the cubic polynomial for trajectory
    % genration for joint-2
    a0_2=theta0_2;
    a1_2=0;
    a2_2=(3/(t_f(i)^2))*(thetaf_2-theta0_2);
    a3_2=(-2/(t_f(i)^3))*(thetaf_2-theta0_2);
    for t=0:0.001:t_f(i)
        % calculating the joints position, velocity, and  accelaration for
        % joint -1
        theta1=a0_1+a1_1*t+a2_1*t^2+a3_1*t^3;
        theta1_dot=a1_1+2*a2_1*t+3*a3_1*t^2;
        theta1_dotdot=2*a2_1+6*a3_1*t;
        % for joint-2
        theta2=a0_2+a1_2*t+a2_2*t^2+a3_2*t^3;
        theta2_dot=a1_2+2*a2_2*t+3*a3_2*t^2;
        theta2_dotdot=2*a2_2+6*a3_2*t;
        
        %% calculating the components of joint torques
        D_term = [(1/3)*m1*l^2+(4/3)*m2*l^2+ m2*cos(theta2)*l^2   (1/3)*m2*l^2*cos(theta2);
                           (1/3)*m2*l^2 + (1/2)*m2*l^2*cos(theta2)   (1/3)*m2*l^2]...
                           *[theta1_dotdot; theta2_dotdot];
          H_term= [ (-1/2)*m2*sin(theta2)*l^2*theta2_dot^2-m2*sin(theta2)*l^2*theta1_dot*theta2_dot;
                                        (1/2)*m2*sin(theta2)*l^2*theta1_dot^2];
         C_term= [(1/2)*m1*g*l*cos(theta1) + (1/2)*m2*g*l*cos(theta1+theta2) + m2*g*l*cos(theta1);
                         (1/2)*m2*g*l*cos(theta1+theta2)];
                     %% calculating the joint torque
tau_current= D_term+H_term+C_term;
%% storing these for plotting later
D=[D,D_term];
H=[H,H_term];
C=[C,C_term];
tau=[tau, tau_current];
time=[time, t];
    

                                    
    end
    % plotting the results for current time of motion
    figure;
    plot (time,tau(1,:),time,tau(2,:), '--','linewidth',1.5);
    grid on;
    xlabel('time (s)');
    ylabel('joint torque (Nm)');
    title(['variation of joint torques (total time of travel = ', num2str(t_f(i)),'s)']);
    legend('joint torque-1','joint torque -2');
    saveas(gcf,sprintf('Joint_torques_Time%f.png',t_f(i)));
    figure;
    plot(time,D(1,:),time,D(2,:),'--',time,H(1,:),time,H(2,:),'--',time,C(1,:),time,C(2,:),':','linewidth',1.5);
    grid on;
    xlabel('time (s)');
    ylabel('joint torque components (Nm)');
    title(['Variation of joint torque components  (Total time of travel = ', num2str(t_f(i)), 's)']);
    legend('Inertia D_1','Inertia D_2', 'Coriolis + centrifugal H_1', 'Coriolis + centrifugal H_2','gravity C_1','gravity C_2');
    saveas (gcf,sprintf('Torque_components_time%f.png', t_f(i)));
    
    %% Reinitializing the storage variables for next run
    D=[];
    H=[];
    C=[];
    tau=[];
    time=[];
    
    
end
