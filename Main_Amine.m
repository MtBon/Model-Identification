%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ANT-X SIMULATOR - MAIN                                                  %
% Authors:  Mattia Giurato (mattia.giurato@polimi.it)                     %
%           Paolo Gattazzo (paolo.gattazzo@polimi.it)                     %
% Date: 13/12/2017                                                        %
% Adapted to ANT-X 2DoF by:  Salvatore Meraglia (salvatore.meraglia@polimi.it)%
% Date: 22/12/2022                                                        %
%
% Further modified to include structure three-state identified longitudinal model
% 06/01/23 ML
% 
% Further modified to pick outputs with measurement error
% 03/01/24 ML
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars;
close all;
addpath('datasets','common','common/simulator-toolbox','common/simulator-toolbox/attitude_library','common/simulator-toolbox/trajectory_library');
clc;

%% Model parameters

% Initial model (state: longitudinal velocity, pitch rate, pitch angle; input: normalised pitching moment; outputs: state and longitudinal acceleration)

Xu=-0.1068;

Xq=0.1192;

Mu=-5.9755;

Mq=-2.6478;

Xd=-10.1647;

Md=450.71;

A=[Xu, Xq, -9.81; Mu, Mq, 0; 0, 1, 0];

B=[Xd; Md; 0];

C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu, Xq, 0]; 

D=[0; 0; 0; Xd];

% Noise

%noise.Enabler = 0;
noise.Enabler = 1;

noise.pos_stand_dev = noise.Enabler * 0.0011;                            	%[m]

noise.vel_stand_dev = noise.Enabler * 0.01;                               %[m/s]

noise.attitude_stand_dev = noise.Enabler * deg2rad(0.33);                 %[rad]
noise.ang_rate_stand_dev = noise.Enabler * deg2rad(1);                   %[rad/s]

% Delays

delay.position_filter = 1;
delay.attitude_filter = 1;
delay.mixer = 1;

%% Load controller parameters

parameters_controller                    

%% M injection example (sweeep: first column time vector, secondo column time history of pitching moment) 

load ExcitationM

SetPoint=[0,0];

%% Values selected

t=ExcitationM(:,1);

simulation_time=t(end)-t(1);

%% Launch SIMULATOR

simout= sim('Simulator_Single_Axis.slx'); 
sim_time = 0:sample_time:simulation_time; % Time vector from Simulink output

%% Prepare the data for identification

% %filtering
% order = 2;
% fs=1/sample_time;
% filtered_Mtot = BWfilter(simout.Mtot,order,fs);
% filtered_q = BWfilter(simout.q,order,fs);
% filtered_ax = BWfilter(simout.ax,order,fs);
% 
% inputs = filtered_Mtot;
% outputs = [filtered_q filtered_ax];

outputs = [simout.q simout.ax];
inputs=simout.Mtot;

data = iddata(outputs, inputs, sample_time, 'Name', 'QuadRotor');
data.InputName = 'Total Pitching Moment';
data.InputUnit = 'N/m';
data.OutputName = {'Pitch Rate', 'Longitudinal Acceleration'};
data.OutputUnit = {'rad/s', 'm/s^2'};
data.Tstart = t(1);
data.TimeUnit = 's';

%% Delete temporary files

if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end
%% Frequency domain identification

data_f = fft(data);
parameters_guess = [0,0,0,0,0,0];
init_sys = idgrey(@Dynamics,parameters_guess,'c');

% Perform the parameter estimation
options = greyestOptions('Display', 'on'); 
estimated_sys = greyest(data_f, init_sys, options);


Xu_est = estimated_sys.A(1,1);
Xq_est = estimated_sys.A(1,2);
Mu_est = estimated_sys.A(2,1);
Mq_est = estimated_sys.A(2,2);
Xd_est = estimated_sys.B(1);
Md_est = estimated_sys.B(2);

paramCovariance = getcov(estimated_sys);
paramVariance = diag(paramCovariance);
totalVariance = sum(paramVariance)

Xu_std = paramVariance(1);
Xq_std = paramVariance(2);
Mu_std = paramVariance(3);
Mq_std = paramVariance(4);
Xd_std = paramVariance(5);
Md_std = paramVariance(6);

A=[Xu_est, Xq_est, -9.81; Mu_est, Mq_est, 0; 0, 1, 0];
B=[Xd_est; Md_est; 0];
C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu_est, Xq_est, 0]; 
D=[0; 0; 0; Xd_est];

%% Plot gaussians 

plot_gaussian(Xd_est,Xd_std,'Xd')
plot_gaussian(Xu_est,Xu_std,'Xu')
plot_gaussian(Xq_est,Xq_std,'Xq')
plot_gaussian(Mq_est,Mq_std,'Mq')
plot_gaussian(Mu_est,Mu_std,'Mu')
plot_gaussian(Md_est,Md_std,'Md')

%% Validate the estimated model
figure()
compare(data_f, estimated_sys);

%% Optimization
disp('Starting the Optimization')
ms = MultiStart();

% Options
opts = optimoptions('fmincon');
opts.Algorithm = 'interior-point';
opts.MaxIterations = 500;
opts.MaxFunctionEvaluations = 500;
opts.ConstraintTolerance = 1e-9;
opts.StepTolerance = 1e-9;
opts.OptimalityTolerance = 1e-9;
opts.Display = 'iter';

% Initial conditions and boundaries
theta0=[Xu_est,Xq_est,Mu_est,Mq_est,Xd_est,Md_est];
eta0=[25,50,90*rand];
lb=[0,0,0];
ub=[inf,inf,90];

%% Problem definition
objective_function = @(eta) perf_metric(eta,sample_time,theta0);
problem = createOptimProblem("fmincon", 'x0', eta0, 'objective', objective_function, 'lb', lb, 'ub', ub', 'nonlcon', @const, 'options', opts);

%% Solution
tic;
[x, fval, exitflag, output, solutions] = run(ms, problem, 10); % Run
time_of_run = toc;
%%

% optimized model
 f1=x(1);f2=x(2);T=x(3);
    t=0:sample_time:T;
    u=0.1*sin(2*pi.*(f1+(f2-f1).*t/T).*t);
    ExcitationM=[t;u]';
   
    simulation_time=T;
    sim_time = 0:sample_time:simulation_time; % Time vector from Simulink output
    simout= sim('Simulator_Single_Axis.slx'); 
    outputs = [simout.q simout.ax];
    inputs=simout.Mtot;

     data_f = fft(data);
    parameters_guess = theta0;
    opt_sys = idgrey(@Dynamics,parameters_guess,'c');
    
    % Perform the parameter estimation
    options = greyestOptions('Display', 'off'); 
    opt_sys = greyest(data_f, opt_sys, options);
    
    
    parCovariance = getcov(opt_sys);
    parVariance = diag(parCovariance);

    figure()
  compare(data_f, opt_sys);

%% MonteCarlo Simulation
disp('Starting the MonteCarlo Simulation');

n=10;
mu=theta0';
eta0=x;
sigma=parVariance;
% Set the random seed for reproducibility 
rng(0);
% Generate the random field, transform the random numbers to the desired distribution 
random_field = mu + sigma .* randn(length(mu), n);
eta_mc=zeros(3,n);
fval_mc=zeros(1,n);

for i=1:n
    theta0=random_field(:,i);
    objective_function = @(eta) perf_metric(eta,sample_time,theta0);
    [eta_mc(:,i), fval_mc(i), exitflag, output, solutions]=fmincon(objective_function,eta0,[],[],[],[],lb,ub,@const,opts);
    
    f1 = eta_mc(1,i);
    f2 = eta_mc(2,i);
    T = eta_mc(3,i);

    t=0:sample_time:T;
    u=0.1*sin(2*pi.*(f1+(f2-f1).*t/T).*t);
    ExcitationM=[t;u]';
   
    simulation_time=T;
    sim_time = 0:sample_time:simulation_time; % Time vector from Simulink output
    simout= sim('Simulator_Single_Axis.slx'); 
    outputs = [simout.q simout.ax];
    inputs=simout.Mtot;

    data_f = fft(data);
    parameters_guess = theta0;
    opt_sys = idgrey(@Dynamics,parameters_guess,'c');
    
    % Perform the parameter estimation
    options = greyestOptions('Display', 'off'); 
    sys = greyest(data_f, opt_sys, options);

    Xu_est = sys.A(1,1);
    Xq_est = sys.A(1,2);
    Mu_est = sys.A(2,1);
    Mq_est = sys.A(2,2);
    Xd_est = sys.B(1);
    Md_est = sys.B(2);

    Theta(i,:)=[Xu_est,Xq_est,Mu_est,Mq_est,Xd_est,Md_est];
    
end

for i = 1 : n

    f1 = eta_mc(1,i);
    f2 = eta_mc(2,i);
    T = eta_mc(3,i);

    for j = 1 : n

        Xu_est = Theta(j,1);
        Xq_est = Theta(j,2);
        Mu_est = Theta(j,3);
        Mq_est = Theta(j,4);
        Xd_est = Theta(j,5);
        Md_est = Theta(j,6);

        A=[Xu_est, Xq_est, -9.81; Mu_est, Mq_est, 0; 0, 1, 0];
        B=[Xd_est; Md_est; 0];
        C=[1, 0, 0; 0, 1, 0; 0, 0, 1; Xu_est, Xq_est, 0]; 
        D=[0; 0; 0; Xd_est];

        t=0:sample_time:T;
        u=0.1*sin(2*pi.*(f1+(f2-f1).*t/T).*t);
        ExcitationM=[t;u]';
       
        simulation_time=T;
        sim_time = 0:sample_time:simulation_time; % Time vector from Simulink output
        simout= sim('Simulator_Single_Axis.slx'); 
        outputs = [simout.q simout.ax];
        inputs=simout.Mtot;
    
        data_f = fft(data);
        parameters_guess = theta0;
        opt_sys = idgrey(@Dynamics,parameters_guess,'c');
        
        % Perform the parameter estimation
        options = greyestOptions('Display', 'off'); 
        opt_sys = greyest(data_f, opt_sys, options);
        
        
        parCovariance = getcov(opt_sys);
        parVariance = diag(parCovariance);
        J(j) = sum(paramVariance);

    end
    avg(i) = mean(J);
    wc(i) = max(J);

end

[best_av, ind_av] = min(avg);
[best_wc, ind_wc] = min(wc);

eta_av = eta_mc(:,ind_av);
eta_wc = eta_mc(:,ind_wc);


 %% Simulation for worst case
 f1_wc = eta_wc(1);
 f2_wc = eta_wc(2);
 T_wc =  eta_wc(3);


    t_wc=0:sample_time:T_wc;
    u_wc=0.1*sin(2*pi.*(f1_wc+(f2_wc-f1_wc).*t_wc/T_wc).*t_wc);
    ExcitationM=[t_wc;u_wc]';

   
    simulation_time=T_wc;
    sim_time = 0:sample_time:simulation_time; % Time vector from Simulink output
    simout= sim('Simulator_Single_Axis.slx'); 
    outputs = [simout.q simout.ax];
    inputs=simout.Mtot;

    data_f = fft(data);
    parameters_guess = theta0;
    opt_sys = idgrey(@Dynamics,parameters_guess,'c');
    
    % Perform the parameter estimation
    options = greyestOptions('Display', 'off'); 
    wc_sys = greyest(data_f, opt_sys, options);


    Xu_wc = wc_sys.A(1,1);
    Xq_wc = wc_sys.A(1,2);
    Mu_wc = wc_sys.A(2,1);
    Mq_wc = wc_sys.A(2,2);
    Xd_wc = wc_sys.B(1);
    Md_wc = wc_sys.B(2);

    Theta_wc = [Xu_wc,Xq_wc,Mu_wc,Mq_wc,Xd_wc,Md_wc];
        
    parCov_wc = getcov(wc_sys);
    parVar_wc = diag(parCov_wc);

    % Validate the estimated model
figure()
compare(data_f, wc_sys);

%% Simulation for Average

f1_av=eta_av(1);
 f2_av=eta_av(2);
 T_av=eta_av(3);

    t_av=0:sample_time:T_av;
    u_av=0.1*sin(2*pi.*(f1_av+(f2_av-f1_av).*t_av/T_av).*t_av);
    ExcitationM=[t_av;u_av]';


    simulation_time=T_av;
    sim_time = 0:sample_time:simulation_time; % Time vector from Simulink output
    simout= sim('Simulator_Single_Axis.slx'); 
    outputs = [simout.q simout.ax];
    inputs=simout.Mtot;

     data_f = fft(data);
    parameters_guess = theta0;
    av_sys = idgrey(@Dynamics,parameters_guess,'c');
    
    % Perform the parameter estimation
    options = greyestOptions('Display', 'off'); 
    av_sys = greyest(data_f, av_sys, options);

    Xu_av = av_sys.A(1,1);
    Xq_av = av_sys.A(1,2);
    Mu_av = av_sys.A(2,1);
    Mq_av = av_sys.A(2,2);
    Xd_av = av_sys.B(1);
    Md_av = av_sys.B(2);
    
   Theta_av = [Xu_av,Xq_av,Mu_av,Mq_av,Xd_av,Md_av];
    parCov_av = getcov(av_sys);
    parVar_av = diag(parCov_av);

    % Validate the estimated model
    figure()
    compare(data_f, av_sys);



%% Functions

function [A, B, C, D] = Dynamics(parameters,varargin)

    Xu = parameters(1);
    Xq = parameters(2);
    Mu = parameters(3);
    Mq = parameters(4);
    Xd = parameters(5);
    Md = parameters(6);
    
    A =     [Xu, Xq, -9.81; 
             Mu, Mq, 0; 
             0, 1, 0];
        
        B = [Xd; 
             Md; 
             0];
        
        C = [0, 1, 0;    
             Xu, Xq, 0]; 
        
        D = [0;          
             Xd];        

end

function plot_gaussian(mu,sigma,var)
    x=linspace(mu-5*sigma,mu+5*sigma,1000);
    y=normpdf(x,mu,sigma);
    figure
    plot(x,y)
    grid on
    xlabel(var)
    ylabel('Probability density')
    title('Normal distribution of', var)
end

function J=perf_metric(eta,sample_time,theta0)

    f1=eta(1);f2=eta(2);T=eta(3);
    t=0:sample_time:T;
    u=0.1*sin(2*pi.*(f1+(f2-f1).*t/T).*t);
    ExcitationM=[t;u]';
    assignin('base', 'ExcitationM', ExcitationM);
    simulation_time=T;
    sim_time = 0:sample_time:simulation_time; % Time vector from Simulink output
    simout= sim('Simulator_Single_Axis.slx'); 
    outputs = [simout.q simout.ax];
    inputs=simout.Mtot;
    
    data = iddata(outputs, inputs, sample_time, 'Name', 'QuadRotor');
    data.InputName = 'Total Pitching Moment';
    data.InputUnit = 'N/m';
    data.OutputName = {'Pitch Rate', 'Longitudinal Acceleration'};
    data.OutputUnit = {'rad/s', 'm/s^2'};
    data.Tstart = t(1);
    data.TimeUnit = 's';
    

    %% Delete temporary files

    if exist('slprj','dir')
        rmdir('slprj', 's')                                                    
    end
    
    
    %% Frequency domain identification
    
    data_f = fft(data);
    parameters_guess = theta0;
    init_sys = idgrey(@Dynamics,parameters_guess,'c');
    
    % Perform the parameter estimation
    options = greyestOptions('Display', 'off'); 
    estimated_sys = greyest(data_f, init_sys, options);
    
    
    paramCovariance = getcov(estimated_sys);
    paramVariance = diag(paramCovariance);
    J = sum(paramVariance);

    Xu_est = estimated_sys.A(1,1);
    Xq_est = estimated_sys.A(1,2);
    Mu_est = estimated_sys.A(2,1);
    Mq_est = estimated_sys.A(2,2);
    Xd_est = estimated_sys.B(1);
    Md_est = estimated_sys.B(2);

    aze=[Xu_est,Xq_est,Mu_est,Mq_est,Xd_est,Md_est];
    
end

function [c,ceq] = const(y)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Inequality and equality constraints
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f1 = y(1);
    f2 = y(2);
    c = f1 - f2;
    ceq = 0;
end

function [filtered_signal] = BWfilter(signal,order,fs)
    % fs = sampling frequency
    % order = order of the filter
    
    fc = 2; % cut-off frequency for all signals
    Wn = fc/(fs/2); % normalized cut-off frequency
    
    [b,a] = butter(order,Wn);
    filtered_signal = filtfilt(b,a,signal); % zero-phase digital filtering

end
%% END OF CODE