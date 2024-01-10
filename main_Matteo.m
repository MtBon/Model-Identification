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

sym = sim('Simulator_Single_Axis');

%% Delete temporary files

if exist('slprj','dir')
    rmdir('slprj', 's')                                                    
end

%%
% Representation of variables
figure(1)
plot(sym.q)
legend('Pitch Rate')
ylabel('[Rad/s]')
grid on;

figure(2)
plot(sym.ax)
legend('Longitudinal Acceleration');
ylabel('[m/s^2]')
grid on;

%%
% Input
u = sym.Mtot;

% Output
y = [sym.q sym.ax];

data = iddata(y, u, sample_time);
options = ssregestOptions('ARXorder',[0 1 3 1; 0 4 3 1])
%x = [sym.vx(1:length(sym.q)) sym.q sym.theta];

sys_closed_loop = ssregest(data,3,options);
A_closed_loop = sys_closed_loop.A;
B_closed_loop = sys_closed_loop.B;
C_closed_loop = sys_closed_loop.C;
D_closed_loop = sys_closed_loop.D;
K = sys_closed_loop.K;
eig(A_closed_loop-K*C_closed_loop)


%% Grey Box Estimation



% data.InputName = 'Total Pitching Moment';
% data.InputUnit = 'N/m';
% data.OutputName = {'Pitch Rate', 'Longitudinal Acceleration'};
% data.OutputUnit = {'rad/s', 'm/s^2'};
% data.Tstart = t(1);
% data.TimeUnit = 's';

order = [2 1 3];
parameters_guess = [Xu; Xq; Mu; Mq; Xd; Md];
init_sys = idgrey(@Dynamics,parameters_guess,'c',order);

% Perform the parameter estimation
options = greyestOptions('Display', 'on'); 
estimated_sys = greyest(data, init_sys, options);

% Validate the estimated model
compare(data, estimated_sys);

%% Functions

function [A, B, C, D] = Dynamics(parameters,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This Function represent the Longitudinal Dynamics of the Quadrotor used
% for the model identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xu = parameters(1);
Xq = parameters(2);
Mu = parameters(3);
Mq = parameters(4);
Xd = parameters(5);
Md = parameters(6);

A = [Xu, Xq, -9.81; 
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


%% END OF CODE