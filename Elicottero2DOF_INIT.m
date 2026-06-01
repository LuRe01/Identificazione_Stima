%% FILE DI INIZIALIZZAZIONE PER STIMA ELICOTTERO 2DOF
% Autori: Edoardo Di Scalzo, Luigi Realacci, Francesco Bonciani
clear 
close all
clc

%-------------------------------------------------------
% Parametri del modello
%-------------------------------------------------------

% Parametri Inerziali e di massa
params.m = 0.2;             % [kg] Massa del sistema 
params.J_alpha = 0.012;     % [kg*m^2] Mom. di Inerzia rispetto all'asse di pitch 
params.J_y = 2.3e-4;        % [kg*m^2] Componente di Inerzia Trasversale 
params.J_z = 36.4e-4;       % [kg*m^2] Componente di Inerzia Verticale 
params.I_b = 2.3e-4;        % [kg*m^2] Componente di Inerzia della base 


% Parametri geometrici e meccanici
params.g = 9.81;            % [m/s^2] Accelerazione di gravità 
params.l = 0.2;             % [m] Distanza tra perno di rotazione e applicazione forze rotori 
params.c_alpha = 0.01;      % [N*m*s/rad] Coefficiente di attrito viscoso asse di pitch 
params.c_beta = 0.01;       % [N*m*s/rad] Coefficiente di attrito viscoso asse di yaw

% Parametri Aerodinamici
params.epsilon_p = 0.1;     % Coefficiente di cross-thrust del rotore di coda sull'asse del pitch
params.epsilon_y = 0.1;     % Coefficiente di cross-thrust del rotore principale sull'asse dello yaw

% Sample time
params.st = 0.00435;        % Pari al sample time del sensore piu veloce (accelerometro)

%-------------------------------------------------------
% Condizioni Iniziali
%-------------------------------------------------------

alpha_0 = 0;            % [rad] Angolo di pitch iniziale
alpha_dot_0 = 0;        % [rad/s] Derivata angolo di pitch iniziale
beta_0 = 0;             % [rad] Angolo di yaw iniziale
beta_dot_0 = 0;         % [rad/s] Derivata angolo di yaw iniziale


% Vettori condizioni iniziali (Vettori colonna!)
q_0 = [alpha_0 beta_0]';                 % Vettore condizioi iniziali per alpha e beta
q_dot_0 = [alpha_dot_0 beta_dot_0]';     % Vettore condizioni iniziali per alpha_dot e beta_dot

%-------------------------------------------------------
% Input del sistema (forze F1 ed F2)
%-------------------------------------------------------
F1 = 0.5;       % [N] Forza del rotore principale
F2 = 0;         % [N] Forza del rotore di coda

%% SENSORI: Accelerometro e Magnetometro (Vectornav vn-100)
%-------------------------------------------------------
% Parametri Sensori
%-------------------------------------------------------
seme_rumore = 12345;        % Seme per rumore additivo dei sensori

% Accelerometro
sigma_acc = 0.0137;         % [m/s^2] Dev. Std. Accelerometro
st_acc = 0.00435;           % Sample time Accelerometro (230 Hz)

% Magnetometro
M_0 = 46.5;                 % [uT] Intensità campo magnetico locale (Es. Livorno)
sigma_magn = 0.14;          % [muT] Dev. Std. Magnetometro
st_magn = 0.005;            % Sample time Magnetometro (200 Hz)

%% Parametri EKF

%-------------------------------------------------------
% Condizioni iniziali
%-------------------------------------------------------
P_0_alpha = deg2rad(5);
P_0_alpha_dot = 0.20;
P_0_beta = deg2rad(8);
P_0_beta_dot = 0.20;

x_0 = [0 0 0 0]';           % Stato iniziale

% Matrice delle covarianze iniziale
P0 = diag([P_0_alpha P_0_alpha_dot P_0_beta P_0_beta_dot]);


% Covarianza rumore di processo
% params.sigma_alpha = 1e-5;                  % Dev. Std. alpha
% params.sigma_alpha_dot = 1e-4;              % Dev. Std. alpha_dot
% params.sigma_beta = 1e-5;                   % Dev. Std. beta
% params.sigma_beta_dot = 1e-4;               % Dev. Std. beta_dot
params.sigma_F = 0.05;                        % [N] Dev.Std. Ingresso

params.R_acc = sigma_acc^2;
params.R_magn = sigma_magn^2;

% Flag per attivare/Disattivare Distanza di Mahalanobis
MAHALANOBIS = false;


% Chiamata al simulatore su Simulink
sim('Elicottero2DOF_Sim.slx');