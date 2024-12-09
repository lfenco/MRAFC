function [] = var_def()

global ti dt tf radtodeg degtorad wid L_anim A_anim xdeseado Ll
global Unc m W B rho Ixx Iyy Izz xpfore xpaft xpforev xpaftv yp zp
global xb yb zb rh lh Ah Cdh Cfh Cbh Vh ma_h span chord Ast Cdst ma_st
global rph lph Aph Cdph Cbph Vph ma_ph xG yG zG
global Xud Yvd Yrd Zwd Zqd Kpd Mwd Mqd Nrd Nvd
global Xuu Yvv Zww Kpp Mqq Mww Nrr Nvv
global Tforce M_RB M_A LT LT6 vel velz zr
global Mb GG2 GGadd m2 W2 B2 Xud2 M_RBm M_Am Mbm GG2m
global angl xp_anim yp_anim zp_anim wp_anim wn dseta Am Bm Cm Dm
global cpath cpd clqr cmrac cmrafc cmodelmrac IAE_pd IAE_lqr IAE_mrac IAE_mracf
ti = 0;   dt = 0.005;   tf = 125;
radtodeg = 180/pi;
degtorad = pi/180;
wid =1.3; % Ancho de la linea en plot

%Parametros animacion
L_anim = 0.6;%1.2;
A_anim = 0.12;%0.6;
xdeseado = 0;
Ll=0.05;

Unc = 1;
%%%%%% set values of AUV hydro coefficients %%%%%%
m= 18.04*1;         % vehicle mass, including the entrained water in kg.
W =m*9.81;          % weight (N) = mg
B = W;              % neutral buoyanct condition
rho = 1000;         % assume density of fresh water [kg/mA3]

Ixx = 0.04933*Unc;      % mass moment of inertia, in kg-mA2
Iyy = 3.97806*Unc;      % see spreadsheet for Ixx, Iyy, Izz estimates
Izz = 3.98686*Unc;

% datos para los trustes alt/fore
xpfore = .327;      % distance in x-axis from origin to foreward thruster (in)
xpaft = .510;       % from origin to aft thruster (in)

xpforev = .04;      % distance in x-axis from origin to foreward thruster (in)
xpaftv = .04;       % from origin to aft thruster (in)
% datos para los trustes port/stbd
yp = .2;            % distance in y-axis from origin to port and starboard thrusters
zp =-0.002;         % position at which port/stbd thrusts are applied w.r.t. origin

xb = 0.0053;
yb = 0;
zb = 0;             % vertical center of buoyancy for a circular cylinder main hull

% sections below input vehicle geometry for computing hydrdynamic coefficients
% the main vehicle hull
rh = 0.0635;        % outer radius of the hull in meter
lh = 1.4235;        % mean length of the hull. Account for one half of hemisphere
                    % origin is at one half of total vehicle length

Ah = 2*rh*lh;       % projected area of hull in heave and sway
Cdh = 0.82;         % cross flow: drag coefficient of cylinder with l/d ratio of 10
Cfh = 0.004;        % axial flow, wetted area: skin friction drag on cylinder with l/d = 10
Cbh = 0.2;          % axial flow: base drag coefficient (see Hoerner - Drag, pg3-12)
Vh = pi*rh^2*lh;    % volume displaced by the hull
ma_h = rho*Vh;      % added mass of circular cylinder in heave and sway

% the struts data. struts are considered as thin plates
span = 0.090;       % distance from hull to propeller housing
chord = 0.063;      % span and chord seen from bow of vehicle,
                    % nomenclature similar to lifting fin
Ast = span*chord;   % projected area of each strut in heave
Cdst = 1.0;         % drag coeffient of plate in heave
ma_st = 0.25*pi*rho*chord^2*span;   % added mass of each strut in heave,
                                    % other directions are zero

% Note: the propeller guards were removed in late Fall 98.
% So, instead of modeling propeller guards (housing), the motor casing is modeled
% port/stbd propeller motor casing
rph = 0.011;        % outer radius of motor casing
lph = 0.095;        % length of port/stbd thruster motor
Aph = 2*rph*lph;    % projected area of each housing in heave and sway
Cdph = 0.64;        % cross flow: drag coefficient of cylinder
Cbph = 0.2;         % axial flow: base drag over cylinder
                    % axial flow: skin friction drag about zero pg3-12

Vph = pi*rph^2*lph; % volume displaced by each motor housing
ma_ph = rho*Vph;    % added mass of each housing in heave and sway

% Mass RB
xG = 0.0053;        % center of mass is fwd of vehicle reference origin in x-axis (m)
yG = 0;             % assume no heeling (m)
zG = 0.0083;        % center of mass is below the origin in +ve z-axis (m)

Tforce =1.75;       % thruster force in N, usage is optional

M_RB = [m       0       0       0       m*zG    -m*yG
         0       m       0      -m*zG     0       m*xG
         0       0       m       m*yG   -m*xG     0
         0      -m*zG    m*yG     Ixx     0       0
         m*zG    0      -m*xG     0       Iyy     0
        -m*yG    m*xG     0       0       0       Izz];
      
%Added Mass
Xud = -0.4 *m;        % X u-dot added mass in surge (5% to 10% of m for streamlined bodies)
Yvd = -29;            % Y v-dot added mass in sway [kg]
Yrd = 1;              % added mass in sway due to yaw. Yrd is positive because the stern dominates
Zwd = -30;            % added mass in heave
Zqd = -0.5;           % added mass in heave due to pitch
Kpd = -1;             % added moment in roll
Mwd = Zqd;            % added moment in pitch due to heave
Mqd = -4;             % added moment in pitch
Nrd = -10;            % added moment in yaw
Nvd = Yrd;            % added moment in yaw due to sway. Nvd is positive

M_A = [ Xud     0       0       0       0       0
        0       Yvd     0       0       0       Yrd
        0       0       Zwd     0       Zqd     0
        0       0       0       Kpd     0       0
        0       0       Mwd     0       Mqd     0
        0       Nvd     0       0       0       Nrd];

%Hydrodynamic Damping
Xuu = -40*Unc;              %account for friction and axial drag
Yvv = -200*Unc;             % viscous damping coefficient in sway 
Zww = -250*Unc;             % viscous damping coefficient in heave
Kpp = -0.07*Unc;            % roll damping
Mqq = -70*Unc;              % damping moment in pitch
Mww = -20*Unc;              % Pitch damping due to heave motion
Nrr = -55*Unc;              % damping moment in yaw
Nvv = 14*Unc;               % yaw damping due to sway motion. Nvv is positive.

%Torque Forces for 4 thruster
LT = [1       1       0       0
      0       0       0       0
      0       0       1       1
      0       0       0       0
      zp      zp  -xpfore    xpaft
      yp     -yp      0       0];

%Torque Forces for 6 thruster 
LT6 = [1       1       0       0        0       0
       0       0       1       1        0       0
       0       0       0       0        1       1
       0       0     xpforev  -xpaftv   0       0
      zp      zp       0       0      -xpfore xpaft
      yp     -yp     xpforev  -xpaftv   0       0];
  
vel = 0.08;
velz = 0.08;
zr = 1;

%============================ MODEL PLANT ================================
% Mass
Mb = (M_RB - M_A);

% Restoring Forces
%-----------------
GG2 = [0   0   0       0              (W-B)        0
       0   0   0    -(W-B)              0          0
       0   0   0       0                0          0
       0   0   0    (zG*W-zb*B)         0          0
       0   0   0       0            (zG*W-zb*B)    0
       0   0   0   -(xG*W-xb*B)    -(yG*W-yb*B)    0];

GGadd = [ 0
          0
       -(W-B)
     -(yG*W-yb*B)
      (xG*W-xb*B)
          0];

%===================== PARAMETERS FOR REFERENCE MODEL ====================      
m2= 18.04;          % vehicle mass, including the entrained water in kg.
W2 =m2*9.81;          % weight (N) = mg
B2 = W2;             % neutral buoyanct condition
Xud2 = -0.4 *m2;

M_RBm = [m2       0       0       0       m2*zG    -m2*yG
         0       m2       0      -m2*zG     0       m2*xG
         0       0       m2       m2*yG   -m2*xG     0
         0      -m2*zG    m2*yG     Ixx     0       0
         m2*zG    0      -m2*xG     0       Iyy     0
        -m2*yG    m2*xG     0       0       0       Izz];

M_Am = [ Xud2     0       0       0       0       0
        0       Yvd     0       0       0       Yrd
        0       0       Zwd     0       Zqd     0
        0       0       0       Kpd     0       0
        0       0       Mwd     0       Mqd     0
        0       Nvd     0       0       0       Nrd];
Mbm = (M_RBm - M_Am);

GG2m = [0   0   0       0              (W2-B2)        0
        0   0   0    -(W2-B2)              0          0
        0   0   0       0                0          0
        0   0   0    (zG*W2-zb*B2)         0          0
        0   0   0       0            (zG*W2-zb*B2)    0
        0   0   0   -(xG*W2-xb*B2)    -(yG*W2-yb*B2)    0];
    
%============================== ANIMATION FRAME ==========================
    angl=degtorad*90*2;
    %Draw frame
    xp_anim = [ -0.5  5.5  5.5  -0.5  -0.5 ]';
    yp_anim = [ -3.0 -3.0  3.0   3.0  -3.0 ]';
    zp_anim = [ -0.5  5.5 ]';
    wp_anim = [ xdeseado   xdeseado ]';
 
%======================= MODEL REFERENCE LTI FOR MRAC ====================
    wn = 5; %natural frequency of the system
    dseta = 1; %damping factor
    
    Am = [zeros(6)  eye(6)
       -(wn^2)*eye(6) -2*dseta*wn*eye(6)];    
    Bm = [zeros(6)
        (wn^2)*eye(6)];
   
    Cm = [eye(6) zeros(6)];
    Dm = zeros(6);
    
%=============================== COLOR ===================================
cpath = [0.9290 0.6940 0.1250];%NARANJA DARK
cpd = [0.4660 0.6740 0.1880];%VERDE DARK
clqr = [0 0.4470 0.7410]; %AZUL DARK
cmrac = [0.6350 0.0780 0.1840];% ROJO DARK
cmrafc = [0 0 0];%[0 0.4470 0.7410]; %AZUL DARK
cmodelmrac = [0.8500 0.3250 0.0980];%[0 1 0];% VERDE

%=========================================================================
IAE_pd = 0;
IAE_lqr = 0;
IAE_mrac = 0;
IAE_mracf = 0;

end