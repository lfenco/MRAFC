clear; close all; clc;

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
global ntdd ruido EnbNoise EnbTimeBreak ltdd tt3 timebreak
global CntPoints ti2 dt2 tf2 t2 nn wxdd wydd wzdd wv waypoints depth velz
global tdd waypoints2 waypoints3 xdd ydd zdd xdd3 ydd3 zdd3
global qeta qnu QQ RR PP KK xp up xm sp sm err Dh tau Tee Txx Tuu K_adap
global AA BB CC EE Ap Bp Cp Ep D Af Bf Cf Df Ke Kx Ku Ta Tb Ki
global flag xmf Kd_mrafc Kp_mrafc KpKd_mracf Ljastt EGV sngl1
global loaddataM AA2 BB2 Plmi2 Kjastt
global degtorad xx6c fdpertBNpsi fdpertApsi fdpertBpsi fdpertCpsi
global kX1 kX2 kX3 kX4 kX5 kX6 kX1p kX2p kX3p kX4p kX5p kX6p pX pXp pXt
global sig1Npsi sig2Npsi sig3Npsi sig4Npsi sig5Npsi sig6Npsi sig0psi
global sig1psi sig2psi sig3psi sig4psi sig5psi sig6psi
global l EGV sngl1 mu_j Bm2 Plmi2 errx xxf xxfo sum_mu_j rr Kjdot Ljdot Kj Lj dt

disp(' ')
disp(' ==== MODEL AND CONTROLLER AUV ===== ')
disp(' == (By: Fenco Bravo Lugui Paolo) == ')
disp(' ===== (Update: 07-Dec-2024) ===== ')
disp(' ')

%====================== DEFINITION OF VARIABLES ==========================
var_def;
EnbNoise = input('input enable noise [Not noise = 0 -- noise = 1] : ');
dist_amp = input('input disturbance amplitude [0.8] : ');
EnbTimeBreak = input('input enable break [Not Break = 0 -- Break = 1]: ');
if (EnbTimeBreak)
    timebreak = input('input time break: ');
end
CntPoints = input('input Cnt Points: [8Points = 0 -- 16Points = 1] : ');
saturation = input('not saturation=0 -- saturation=1: '); 
loaddataM = input('model reference: new calculation = 0 -- load data = 1: '); 
regulation = input('input [regulation  = 0 -- tracking = 1] : ');
Figure1 = 0;
NewFuncMemb = 1;

%============================== PATH TRACKING  ===========================
trackpath;
if (Figure1 == 1)
    figure(1);
    subplot(3,1,1); plot(tdd,xr,'m'); xlabel('time'); ylabel('x(m)');
    subplot(3,1,2); plot(tdd,yr,'b'); xlabel('time'); ylabel('y(m)');
    subplot(3,1,3); plot(xr,yr,'r'); axis('equal'); xlabel('x(m)'); ylabel('y(m)');
    suptitle('nonparametric desired path');
end

    %=====================================================================
    % Calculate transformation matrix T and use it to calculate xdot, ydot, and zdot
    nsl=4;% number subspace linear
    for ii = 1:nsl
        if (ii == 1)
            cphi = cos(0) ;
            ctheta = cos(0);
            cpsi = cos(pi) ;
            ttheta = tan(0);
            sectheta = sec(0);
            sphi = sin(0) ;
            stheta = sin(0);
            spsi = sin(pi);
        elseif(ii == 2)
            cphi = cos(0) ;
            ctheta = cos(0);
            cpsi = cos(pi/2);
            ttheta = tan(0);
            sectheta = sec(0);
            sphi = sin(0) ;
            stheta = sin(0);
            spsi = sin(pi/2);
        elseif(ii == 3)
            cphi = cos(0) ;
            ctheta = cos(0);
            cpsi = cos(0);
            ttheta = tan(0);
            sectheta = sec(0);
            sphi = sin(0) ;
            stheta = sin(0);
            spsi = sin(0);
        elseif(ii == 4)
            cphi = cos(0) ;
            ctheta = cos(0);
            cpsi = cos(-pi/2);
            ttheta = tan(0);
            sectheta = sec(0);
            sphi = sin(0) ;
            stheta = sin(0);
            spsi = sin(-pi/2);
        end
        
        R_Theta = zeros(3,3);
        R_Theta(1,1) = cpsi*ctheta;
        R_Theta(1,2) = cpsi*stheta*sphi - spsi*cphi;
        R_Theta(1,3) = spsi*sphi + cpsi*stheta*cphi;
        R_Theta(2,1) = spsi*ctheta;
        R_Theta(2,2) = spsi*stheta*sphi + cpsi*cphi;
        R_Theta(2,3) = -cpsi*sphi + spsi*stheta*cphi;
        R_Theta(3,1) = -stheta;
        R_Theta(3,2) = ctheta*sphi;
        R_Theta(3,3) = ctheta*cphi; 

        T_Theta = zeros(3,3);
        T_Theta(1,1) = 1;
        T_Theta(1,2) = sphi*ttheta;
        T_Theta(1,3) = cphi*ttheta;
        T_Theta(2,1) = 0;
        T_Theta(2,2) = cphi;
        T_Theta(2,3) = -sphi;
        T_Theta(3,1) = 0;
        T_Theta(3,2) = sphi*sectheta;
        T_Theta(3,3) = cphi*sectheta;

        J_Theta(:,:,ii) = [R_Theta zeros(3)
                    zeros(3) T_Theta];
    end
    
    %=====================================================================      
    for ii = 1:nsl
        %======= MODEL LINEAL AUV ==========
        AA(:,:,ii) = [ zeros(6)  J_Theta(:,:,ii)
                    -(Mb)\GG2 zeros(6) ];
        %======= MODEL FOR REFERENCE =======        
        AA2(:,:,ii) = [ zeros(6)  J_Theta(:,:,ii)
                        -(Mbm)\GG2m zeros(6) ];
    end               
    
    BB = [ zeros(6)
           inv(Mb) ];
    BB2 = [ zeros(6)
           inv(Mbm) ];
    EE = [ zeros(6)
           inv(Mb) ];
    
    CC = [eye(6) zeros(6)];
    DD = zeros(6);
    
    initial_conditions; % initial conditions   
    lmi_calculation; %lmi calculation
    
    %=========================== MODEL LINEAL T-S =======================
    for ii = 1:nsl
        for jj = 1:nsl            
            Am2(:,:,ii,jj) = AA2(:,:,ii)-BB2*Kjastt(:,:,jj);
        end            
    end
    Bm2 = BB2*Ljastt;
    
    %==================== PLOT MEMBERSHIP FUNCTION ===================
    membership_function;
    Kj = zeros(6,12,pXt); %Inicializa Matriz Kj(i,j,h)
    Lj = zeros(6,6,pXt); %Inicializa Matriz Lj(i,j,h)
    
    for gg = 1:pXt
      Kj(:,:,gg) =  Kjastt(:,:,gg);
      Lj(:,:,gg) =  Ljastt;    
    end

%========================= CONDITION OF EXTERNAL FORCES ==================
Tau_ext = zeros(6,1);
fext = dist_amp.*ones(6,1);
wext = degtorad*20;

xini = [0; 0; 0; 0; 0; 0; 2.5; -2.5; 0; 0; 0; angl];
T = zeros(6) ;
TT = zeros(6);

phi = xini(10) ;
the = xini(11) ;
psi = xini(12) ;
cphi = cos(phi);
cthe = cos(the);
cpsi = cos(psi);
sphi = sin(phi);
sthe = sin(the);
spsi = sin(psi);

T(1,1) = cpsi*cthe;
T(1,2) = cpsi*sthe*sphi - spsi*cphi;
T(1,3) = spsi*sphi + cpsi*sthe*cphi;
T(2,1) = spsi*cthe ;
T(2,2) = spsi*sthe*sphi + cpsi*cphi;
T(2,3) = -cpsi*cphi + spsi*sthe*cphi;
T(3,1) = -sthe ;
T(3,2) = cthe*sphi ;
T(3,3) = cthe*cphi ;

TT(1:3,1:3) =T(1:3,1:3)';
TT(4,4) = 1;
TT(4,5) = 0;
TT(4,6) = -sthe;
TT(5,4) = 0;
TT(5,5) = cphi;
TT(5,6) = sphi*cthe ;
TT(6,4) = 0 ;
TT(6,5) = -sphi ;
TT(6,6) = cphi*cthe ;

%============= VARIABLE DE ESTADO LQR ===========
%nuLqr = [u; v; w; p; q; r];
nuLqr = TT*xini(1:6);
%etaLqr = [x; y; z; phi; theta; psi];
etaLqr = [xini(7); xini(8); xini(9); xini(10); xini(11); xini(12)];
%tauLqr = [Tport; Tstbd; Tfore; Taft];
tauLqr4thr = [0; 0; 0; 0];
%tauLqr = [Tport; Tstbd; Tfore_h; Taft_h; Tfore_v; Taft_v];
tauLqr6thr = [0; 0; 0; 0; 0; 0];

%============= VARIABLE DE ESTADO MRAC ===========
%tauMrac = [Tport; Tstbd; Tfore; Taft];
tauAdp = [0; 0; 0; 0];
%============= VARIABLE DE ESTADO MARCF ===========
%nuMrafc = [u; v; w; p; q; r];
nuMrafc = TT*xini(1:6);
%etaMrafc = [x; y; z; phi; theta; psi];
etaMrafc = [xini(7); xini(8); xini(9); xini(10); xini(11); xini(12)];
%tauMrafc4thr = [Tport; Tstbd; Tfore; Taft];
tauMrafc4thr = [0; 0; 0; 0];
%tauMrafc6thr = [Tport; Tstbd; Tfore_h; Taft_h; Tfore_v; Taft_v];
tauMrafc6thr = [0; 0; 0; 0; 0; 0];

disp('pause, press any key to continue.. ')
pause

% referencias
r11 = 2.0; r21 = -2.0; r31 = 1.0; r41 = degtorad*90*2;
r12 = 1.2; r22 = -1.2; r32 = 2.0; r42 = degtorad*90*2;
r13 = 0.6; r23 = 0.8; r33 = 3.0; r43 = degtorad*55*2;
%============================= MAIN LOOP =================================
xr2 = zeros(ltdd,1);
yr2 = zeros(ltdd,1);
xdd2 = zeros(ltdd,1);
ydd2 = zeros(ltdd,1);
zdd2 = zeros(ltdd,1);
   
k = 1;
for t = ti:dt:tf

    if(t>=0 && t<=tf/4)
     ref1 = r11; ref2 = r21; ref3 = r31; ref4 = r41;  
    elseif(t>tf/3 && t<=2*tf/3)
     ref1 = r12; ref2 = r22; ref3 = r32; ref4 = r42;
    elseif(t>2*tf/3 && t<=tf)
     ref1 = r13; ref2 = r23; ref3 = r33; ref4 = r43;
    end
    
    xdd2(k,1)=xdd(k,1);
    ydd2(k,1)=ydd(k,1);
    zdd2(k,1)=zdd(k,1);
 
    if ((t > 11.2) &  (t < 20.2))
        %External Forces
        Tau_ext = fext*cos(2*wext*t);
    else
        Tau_ext = zeros(6,1)*cos(2*wext*t);       
    end
    
%============================= MODEL LQR ================================  
    xxlqr = [etaLqr; nuLqr];

    x1 = etaLqr(1,1) + ruido(k,1);%==== x
    y1 = etaLqr(2,1) + ruido(k,1);%==== y
    z1 = etaLqr(3,1) + ruido(k,1);%==== z
    phi1 = etaLqr(4,1) + ruido(k,1);%==== phi
    theta1 = etaLqr(5,1) + ruido(k,1);%==== theta
    psi1 = etaLqr(6,1) + ruido(k,1);%==== psi
    u1 = nuLqr(1,1);
    v1 = nuLqr(2,1);
    w1 = nuLqr(3,1);
    p1 = nuLqr(4,1);
    q1 = nuLqr(5,1);
    r1 = nuLqr(6,1);

    psiastLqr = atan2(2.5-x1,y1);
    psiast2Lqr(k,1)=radtodeg*psiastLqr;
    if (regulation == 0)
        xxlqro = [ref1;ref2;ref3;0;0;ref4;zeros(6,1)];
    else
        xxlqro = [ xdd(k,1); ydd(k,1); zdd(k,1); 0; 0; psiastLqr; zeros(6,1)];
    end
    errxLqr = xxlqr - xxlqro;
    TorqLqr = -KK*errxLqr;
    
    %Performance index
    error_lqr(k,:) = abs(errxLqr);
    IAE_lqr = IAE_lqr + error_lqr(k,:);

    tauLqr6thr = LT6\TorqLqr;
    if (saturation == 1)
       %Limit for tauLqr6
        if(tauLqr6thr(1,1) > Tforce)
            tauLqr6thr(1,1) = Tforce;
        elseif(tauLqr6thr(1,1) < -Tforce)
            tauLqr6thr(1,1) = -Tforce;
        end   
        if(tauLqr6thr(2,1) > Tforce)
            tauLqr6thr(2,1) = Tforce;
        elseif(tauLqr6thr(2,1) < -Tforce)
            tauLqr6thr(2,1) = -Tforce;
        end    
        if(tauLqr6thr(3,1) > Tforce)
            tauLqr6thr(3,1) = Tforce;
        elseif(tauLqr6thr(3,1) < -Tforce)
            tauLqr6thr(3,1) = -Tforce;
        end

        if(tauLqr6thr(4,1) > Tforce)
            tauLqr6thr(4,1) = Tforce;
        elseif(tauLqr6thr(4,1) < -Tforce)
            tauLqr6thr(4,1) = -Tforce;
        end
        if(tauLqr6thr(5,1) > Tforce)
            tauLqr6thr(5,1) = Tforce;
        elseif(tauLqr6thr(5,1) < -Tforce)
            tauLqr6thr(5,1) = -Tforce;
        end
        if(tauLqr6thr(6,1) > Tforce)
            tauLqr6thr(6,1) = Tforce;
        elseif(tauLqr6thr(6,1) < -Tforce)
            tauLqr6thr(6,1) = -Tforce;
        end
    end
    
    % Calculate transformation matrix T and use it to calculate xdot, ydot, and zdot 
    cphi1 = cos(phi1) ;
    ctheta1 = cos(theta1);
    cpsi1 = cos(psi1) ;
    ttheta1 = tan(theta1);
    sectheta1 = sec(theta1);
    sphi1 = sin(phi1) ;
    stheta1 = sin(theta1);
    spsi1 = sin(psi1);
       
    R_Theta = zeros(3,3);
    R_Theta(1,1) = cpsi1*ctheta1;
    R_Theta(1,2) = cpsi1*stheta1*sphi1 - spsi1*cphi1;
    R_Theta(1,3) = spsi1*sphi1 + cpsi1*stheta1*cphi1;
    R_Theta(2,1) = spsi1*ctheta1;
    R_Theta(2,2) = spsi1*stheta1*sphi1 + cpsi1*cphi1;
    R_Theta(2,3) = -cpsi1*sphi1 + spsi1*stheta1*cphi1;
    R_Theta(3,1) = -stheta1;
    R_Theta(3,2) = ctheta1*sphi1;
    R_Theta(3,3) = ctheta1*cphi1; 

    T_Theta = zeros(3,3);
    T_Theta(1,1) = 1;
    T_Theta(1,2) = sphi1*ttheta1;
    T_Theta(1,3) = cphi1*ttheta1;
    T_Theta(2,1) = 0;
    T_Theta(2,2) = cphi1;
    T_Theta(2,3) = -sphi1;
    T_Theta(3,1) = 0;
    T_Theta(3,2) = sphi1*sectheta1;
    T_Theta(3,3) = cphi1*sectheta1;

    J_Theta = [R_Theta zeros(3)
                zeros(3) T_Theta];

    %========================== MATRIX CALCULATION =======================
    %Coriolis
    C_RB = [ 0          0               0        m*zG*r1        -m*(xG*q1-w1)     -m*(xG*r1+v1)
             0          0               0       -m*w1            m*(zG*r1+xG*p1)   m*u1
             0          0               0       -m*(zG*p1-v1)    -m*(zG*q1+u1)      m*xG*p1
            -m*zG*r1    0               m*zG*p1     0               Izz*r1        -Iyy*q1
             m*xG*q1   -m*(zG*r1+xG*p1) m*zG*q1    -Izz*r1           0             Ixx*p1 
             m*xG*r1    0              -m*xG*p1     Iyy*q1          -Ixx*p1          0];

    %Coriolis Mass Add
    C_A = [  0                  0               0               0             -(Zwd*w1+Zqd*q1)  (Yvd*v1+Yrd*r1)
             0                  0               0              (Zwd*w1+Zqd*q1)       0             -Xud*u1
             0                  0               0             -(Yvd*v1+Yrd*r1)      Xud*u1            0
             0            -(Zwd*w1+Zqd*q1) (Yvd*v1+Yrd*r1)   0                   -(Nvd*v1+Nrd*r1)   (Mwd*w1+Mqd*q1)
            (Zwd*w1+Zqd*q1)       0              -Xud*u1          (Nvd*v1+Nrd*r1)       0              -Kpd*p1 
           -(Yvd*v1+Yrd*r1)     Xud*u1             0             -(Mwd*w1+Mqd*q1)      Kpd*p1            0];   

    %Hydrodynamic Damping Dn(Vr) (Parte no lineal)    
    Dn = [ -Xuu*abs(u1)  0           0           0           0           0
            0           -Yvv*abs(v1) 0           0           0           0
            0           0           -Zww*abs(w1) 0           0           0
            0           0           0           -Kpp*abs(p1) 0           0
            0           0           -Mww*abs(w1) 0           -Mqq*abs(q1) 0
            0           -Nvv*abs(v1) 0           0           0           -Nrr*abs(r1)];     

    %Restoring Forces
    gn = [(W-B)*sin(theta1)
          -(W-B)*cos(theta1)*sin(phi1)
          -(W-B)*cos(theta1)*cos(phi1)
          (zG*W-zb*B)*cos(theta1)*sin(phi1)
          (zG*W-zb*B)*sin(theta1)-xG*(W-B)*cos(theta1)*cos(phi1)
          -xG*(W-B)*cos(theta1)*sin(phi1)];

    %========================= Non-linear equation =======================
    etapLqr = J_Theta*nuLqr;
    nupLqr = Mb\(-(C_RB + C_A + Dn)*nuLqr - gn + LT6*tauLqr6thr + Tau_ext);

    %=====================================================================
    etaLqr = etaLqr + etapLqr*dt;
    nuLqr = nuLqr + nupLqr*dt;

    eta2Lqr(1,1) = etaLqr(1,1);
    eta2Lqr(2,1) = etaLqr(2,1);
    eta2Lqr(3,1) = etaLqr(3,1);
    eta2Lqr(4,1) = radtodeg*etaLqr(4,1);
    eta2Lqr(5,1) = radtodeg*etaLqr(5,1);
    eta2Lqr(6,1) = radtodeg*etaLqr(6,1);

    % Save for Plot
    Tau_ext2(k,:)= Tau_ext';
    xEuLqrNu(k,:) = nuLqr';
    xEuLqrEta(k,:) = etaLqr';
    xEuLqrEta2(k,:) = eta2Lqr';
    uLqr2(k,:) = tauLqr6thr';
    errxEuLqr(k,:) = errxLqr';

% ============================ MRAC adaptive model =======================   
    x7 = xm(1,1);%==== x
    x8 = xm(2,1);%==== y 
    
    psiast = atan2(2.5-x7,x8);
    psiast2mrac(k,1)=radtodeg*psiast;
    if (regulation == 0)
        um = [ref1;ref2;ref3;0;0;ref4];
        xxmarco = [ref1;ref2;ref3;0;0;ref4;zeros(6,1)];
    else
        um = [xdd(k,1); ydd(k,1);zdd(k,1);0;0;psiast];
        xxmarco = [xdd(k,1); ydd(k,1); zdd(k,1); 0; 0; psiast;zeros(6,1)];
    end
    %Performance index
    error_mrac(k,:) = abs(err);
    IAE_mrac = IAE_mrac + error_mrac(k,:);
    
    %PLANT
    xdot = Ap*xp + Bp*up + Ep*Tau_ext;
    yp = Cp*xp + ruido(k,1).*ones(6,1);

    %FEEDFORWARD
    sdot = Af*sp + Bf*up;
    rp = Cf*sp + Df*up;

    %Output
    zp = yp + rp; 

    %Model
    xmdot = Am*xm + Bm*um;
    ym = Cm*xm + Dm*um;

    if (t==0)
        disp('argumento solo planta');
    end
    %Control Error
    err = ym-zp;

    %Control Law
    rT = [err' xm' um'];

    Ki = Ki + err*rT*Ta*dt;
    Kpe = err*rT*Tb;
    Kt = Kpe + Ki;
    up = Kt*rT';

    tauAdp6 = LT6\up ;

    %Limit for tauLqr6
    if (saturation == 1)
        if(tauAdp6(1,1) > Tforce)
            tauAdp6(1,1) = Tforce;
        elseif(tauAdp6(1,1) < -Tforce)
            tauAdp6(1,1) = -Tforce;
        end   
        if(tauAdp6(2,1) > Tforce)
            tauAdp6(2,1) = Tforce;
        elseif(tauAdp6(2,1) < -Tforce)
            tauAdp6(2,1) = -Tforce;
        end

        if(tauAdp6(3,1) > Tforce)
            tauAdp6(3,1) = Tforce;
        elseif(tauAdp6(3,1) < -Tforce)
            tauAdp6(3,1) = -Tforce;
        end

        if(tauAdp6(4,1) > Tforce)
            tauAdp6(4,1) = Tforce;
        elseif(tauAdp6(4,1) < -Tforce)
            tauAdp6(4,1) = -Tforce;
        end
        if(tauAdp6(5,1) > Tforce)
            tauAdp6(5,1) = Tforce;
        elseif(tauAdp6(5,1) < -Tforce)
            tauAdp6(5,1) = -Tforce;
        end
        if(tauAdp6(6,1) > Tforce)
            tauAdp6(6,1) = Tforce;
        elseif(tauAdp6(6,1) < -Tforce)
            tauAdp6(6,1) = -Tforce;
        end

        up = LT6*tauAdp6;
    end
    
    %Euler
    xp = xp + xdot*dt;
    sp = sp + sdot*dt;
    xm = xm + xmdot*dt;

    % Save fot Plot
    Tau_ext2(k,:)= Tau_ext';
    %==============
    xmAdp(k,:) = xm;
    xpAdp(k,:) = xp;
       
    ymAdp(k,:) = ym;
    zpAdp(k,:) = zp;
    upAdp2(k,:) = tauAdp6';    
    errAdp(k,:) = err;

% ========================== MRAFC adaptive model ========================
    if (t==0)
     length(t)
    end

    xxf = [etaMrafc; nuMrafc];

    u3 = nuMrafc(1,1);
    v3 = nuMrafc(2,1);
    w3 = nuMrafc(3,1);
    p3 = nuMrafc(4,1);
    q3 = nuMrafc(5,1);
    r3 = nuMrafc(6,1);
    x3 = etaMrafc(1,1) + ruido(k,1);%==== x
    y3 = etaMrafc(2,1) + ruido(k,1);%==== y
    z3 = etaMrafc(3,1) + ruido(k,1);%==== z
    phi3 = etaMrafc(4,1) + ruido(k,1);%==== phi
    theta3 = etaMrafc(5,1) + ruido(k,1);%==== theta
    psi3 = etaMrafc(6,1) + ruido(k,1);%==== psi
    
    u3m = xmf(7,1);%==== u
    v3m = xmf(8,1);%==== v
    w3m = xmf(9,1);%==== w
    p3m = xmf(10,1);%==== p
    q3m = xmf(11,1);%==== q
    r3m = xmf(12,1);%==== r
    x3m = xmf(1,1);%==== x
    y3m = xmf(2,1);%==== y
    z3m = xmf(3,1);%==== z
    phi3m = xmf(4,1);%==== phi
    theta3m = xmf(5,1);%==== theta
    psi3m = xmf(6,1);%==== psi
    
    psiast = atan2(2.5-x3m,y3m);
    if (regulation == 0)
        xxfo = [ref1;ref2;ref3;0;0;ref4;zeros(6,1)];
    else
        xxfo = [xdd(k,1); ydd(k,1);zdd(k,1);0;0;psiast;zeros(6,1)];
    end
    errd = xmf - xxfo;
    rr = -KpKd_mracf*errd; % PD
    
    % =========================== MODEL =========================         
    errx = xxf-xmf; % CONTROL ERROR (X_planta - X_modelo)

    % ====================== Performance index error =====================  
    err_iae = xxf -[xdd(k,1); ydd(k,1);zdd(k,1);0;0;psiast;zeros(6,1)];
        
    % Performance index
    error_mrafc(k,:) = abs(err_iae);
    IAE_mracf = IAE_mracf + error_mrafc(k,:);   
    
    if (NewFuncMemb == 1)
        psi3abs = psi3;
        
       %============== RANGE <-180 & >180 ===============
       if(psi3abs <= sig6Npsi)
          fpx6(1,1) = 1;
       elseif((psi3abs > sig6Npsi) && (psi3abs<= sig4Npsi))
          fpx6(1,1) = (psi3abs - sig4Npsi)/(sig6Npsi - sig4Npsi);
       elseif((psi3abs > sig4psi) && (psi3abs<= sig6psi))
          fpx6(1,1) = (psi3abs - sig4psi)/(sig6psi - sig4psi);
       elseif(psi3abs > sig6psi)
          fpx6(1,1) = 1;
       else
          fpx6(1,1) = 0;
       end
       
    
       fdpertDDpsi(k,1) = fpx6(1,1);
       
       %============== RANGE 10 A 170 ===============
       if(psi3abs <= sig1psi)
          fpx6(2,1) = 0;
       elseif((psi3abs > sig1psi) && (psi3abs<= sig3psi))
          fpx6(2,1) = (psi3abs - sig1psi)/(sig3psi - sig1psi);
       elseif((psi3abs > sig3psi) && (psi3abs <= sig5psi))
          fpx6(2,1) = (psi3abs - sig5psi)/(sig3psi - sig5psi);
       elseif(psi3abs > sig5psi)
          fpx6(2,1) = 0;
       end 
       fdpertCCpsi(k,1) = fpx6(2,1);
       
       %============== RANGE -80 A 80 ===============

       if(psi3abs <= sig2Npsi)
          fpx6(3,1) = 0;
       elseif((psi3abs > sig2Npsi) && (psi3abs<= sig0psi))
          fpx6(3,1) = (psi3abs-sig2Npsi)/(sig0psi - sig2Npsi);
       elseif((psi3abs > sig0psi) && (psi3abs<= sig2psi))
          fpx6(3,1) = (psi3abs - sig2psi)/(sig0psi - sig2psi);
       elseif(psi3abs > sig2psi)
          fpx6(3,1) = 0;
       end
       fdpertBBpsi(k,1) = fpx6(3,1);
       
       %============== RANGE -170 A 10 ===============
       if(psi3abs <= sig5Npsi)
          fpx6(4,1) = 0;
       elseif((psi3abs > sig5Npsi) && (psi3abs<= sig3Npsi))
          fpx6(4,1) = (psi3abs - sig5Npsi)/(sig3Npsi - sig5Npsi);
       elseif((psi3abs > sig3Npsi) && (psi3abs <= sig1Npsi))
          fpx6(4,1) = (psi3abs - sig1Npsi)/(sig3Npsi - sig1Npsi);
       elseif(psi3abs > sig1Npsi)
          fpx6(4,1) = 0;
       end 
       fdpertAApsi(k,1) = fpx6(4,1);
       
    end
    
    l = 1;
    for h2 = 1:kX6        
        mu_j(l,1) = fpx6(h2,1);
        l = l + 1;
    end
    
    sum_mu_j = sum(mu_j);

    %************************ functions for model **********************    
    if(NewFuncMemb == 1)

       psi3mabs = psi3m;
       %============== RANGE <-180 & >180 ===============
       if(psi3mabs <= sig6Npsi)
          fpx6m(1,1) = 1;
       elseif((psi3mabs > sig6Npsi) && (psi3mabs<= sig4Npsi))
          fpx6m(1,1) = (psi3mabs - sig4Npsi)/(sig6Npsi - sig4Npsi);
       elseif((psi3mabs > sig4psi) && (psi3mabs<= sig6psi))
          fpx6m(1,1) = (psi3mabs - sig4psi)/(sig6psi - sig4psi);
       elseif(psi3mabs > sig6psi)
          fpx6m(1,1) = 1;
       else
          fpx6m(1,1) = 0;
       end
       fdpertDDpsim(k,1) = fpx6m(1,1);
       
       %============== RANGE 10 A 170 ===============
       if(psi3mabs <= sig1psi)
          fpx6m(2,1) = 0;
       elseif((psi3mabs > sig1psi) && (psi3mabs<= sig3psi))
          fpx6m(2,1) = (psi3mabs - sig1psi)/(sig3psi - sig1psi);
       elseif((psi3mabs > sig3psi) && (psi3mabs <= sig5psi))
          fpx6m(2,1) = (psi3mabs - sig5psi)/(sig3psi - sig5psi);
       elseif(psi3mabs > sig5psi)
          fpx6m(2,1) = 0;
       end 
       fdpertCCpsim(k,1) = fpx6m(2,1);
       
       %============== RANGE -80 A 80 ===============
       if(psi3mabs <= sig2Npsi)
          fpx6m(3,1) = 0;
       elseif((psi3mabs > sig2Npsi) && (psi3mabs<= sig0psi))
          fpx6m(3,1) = (psi3mabs-sig2Npsi)/(sig0psi - sig2Npsi);
       elseif((psi3mabs > sig0psi) && (psi3mabs<= sig2psi))
          fpx6m(3,1) = (psi3mabs - sig2psi)/(sig0psi - sig2psi);
       elseif(psi3mabs > sig2psi)
          fpx6m(3,1) = 0;
       end
       fdpertBBpsim(k,1) = fpx6m(3,1);
       
       %============== RANGE -170 A 10 ===============
       if(psi3mabs <= sig5Npsi)
          fpx6m(4,1) = 0;
       elseif((psi3mabs > sig5Npsi) && (psi3mabs<= sig3Npsi))
          fpx6m(4,1) = (psi3mabs - sig5Npsi)/(sig3Npsi - sig5Npsi);
       elseif((psi3mabs > sig3Npsi) && (psi3mabs <= sig1Npsi))
          fpx6m(4,1) = (psi3mabs - sig1Npsi)/(sig3Npsi - sig1Npsi);
       elseif(psi3mabs > sig1Npsi)
          fpx6m(4,1) = 0;
       end 
       fdpertAApsim(k,1) = fpx6m(4,1);

    end
    
    lm = 1;
    for h2 = 1:kX6                                    
        mu_jm(lm,1) = fpx6m(h2,1);
        lm = lm + 1;
    end       
    sum_mu_jm = sum(mu_jm);
    
    % =========================  LEY FUZZY 1 MRAFC =======================
    for j = 1:l-1
        uu(j,:) = ((mu_jm(j,1))*(-Kj(:,:,j)*(xxf-xxfo) + Lj(:,:,j)*rr))/sum_mu_jm;
    end
    uf = sum(uu);
    uf= uf';        
        
    tauMrafc6thr = LT6\uf;
    if (saturation == 1)
        %Limit for tauMracf6
        if(tauMrafc6thr(1,1) > Tforce)
            tauMrafc6thr(1,1) = Tforce;
        elseif(tauMrafc6thr(1,1) < -Tforce)
            tauMrafc6thr(1,1) = -Tforce;
        end

        if(tauMrafc6thr(2,1) > Tforce)
            tauMrafc6thr(2,1) = Tforce;
        elseif(tauMrafc6thr(2,1) < -Tforce)
            tauMrafc6thr(2,1) = -Tforce;
        end

        if(tauMrafc6thr(3,1) > Tforce)
            tauMrafc6thr(3,1) = Tforce;
        elseif(tauMrafc6thr(3,1) < -Tforce)
            tauMrafc6thr(3,1) = -Tforce;
        end

        if(tauMrafc6thr(4,1) > Tforce)
            tauMrafc6thr(4,1) = Tforce;
        elseif(tauMrafc6thr(4,1) < -Tforce)
            tauMrafc6thr(4,1) = -Tforce;
        end
        if(tauMrafc6thr(5,1) > Tforce)
            tauMrafc6thr(5,1) = Tforce;
        elseif(tauMrafc6thr(5,1) < -Tforce)
            tauMrafc6thr(5,1) = -Tforce;
        end
        if(tauMrafc6thr(6,1) > Tforce)
            tauMrafc6thr(6,1) = Tforce;
        elseif(tauMrafc6thr(6,1) < -Tforce)
            tauMrafc6thr(6,1) = -Tforce;
        end
    end
      
    % Calculate transformation matrix T and use it to calculate xdot, ydot, and zdot 
    cphi3 = cos(phi3) ;
    ctheta3 = cos(theta3);
    cpsi3 = cos(psi3) ;
    ttheta3 = tan(theta3);
    sectheta3 = sec(theta3);
    sphi3 = sin(phi3) ;
    stheta3 = sin(theta3);
    spsi3 = sin(psi3);
           
    R_Theta = zeros(3,3);
    R_Theta(1,1) = cpsi3*ctheta3;
    R_Theta(1,2) = cpsi3*stheta3*sphi3 - spsi3*cphi3;
    R_Theta(1,3) = spsi3*sphi3 + cpsi3*stheta3*cphi3;
    R_Theta(2,1) = spsi3*ctheta3;
    R_Theta(2,2) = spsi3*stheta3*sphi3 + cpsi3*cphi3;
    R_Theta(2,3) = -cpsi3*sphi3 + spsi3*stheta3*cphi3;
    R_Theta(3,1) = -stheta3;
    R_Theta(3,2) = ctheta3*sphi3;
    R_Theta(3,3) = ctheta3*cphi3; 

    T_Theta = zeros(3,3);
    T_Theta(1,1) = 1;
    T_Theta(1,2) = sphi3*ttheta3;
    T_Theta(1,3) = cphi3*ttheta3;
    T_Theta(2,1) = 0;
    T_Theta(2,2) = cphi3;
    T_Theta(2,3) = -sphi3;
    T_Theta(3,1) = 0;
    T_Theta(3,2) = sphi3*sectheta3;
    T_Theta(3,3) = cphi3*sectheta3;

    J_Theta = [R_Theta zeros(3)
                zeros(3) T_Theta];

    %========================== MATRIX CALCULATION =======================
    %Coriolis
    C_RB = [ 0          0               0        m*zG*r3        -m*(xG*q3-w3)     -m*(xG*r3+v3)
             0          0               0       -m*w3            m*(zG*r3+xG*p3)   m*u3
             0          0               0       -m*(zG*p3-v3)    -m*(zG*q3+u3)      m*xG*p3
            -m*zG*r3    0               m*zG*p3     0               Izz*r3        -Iyy*q3
             m*xG*q3   -m*(zG*r3+xG*p3) m*zG*q3    -Izz*r3           0             Ixx*p3 
             m*xG*r3    0              -m*xG*p3     Iyy*q3          -Ixx*p3          0];

    %Coriolis Mass Add
    C_A = [  0                  0               0               0             -(Zwd*w3+Zqd*q3)  (Yvd*v3+Yrd*r3)
             0                  0               0              (Zwd*w3+Zqd*q3)       0             -Xud*u3
             0                  0               0             -(Yvd*v3+Yrd*r3)      Xud*u3            0
             0            -(Zwd*w3+Zqd*q3) (Yvd*v3+Yrd*r3)   0                   -(Nvd*v3+Nrd*r3)   (Mwd*w3+Mqd*q3)
            (Zwd*w3+Zqd*q3)       0              -Xud*u3          (Nvd*v3+Nrd*r3)       0              -Kpd*p3 
           -(Yvd*v3+Yrd*r3)     Xud*u3             0             -(Mwd*w3+Mqd*q3)      Kpd*p3            0];   

    %Hydrodynamic Damping Dn(Vr) (Parte no lineal)    
    Dn = [ -Xuu*abs(u3)  0           0           0           0           0
            0           -Yvv*abs(v3) 0           0           0           0
            0           0           -Zww*abs(w3) 0           0           0
            0           0           0           -Kpp*abs(p3) 0           0
            0           0           -Mww*abs(w3) 0           -Mqq*abs(q3) 0
            0           -Nvv*abs(v3) 0           0           0           -Nrr*abs(r3)];     

    %Restoring Forces
    gn = [(W-B)*sin(theta3)
          -(W-B)*cos(theta3)*sin(phi3)
          -(W-B)*cos(theta3)*cos(phi3)
          (zG*W-zb*B)*cos(theta3)*sin(phi3)
          (zG*W-zb*B)*sin(theta3)-xG*(W-B)*cos(theta3)*cos(phi3)
          -xG*(W-B)*cos(theta3)*sin(phi3)];

    %========================= Non-linear equation =======================
    etapMrafc = J_Theta*nuMrafc;
    nupMracf = Mb\(-(C_RB + C_A + Dn)*nuMrafc - gn + LT6*tauMrafc6thr + Tau_ext);

    %=====================================================================
    ss = 1;        
    for ii = 1:l-1
        for jj = 1:l-1      
        xmdf(ss,:) = (mu_jm(ii,1)*mu_jm(jj,1)*(Am2(:,:,ii,jj)*(xmf-xxfo) + Bm2*rr))/(sum_mu_jm*sum_mu_jm);
        ss = ss + 1;
        end
    end

    xmdotf = sum(xmdf);
    xmdotf= xmdotf';
    xmf = xmf + xmdotf*dt;       

    etaMrafc = etaMrafc + etapMrafc*dt;
    nuMrafc = nuMrafc + nupMracf*dt;

    etaMrafc2(1,1) = etaMrafc(1,1);
    etaMrafc2(2,1) = etaMrafc(2,1);
    etaMrafc2(3,1) = etaMrafc(3,1);
    etaMrafc2(4,1) = radtodeg*etaMrafc(4,1);
    etaMrafc2(5,1) = radtodeg*etaMrafc(5,1);
    etaMrafc2(6,1) = radtodeg*etaMrafc(6,1);

    % Save for Plot
    Tau_ext2(k,:)= Tau_ext';
    xmMracf(k,:) = xmf;

    xEuMracfNu(k,:) = nuMrafc';
    xEuMracfEta(k,:) = etaMrafc';
    xEuMracfEta2(k,:) = etaMrafc2';
    uMracf2(k,:) = tauMrafc6thr'; % Control MRACF

    errxk(k,:) = errx;
    errxMracf(k,:) = errd;

    %====== TO CALCULATE AVERAGE SPEED OF NONLINEAR SYSTEM ======
    SpeedX(k,1) = etapMrafc(1,1);
    SpeedX(k,2) = SpeedX(k,1)*2/pi;
    SpeedY(k,1) = etapMrafc(2,1);
    SpeedY(k,2) = SpeedY(k,1)*2/pi;
    
    %-------------------------- update data ------------------------------
    update_data;
    tt(k,1) = t; 
    rzk(k,1) = zr;
    r_step1(k,1) = ref1;
    r_step2(k,1) = ref2;
    r_step3(k,1) = ref3;
    r_step4(k,1) = ref4;
    k = k + 1;

    if(mod(k,1000) == 0)
       kcount=floor(k/1000)
    end
    if (EnbTimeBreak)      
       if(t == timebreak)
           break;
       end
    end
end
%%
%=============================    GRAPHYS   =============================
if(vel >= 0.5)
    t_anim = 20;
elseif((vel >= 0.1) && (vel < 0.5) )
    t_anim = 50;
elseif((vel >= 0.01) && (vel < 0.1) )
    t_anim = 200;
end

if(dist_amp > 0)
    figure(13);
    plot(tt,Tau_ext2(:,1),'m','LineWidth',wid);    grid; title('\tau_{ext}');
    axis([0 104 -0.5 0.5]);
end

%=============    DRAW PATH AND LQR, MRAC AND MRACF CONTROLLER   =========
figure(4);
subplot(2,1,1);   plot(xEuLqrEta(1,1),xEuLqrEta(1,2),'o',xEuMracfEta(1,1),xEuMracfEta(1,2),'o',xdd2(1,1),ydd2(1,1),'o',xmAdp(1,1),xmAdp(1,2),'o',xpAdp(1,1),xpAdp(1,2),'o','LineWidth',wid);
hold on
subplot(2,1,1);   pp211 = plot(xEuLqrEta(:,1),xEuLqrEta(:,2),xpAdp(:,1),xpAdp(:,2),xEuMracfEta(:,1),xEuMracfEta(:,2),xdd2,ydd2,'LineWidth',wid);
pp211(1).Color = clqr; 
pp211(1).DisplayName = 'lqr';
pp211(2).Color = cmrac; 
pp211(2).DisplayName = 'mrac';
pp211(3).Color = cmrafc; 
pp211(3).DisplayName = 'mrafc';
pp211(4).Color = cpath; 
pp211(4).DisplayName = 'trajectory';
legend([pp211(1) pp211(2) pp211(3) pp211(4)],'location','northeast');
title('Top view (inertial-y vs. inertial-x)');
axis('equal');grid on; zoom on;
ylabel('y[m]','FontWeight','bold');
xlabel('x[m]','FontWeight','bold');

subplot(2,1,2);   plot(xEuLqrEta(1,1),xEuLqrEta(1,3),'o',xEuMracfEta(1,1),xEuMracfEta(1,3),'o',xdd2(1,1),zdd2(1,1),'o',xmAdp(1,1),xmAdp(1,3),'o',xpAdp(1,1),xpAdp(1,3),'o','LineWidth',wid);
hold on
subplot(2,1,2);   pp212 = plot(xEuLqrEta(:,1),xEuLqrEta(:,3),xpAdp(:,1),xpAdp(:,3),xEuMracfEta(:,1),xEuMracfEta(:,3),xdd2,zdd2,'LineWidth',wid);
pp212(1).Color = clqr; 
pp212(1).DisplayName = 'lqr';
pp212(2).Color = cmrac; 
pp212(2).DisplayName = 'mrac';
pp212(3).Color = cmrafc; 
pp212(3).DisplayName = 'mrafc';
pp212(4).Color = cpath; 
pp212(4).DisplayName = 'trajectory';
legend([pp212(1) pp212(2) pp212(3) pp212(4)],'location','northeast');
title('Side view (inertial-z vs. inertial-x)');set(gca, 'Ydir', 'reverse');
axis('equal'); grid on; zoom on;
ylabel('y[m]','FontWeight','bold');
xlabel('z[m]','FontWeight','bold');

%============================ FIGURE 3D 2 ============================    
figure(66); 
p1 = plot3(yr2,xr2,rzk,'Color',cpath,'DisplayName','trajectory','LineWidth',wid); % Grafico Senal Deseada
hold on
p2 = plot3(xEuLqrEta(:,2),xEuLqrEta(:,1),xEuLqrEta(:,3),'Color',clqr,'DisplayName','lqr','LineWidth',wid); % Grafico Control LQR
p3 = plot3(zpAdp(:,2),zpAdp(:,1),zpAdp(:,3),'Color',cmrac,'DisplayName','mrac','LineWidth',wid); % Grafico Control Adaptativo
p4 = plot3(xEuMracfEta(:,2),xEuMracfEta(:,1),xEuMracfEta(:,3),'Color',cmrafc,'DisplayName','mracf','LineWidth',wid); % Grafico Control MRACF

hAxis = gca;
hAxis.XRuler.FirstCrossoverValue  = 0; % X crossover with Y axis
hAxis.YRuler.FirstCrossoverValue  = 0; % Y crossover with X axis
hAxis.ZRuler.FirstCrossoverValue  = 0; % Z crossover with X axis
hAxis.XRuler.SecondCrossoverValue = 0; % X crossover with Z axis
hAxis.YRuler.SecondCrossoverValue = 0; % Y crossover with Z axis
hAxis.ZRuler.SecondCrossoverValue = 0; % Z crossover with Y axis
legend([p1 p2 p3 p4],'location','northeast');
hold off
title(['3D trajectory']);
xlabel('y(m)'); ylabel('x(m)'); zlabel('z(m)'); % etiqueta de ejes
zlim([-2 10]);
grid; % Activa Regilla
set(gca, 'Xdir', 'reverse'); 
set(gca, 'Ydir', 'reverse'); 
set(gca, 'Zdir', 'reverse');   

%============================== GRAPHYS ERROR ============================
figure(70);
subplot(3,2,1);   pp321 = plot(tt,errxEuLqr(:,1),tt,errAdp(:,1),tt,errxMracf(:,1),'LineWidth',wid);  grid; title('Error x(m)');
pp321(1).Color = clqr; 
pp321(1).DisplayName = 'lqr';
pp321(2).Color = cmrac; 
pp321(2).DisplayName = 'mrac';
pp321(3).Color = cmrafc; 
pp321(3).DisplayName = 'mrafc';
legend([pp321(1) pp321(2) pp321(3)],'location','northeast','Orientation','horizontal');
ylabel('x[m]','FontWeight','bold');

subplot(3,2,3);   pp323 = plot(tt,errxEuLqr(:,2),tt,errAdp(:,2),tt,errxMracf(:,2),'LineWidth',wid);  grid; title('Error y(m)');
pp323(1).Color = clqr; 
pp323(1).DisplayName = 'lqr';
pp323(2).Color = cmrac; 
pp323(2).DisplayName = 'mrac';
pp323(3).Color = cmrafc; 
pp323(3).DisplayName = 'mrafc';
ylabel('y[m]','FontWeight','bold');

subplot(3,2,5);   pp325 = plot(tt,errxEuLqr(:,3),tt,errAdp(:,3),tt,errxMracf(:,3),'LineWidth',wid);  grid; title('Error z(m)'); set(gca, 'Ydir', 'reverse');
pp325(1).Color = clqr; 
pp325(1).DisplayName = 'lqr';
pp325(2).Color = cmrac; 
pp325(2).DisplayName = 'mrac';
pp325(3).Color = cmrafc; 
pp325(3).DisplayName = 'mrafc';
ylabel('z[m]','FontWeight','bold');
xlabel('time[s]','FontWeight','bold');

subplot(3,2,2);   pp322 = plot(tt,errxEuLqr(:,4),tt,errAdp(:,4),tt,errxMracf(:,4),'LineWidth',wid);    grid; title('Error Roll (degrees)');  
pp322(1).Color = clqr; 
pp322(1).DisplayName = 'mlqr';
pp322(2).Color = cmrac; 
pp322(2).DisplayName = 'mrac';
pp322(3).Color = cmrafc; 
pp322(3).DisplayName = 'mrafc';
ylabel('\Phi[deg]','FontWeight','bold');

subplot(3,2,4);   pp324 = plot(tt,errxEuLqr(:,5),tt,errAdp(:,5),tt,errxMracf(:,5),'LineWidth',wid);    grid; title('Error Pitch (degrees)');   
pp324(1).Color = clqr; 
pp324(1).DisplayName = 'lqr';
pp324(2).Color = cmrac; 
pp324(2).DisplayName = 'mrac';
pp324(3).Color = cmrafc; 
pp324(3).DisplayName = 'mrafc';
ylabel('\Theta[deg]','FontWeight','bold');

subplot(3,2,6);   pp326 = plot(tt,errxEuLqr(:,6),tt,errAdp(:,6),tt,errxMracf(:,6),'LineWidth',wid);    grid; title('Error Yaw (degrees)');
pp326(1).Color = clqr; 
pp326(1).DisplayName = 'lqr';
pp326(2).Color = cmrac; 
pp326(2).DisplayName = 'mrac';
pp326(3).Color = cmrafc; 
pp326(3).DisplayName = 'mrafc';
ylabel('\Psi[deg]','FontWeight','bold');
xlabel('time[s]','FontWeight','bold');

%============================= GRAPHYS REGULATION ======================== 
figure(95);
subplot(3,2,1);   pp321 = plot(tt,uLqr2(:,1),tt,upAdp2(:,1),tt,uMracf2(:,1),'LineWidth',wid);    grid; title('Tport - Left');
pp321(1).Color = clqr; 
pp321(1).DisplayName = 'lqr';
pp321(2).Color = cmrac; 
pp321(2).DisplayName = 'mrac';
pp321(3).Color = cmrafc; 
pp321(3).DisplayName = 'mrafc';
legend([pp321(1) pp321(2) pp321(3)],'location','northeast','Orientation','horizontal');
ylabel('u_1[N]','FontWeight','bold');

subplot(3,2,3);   pp323 = plot(tt,uLqr2(:,3),tt,upAdp2(:,3),tt,uMracf2(:,3),'LineWidth',wid);    grid; title('Taft_h - Back');
pp323(1).Color = clqr; 
pp323(1).DisplayName = 'lqr';
pp323(2).Color = cmrac; 
pp323(2).DisplayName = 'mrac';
pp323(3).Color = cmrafc; 
pp323(3).DisplayName = 'mrafc';
ylabel('u_3[N]','FontWeight','bold');

subplot(3,2,5);   pp325 = plot(tt,uLqr2(:,5),tt,upAdp2(:,5),tt,uMracf2(:,5),'LineWidth',wid);    grid; title('Taft_v - Back');
pp325(1).Color = clqr; 
pp325(1).DisplayName = 'lqr';
pp325(2).Color = cmrac; 
pp325(2).DisplayName = 'mrac';
pp325(3).Color = cmrafc; 
pp325(3).DisplayName = 'mrafc';
ylabel('u_5[N]','FontWeight','bold');
xlabel('time[s]','FontWeight','bold');

subplot(3,2,2);   pp322 = plot(tt,uLqr2(:,2),tt,upAdp2(:,2),tt,uMracf2(:,2),'LineWidth',wid);    grid; title('Tstbd - Right');
pp322(1).Color = clqr; 
pp322(1).DisplayName = 'lqr';
pp322(2).Color = cmrac; 
pp322(2).DisplayName = 'mrac';
pp322(3).Color = cmrafc; 
pp322(3).DisplayName = 'mrafc';
ylabel('u_2[N]','FontWeight','bold');

subplot(3,2,4);   pp324 = plot(tt,uLqr2(:,4),tt,upAdp2(:,4),tt,uMracf2(:,4),'LineWidth',wid);    grid; title('Tfore_h - Forward');
pp324(1).Color = clqr; 
pp324(1).DisplayName = 'lqr';
pp324(2).Color = cmrac; 
pp324(2).DisplayName = 'mrac';
pp324(3).Color = cmrafc; 
pp324(3).DisplayName = 'mrafc';
ylabel('u_4[N]','FontWeight','bold');

subplot(3,2,6);   pp326 = plot(tt,uLqr2(:,6),tt,upAdp2(:,6),tt,uMracf2(:,6),'LineWidth',wid);    grid; title('Tfore_v - Forward');
pp326(1).Color = clqr; 
pp326(1).DisplayName = 'lqr';
pp326(2).Color = cmrac; 
pp326(2).DisplayName = 'mrac';
pp326(3).Color = cmrafc; 
pp326(3).DisplayName = 'mrafc';
ylabel('u_6[N]','FontWeight','bold');
xlabel('time[s]','FontWeight','bold');

%======================== GRAPHYS REGULATION COMPACT =====================
figure(96);
if (regulation == 0)
    subplot(3,2,1);   pp321 = plot(tt,xEuLqrEta(:,1),tt,zpAdp(:,1),tt,xEuMracfEta(:,1), tt,r_step1,'LineWidth',wid);  grid; title('inertial frame x(m)');
else
    subplot(3,2,1);   pp321 = plot(tt,xEuLqrEta(:,1),tt,zpAdp(:,1),tt,xEuMracfEta(:,1), tt,xdd2,'LineWidth',wid);  grid; title('inertial frame x(m)');
end
pp321(1).Color = clqr; 
pp321(1).DisplayName = 'lqr';
pp321(2).Color = cmrac; 
pp321(2).DisplayName = 'marc';
pp321(3).Color = cmrafc; 
pp321(3).DisplayName = 'mrafc';
pp321(4).Color = cpath; 
pp321(4).DisplayName = 'reference';
legend([pp321(1) pp321(2) pp321(3) pp321(4)],'location','southwest','Orientation','horizontal');
ylabel('x[m]','FontWeight','bold');

if (regulation == 0)
    subplot(3,2,3);   pp323 = plot(tt,xEuLqrEta(:,2),tt,zpAdp(:,2),tt,xEuMracfEta(:,2), tt,r_step2,'LineWidth',wid);  grid; title('inertial frame y(m)');
else
    subplot(3,2,3);   pp323 = plot(tt,xEuLqrEta(:,2),tt,zpAdp(:,2),tt,xEuMracfEta(:,2), tt,ydd2,'LineWidth',wid);  grid; title('inertial frame y(m)');
end
pp323(1).Color = clqr; 
pp323(1).DisplayName = 'lrq';
pp323(2).Color = cmrac; 
pp323(2).DisplayName = 'marc';
pp323(3).Color = cmrafc; 
pp323(3).DisplayName = 'mrafc';
pp323(4).Color = cpath; 
pp323(4).DisplayName = 'reference';
ylabel('y[m]','FontWeight','bold');

if (regulation == 0)
    subplot(3,2,5);   pp325 = plot(tt,xEuLqrEta(:,3),tt,zpAdp(:,3),tt,xEuMracfEta(:,3), tt,r_step3,'LineWidth',wid);  grid; title('inertial frame z(m)'); set(gca, 'Ydir', 'reverse');
else
    subplot(3,2,5);   pp325 = plot(tt,xEuLqrEta(:,3),tt,zpAdp(:,3),tt,xEuMracfEta(:,3), tt,zdd2,'LineWidth',wid);  grid; title('inertial frame z(m)'); set(gca, 'Ydir', 'reverse');
end
pp325(1).Color = clqr; 
pp325(1).DisplayName = 'lqr';
pp325(2).Color = cmrac; 
pp325(2).DisplayName = 'mrac';
pp325(3).Color = cmrafc; 
pp325(3).DisplayName = 'mrafc';    
pp325(4).Color = cpath; 
pp325(4).DisplayName = 'reference';
xlabel('time[s]','FontWeight','bold');
ylabel('z[m]','FontWeight','bold');

subplot(3,2,2);   pp322 = plot(tt,radtodeg*xEuLqrEta(:,4),tt,radtodeg*zpAdp(:,4),tt,radtodeg*xEuMracfEta(:,4),'LineWidth',wid);    grid; title('Roll phi (degrees)');
pp322(1).Color = clqr; 
pp322(1).DisplayName = 'lqr';
pp322(2).Color = cmrac; 
pp322(2).DisplayName = 'marc';
pp322(3).Color = cmrafc; 
pp322(3).DisplayName = 'mrafc';
ylabel('\Phi[deg]','FontWeight','bold');

subplot(3,2,4);   pp324 = plot(tt,radtodeg*xEuLqrEta(:,5),tt,radtodeg*zpAdp(:,5),tt,radtodeg*xEuMracfEta(:,5),'LineWidth',wid);    grid; title('Pitch theta (degrees)');
pp324(1).Color = clqr; 
pp324(1).DisplayName = 'lqr';
pp324(2).Color = cmrac; 
pp324(2).DisplayName = 'marc';
pp324(3).Color = cmrafc; 
pp324(3).DisplayName = 'mrafc';
ylabel('\Theta[deg]','FontWeight','bold');

subplot(3,2,6);   pp326 = plot(tt,radtodeg*xEuLqrEta(:,6),tt,radtodeg*zpAdp(:,6),tt,radtodeg*xEuMracfEta(:,6),'LineWidth',wid);    grid; title('Yaw psi (degrees)');
pp326(1).Color = clqr; 
pp326(1).DisplayName = 'lqr';
pp326(2).Color = cmrac; 
pp326(2).DisplayName = 'marc';
pp326(3).Color = cmrafc; 
pp326(3).DisplayName = 'mrafc';
ylabel('\Psi[deg]','FontWeight','bold');
xlabel('time[s]','FontWeight','bold');

%performance index   
rmse_lqr = rms(error_lqr);
rmse_mrac = rms(error_mrac);
rmse_mrafc = rms(error_mrafc);