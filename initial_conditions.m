function [] = initial_conditions()
global qeta qnu QQ RR PP KK xp up xm sp sm err Dh tau Tee Txx Tuu K_adap
global angl AA BB CC EE Ap Bp Cp Ep D Af Bf Cf Df Ke Kx Ku Ta Tb Ki
global flag xmf Kd_mrafc Kp_mrafc KpKd_mracf Ljastt EGV sngl1
%======================= LQR initial conditions  =========================  
    % pesos
    qeta = 1e8;
    qnu = 1e4;

    QQ = [qeta*eye(6)  zeros(6)
        zeros(6)  qnu*eye(6)];
    RR = [ 1 ];
    PP = are(AA(:,:,1),BB*inv(RR)*BB',QQ);% original
    KK = inv(RR)*BB'*PP;    
 

%======================== MARC initial conditions  =======================
    %Condicion inicial planta
    xp = [2.5 -2.5 0 0 0 angl 0 0 0 0 0 0]';
    up = [0 0 0 0 0 0]';
   
    xm = [2.5 -2.5 0 0 0 angl 0 0 0 0 0 0]'; %Condicion inicial modelo   
    sp = [0 0 0 0 0 0]'; %Condicion inicial feedforward

    %Condicion inicial feedforward model
    sm = [0 0 0 0 0 0]';
    err = [0 0 0 0 0 0]';
    
    Dh= 0.27;    
    tau =2.7;
    Tee = 400;
    Txx =0.75;
    Tuu = 17;
    K_adap = 50*eye(6);
    
    % State space form of Plant
    % ****************************
    Ap = AA(:,:,3); Bp = BB; Cp = CC; Ep = EE;
       
    % State space form of a FF Compensator
    % *****************************
    D = Dh*eye(6);
    
    Af = -1/tau;
    Bf = 1/tau;
    Cf = D;
    Df = 0;
    
    Ke = Tee*eye(6);    
    Kx = Txx*eye(12);
    Ku = Tuu*eye(6);    

    Ta = [Ke         zeros(6,12)  zeros(6,6)
        zeros(12,6)     Kx       zeros(12,6)
        zeros(6,6)  zeros(6,12)     Ku    ];

    Tb = Ta;
    Ki = zeros(6,24);
    
%========================= MRAFC initial conditions  =====================
    %Condicion inicial modelo                                            
    if (flag==0)                                        
        xmf = [0 0 0 0 0 angl 0 0 0 0 0 0]'; % ESTADO PARA INICAR PARA                                             %TRIANGULO
    else
        xmf = [2.5 -2.5 0 0 0 angl 0 0 0 0 0 0]'; %ESTADO PARA INICIAR DESDE
                                                %UN PUNTO TANGENTE                                    
    end
    
    Kd_mrafc = 0.05*eye(6);
    Kp_mrafc = 1*Kd_mrafc;
    KpKd_mracf = [Kp_mrafc Kd_mrafc];
    
    Ljastt = eye(6);
    EGV = eig(Ljastt);
    if(all(EGV > 0))
        sngl1 = 1;
    elseif(all(EGV < 0))   
        sngl1 = -1;
    end
end