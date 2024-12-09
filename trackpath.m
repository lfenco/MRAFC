function [] = trackpath()
global ti dt tf ntdd ruido EnbNoise EnbTimeBreak ltdd tt3 timebreak
global CntPoints ti2 dt2 tf2 t2 nn wxdd wydd wzdd wv waypoints depth velz
global tdd waypoints2 waypoints3 xdd ydd zdd xdd3 ydd3 zdd3

tdd = ti:dt:tf;
ntdd = length(tdd);
ruido = EnbNoise*0.15*randn(ntdd,1);
depth=velz.*tdd;

if (EnbTimeBreak == 0)
    ltdd=length(tdd);
else
    tt3 = ti:dt:timebreak; tt3 = tt3';
    ltdd=length(tt3);
end

if(CntPoints == 0)
    %trayectoria circular 8 puntos
    ti2 = 0;   dt2 = 15.625;     tf2 = tf;
    t2 = ti2:dt2:tf2;
    nn=1;
    wxdd = 2.5*cos(2*nn*pi*t2/tf2 + pi/2)+2.5;
    wydd = 2.5*sin(2*nn*pi*t2/tf2 - pi/2);
    wzdd = velz.*t2;
    wv = 0.122;
    waypoints = [ wxdd(1)  wydd(1)   wzdd(1)    0   wv
                  wxdd(2)  wydd(2)   wzdd(2)    0   wv
                  wxdd(3)  wydd(3)   wzdd(3)    0   wv
                  wxdd(4)  wydd(4)   wzdd(4)    0   wv
                  wxdd(5)  wydd(5)   wzdd(5)    0   wv
                  wxdd(6)  wydd(6)   wzdd(6)    0   wv
                  wxdd(7)  wydd(7)   wzdd(7)    0   wv
                  wxdd(8)  wydd(8)   wzdd(8)    0   wv
                  wxdd(9)  wydd(9)   wzdd(9)    0   wv];    
elseif(CntPoints == 1)   
    %trayectoria circular 16 puntos
    ti2 = 0;   dt2 = 7.8125;     tf2 = tf;
    t2 = ti2:dt2:tf2;
    nn=1;
    depth=velz.*t2;
    wxdd = 2.5*cos(2*nn*pi*t2/tf2 + pi/2)+2.5;
    wydd = 2.5*sin(2*nn*pi*t2/tf2 - pi/2);
    wzdd = velz.*t2;
    wv = 0.1248;

    waypoints = [ wxdd(1)  wydd(1)   wzdd(1)    0   wv
                  wxdd(2)  wydd(2)   wzdd(2)    0   wv
                  wxdd(3)  wydd(3)   wzdd(3)    0   wv
                  wxdd(4)  wydd(4)   wzdd(4)    0   wv
                  wxdd(5)  wydd(5)   wzdd(5)    0   wv
                  wxdd(6)  wydd(6)   wzdd(6)    0   wv
                  wxdd(7)  wydd(7)   wzdd(7)    0   wv
                  wxdd(8)  wydd(8)   wzdd(8)    0   wv
                  wxdd(9)  wydd(9)   wzdd(9)    0   wv
                  wxdd(10)  wydd(10)   wzdd(10)    0   wv
                  wxdd(11)  wydd(11)   wzdd(11)    0   wv
                  wxdd(12)  wydd(12)   wzdd(12)    0   wv
                  wxdd(13)  wydd(13)   wzdd(13)    0   wv
                  wxdd(14)  wydd(14)   wzdd(14)    0   wv
                  wxdd(15)  wydd(15)   wzdd(15)    0   wv
                  wxdd(16)  wydd(16)   wzdd(16)    0   wv
                  wxdd(17)  wydd(17)   wzdd(17)    0   wv];
end

%trayectory_generation(waypoints,type_int,type_yaw,r,t)
waypoints2 = trayectory_generation(waypoints,0,1,0.2,tdd);
xdd = waypoints2(:,1);
ydd = waypoints2(:,2);
zdd = waypoints2(:,3);

%sirve para solo para grafico
waypoints3 = trayectory_generation(waypoints,0,1,0.2,tdd);
xdd3 = waypoints3(:,1);
ydd3 = waypoints3(:,2);
zdd3 = waypoints3(:,3);
end