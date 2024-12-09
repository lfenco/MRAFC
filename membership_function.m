function [] = membership_function()
global degtorad xx6c fdpertBNpsi fdpertApsi fdpertBpsi fdpertCpsi
global kX1 kX2 kX3 kX4 kX5 kX6 kX1p kX2p kX3p kX4p kX5p kX6p pX pXp pXt
global sig1Npsi sig2Npsi sig3Npsi sig4Npsi sig5Npsi sig6Npsi sig0psi
global sig1psi sig2psi sig3psi sig4psi sig5psi sig6psi wid

disp('-> [----- Grafica nuevas funciones de membresia -----]')
sig1Npsi = -1*degtorad; sig2Npsi = -89*degtorad; 
sig3Npsi = -90*degtorad; sig4Npsi = -91*degtorad;
sig5Npsi = -179*degtorad; sig6Npsi = -180*degtorad;
sig0psi = 0*degtorad; 
sig1psi = 1*degtorad; sig2psi = 89*degtorad;   
sig3psi = 90*degtorad; sig4psi = 91*degtorad; 
sig5psi = 179*degtorad; sig6psi = 180*degtorad;  

cc = 1;
for xx6 = -pi:0.001:pi
   xx6c(cc,1) = xx6;
   xx6abs = xx6;
   %============== RANGO -170 A 10 ===============
   if(xx6abs <= sig5Npsi)
      fdpBNpsi = 0;
   elseif((xx6abs > sig5Npsi) & (xx6abs<= sig3Npsi))
      fdpBNpsi = (xx6abs - sig5Npsi)/(sig3Npsi - sig5Npsi);
   elseif((xx6abs > sig3Npsi) & (xx6abs <= sig1Npsi))
      fdpBNpsi = (xx6abs - sig1Npsi)/(sig3Npsi - sig1Npsi);
   elseif(xx6abs > sig1Npsi)
      fdpBNpsi = 0;
   end 
   fdpertBNpsi(cc,1) = fdpBNpsi;

   %============== RANGO -80 A 80 ===============

   if(xx6abs <= -sig2psi)
      fdpApsi = 0;
   elseif((xx6abs > sig2Npsi) & (xx6abs<= sig0psi))
      fdpApsi = (xx6abs-sig2Npsi)/(sig0psi - sig2Npsi);
   elseif((xx6abs > sig0psi) & (xx6abs<= sig2psi))
      fdpApsi = (xx6abs - sig2psi)/(sig0psi - sig2psi);
   elseif(xx6abs > sig2psi)
      fdpApsi = 0;
   end
   fdpertApsi(cc,1) = fdpApsi;

   %============== RANGO 10 A 170 ===============
   if(xx6abs <= sig1psi)
      fdpBpsi = 0;
   elseif((xx6abs > sig1psi) & (xx6abs<= sig3psi))
      fdpBpsi = (xx6abs - sig1psi)/(sig3psi - sig1psi);
   elseif((xx6abs > sig3psi) & (xx6abs <= sig5psi))
      fdpBpsi = (xx6abs - sig5psi)/(sig3psi - sig5psi);
   elseif(xx6abs > sig5psi)
      fdpBpsi = 0;
   end 
   fdpertBpsi(cc,1) = fdpBpsi;

   %============== RANGO <-180 & >180 ===============
   if(xx6abs <= sig6Npsi)
      fdpCpsi = 1;
   elseif((xx6abs > sig6Npsi) & (xx6abs<= sig4Npsi))
      fdpCpsi = (xx6abs - sig4Npsi)/(sig6Npsi - sig4Npsi);
   elseif((xx6abs > sig4psi) & (xx6abs<= sig6psi))
      fdpCpsi = (xx6abs - sig4psi)/(sig6psi - sig4psi);
   elseif(xx6abs > sig6psi)
      fdpCpsi = 1;
   end
   fdpertCpsi(cc,1) = fdpCpsi;   
   cc = cc + 1;
end

figure(9);
plot(xx6c,fdpertBNpsi,xx6c,fdpertApsi,xx6c,fdpertBpsi,xx6c,fdpertCpsi,'LineWidth',2);
axis([-3.14 3.14 0 1.0]);
grid; % Activa Regilla
yticks([0 0.2 0.4 0.6 0.8 1.0])
%legend('M1','M2','M3','M4','M5','Orientation','horizontal')

kX1 = 1;kX2 = 1;kX3 = 1;kX4 = 1;kX5 = 1;kX6 = 4;
kX1p = 1;kX2p = 1;kX3p = 1;kX4p = 1;kX5p = 1;kX6p = 1;
pX = kX1*kX2*kX3*kX4*kX5*kX6;
pXp = kX1p*kX2p*kX3p*kX4p*kX5p*kX6p;
pXt = pX*pXp;

end