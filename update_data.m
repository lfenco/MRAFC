function [] = update_data()
global l EGV sngl1 mu_j Bm2 Plmi2 errx xxf xxfo sum_mu_j rr Kjdot Ljdot Kj Lj dt

for j = 1:l-1
    EGV = eig(Lj(:,:,j));
    if(all(EGV > 0))
        sngl1 = 1;
    elseif(all(EGV < 0))
        sngl1 = -1;
    end    
    Kjdot(:,:,j) = (mu_j(j,1)*sngl1*Bm2'*Plmi2*errx*(xxf-xxfo)')/sum_mu_j;
    Ljdot(:,:,j) = -(mu_j(j,1)*sngl1*Bm2'*Plmi2*errx*rr')/sum_mu_j;
    Kj(:,:,j) = Kj(:,:,j) + Kjdot(:,:,j)*dt;
    Lj(:,:,j) = Lj(:,:,j) + Ljdot(:,:,j)*dt;            
end
end