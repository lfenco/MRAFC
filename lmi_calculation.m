function [] = lmi_calculation()
global loaddataM AA2 BB2 Plmi2 Kjastt
    
    %===================== FOR REFERENCE MODEL ========================
    if (loaddataM == 1)
        load lmiModel09102023
    else  
        setlmis([]);
        Qlmi2=lmivar(1,[size(AA2(:,:,1),1) 1]);
        M11=lmivar(2,[6 12]);
        M21=lmivar(2,[6 12]);
        M31=lmivar(2,[6 12]);
        M41=lmivar(2,[6 12]);

        xcond=newlmi;
        lmiterm([-xcond 1 1 Qlmi2],1,1);
        stcond1=newlmi;
        lmiterm([stcond1 1 1 Qlmi2],AA2(:,:,1) ,1,'s');
        lmiterm([stcond1 1 1 M11],BB2,-1,'s');
        stcond2=newlmi;
        lmiterm([stcond2 1 1 Qlmi2],AA2(:,:,2) ,1,'s');
        lmiterm([stcond2 1 1 M21],BB2,-1,'s');%B2
        stcond3=newlmi;
        lmiterm([stcond3 1 1 Qlmi2],AA2(:,:,3) ,1,'s');
        lmiterm([stcond3 1 1 M31],BB2,-1,'s');
        stcond4=newlmi;
        lmiterm([stcond4 1 1 Qlmi2],AA2(:,:,4) ,1,'s');
        lmiterm([stcond4 1 1 M41],BB2,-1,'s');%B2
        stcond5=newlmi;
        lmiterm([stcond5 1 1 Qlmi2],AA2(:,:,1) ,1,'s');
        lmiterm([stcond5 1 1 Qlmi2],AA2(:,:,2) ,1,'s');
        lmiterm([stcond5 1 1 M21],BB2,-1,'s');
        lmiterm([stcond5 1 1 M11],BB2,-1,'s');%B2

        stcond6=newlmi;    
        lmiterm([stcond6 1 1 Qlmi2],AA2(:,:,3) ,1,'s');
        lmiterm([stcond6 1 1 Qlmi2],AA2(:,:,4) ,1,'s');
        lmiterm([stcond6 1 1 M41],BB2,-1,'s');
        lmiterm([stcond6 1 1 M31],BB2,-1,'s');%B2
        %------------------------------------
        lmisys=getlmis;
        [tmin,xfeas]=feasp(lmisys);
        Xopt2=dec2mat(lmisys,xfeas,Qlmi2);
        M11opt=dec2mat(lmisys,xfeas,M11);
        M21opt=dec2mat(lmisys,xfeas,M21);
        M31opt=dec2mat(lmisys,xfeas,M31);
        M41opt=dec2mat(lmisys,xfeas,M41);
        Plmi2=inv(Xopt2);
        Kjastt(:,:,1) = M11opt*Plmi2;
        Kjastt(:,:,2) = M21opt*Plmi2;
        Kjastt(:,:,3) = M31opt*Plmi2;
        Kjastt(:,:,4) = M41opt*Plmi2;

        %*************************************
        % Saving data model
        %*************************************
        save lmiModel09102023 Plmi2 Kjastt;
    end
end