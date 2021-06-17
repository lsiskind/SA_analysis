function[abkoz2, abkoz4, tabkoz2, tabkoz4] = alphabetaIntrusions(PSARJ_flt, P42A_BGA_flt, P44A_BGA_flt, timevec_flt);

% alpha angles at which bga intrusions occur


tabkoz2 = timevec_flt(((PSARJ_flt>=255 & PSARJ_flt <= 289) & (P42A_BGA_flt >=24 & P42A_BGA_flt <= 150)) | ((PSARJ_flt>=255 & PSARJ_flt <= 289) & (P42A_BGA_flt >=229 & P42A_BGA_flt <= 337)));
abkoz2 = PSARJ_flt(((PSARJ_flt>=255 & PSARJ_flt <= 289) & (P42A_BGA_flt >=24 & P42A_BGA_flt <= 150)) | ((PSARJ_flt>=255 & PSARJ_flt <= 289) & (P42A_BGA_flt >=229 & P42A_BGA_flt <= 337)));

tabkoz4 = timevec_flt(((PSARJ_flt>=75 & PSARJ_flt <= 110) & (P44A_BGA_flt >=41 & P44A_BGA_flt<=149)) | ((PSARJ_flt>=75 & PSARJ_flt <= 110) & (P44A_BGA_flt>=230 & P44A_BGA_flt<=338)));
abkoz4 = PSARJ_flt(((PSARJ_flt>=75 & PSARJ_flt <= 110) & (P44A_BGA_flt >=41 & P44A_BGA_flt<=149)) | ((PSARJ_flt>=75 & PSARJ_flt <= 110) & (P44A_BGA_flt>=230 & P44A_BGA_flt<=338)));

end