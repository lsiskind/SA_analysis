function [afov2, tafov2, akoz2, takoz2, afov4, tafov4, akoz4, takoz4] = alphaIntrusions(PSARJ_flt, timevec_flt)

% time instances and alpha angles of dynamic intrusions into stray light FOV
tafov2 = timevec_flt(PSARJ_flt >= 268 & PSARJ_flt <= 276);
afov2 = PSARJ_flt(PSARJ_flt >= 268 & PSARJ_flt <= 276);

tafov4 = timevec_flt(PSARJ_flt >= 89 & PSARJ_flt <= 96);
afov4 = PSARJ_flt(PSARJ_flt >= 89 & PSARJ_flt <= 96);

% time instances and alpha angles of dyncamic intrusions into stray light KOZ envelope
takoz2 = timevec_flt(PSARJ_flt>=255 & PSARJ_flt<=289);
akoz2 = PSARJ_flt(PSARJ_flt>=255 & PSARJ_flt<=289);

takoz4 = timevec_flt(PSARJ_flt>=75 & PSARJ_flt<=110); 
akoz4 = PSARJ_flt(PSARJ_flt>=75 & PSARJ_flt<=110);


end