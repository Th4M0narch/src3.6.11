function T = metrics_to_table(label, M)
T = table({label},M.N,M.FAC2,M.FAC5,M.FAC10,M.MG,M.VG,M.FB,M.NMSE,M.NAD,M.CC,...
    'VariableNames',{'Case','N','FAC2','FAC5','FAC10','MG','VG','FB','NMSE','NAD','CC'});
end
