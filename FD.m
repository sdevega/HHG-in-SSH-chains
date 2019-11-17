function FermiDirac = FD(EE)
% Fermi-Dirac distribution at zero temperature

if(EE < 0) 
   FermiDirac = 1.0;
elseif(EE == 0)
   FermiDirac = 0.5;
else
   FermiDirac = 0.0;
end

end