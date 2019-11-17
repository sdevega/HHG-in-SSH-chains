function drhodt = drho(t,rhot)
% EOM: drho = -i [H,rho] - gamma (rho - rho0) 

global N gamma X H0 rho0 vC  E0

global Delta omega t0

rhot = reshape(rhot,[N,N]);         % 1) reshape matrix from 1xN*N to NxN
%rhot = sparse(rhot);               % 2) sparse matrix to save memory and speed up -- no good

Eext    = 2.0 * E0 * sin(omega*(t-t0)) * exp(-(t-t0)^2/Delta^2);

phi_ext = -X*Eext;                % diagonal matrix
phi_ind = -2*vC*diag(rhot-rho0);  % vector 
phi = phi_ext + diag(phi_ind);    % diagonal matrix

H0rhot  = H0*rhot;
phirhot = phi*rhot;
drhodt0 = -1j*(H0rhot-H0rhot')+1j*(phirhot-phirhot')-gamma*(rhot-rho0);

drhodt  = drhodt0(:);    % flatten the matrix again
end



