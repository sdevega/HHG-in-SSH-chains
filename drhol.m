function drhodt = drhol(t,rhot)
% EOM: drho = -i [H,rho] - gamma (rho - rho0) 

global N gamma X H0 rho0 vC  E0

global Delta omega

rhot = reshape(rhot,[N,N]);         % 1) reshape matrix from 1xN*N to NxN
%rhot = sparse(rhot);               % 2) sparse matrix to save memory and speed up -- no good

if t<Delta
   Eext    = 2.0*E0*cos(omega*t); % delta pulse
else
   Eext    = 0;
end

phi_ext = -X*Eext;                % diagonal matrix
phi_ind = -2*vC*diag(rhot-rho0);  % vector 
phi = phi_ext + diag(phi_ind);    % diagonal matrix

drhodt0 = -1j*(H0*rhot - rhot*H0)+1j*(phi*rhot-rhot*phi)-gamma*(rhot-rho0);

drhodt  = drhodt0(:);    % flatten the matrix again
end

