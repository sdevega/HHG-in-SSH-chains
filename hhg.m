%---------------------------------------------------------------%
%                                                               %  
% HHG in SSH chains  (TB with 2 hoppings)                       %
%                                                               %
%---------------------------------------------------------------%
%                                                               %
% Numerical solution of the EOM:                                %
%                                                               %
% drho = -i [H,rho] - gamma (rho - rho0)                        %
%                                                               %
%---------------------------------------------------------------%
%                                                               %
% H = H0 - e phi = H0 - e (phi_ind + phi_ext)                   %
%                                                               %
% phi_ind = v (rho - rho0)                                      %
% phi_ext = x*E                                                 %
% E = 2 E0 sin(omega*(t-t0)) * exp(-(t-t0)^2/Delta^2)           %
%    E0 = sqrt(2 pi I0/c)                                       %
% v_{ll'} = 1/a/|l-l'| for l neq l' and v0=0.577 otherwise      %
%                                                               %
%                                                               %
% H0      TB Hamiltonian                                        %
% I0      pulse intensity                                       %
% omega   input frequency                                       %
% Delta   pulse width                                           %
% t0      pulse center                                          %
% v0      onsite Coulomb interaction                            %
% a       atom spacing                                          %
% l       atom site                                             %
%                                                               %
%---------------------------------------------------------------%

%clear all
close all

%% --- General constants

global c 

eV = 27.2113834;         % 1au=27eV
a0 = 0.529177208;        % a0 in nm
c  = 137.03599971;       % c in au
nm = 10/a0;              % 1nm=10/a0 au
vf = c/300;              % Fermi velocity graphene
hbar  = 6.582119514e-16;       % eV*s
fs    = 0.024188843265;        % 1 au of time in fs   
au_s  = 2.4188843e-17;         % 1 au of time in sec   
kT    = 0.0001*8.61733e-5/eV;  % Boltzman constant times temperature (T). Here: T=0
au_kg = 9.1093819e-31;         % 1 au of mass in kg 


%% --- Geometry and input constants

global N gamma

N = 25;  % number of atoms in the SSH chain

t1    = -2.3/eV;   % hopping 1
t2    = -2.8/eV;   % hopping 2
a     = 0.1421*nm; % atomic spacing
gamma = 0.5*50e-3/eV;  % damping

global I0 Delta t0 omega E0

I0    = 1.0e14*au_s^3/au_kg;   % delta pulse "intensity"
Delta = 100/fs/1.6651092;      % delta pulse width
t0    = 300/fs;                % pulse center in time
v0    = 0.577761;              % onsite Coulomb interaction
E0    = sqrt(2*pi*I0/c);

global H0 X vC rho0

omega = 1.22652/eV;            % input freq

close all
%% --- Unpertutbed Hamiltonian: H0

H0 = zeros(N,N);      % Hamiltonian
for ii=1:N
   if rem(ii,2)~=0    % odd
      if ii<N
         H0(ii,ii+1)=t1;
      end
      if ii>1
         H0(ii,ii-1)=t2;
      end
   else
      H0(ii,ii-1)=t1;
      if ii<N
         H0(ii,ii+1)=t2;
      end
   end
end

%% --- Position matrix

X = zeros(N,N);
for ii=1:N
   X(ii,ii)=a*(ii-1);
end

%X = sparse(X);

%% --- Coulomb interaction

vC = zeros(N,N);
for ii=1:N
   for jj=1:N
      if(jj==ii)
         vC(ii,jj)=v0;
      else
         vC(ii,jj)=1/a/abs(ii-jj);
      end
   end
end

%% --- Diagonalization of the unperturbed Hamiltonian

% [V,D] = eig(A) returns diagonal matrix D of eigenvalues and matrix V whose columns are
% the corresponding right eigenvectors, so that A*V = V*D

[V,D] = eig(H0); 


%% --- Equilibrium density matrix (zero temperature)

% in the state representation: 
%   rho0_{jj'} = Kronecker_delta_{jj'} * Fermi-Dirac_{jj'}  for  j neq j'
%              = 0                                          for  j = j'
% change of basis to go to site representation: rho0_site = V * rho_state * V'

rho0j = zeros(N,N); 
for ii=1:N
   rho0j(ii,ii)=FD(D(ii,ii));  % in the state representation
end

rho0 = V*(rho0j*V');           % in the site representation

%rho0 = sparse(rho0);

%% --- ODE solver

% general info:       https://es.mathworks.com/help/matlab/ordinary-differential-equations.html
% specific for ode45: https://es.mathworks.com/help/matlab/ref/ode45.html

% [t,y] = ode45(odefun,tspan,y0), where tspan = [t0 tf], integrates the system of
% differential equations y'=f(t,y) from t0 to tf with initial conditions y0.
% Nonstiff ODE
% Each row in the solution array y corresponds to a value returned in column vector t.
% [t,y] = ode15s(odefun,tspan,y0);  % stiff ODEs

tspan   = [0 : 0.01/fs : 800/fs];
opts = odeset('RelTol',1e-15,'AbsTol',1e-15);
tic
[t,rho] = ode45(@drho,tspan,rho0(:),opts);
toc


fprintf('ODE solver time steps: %e \n',length(t));
fprintf('\n\n');

%% --- Induced dipole in time:  pind = -2e tr[ x (rho - rho0) ]

pind = zeros(length(t),1);
for ii=1:length(t)
   pind(ii,1) = -2*trace(X*(reshape(rho(ii,:),[N,N]) - rho0));
end

%close all
%% test
figure
plot(t.*fs,real(pind),'-k')
hold on
plot(t.*fs,imag(pind),'-m')
xlabel('$t$ (fs)', 'Interpreter', 'latex');
ylabel('$p_{ind}$ (a.u.)', 'Interpreter', 'latex');
title('matlab')
set(gca,'FontSize',15);
set(gca,'TickDir','out');


% model: from charpa 
%nin3   = sprintf('cchain_N%i_Q0_gam50_I14_gauss_100fs_a0.1421_t1%g_t2%g_HF1_w%g.dat.td0',N,t1*eV,t2*eV,omega*eV);
%dataw3 = dlmread(nin3,'',0,0);
%tm     = dataw3(:,1);  % fs
%pretm  = dataw3(:,2);  % a.u.
%pimtm  = dataw3(:,3);  % a.u.

%figure
%plot(tm,pretm,'-b')
%hold on
%plot(tm,pimtm,'-r')
%xlabel('$t$ (fs)', 'Interpreter', 'latex');
%ylabel('$p_{ind}$ (a.u.)', 'Interpreter', 'latex');
%title('model');
%set(gca,'FontSize',15);
%set(gca,'TickDir','out');


%% --- Induced dipole in frequency. Fourier transform of pind

% Fourier transform definition:
%   p_ind(w) = int_tmin^tmax dt p_ind(t) exp(i * w * t) 


nw = 1000; 
w  = linspace(0.001/eV,7/eV,nw);

nout2 = sprintf('chain_ode45_hhg_N%i_t1%g_t2%g_I14_w%g.dat',N,t1*eV,t2*eV,omega*eV);
fout2 = fopen(nout2,'w');

pindw = zeros(nw,1);

fprintf('FT time for %d freqs.: \n',nw);


% -- direct integral of pind(t)*exp(i w t)
% reference for trapz (integral): https://es.mathworks.com/help/matlab/ref/trapz.html
% trapz(X,Y) integrates Y with respect to the coordinates or scalar spacing specified by X.
tic
for ii=1:nw
   pinter = pind.*exp(1j*w(ii).*t);
   pindw(ii,1) = trapz(t,pinter);
   fprintf(fout2,'%g %g %g \n',w(ii)*eV,real(pindw(ii,1)),imag(pindw(ii,1)));
end
toc


% approximating to a line: f(t) = aa*t + bb
% analytic integral: int dt (aa*t + bb)*exp(cc*t) = (aa*t + bb - aa/cc)*exp(cc*t)/cc
%{
% -- matrix treatment
pind1 = pind(2:length(pind),1);
pind0 = pind(1:(length(pind)-1),1);
time1 = t(2:length(t),1);
time0 = t(1:(length(t)-1),1);

tic
for jj=1:nw
   cc  = 1j*w(jj);
   aa  = (pind1-pind0)./(time1-time0);
   bb  = pind0 - aa.*time0;
   ff0 = (aa.*time0 + bb - aa./cc).*exp(cc.*time0)./cc;
   ff1 = (aa.*time1 + bb - aa./cc).*exp(cc.*time1)./cc;
   
   pindw(jj) = sum(ff1-ff0);
   
   fprintf(fout2,'%g %g %g \n',w(jj)*eV,real(pindw(jj)),imag(pindw(jj)));
end
toc
fprintf('\n\n');
%}
fclose(fout2);



% --- test
figure
%plot(w.*eV,real(pindw),'-k','Linewidth',2)
hold on
semilogy(w.*eV,imag(pindw),'-m','Linewidth',2)
xlabel('$\hbar\omega_{out}$ (eV)', 'Interpreter', 'latex');
ylabel('Im \{$p_{ind}$\} (a.u.)', 'Interpreter', 'latex');
title('matlab');
legend('real','imag')
set(gca,'FontSize',15);
set(gca,'TickDir','out');




%% --- Linear plots 

Nc  =  6;      % number of atoms in the SSH chain
t1c = -2.3;    % hopping 1
t2c = -2.8;    % hopping 2
wc  =  1.22652;

% https://es.mathworks.com/help/matlab/math/choose-an-ode-solver.html


% retrieved in Matlab with RK45 - ode45
nin1    = sprintf('chain_hhg_N%i_t1%g_t2%g_I14_w%g.dat',Nc,t1c,t2c,wc);
dataw1  = dlmread(nin1,'',0,0);
wf1     = dataw1(:,1);  % eV
pre1    = dataw1(:,2);  % a.u.
pim1    = dataw1(:,3);  % a.u.


p01 = sqrt(pre1.^2 + pim1.^2);


figure
semilogy(wf1,p01,'xk','Linewidth',2)

xlabel('$\hbar\omega_{out}$ (eV)', 'Interpreter', 'latex');
ylabel('$|p_{ind}|$ (a.u.)', 'Interpreter', 'latex');
xlim([0 7])
legend('ode45')
set(gca,'FontSize',15);
set(gca,'TickDir','out');



%% --- Other FT

%{
% -- brute force FT 
% C-like Fourier transform approximating to a line -- same as matrix
for jj=1:nw
   ff = 0; 
   for ii=1:length(t)-1
      aa  = (pind(ii+1)-pind(ii)) / (t(ii+1)-t(ii));
      bb  = pind(ii) - aa*t(ii);
      cc  = 1j*w(jj);
      ff1 = (aa*t(ii)   + bb - aa/cc)*exp(cc*t(ii))  /cc;
      ff2 = (aa*t(ii+1) + bb - aa/cc)*exp(cc*t(ii+1))/cc;
      ff  = ff + ff2 - ff1;
   end
   pindw(jj) = ff;%*0.5 / Delta / E0;
   fprintf(fout2,'%g %g %g \n',w(jj)*eV,real(ff),imag(ff));
end
%}
