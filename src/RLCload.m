
%  1D-FDTD code for the solution of the transmission line equations 
%  with parallel RCL elements added to the load.
%  Reflection patterns of voltage and current waves at the load end 
%  are displayed. 
%  Calculation of S11 term is performed.

%  Ana C. C. Couto, CNPq Science Without Borders,
%  supervised by Costas D. Sarris,   
%  Department of Electrical and Computer Engineering,
%  University of Toronto 
%
%  ECE 1252H, Summer 2016

close all;
clear ; 
clf   ;  % initialization
clc   ;

% UNITS
meters = 1;
centimeters = 1e-2 * meters;
millimeters = 1e-3 * meters;
inches      = 2.54 * centimeters;
feet        = 12 * inches;
seconds     = 1;
hertz       = 1/seconds;
kilohertz   = 1e3 * hertz;
megahertz   = 1e6 * hertz;
gigahertz   = 1e9 * hertz;

%%%% FDTD PARAMETERS AND PROGRAM CONSTANTS
c0 = 299792458 * meters/seconds;   % speed of light in free space
e0 = 8.8541878176e-12 * 1/meters;  % permittivity of free space 
u0 = 1.2566370614e-6 * 1/meters;   % permeability of free space

fmax  = 20*gigahertz      ;        % gaussian excitation parameters
Ts    = 1./(2.*fmax) ;
t0    = 3.*Ts        ; 

time_steps = 1658;                  % time steps to be run
Ncells     = 800;                   % number of cells 
dx         = c0/fmax/40. * meters;  % cell size is chosen lambda/40
s          = 0.99;                   % CFL number
dt         = s*dx/c0;               % time step duration
tm         = dt*(1:time_steps) ;    
N_FFT      = time_steps;
FREQ       = linspace(0,fmax,N_FFT);
freq       = (1/dt)*(0:N_FFT/2)/N_FFT ;
w_array    = 2*pi*FREQ;             % angular frequency array
nsamp      = Ncells-10 ;

L_dis   = ((pi*1e-7)/2)*ones(1,Ncells) ;  % distributed inductance
C_dis   = (1.0./(L_dis*(c0*c0)));    % distributed capacitance
R_dis   = zeros(1,Ncells) ;          % distributed resistance
G_dis   = zeros(1,Ncells) ;          % distributed conductance
Zc      = sqrt(L_dis/C_dis)  ;       % impedance of the line

R_load  = 70e100;                    % resistance of the load
C_load  = 1e-12;                     % capacitance of the load
L_load  = 2e119;                     % inductance of the load

YL = (1./R_load) + (1./(1i.*w_array.*L_load)) + 1i.*w_array.*C_load;
ZL = 1./YL;

% ARRAYS OF VOLTAGE AND CURRENT
vn(1:Ncells+1) = 0 ;  
in(1:Ncells)   = 0 ; 
inind          = 0 ;
sn             = 0 ;

% TRANSMISSION LINE EQS' AND BOUNDARY CONDITION CONSTANTS
vrel = (G_dis*dt)./(2.*C_dis) ;
irel = (R_dis*dt)./(2.*L_dis) ;

cin  = (1 - irel)./(1 + irel)       ;
ain  = dt./( L_dis.*dx.*(1.+irel) ) ; 

cvn  = (1 - vrel)./(1 + vrel)       ;
avn  = dt./( C_dis.*dx.*(1.+vrel) ) ; 

mur_coef = (s - 1)/(s + 1) ; 

cint = dt./L_load;

% constants for phase equation of 3 elements together
cat = 1 + dt.*(G_dis.*dx + 1./R_load)./(2.*C_load + 2.*C_dis.*dx);
cbt = 1 - dt.*(G_dis.*dx + 1./R_load)./(2.*C_load + 2.*C_dis.*dx);
cct = dt./(C_load + C_dis.*dx);
cdt = (dt.*dt)./(4.*L_load.*(C_load + C_dis.*dx));
cet = cdt;
cft = dt./(C_load + C_dis.*dx);
cgt = dt./(2.*L_load);

% == MOVIE INITIALIZATION ==

x=linspace(dx, Ncells*dx, Ncells) ; 
 
subplot(2,1,1),plot(x,vn(1:Ncells),'r'),axis([0 .5 -1 1]);
ylabel('Voltage');

subplot(2,1,2),plot(x,in,'b'),axis([0 .5 -3.0e-3 3.0e-3]);
xlabel('x (meters)');ylabel('Current');

rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(time_steps/20,gcf,rect);

vrec = zeros(1, time_steps) ;

%
% Marching-in-time loop
%
 
tn = 0;
for n = 1:time_steps 

  % ================
  %  Voltage update
  % ================ 

   vn(2:Ncells-1) = cvn(2:Ncells-1).*vn(2:Ncells-1) -...
                  avn(2:Ncells-1).*( in(2:Ncells-1) - in(1:Ncells-2) )  ;  

   % -- boundary conditions
   tn=n*dt ; 
   us = exp(-((tn-t0)/Ts)^2);  % source term
   us_array = exp(-((tm-t0)./Ts).*((tm-t0)./Ts));
   us1 = exp(-((tn+dt-t0)./Ts).*((tn+dt-t0)./Ts)); % source term on time n+1
   cus = 0.5*(us+us1); % avarage term that gives source on time n+1/2
   
   % boundary equations at last cell
   
   % calculation of inductor current via voltage integral
   sn = sn + vn(Ncells); 
   inind = (dt./L_load).*sn;
               
   % using phase function and accounting for time difference for i and v             
    vn(Ncells) = ((cbt(Ncells)-cet(Ncells))./(cat(Ncells)+cdt(Ncells))).*vn(Ncells)-...
                 inind.*cct(Ncells)./(cat(Ncells)+cdt(Ncells)) + ...
                 in(Ncells-1).*cft(Ncells)./(cat(Ncells)+cdt(Ncells));
   
   % -- source (transparent)

   %tn = tn + dt;
   vn(1) = us;
   
   % recorded field

   vrec(n) = vn(nsamp);

   % current update

   in(1:Ncells) = cin(1:Ncells).*in(1:Ncells) -... 
                  ain(1:Ncells).*( vn(2:Ncells+1) - vn(1:Ncells) ) ;                     
%
% VISUALISATION
%

 if mod(n,20)==0;

   rtime=num2str(round(n*dt/1.0e-12));
   cstep=num2str(n) ;

   subplot(2,1,1),plot(x,vn(1:Ncells),'r'),axis([0 Ncells*dx -2 2]);
   title(['time = ',rtime,' ps',' step=',cstep]);
   ylabel('Voltage');

   subplot(2,1,2),plot(x,in,'b'),axis([0 Ncells*dx -.08 .08]);
   xlabel('x (meters)');ylabel('Current');

  M(:,n/20)=getframe(gcf,rect);
 
 end

end

movie(gcf,M,0,20,rect) 


%% Conversion to frequency domain

% fft using source
tinc = tm - (nsamp-1)*dx/c0 ;
vinc = exp(-((tinc-t0)./Ts).*((tinc-t0)./Ts))  ;
Fvref = fft(vrec-vinc, N_FFT) ;
Fvinc = fft(vinc , N_FFT) ;
S11_fft1  = Fvref./Fvinc ;

% THEORY

S11_th = (ZL - Zc)./(ZL + Zc) ;

figure('units','centimeters','position',[.1 .1 30 15])
plot(1.e-9.*FREQ, 20.*log10(abs(S11_th(1:N_FFT))), ...
    '--', 'Linewidth', 2)
axis([0 20 -30 10])
xlabel('Frequency [GHz]', 'fontsize', 16) ;
ylabel('S_{11} [dBs]', 'fontsize', 16)
grid on
hold on
plot(1.e-9*freq,  20.*log10(abs(S11_fft1(1:N_FFT/2+1)))) ;
axis([0 20 -30 20])
xlabel('Frequency [GHz]', 'fontsize', 16) ;
ylabel('fft source S_{11} [dBs]', 'fontsize', 16)
legend('TL-theory', 'FDTD') ;
