
%  1D-FDTD code for the solution of the transmission line equations. 
%  parallel RC elements are added to the load end.
%  a matching network (at f = 100MHz) is inserted, with a capacitor placed
%  at a distance d = 0.207*lambda from the load, in parallel with the source.
%
%  Calculation of S11 term is performed.
%
%                            <---------------- d ---------------->
%       ----------------------------------------------------------
%                            |                                    |
%                            |                                   ---  
%                            |                                  |   |
%                            -  C_mid                    C_load -   > R_load
%                            -                                  -   <
%                            |                                  |___|
%                            |                                    |
%       ---------------------------------------------------------- 
%
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

fmax  = 1e-1*gigahertz      ;        % gaussian excitation parameters
Ts    = 1./(2.*fmax) ;
t0    = 3.*Ts        ; 

fsin  = fmax;                       % fixed frequency for sin wave
w     = 2*pi*fsin;                  % omega frequency
k     = w./c0;                      % wave number

time_steps = 1998;                  % time steps to be run
Ncells     = 800;                   % number of cells 
dx         = c0/fmax/40. * meters;  % cell size is chosen lambda/40
s          = 0.8;                   % CFL number
dt         = s*dx/c0;               % time step duration
tm         = dt*(1:time_steps) ;    
N_FFT      = time_steps;
FREQ       = linspace(0,fmax,N_FFT);
freq       = (1/dt)*(0:N_FFT/2)/N_FFT ;
w_array    = 2*pi*FREQ;             % angular frequency array
K          = exp(-1i*2*pi*dt*FREQ);
REF        = zeros(1,N_FFT);        % initializing fourier transforms
REC        = zeros(1,N_FFT);
SRC        = zeros(1,N_FFT);
nsamp      = Ncells-50 ;

L_dis   = ((pi*1e-7)/1.88365)*ones(1,Ncells) ;  % distributed inductance
C_dis   = (1.0./(L_dis*(c0*c0)));    % distributed capacitance
R_dis   = zeros(1,Ncells) ;          % distributed resistance
G_dis   = zeros(1,Ncells) ;          % distributed conductance
Zc      = sqrt(L_dis/C_dis)  ;       % impedance of the line
Rs      = 0.01;                       % resistance of the source

R_load  = 125;                        % resistance of the load
C_load  = 2.548e-11;                 % capacitance of the load
L_load  = 2e119;                     % inductance of the load
C_mid   = 50e-12;                    % capacitor inserted in the middle

YL = (1./R_load) + 1i.*w_array.*C_load; % admittance of the load
ZL = 1./YL;                             % impedance of the load

% LOAD ELEMENTS LOCATION AND UPDATE
  ncap         = Ncells-8;         % cell of load capacitor
  C_dis(ncap)  = C_mid./dx;        % updated capacitance at the middle


% ARRAYS OF VOLTAGE AND CURRENT
vn(1:Ncells+1) = 0 ;  
in(1:Ncells)   = 0 ; 
inind          = 0 ;
sn             = 0 ;

v1m_prev  = 0 ;  % auxiliary variable for Boundary C.


% TRANSMISSION LINE EQS' AND BOUNDARY CONDITION CONSTANTS
vrel = (G_dis*dt)./(2.*C_dis) ;
irel = (R_dis*dt)./(2.*L_dis) ;

cin  = (1 - irel)./(1 + irel)       ;
ain  = dt./( L_dis.*dx.*(1.+irel) ) ; 

cvn  = (1 - vrel)./(1 + vrel)       ;
avn  = dt./( C_dis.*dx.*(1.+vrel) ) ; 

mur_coef = (s - 1)/(s + 1) ; 

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

M=moviein(time_steps/80,gcf,rect);

vrec = zeros(1, time_steps) ;

%
% Marching-in-time loop
%
 
tn = 0;
xn = 0;
for n = 1:time_steps 

  % ================
  %  Voltage update
  % ================ 

   vn(2:Ncells-1) = cvn(2:Ncells-1).*vn(2:Ncells-1) -...
                  avn(2:Ncells-1).*( in(2:Ncells-1) - in(1:Ncells-2) )  ;  

   % -- boundary conditions
   us = exp(-((tn-t0)/Ts)^2);  % source term
   us_array = exp(-((tm-t0)./Ts).*((tm-t0)./Ts));
   us1 = exp(-((tn+dt-t0)./Ts).*((tn+dt-t0)./Ts)); % source term on time n+1
   cus = 0.5*(us+us1); % avarage term that gives source on time n+1/2
   
   vn(1) = v1m_prev + mur_coef*(vn(2) - vn(1))               ; 
   
   v1m_prev = vn(2);
   
   % boundary equations at last cell
   
   % calculation of inductor current via voltage integral
   sn = sn + vn(Ncells); 
   inind = (dt./L_load).*sn;
            
   % using phase function and accounting for time difference for i and v             
    vn(Ncells) = ((cbt(Ncells)-cet(Ncells))./(cat(Ncells)+cdt(Ncells))).*vn(Ncells)-...
                 inind.*cct(Ncells)./(cat(Ncells)+cdt(Ncells)) + ...
                 in(Ncells-1).*cft(Ncells)./(cat(Ncells)+cdt(Ncells));
               
   % -- source (transparent)
   tn = tn + dt;
   vn(1) = us;                         % gaussian source
   
   % recorded field

   vrec(n) = vn(nsamp);

   % current update

   in(1:Ncells) = cin(1:Ncells).*in(1:Ncells) -... 
                  ain(1:Ncells).*( vn(2:Ncells+1) - vn(1:Ncells) ) ; 
                
%
% VISUALISATION
%

 if mod(n,80)==0;

   rtime=num2str(round(n*dt/1.0e-12));
   cstep=num2str(n) ;

   subplot(2,1,1),plot(x,vn(1:Ncells),'r'),axis([0 Ncells*dx -2 2]);
   title(['time = ',rtime,' ps',' step=',cstep]);
   ylabel('Voltage');

   subplot(2,1,2),plot(x,in,'b'),axis([0 Ncells*dx -.08 .08]);
   xlabel('x (meters)');ylabel('Current');

  M(:,n/80)=getframe(gcf,rect);
 
 end

end

movie(gcf,M,0,80,rect) 


%% Conversion to frequency domain

% fft using source
tinc = tm - (nsamp-1)*dx/c0 ;
vinc = exp(-((tinc-t0)./Ts).*((tinc-t0)./Ts))  ;
Fvref = fft(vrec-vinc, N_FFT) ;
Fvinc = fft(vinc , N_FFT) ;
S11_fft1  = Fvref./Fvinc ;

figure(3)
plot(1.e-7*freq,  20.*log10(abs(S11_fft1(1:N_FFT/2+1))),...
    '--', 'Linewidth', 2) ;
axis([0 11 -40 5])
xlabel('Frequency [10 MHz]', 'fontsize', 16) ;
ylabel('fft source S_{11} [dBs]', 'fontsize', 16)
grid on
