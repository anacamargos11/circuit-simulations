
%  1D-FDTD code for the solution of  a matched transmission line  
%  with a resistor added to the load end.
%
%  This code simulates the incident wave for the 'oneresistor' code.

%  Ana C. C. Couto, CNPq Science Without Borders,
%  supervised by Costas D. Sarris,   
%  Department of Electrical and Computer Engineering,
%  University of Toronto 
%
%  ECE 1252H, Summer 2016

close all;
clear ; % initialization  
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


time_steps = 1798;                  % time steps to be run
Ncells     = 800;                   % number of cells 
dx         = c0/fmax/20. * meters;  % cell size is chosen lambda/20
s          = 0.9;                   % CFL number
dt         = s*dx/c0;               % time step duration
tm         = dt*(1:time_steps) ;    
N_FFT      = time_steps;
FREQ       = linspace(0,fmax,N_FFT);
K          = exp(-1i*2*pi*dt*FREQ);
RECM       = zeros(1,N_FFT);        % initializing fourier transforms
SRC        = zeros(1,N_FFT);
nsamp      = Ncells-10 ;

L_dis = (pi*1e-7)/2 ;              % distributed inductance
C_dis = (1.0/(L_dis*(c0*c0)));     % distributed capacitance
R_dism = zeros(1,Ncells) ;          % distributed resistance
G_dism = zeros(1,Ncells) ;          % distributed conductance
Zc    = sqrt(L_dis/C_dis)  ;       % impedance of the line
Z_lim = Zc/dx;                     % limit for comparison of impedances
RLm   = Zc;                        % resistance of the load
Rs    = 0.1;                       % resistance of the source


% RESISTOR LOCATION AND CONDUCTANCE UPDATE
nres        = Ncells;                 % cell of load resistor
G_dism(nres) = 1/RLm/dx ;            % updated conductance at the load
R_dism(nres) = RLm*dx;               % updated resistance at the load

% ARRAYS OF VOLTAGE AND CURRENT
vnm(1:Ncells+1) = 0 ;  
inm(1:Ncells)   = 0 ; 

v1_prev  = 0 ;
vN_prev  = 0 ;                     % auxiliary variables for Boundary C.
vN_prev1 = 0 ; 
v1_prev1 = 0 ;

% TRANSMISSION LINE EQS' AND MUR'S CONSTANTS
vrelm = (G_dism*dt)./(2.*C_dis) ;
irelm = (R_dism*dt)./(2.*L_dis) ;

cinm  = (1 - irelm)./(1 + irelm)       ;
ainm  = dt./( L_dis.*dx.*(1.+irelm) ) ; 

cvnm  = (1 - vrelm)./(1 + vrelm)       ;
avnm  = dt./( C_dis.*dx.*(1.+vrelm) ) ; 

mur_coef = (s - 1)/(s + 1) ; 


% == MOVIE INITIALIZATION ==

x=linspace(dx, Ncells*dx, Ncells) ; 
 
%subplot(2,1,1),plot(x,vnm(1:Ncells),'r'),axis([0 .5 -1 1]);
%ylabel('Voltage');

%subplot(2,1,2),plot(x,inm,'b'),axis([0 .5 -3.0e-3 3.0e-3]);
%xlabel('x (meters)');ylabel('Current');

%rect=get(gcf,'Position');
%rect(1:2)=[0 0];

%M=moviein(time_steps/20,gcf,rect);

%vrefm = zeros(1, time_steps) ;
vincm = zeros(1, time_steps) ;


%
% Marching-in-time loop
%
 
tn = 0;
for n = 1:time_steps 

  % ================
  %  Voltage update
  % ================ 

   vnm(2:Ncells) = cvnm(2:Ncells).*vnm(2:Ncells) -...
                  avnm(2:Ncells).*( inm(2:Ncells) - inm(1:Ncells-1) )  ;  
  

   % -- boundary conditions
   us = exp(-((tn-t0)./Ts).*((tn-t0)./Ts));  % source term
   us1 = exp(-((tn+dt-t0)./Ts).*((tn+dt-t0)./Ts)); % source term on time n+1
   cus = 0.5*(us+us1); % avarage term that gives source on time n+1/2
   
   vnm(1) = v1_prev + 2*(cus-inm(1)*Rs)  ;     % boundary equation at 1st cell
   
   vnm(Ncells) = RLm*inm(Ncells-1) ;         % corrected boundary equation
   
   v1_prev  = -vnm(1); 
     
   % -- source (transparent)

   tn = tn + dt;
   vnm(1) = vnm(1) + us;
   
   % recorded field

   vincm(n) = vnm(nsamp) ;

   % current update

   inm(1:Ncells) = cinm(1:Ncells).*inm(1:Ncells) -... 
                  ainm(1:Ncells).*( vnm(2:Ncells+1) - vnm(1:Ncells) ) ; 
              
   for nf = 1:N_FFT       
       RECM(nf) = RECM(nf) + (K(nf)^time_steps)*vnm(nsamp);     
   end
  
   
%
% VISUALISATION
%

 %if mod(n,20)==0;

   %rtime=num2str(round(n*dt/1.0e-12));
   %cstep=num2str(n) ;

  % subplot(2,1,1),plot(x,vnm(1:Ncells),'r'),axis([0 Ncells*dx -2 2]);
  % title(['time = ',rtime,' ps',' step=',cstep]);
  % ylabel('Voltage');%

  % subplot(2,1,2),plot(x,inm,'b'),axis([0 Ncells*dx -.1 .1]);
  % xlabel('x (meters)');ylabel('Current');

 % M(:,n/20)=getframe(gcf,rect);
 
 %end

end

%movie(gcf,M,0,20,rect)

Fvincm = fft(vincm , N_FFT) ;  % purely incident source transform

