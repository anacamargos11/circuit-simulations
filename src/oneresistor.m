
%  1D-FDTD code for the solution of the transmission line equations 
%  with a resistor added to the load end.
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
run match1res;

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
REF        = zeros(1,N_FFT);        % initializing fourier transforms
REC        = zeros(1,N_FFT);
SRC        = zeros(1,N_FFT);
nsamp      = Ncells-10 ;

L_dis = (pi*1e-7)/2 ;              % distributed inductance
C_dis = (1.0/(L_dis*(c0*c0)));     % distributed capacitance
R_dis = zeros(1,Ncells) ;          % distributed resistance
G_dis = zeros(1,Ncells) ;          % distributed conductance
Zc    = sqrt(L_dis/C_dis)  ;       % impedance of the line
Z_lim = Zc/dx;                     % limit for comparison of impedances
RL    = 80;                        % resistance of the load
Rs    = 0.1;                       % resistance of the source


% RESISTOR LOCATION AND CONDUCTANCE UPDATE
nres        = Ncells;                 % cell of load resistor
G_dis(nres) = 1/RL/dx ;            % updated conductance at the load
R_dis(nres) = RL*dx;               % updated resistance at the load

% ARRAYS OF VOLTAGE AND CURRENT
vn(1:Ncells+1) = 0 ;  
in(1:Ncells)   = 0 ; 

v1_prev  = 0 ; % auxiliary variable for Boundary C.

% TRANSMISSION LINE EQS' AND MUR'S CONSTANTS
vrel = (G_dis*dt)./(2.*C_dis) ;
irel = (R_dis*dt)./(2.*L_dis) ;

cin  = (1 - irel)./(1 + irel)       ;
ain  = dt./( L_dis.*dx.*(1.+irel) ) ; 

cvn  = (1 - vrel)./(1 + vrel)       ;
avn  = dt./( C_dis.*dx.*(1.+vrel) ) ; 

mur_coef = (s - 1)/(s + 1) ; 


% == MOVIE INITIALIZATION ==

x=linspace(dx, Ncells*dx, Ncells) ; 
 
subplot(2,1,1),plot(x,vn(1:Ncells),'r'),axis([0 .5 -1 1]);
ylabel('Voltage');

subplot(2,1,2),plot(x,in,'b'),axis([0 .5 -3.0e-3 3.0e-3]);
xlabel('x (meters)');ylabel('Current');

rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(time_steps/20,gcf,rect);

%vref = zeros(1, time_steps) ;

%
% Marching-in-time loop
%
 
tn = 0;
for n = 1:time_steps 

  % ================
  %  Voltage update
  % ================ 

   vn(2:Ncells) = cvn(2:Ncells).*vn(2:Ncells) -...
                  avn(2:Ncells).*( in(2:Ncells) - in(1:Ncells-1) )  ;    

   % -- boundary conditions
   us = exp(-((tn-t0)./Ts).*((tn-t0)./Ts));  % source term
   us1 = exp(-((tn+dt-t0)./Ts).*((tn+dt-t0)./Ts)); % source term on time n+1
   cus = 0.5*(us+us1); % avarage term that gives source on time n+1/2
   
   vn(1) = v1_prev + 2*(cus-in(1)*Rs)  ;     % boundary equation at 1st cell

   vn(Ncells) = RL*in(Ncells-1) ;           % boundary equation at last cell
   
   v1_prev  = -vn(1); 
   iN_prev  = -in(Ncells-1);
   
   % -- source (transparent)

   tn = tn + dt;
   vn(1) = vn(1) + us;
   
   % recorded field

   vrec(n) = vn(nsamp);

   % current update

   in(1:Ncells) = cin(1:Ncells).*in(1:Ncells) -... 
                  ain(1:Ncells).*( vn(2:Ncells+1) - vn(1:Ncells) ) ; 
              
   % CURRENT boundary condition 
   if (RL > 75)
        in(Ncells-1) = iN_prev + (vn(Ncells)+vn(Ncells-1))/RL;
   end
   
   % Update Fourier Transforms  
   for nf = 1:N_FFT       
       REC(nf) = REC(nf) + (K(nf)^time_steps)*vn(Ncells-10);           
   end
    
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

freq = (1/dt)*(0:N_FFT/2)/N_FFT ;

REF = REC - RECM;
S11 = REF./RECM;

figure(2)

plot(1.e-9*freq,  20.*log10(abs(S11(1:N_FFT/2+1)))) ;
axis([0 250 -21 1])
xlabel('Frequency [GHz]', 'fontsize', 16) ;
ylabel('S_{11} [dBs]', 'fontsize', 16)
grid on
hold on

% THEORY

RL_array = RL*ones(1,N_FFT);
Zc_array = Zc*ones(1,N_FFT);

S11_th = (RL_array - Zc_array)./(RL_array + Zc_array) ;

plot(1.e-9.*freq, 20.*log10(abs(S11_th(1:N_FFT/2+1))), '--', 'Linewidth', 2)

legend('FDTD', 'TL-theory') ;

axis([0 250 -21 1])
xlabel('Frequency [GHz]', 'fontsize', 16) ;
ylabel('S_{11} [dBs]', 'fontsize', 16)
grid on

