
%  1D-FDTD code for the solution of the transmission line equations. 
%  Sine-shaped soft sorce is introduced. 
%  Lumped resistors are added at the load end.
%  Reflected wave interferes with incident wave, resulting in
%  a standing wave pattern.
%  V and Zc*I waves are displayed.

%  Ana C. C. Couto, CNPq Science Without Borders,
%  supervised by Costas D. Sarris,   
%  Department of Electrical and Computer Engineering,
%  University of Toronto 
%
%  ECE 1252H, Summer 2016


close all;
clear ; 
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

fmax = 20.e9      ;                % gaussian excitation parameters
w    = 2*pi*fmax;                  % omega frequency
k    = w./c0;                      % wave number

time_steps = 3000;                  % time steps to be run
Ncells     = 800;                   % number of cells 
dx         = c0/fmax/20. * meters;  % cell size is chosen lambda/20
s          = 0.915;                   % CFL number
dt         = s*dx/c0;               % time step duration

L_dis = 4*pi*1e-7 ;                % distributed inductance
C_dis = 1.0/( L_dis*(c0*c0) ) ;    % distributed capacitance
R_dis = zeros(1,Ncells) ;          % distributed resistance
G_dis = zeros(1,Ncells) ;          % distributed conductance
Zc = sqrt(L_dis/C_dis)  ;          % impedance of the line
RL = 20;                           % resistance of the load

% RESISTOR LOCATION AND CONDUCTANCE UPDATE
nres  = 800;             % cell of load resistor
G_dis(nres) = 1/RL/dx ;  % updated conductance at the load

% ARRAYS OF VOLTAGE AND CURRENT
vn(1:Ncells+1) = 0 ;  
in(1:Ncells)   = 0 ; 

v1_prev = 0 ;
vN_prev = 0 ;   % auxiliary variables for Mur's ABC

% TRANSMISSION LINE EQS' AND MUR'S CONSTANTS
vrel = (G_dis*dt)./(2.*C_dis) ;
irel = (R_dis*dt)./(2.*L_dis) ;

cin = (1 - irel)./(1 + irel)       ;
ain = dt./( L_dis.*dx.*(1.+irel) ) ; 

cvn = (1 - vrel)./(1 + vrel)       ;
avn = dt./( C_dis.*dx.*(1.+vrel) ) ; 

mur_coef = (s - 1)/(s + 1) ; 


% == MOVIE INITIALIZATION ==

x=linspace(dx, Ncells*dx, Ncells) ; 

figure('units','centimeters','position',[.1 .1 60 15])
axes('position', [.03 .06 0.95 0.9])
plot(x,vn(1:Ncells),'r',x,Zc*in,'b'),axis([0 .5 -1 1]);
xlabel('x (meters)');ylabel('V, Zc*I');
legend('Voltage','Zc*Current');

rect=get(gcf,'Position');
rect(1:2)=[0 0];

M=moviein(time_steps/15,gcf,rect);

vref = zeros(1, time_steps) ;


%
% Marching-in-time loop
%
 
tn = 0; % time
xn = 0; % space

for n = 1:time_steps 

  % ================
  %  Voltage update
  % ================ 

   vn(2:Ncells) = cvn(2:Ncells).*vn(2:Ncells) -...
                  avn(2:Ncells).*( in(2:Ncells) - in(1:Ncells-1) )  ;  
  

   % -- boundary conditions

   vn(1)        = v1_prev + mur_coef*(vn(2) - vn(1))               ; 
   vn(Ncells+1) = vN_prev + mur_coef*( vn(Ncells) - vn(Ncells+1) ) ; 
  
   v1_prev = vn(2)        ;
   vN_prev = vn(Ncells)   ; 

   % -- source (transparent)

   tn = tn + dt ;
   xn = xn + dx ;
   vn(5) = vn(5) + sin(k*xn - w*tn);

   % current update

   in(1:Ncells) = cin(1:Ncells).*in(1:Ncells) -... 
                  ain(1:Ncells).*( vn(2:Ncells+1) - vn(1:Ncells) ) ; 
   
%
% VISUALISATION
%

 if mod(n,15)==0;

   rtime=num2str(round(n*dt/1.0e-12));
   cstep=num2str(n) ;

   plot(x,vn(1:Ncells),'r',x,Zc*in,'b'),axis([0 Ncells*dx -1 1]);
   title(['time = ',rtime,' ps',' step=',cstep]);
   xlabel('x (meters)');ylabel('V, Zc*I');
   legend('Voltage','Zc*Current');

  M(:,n/15)=getframe(gcf,rect);
 
 end

end

movie(gcf,M,0,20,rect) 
