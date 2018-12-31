%
% 3-D FDTD code for the solution and field derivatives of the microstrip problem 
% in the paper "Application of the 3-D FDTD Method to the Analysis
% of Planar Microstrip Circuits", by Sheen at al.
%
% The original FDTD code author: Costas D. Sarris 
% The FDTD code is augmented with bi-complex number class developed by: K.A Liu
%          
clear all;

% ===================
%  program constants 
% ===================

light_speed =  2.99792458e8        ;  % speed of light in free space
mu_0  = 4.0*pi*1.0e-7              ;  % free space magnetic permeability
eps_0 = 1.0/((light_speed^2)*mu_0) ;  % free space dielectric permittivity

% =================
%  Mesh Definition
% =================

dx = 0.4064e-3  ;
dy = 0.4233e-3  ;
dz = 0.265e-3  ;

dt = 0.441e-12 ;         % corresponds to about 0.7 of the CFL limit

nxi = 80     ; 
nyi = 100    ;
nzi = 16     ;           % interior cells

nx = nxi + 1 ; 
ny = nyi + 1 ; 
nz = nzi + 1 ;           % add one cell in each direction for ABCs


% =========================================================
% Initialization for Bi-complex environment 
% =========================================================
[j1 j2]=initialize;
h = 1e-10;
Z = zeros(nx,ny,nz); O = ones(nx,ny,nz);

 
% =========================================================
% Gaussian excitation : exp( - ((t-t0)/Ts)^2 ) 
% =========================================================

Ts = 15e-12 ; 
t0 = 3*Ts   ; 
nmax = round(6*Ts/dt) ;   % max number of time steps for excitation

steps = 2000 ;            % time steps to be run


% =========================================================
% Dielectric permittivities
% =========================================================

epssubstrate = 2.2 + h*j2; 
epsr = bi(Z,Z);
epsr(1:nx, 1:ny, 1:3)  = epssubstrate; 
epsr(1:nx, 1:ny, 4:nz) = bi(1.0 , 0) ; 

% Define dielectric permittivities for each field node (generic loop). 

epsxx = epsr;
epsyy = epsr;
epszz = epsr;

% =========================================================
% Geometric properties
% =========================================================

DX = bi(dx*O,Z);
DY = bi(dy*O,Z);
DZ = bi(dz*O,Z);

DX(66 , 51:57 , 4 ) = dy+h*j1;

% ========================================
%  Multipliers for field update equations 
%         (optimized formulation)
% ========================================
aex = bi(Z,Z); 
aey = bi(Z,Z); 
aez = bi(Z,Z); 

aex(1:nx, 2:ny, 2:nz) = bi(dt,0)*bi(dx,0)/(bi(eps_0,0)*epsxx(1:nx, 2:ny, 2:nz)*DY(1:nx, 2:ny, 2:nz)*DZ(1:nx, 2:ny, 2:nz)) ;
aey(2:nx, 1:ny, 2:nz) = bi(dt,0)*bi(dy,0)/(bi(eps_0,0)*epsyy(2:nx, 1:ny, 2:nz)*DX(2:nx, 1:ny, 2:nz)*DZ(2:nx, 1:ny, 2:nz)) ;
aez(2:nx, 2:ny, 1:nz) = bi(dt,0)*bi(dz,0)/(bi(eps_0,0)*epszz(2:nx, 2:ny, 1:nz)*DX(2:nx, 2:ny, 1:nz)*DY(2:nx, 2:ny, 1:nz)) ;

ahx = bi(dt,0)*bi(dx,0)/(bi(mu_0,0) * DY *DZ );
ahy = bi(dt,0)*bi(dy,0)/(bi(mu_0,0) * DX *DZ );
ahz = bi(dt,0)*bi(dz,0)/(bi(mu_0,0) * DX *DY );

% ============================
%  Multipliers for Mur's ABCs
% ============================
mur_coefyx = bi(Z,Z); 
mur_coefyz = bi(Z,Z); 
mur_coefxy = bi(Z,Z); 
mur_coefxz = bi(Z,Z);
mur_coefz = bi(Z,Z);

mur_coefyx(1:nxi,1,2:nzi) = (bi(dt,0)*bi(light_speed,0)/sqrtbi(epsxx(1:nxi,2,2:nzi)) - bi(dy,0))/...
                            (bi(dt,0)*bi(light_speed,0)/sqrtbi(epsxx(1:nxi,2,2:nzi)) + bi(dy,0))  ;  


mur_coefyz(2:nxi,1,1:nzi) = (bi(dt,0)*bi(light_speed,0)/sqrtbi(epszz(2:nxi,2,1:nzi)) - bi(dy,0))/...
                            (bi(dt,0)*bi(light_speed,0)/sqrtbi(epszz(2:nxi,2,1:nzi)) + bi(dy,0))   ;

mur_coefxy(1,1:nyi,2:nzi) = (bi(dt,0)*bi(light_speed,0)/sqrtbi(epsyy(2,1:nyi,2:nzi)) - bi(dx,0))/...
                            (bi(dt,0)*bi(light_speed,0)/sqrtbi(epsyy(2,1:nyi,2:nzi)) + bi(dx,0))   ;

mur_coefxz(1,2:nyi,1:nzi) = (bi(dt,0)*bi(light_speed,0)/sqrtbi(epszz(2,2:nyi,1:nzi)) - bi(dx,0))/...
                            (bi(dt,0)*bi(light_speed,0)/sqrtbi(epszz(2,2:nyi,1:nzi)) + bi(dx,0))   ;

mur_coefz = bi((dt*light_speed-dz)./(dt*light_speed+dz)*O,Z)  ;

% FIELD ARRAYS
ex = bi(Z,Z); 
ey = bi(Z,Z); 
ez = bi(Z,Z); 

hx = bi(Z,Z); 
hy = bi(Z,Z); 
hz = bi(Z,Z); 

% auxiliary arrays for Mur's boundary conditions

exy0 = bi(Z,Z); 
ezy0 = bi(Z,Z); 

exym = bi(Z,Z); 
ezym = bi(Z,Z); 

eyx0 = bi(Z,Z); 
ezx0 = bi(Z,Z); 

eyxm = bi(Z,Z); 
ezxm = bi(Z,Z); 

exzm = bi(Z,Z); 
eyzm = bi(Z,Z); 


%
% === TIME STEPPING LOOP ===
%

for n =  1 : steps 

  load = n/steps*100
  % ELECTRIC FIELD UPDATE EQUATIONS
    
  ex(1:nxi,2:nyi,2:nzi) = ex(1:nxi, 2:nyi, 2:nzi)  +...
  aex(1:nxi,2:nyi, 2:nzi)*(hz(1:nxi,2:nyi,2:nzi)-hz(1:nxi,1:nyi-1,2:nzi)+...
                            hy(1:nxi,2:nyi,1:nzi-1)-hy(1:nxi,2:nyi,2:nzi) );

  ey(2:nxi,1:nyi,2:nzi) = ey(2:nxi,1:nyi,2:nzi)+...
  aey(2:nxi,1:nyi,2:nzi)*(hx(2:nxi,1:nyi,2:nzi)-hx(2:nxi,1:nyi,1:nzi-1)+...
                           hz(1:nxi-1,1:nyi,2:nzi)-hz(2:nxi,1:nyi,2:nzi));
                    
  ez(2:nxi,2:nyi,1:nzi) = ez(2:nxi,2:nyi,1:nzi)+...
  aez(2:nxi,2:nyi,1:nzi)*(hx(2:nxi,1:nyi-1,1:nzi)-hx(2:nxi,2:nyi,1:nzi)+...
                          hy(2:nxi,2:nyi,1:nzi)-hy(1:nxi-1,2:nyi,1:nzi));
                                              
  % Perfect electric conductors 

  ex(30:36,  1:51, 4) = bi(0,0) ;   % ex-nodes on patch
  ex(16:66,  51:57, 4) = bi(0,0) ;   % ex-nodes on microstrip
  ex(46:52,  57:100, 4) = bi(0,0);  
  
  ey(30:36,  1:51, 4) = bi(0,0) ;   % ex-nodes on patch
  ey(16:66,  51:57, 4) = bi(0,0) ;   % ex-nodes on microstrip
  ey(46:52,  57:100, 4) = bi(0,0);   % ey-nodes on microstrip


  % =============================================
  %  Boundary conditions at the excitation plane 
  % =============================================

  if n <= nmax    % source is still active, apply magnetic wall 

   ex(1:nxi,5,2:nzi) = ex(1:nxi, 5, 2:nzi)  +...
   aex(1:nxi,5, 2:nzi)*(bi(2,0)*hz(1:nxi,5,2:nzi)+...
                            hy(1:nxi,5,1:nzi-1)-hy(1:nxi,5,2:nzi) );

  ez(2:nxi,5,1:nzi) = ez(2:nxi,5,1:nzi)+...
   aez(2:nxi,5,1:nzi)*(bi(-2,0)*hx(2:nxi,5,1:nzi)+...
                           hy(2:nxi,5,1:nzi)-hy(1:nxi-1,5,1:nzi)); 

   % Gaussian pulse application

   ez(30:35, 5, 1:3 ) = bi( exp( - ((n*dt-t0)/Ts).^2 ), 0) ;

  end

   
  % Mur's ABC for the open boundary y=0

  ex(1:nxi,1,2:nzi) = exy0(1:nxi,1,2:nzi) +...
                      mur_coefyx(1:nxi,1,2:nzi)*(ex(1:nxi,2,2:nzi)-ex(1:nxi,1,2:nzi));

  ez(2:nxi,1,1:nzi) = ezy0(2:nxi,1,1:nzi) +...
                     mur_coefyz(2:nxi,1,1:nzi)*(ez(2:nxi,2,1:nzi)-ez(2:nxi,1,1:nzi)) ;

  % update auxiliary arrays for Mur's boundary conditions

  exy0(1:nxi,1,2:nzi) = ex(1:nxi,2,2:nzi) ;
  ezy0(2:nxi,1,1:nzi) = ez(2:nxi,2,1:nzi) ;

   
  % Mur's ABC for the open boundary y=ymax
  
  ex(1:nxi,ny,2:nzi) = exym(1:nxi,1,2:nzi) +...
  mur_coefyx(1:nxi,1,2:nzi)*(ex(1:nxi,nyi,2:nzi)-ex(1:nxi,ny,2:nzi));

  ez(2:nxi,ny,1:nzi) = ezym(2:nxi,1, 1:nzi) +...
  mur_coefyz(2:nxi,1,1:nzi)*(ez(2:nxi,nyi,1:nzi)-ez(2:nxi,ny,1:nzi));

  exym(1:nxi,1,2:nzi) = ex(1:nxi,nyi,2:nzi) ;
  ezym(2:nxi,1,1:nzi) = ez(2:nxi,nyi,1:nzi) ;


  % Mur's ABC for the open boundary x=0 

  ey(1,1:nyi,2:nzi) = eyx0(1,1:nyi,2:nzi) +...
  mur_coefxy(1,1:nyi,2:nzi)*(ey(2,1:nyi,2:nzi)-ey(1,1:nyi,2:nzi)) ;

  ez(1,2:nyi,1:nzi) = ezx0(1,2:nyi,1:nzi) +...
  mur_coefxz(1,2:nyi,1:nzi)*(ez(2,2:nyi,1:nzi)-ez(1,2:nyi,1:nzi)) ;  

  eyx0(1,1:nyi,2:nzi) = ey(2,1:nyi,2:nzi) ;  
  ezx0(1,2:nyi,1:nzi) = ez(2,2:nyi,1:nzi) ;


  % Mur's ABC for the open boundary x=xmax

  ey(nx,1:nyi,2:nzi) = eyxm(1,1:nyi,2:nzi) +...
  mur_coefxy(1,1:nyi,2:nzi)*(ey(nxi,1:nyi,2:nzi)-ey(nx,1:nyi,2:nzi)) ;

  ez(nx,2:nyi,1:nzi) = ezxm(1,2:nyi,1:nzi) +...
  mur_coefxz(1,2:nyi,1:nzi)*(ez(nxi,2:nyi,1:nzi)-ez(nx,2:nyi,1:nzi)) ;  

  eyxm(1,1:nyi,2:nzi) = ey(nxi,1:nyi,2:nzi) ;  
  ezxm(1,2:nyi,1:nzi) = ez(nxi,2:nyi,1:nzi) ;


  % Mur's ABC for the open boundary z=zmax
   
  ex(1:nxi,2:nyi,nz) = exzm(1:nxi,2:nyi,1) +...
  mur_coefz(1:nxi,2:nyi,1)*(ex(1:nxi,2:nyi,nzi)-ex(1:nxi,2:nyi,nz));  

  ey(2:nxi,1:nyi,nz) = eyzm(2:nxi,1:nyi,1) +...
  mur_coefz(2:nxi,1:nyi,1)*(ey(2:nxi,1:nyi,nzi)-ey(2:nxi,1:nyi,nz)); 
   
  exzm(1:nxi,2:nyi,1) = ex(1:nxi,2:nyi,nzi) ;
  eyzm(2:nxi,1:nyi,1) = ey(2:nxi,1:nyi,nzi) ;
  
  % MAGNETIC FIELD UPDATE EQUATIONS
  
  hx(1:nxi,1:nyi,1:nzi)=hx(1:nxi,1:nyi,1:nzi)+...
                   ahx(1:nxi,1:nyi,1:nzi)*(ey(1:nxi,1:nyi,2:nz)-ey(1:nxi,1:nyi,1:nzi)+...
                         ez(1:nxi,1:nyi,1:nzi)-ez(1:nxi,2:ny,1:nzi) );
                
  hy(1:nxi,1:nyi,1:nzi)=hy(1:nxi,1:nyi,1:nzi)+...
                   ahy(1:nxi,1:nyi,1:nzi)*(ex(1:nxi,1:nyi,1:nzi)-ex(1:nxi,1:nyi,2:nz)+...
                         ez(2:nx,1:nyi,1:nzi)-ez(1:nxi,1:nyi,1:nzi) );
                
  hz(1:nxi,1:nyi,1:nzi)=hz(1:nxi,1:nyi,1:nzi)+...
                   ahz(1:nxi,1:nyi,1:nzi)*(ex(1:nxi,2:ny,1:nzi) - ex(1:nxi,1:nyi,1:nzi)+...
                         ey(1:nxi,1:nyi,1:nzi)- ey(2:nx,1:nyi,1:nzi));
                     

    if~mod(n,1)           
        show(nx,ny,real(ez));
    end
    
   Re_E(n) = real(ez(49,90,3));
   Im1_E(n) = im1(ez(33,5,3));
   Im2_E(n) = im2(ez(49,90,3));
   Im12_E(n) = im12(ez(49,90,3));
end % end of time stepping loop

fname1 = sprintf('MCSD_re');
fname2 = sprintf('MCSD_i1');
fname3 = sprintf('MCSD_i2');
fname4 = sprintf('MCSD_i12');

csvwrite(fname1, Re_E);
csvwrite(fname2, Im1_E);
csvwrite(fname3, Im2_E);
csvwrite(fname4, Im2_E);

















