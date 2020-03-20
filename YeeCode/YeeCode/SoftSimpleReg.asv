
%This program simulates the Ez, Hx, Hy fields
%in a region which is probably a waveguide.




winstyle = 'docked';
% winstyle = 'normal';

set(0,'DefaultFigureWindowStyle',winstyle)
set(0,'defaultaxesfontsize',18)
set(0,'defaultaxesfontname','Times New Roman')
% set(0,'defaultfigurecolor',[1 1 1])

% clear VARIABLES;
clear
global spatialFactor;
global c_eps_0 c_mu_0 c_c c_eta_0
global simulationStopTimes;
global AsymForcing
global dels
global SurfHxLeft SurfHyLeft SurfEzLeft SurfHxRight SurfHyRight SurfEzRight



dels = 0.75;
spatialFactor = 1;

c_c = 299792458;                  % speed of light
c_eps_0 = 8.8542149e-12;          % vacuum permittivity
c_mu_0 = 1.2566370614e-6;         % vacuum permeability
c_eta_0 = sqrt(c_mu_0/c_eps_0);


tSim = 200e-15 %simulation time
% f = 400e12;
f = 230e12; %Frequency - lower frequencies (10e12) have more continuous peaks (no
% visible waves/gratings - High frequencies(400e12) propagate too fast -> small amplitude -> "travels right through")
lambda = c_c/f;

xMax{1} = 30e-6; %max X axis
nx{1} = 200; %200 points along x
ny{1} = 0.75*nx{1}; %3/4 of x

% xMax{2} = 30e-6; %max X axis
% nx{2} = 200; %200 points along x
% ny{21} = 0.5*nx{2}; %3/4 of x


Reg.n = 1;

mu{1} = ones(nx{1},ny{1})*c_mu_0; %magnetic permeability

epi{1} = ones(nx{1},ny{1})*c_eps_0; %electric permittivity
epi{1}(50:150,55:95)= c_eps_0*11.3;
% epi{1}(125:150,55:95)= c_eps_0*11.3;%"Inclusions" - Reflection waves and width/number of
% waves

sigma{1} = zeros(nx{1},ny{1});
sigmaH{1} = zeros(nx{1},ny{1});


% Reg.n = 2;
% mu{2} = ones(nx{2},ny{2})*c_mu_0; %magnetic permeability
% 
% epi{2} = ones(nx{2},ny{2})*c_eps_0; %electric permittivity
% epi{2}(50:150,55:95)= c_eps_0*11.3;
% 
% sigma{2} = zeros(nx{2},ny{2});
% sigmaH{2} = zeros(nx{2},ny{2});




dx = xMax{1}/nx{1}; %discretization of space
dt = 0.25*dx/c_c; %timestep
nSteps = round(tSim/dt*2); 
yMax = ny{1}*dx;
nsteps_lamda = lambda/dx %steps along frequency

movie = 1; %run as a "movie"
Plot.off = 0;
Plot.pl = 0;
Plot.ori = '13';
Plot.N = 100;
Plot.MaxEz = 1.1;
Plot.MaxH = Plot.MaxEz/c_eta_0;
Plot.pv = [0 0 90];
Plot.reglim = [0 xMax{1} 0 yMax];

%bc structure is used to track boundary conditions
%and keep track of positions for parameters to plot.
%
%bc{1}.s(1)is setting up the wave in the waveguide. 
bc{1}.NumS = 1; %Unsure of effect
bc{1}.s(1).xpos = nx{1}/(4) + 1; %modifies source x position
bc{1}.s(1).type = 'ss'; %changed to 's' no noticable effect
bc{1}.s(1).fct = @PlaneWaveBC;%turns 'paras' into a plane wave

% bc{2}.NumS = 1; %Unsure of effect
% bc{2}.s(2).xpos = nx{1}/(4) + 1; %modifies source x position
% bc{2}.s(2).type = 'ss'; %changed to 's' no noticable effect
% bc{2}.s(2).fct = @PlaneWaveBC;%turns 'paras' into a plane wave


% mag = -1/c_eta_0;
mag = 3; %magnitude of wave
phi = 0; %Phase shift - effects magnitudes upon interference
omega = f*2*pi;
betap = 0;
t0 = 30e-15;
st = -0.05; %waves stay in plane longer - not much loss
% st = 15e-15; %lets waves propagate more freely. More interferenec. Longer interference time
s = 0;
y0 = yMax/2;
sty = 1.5*lambda;
bc{1}.s(1).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};

Plot.y0 = round(y0/dx);

% bc{2}.s(2).paras = {mag,phi,omega,betap,t0,st,s,y0,sty,'s'};
% Plot.y0 = round(y0/dx);

%changing all these to 'e' seems to create some kind of 
%resosnant modes within the material (more reflections/interference zones?)
bc{1}.xm.type = 'a';
bc{1}.xp.type = 'a';
bc{1}.ym.type = 'a';
bc{1}.yp.type = 'a';

% bc{2}.xm.type = 'a';
% bc{2}.xp.type = 'a';
% bc{2}.ym.type = 'a';
% bc{2}.yp.type = 'a';

pml.width = 20 * spatialFactor;
pml.m = 3.5;

Reg.n  = 1;
Reg.xoff{1} = 0;
Reg.yoff{1} = 0;

% Reg.n  = 1;
% Reg.xoff{2} = 0;
% Reg.yoff{2} = 0;

RunYeeReg






