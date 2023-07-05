
%%****************************************************************************************************
%*****************************************************************************************************
%                         This code impelements a dipole antenna
%                         imput impedance simulation using a vectorized 3D FDTD
%                         truncated with PMLs. The output is a figure of
%                         the imput impedance versus frequency and relative
%                         lenght, as well as plots of the current
%                         distribution around the dipole.
%                         
%                         The simulation parameters are:
%                         - dipole lenght: 250 mm
%                         - dx = dy = 1 mm ; dz = 2 mm 
%                         - Volume of simulation: 29 mm X 29 mm X 306 mm
%                         - dt = 1.92 Ps
%                        
%                         Author: Felix Munoz, 2023
%
%****************************************************************************************************
%****************************************************************************************************
%clear
tic
%% Fundamental constants:
e_0 = 8.854e-12;
u_0 = 4e-7*pi;
c_0 = 1/sqrt(e_0*u_0);

%****************************************************************************************************
%                                 Simulation Parameters
%****************************************************************************************************

% Step size. 1 mm in the x and y direction, and 2 mm in the z direction.
[dx,dy] = deal(1e-3); 
dz = 2e-3;

% Time step and number of time steps using Courant criterion:
dt = dx/(sqrt(3)*c_0); 
t_steps = 25e3; % [Can reduce to 20e3, but more course freq resolution]

% Antenna parameters:
length = 250e-3; % Lenght of antenna is 250 mm.
len = floor((length/dz)/2); % antenna lenght in number of cells is 2*len + 1

% Simulation size. 
% In all dimensions, allow for 8-layer PML on each side and six cells between the antenna and the PML boundary.
[Nx,Ny] = deal(1 + 16 + 12); 
Nz = 2*len + 1 + 16 + 12; 


% Voltage source position:
[Sx,Sy,Sz] = deal(round(Nx/2),round(Ny/2),round(Nz/2)); % Feed point center. 

% If visualize in on, a movie of the fields will be shown. 
visualize = 0;
%% Voltage Source
%****************************************************************************************************
%           Define and plot source in time and frequency domain.
%****************************************************************************************************

% Center wavelenght of the excitation. Want frequency content around
% 400-800 MHz.
l_0 = 500e-3;    % center
l_u = 1000.9e-3; % upper
l_l = 100e-3;    % lower 

% Center frequency and spread of the Gaussian
w_0 = (2*pi*c_0)/l_0; %
sigma = (2/w_0)*(l_0/(l_u - l_l)); 

% Modulated Gaussian voltage source:
source = @(t) exp( - ( (t - 4*sigma)./sigma ).^2   ) .* sin(w_0*(t - 4.*sigma));


% Time array
t = (1:t_steps)*dt; 

% Plot source in time:
figure;
subplot(1,2,1)
plot(t/1e-9,source(t)); % Plot source in nano-sec:
xlim( [0 15])
xlabel('Time (ns)')
ylabel('Pulse Amplitude (V)')
title('Modulated Gaussian Pulse vs Time')
grid on 


% Frequeny content of source: 
source_fft = abs(fft(source(t)));
% sampling frequency:
fs = 1/dt; 
N = size(source_fft,2); % size of source 
freq = (0:N-1).* fs/N; % frequency bins

% Plot Voltage source in frequency domain:
subplot(1,2,2)
plot(freq/1e6,source_fft) 
xlim([0 3e3])
xlabel('Freq MHz')
ylabel('|f(w)| (V)')
title('Modulated Gaussian Pulse vs Frequency')
grid on 


%% Initialize Fields
%****************************************************************************************************
%                   Initialize All Fields and Current 
%****************************************************************************************************


% Initialize Fields
[Ex,Dx,Dxp] = deal(zeros(Nx-1,Ny,Nz));
[Ey,Dy,Dyp] = deal(zeros(Nx,Ny-1,Nz));
[Ez,Dz,Dzp] = deal(zeros(Nx,Ny,Nz-1));
[Hx,Bx,Bxp] = deal(zeros(Nx,Ny-1,Nz-1));
[Hy,By,Byp] = deal(zeros(Nx-1,Ny,Nz-1));
[Hz,Bz,Bzp] = deal(zeros(Nx-1,Ny-1,Nz));

% Initialize epsilons and mus
[e_Ex] = deal(ones(Nx-1,Ny,Nz))*e_0;
[e_Ey] = deal(ones(Nx,Ny-1,Nz))*e_0;
[e_Ez] = deal(ones(Nx,Ny,Nz-1))*e_0;
[u_Hx] = deal(ones(Nx,Ny-1,Nz-1)*u_0);
[u_Hy] = deal(ones(Nx-1,Ny,Nz-1)*u_0);
[u_Hz] = deal(ones(Nx-1,Ny-1,Nz)*u_0);

% Initialize current and voltage arrays:
I_z = zeros(1,t_steps);
V_z = zeros(1,t_steps);
I_all_z= ones(2*len + 1,t_steps);

%% PML Coefficients
%****************************************************************************************************
%                               Calculte PML Coefficients:
%%***************************************************************************************************
unshifted = zeros(1,Nx);     % Store shifted sigmas for x and y
shifted   = zeros(1,Nx-1);   % Store unshifted sigmas for x and y

unshifted_z = zeros(1,Nz);     % Store shifted sigmas  for z
shifted_z   = zeros(1,Nz-1);   % Store unshifted sigmas for z

L = 8; % Size of PML 
R = 1e-8; % Reflection coefficient
p_max = -((3+1)/4)* ((c_0)/(L*dx))* log(R); % max sigma 

%****************************************************************************************************
%           Calculate shifted and unshifted sigmas. Since the simulation has dx = dy,
%           sigma_x and sigma_y are the same size. sigma_z is of different size. 
%****************************************************************************************************

for i=1:L
	xxn=(L+1-i)/(L); 
	xn = p_max*(xxn^3);    % Start at (L/L)^3 and end at (1/L)^3

    unshifted(i) = xn;
	unshifted(Nx-i+1) = xn;

    unshifted_z(i) = xn;
	unshifted_z(Nz-i+1) = xn;

    
    % Shifted by 0.5
    xxn=(L+1-i-0.5)/(L); % Start at [(L - 0.5)/L]^3 and end at (0.5/L)^3
	xn=p_max*(xxn^3);

	shifted(i)=xn;
	shifted(Nx-i)=xn;

    shifted_z(i)=xn;
	shifted_z(Nz-i)=xn;
end

% Uncomment to remove PML:
%unshifted = zeros(1,Nx); 
%shifted   = zeros(1,Nx-1); 

%****************************************************************************************************
%                    Pre-calculate update equations constants
%****************************************************************************************************

% Unshifted for x and y:
u_p = 1 + dt.*unshifted./2;
u_m = 1 - dt.*unshifted./2;
% Shifted for x and y:
s_p = 1 + dt.*shifted./2;
s_m = 1 - dt.*shifted./2;

% Unshifted for z:
z_u_p = 1 + dt.*unshifted_z./2;
z_u_m = 1 - dt.*unshifted_z./2;
% Shifted for z:
z_s_p = 1 + dt.*shifted_z./2;
z_s_m = 1 - dt.*shifted_z./2;


% Reshape to appropriate dimensions for element-wise multiplication:
% sigma_x constants:
x_u_p = u_p';
x_u_m = u_m';
x_s_p = s_p';
x_s_m = s_m';
% sigma_y constants:
y_u_p = u_p;
y_u_m = u_m;
y_s_p = s_p;
y_s_m = s_m;
% sigma z constants:
z_u_p = reshape(z_u_p,1,1,size(z_u_p,2));
z_u_m = reshape(z_u_m,1,1,size(z_u_m,2));
z_s_p = reshape(z_s_p,1,1,size(z_s_p,2));
z_s_m = reshape(z_s_m,1,1,size(z_s_m,2));

%% Initialize plots for visualization:
if visualize == 1
figure
end


%% Begin Time stepping.
% Initalize waitbar for progress bar.
h = waitbar(0,'Time stepping...');
%****************************************************************************************************
%                                KING Update Equations 
%****************************************************************************************************
for idx = 1:t_steps
    % Update Hx: [shifted y,z] [N,N-1,N-1]
    Bxp = Bx;
    Bx = Bx.*(y_s_m./y_s_p) + (dt./y_s_p).* (diff(Ey,1,3)/dz - diff(Ez,1,2)/dy);
    Hx = Hx.*(z_s_m./z_s_p) + (1./z_s_p).*( x_u_p.*Bx - x_u_m.*Bxp).*(1./u_Hx);
    
    % Update Hy: [shifted x,z] [N-1,N,N-1]
    Byp = By;
    By = By.*(z_s_m./z_s_p) + (dt./z_s_p).*(diff(Ez,1,1)/dx - diff(Ex,1,3)/dz);
    Hy = Hy.*(x_s_m./x_s_p) + (1./x_s_p).*( y_u_p.*By - y_u_m.*Byp).*(1./u_Hy);
    
    % Update Hz: [Shifted x,y] [N-1,N-1,N]
    Bzp = Bz;
    Bz = Bz.*(x_s_m./x_s_p) + (dt./x_s_p).*(diff(Ex,1,2)/dy - diff(Ey,1,1)/dx);
    Hz = Hz.*(y_s_m./y_s_p) + (1./y_s_p).*( z_u_p.*Bz - z_u_m.*Bzp).*(1./u_Hz);
    
    %*******************************************************************************
    %               Update Ex: [Shifted x] [N-1,N,N] [Cut:s_y,s_z] 
    %*******************************************************************************
    Dxp = Dx(:,2:end-1,2:end-1);
    Dx(:,2:end-1,2:end-1) = Dx(:,2:end-1,2:end-1).*(y_u_m(2:end-1)./y_u_p(2:end-1)) + ...
        (dt./y_u_p(2:end-1)).*( diff(Hz(:,:,2:end-1),1,2)/dy - diff(Hy(:,2:end-1,:),1,3)/dz);

    
    Ex(:,2:end-1,2:end-1) = Ex(:,2:end-1,2:end-1).*(z_u_m(2:end-1)./z_u_p(2:end-1)) + ...
        (1./z_u_p(2:end-1)).*( x_s_p.*Dx(:,2:end-1,2:end-1) - x_s_m.*Dxp ).*(1./e_Ex(:,2:end-1,2:end-1));
    
    %*****************************************************************************************************
    %           Update Ey: [Shifted y] [N,N-1,N] [Cut:s_x,s_z]
    %*****************************************************************************************************
    Dyp = Dy(2:end-1,:,2:end-1);
    Dy(2:end-1,:,2:end-1) = Dy(2:end-1,:,2:end-1).*(z_u_m(2:end-1)./z_u_p(2:end-1)) + ...
        (dt./z_u_p(2:end-1)).*( diff(Hx(2:end-1,:,:),1,3)/dz - diff(Hz(:,:,2:end-1),1,1)/dx);

    
    Ey(2:end-1,:,2:end-1) = Ey(2:end-1,:,2:end-1).*(x_u_m(2:end-1)./x_u_p(2:end-1)) + ...
        (1./x_u_p(2:end-1)).*( y_s_p.*Dy(2:end-1,:,2:end-1) - y_s_m.*Dyp ).*(1/e_Ey(2:end-1,:,2:end-1));
    
    
    %*****************************************************************************************************
    %               Update Ez: [Shifted z] [N,N,N-1] [Cut:s_x,s_y]
    %*****************************************************************************************************
    Dzp = Dz(2:end-1,2:end-1,:);
    Dz(2:end-1,2:end-1,:) = Dz(2:end-1,2:end-1,:).*(x_u_m(2:end-1)./x_u_p(2:end-1)) + ...
        (dt./x_u_p(2:end-1)).*( diff(Hy(:,2:end-1,:),1,1)/dx - diff(Hx(2:end-1,:,:),1,2)/dy); 


    % Force dipole antenna:
    Dz(Sx,Sy,Sz+1:Sz+len) = 0;
    Dz(Sx,Sy,Sz-len:Sz-1) = 0;

    
    Ez(2:end-1,2:end-1,:) = Ez(2:end-1,2:end-1,:).*(y_u_m(2:end-1)./y_u_p(2:end-1)) + ...
        (1./y_u_p(2:end-1)).*( z_s_p.*Dz(2:end-1,2:end-1,:) - z_s_m.*Dzp ).*(1./e_Ez(2:end-1,2:end-1,:));


    %*****************************************************************************************************
    %                       Voltage feed and current probing:
    %*****************************************************************************************************
    Ez(Sx,Sy,Sz) = -source((idx)*dt)/dz; 
  

    % Contour to find the current just above the feed gap using Amperes's
    % countour law:
    I_z(idx) = Hy(Sx,Sy,Sz+1)*dy - Hx(Sx,Sy,Sz+1)*dx - Hy(Sx-1,Sy,Sz+1)*dy + Hx(Sx,Sy-1,Sz+1)*dx;
    V_z(idx) = source((idx)*dt); 

    % Finding the current all along the lenght of the antenna:
    index = 1;
    for l = Sz-len:Sz+len
        I_all_z(index,idx) = Hy(Sx,Sy,l)*dy - Hx(Sx,Sy,l)*dx - Hy(Sx-1,Sy,l)*dy + Hx(Sx,Sy-1,l)*dx;
        index = index + 1;
    end
    
   
    
   %*****************************************************************************************************
   %                            Vizualize all planes of Ex:
   %*****************************************************************************************************
   % 
   if visualize == 1

       subplot(1,3,1)
       imshow(transpose(squeeze((Ez(Sx+1,:,:)))),[],'InitialMagnification',4000);
       colormap('parula') 
       title("Above X") 

       subplot(1,3,2)
       imshow(transpose(squeeze((Ez(Sx,:,:)))),[],'InitialMagnification',4000);
       colormap('parula') 
       title("Through X") 
    
    
       subplot(1,3,3)
       imshow((squeeze(Ez(:,:,Sz))),[],'InitialMagnification',4000);
       colormap('parula')
       title("Through Z")

   end


    drawnow % Need, otherwise images won't update as it loops. 
    waitbar(idx/t_steps,h);
end
close(h)
toc

% *****************************************************************************************************
%                                       POST-PROCESSING
%*******************************************************************************************************


%% Post-processing: Imput Impedance
% *****************************************************************************************************
%                               Calculate imput impedance:
%*******************************************************************************************************
% Calculate frequency bins:
fs = 1/dt; % sampling freq
N = size(I_z,2); % size of source array
freq = (0:N-1).* fs/N; % frequencies


% Calculate I(w) and V(w)
I_z_w = fft(I_z);
V_z_w = fft(V_z);


% Calculate and plot Plot Z(w)
Z_w = V_z_w./I_z_w;


% Plot real and imaginary parts of Z(w) versus frequency
fig = figure;
fig.Position = [440 305 751 393];
plot(freq(2:end)/1e6,imag(Z_w(2:end)));
hold on 
plot(freq(2:end)/1e6,real(Z_w(2:end)));
xlim([370 4000])
title('Dipole Imput Impedance (l=250 mm, a = 1 mm)')
% zero line
plot(xlim,[0,0])
grid on 
grid minor
xlabel('Frequency (MHz)')
ylabel('Input Impedance (Ohms)')
legend('Im(Z)','R(Z)')

%interpolate: Use this only to visualize where resonance occurs
% freq2 = (0:0.03125/2:N-1).* fs/N; % for interpolation
% imag_i = interp1(freq(2:end),imag(Z_w(2:end)),freq2(2:end),'spline');
% real_i = interp1(freq(2:end),real(Z_w(2:end)),freq2(2:end),'spline');
% plot(freq2(2:end)/1e6,imag_i,':.');
% plot(freq2(2:end)/1e6,real_i,':.');



% Plot real and imaginary parts of Z(w) versus l/lamda 
fig = figure;
fig.Position = [440 305 751 393];
plot(freq(2:end).*(length/c_0),imag(Z_w(2:end)));
hold on 
plot(freq(2:end).*(length/c_0),real(Z_w(2:end)));
xlim([.3 3])
title('Dipole Imput Impedance (l=250 mm, a = 1 mm)')
% zero line
plot(xlim,[0,0])
grid on 
grid minor
ylabel('Input Impedance (Ohms)')
xlabel('$\frac{l}{\lambda}$','Interpreter','latex')
legend('Im(Z)','R(Z)')

% interpolate:
% plot(freq2(2:end).*(length/c_0),imag_i,':.');
% plot(freq2(2:end).*(length/c_0),real_i,':.');



%% Post-processing: Current Distribution along dipole

% Frequency for which the antenna is half the waveleght:
lamda = 2*length;
f_1 = (c_0)/lamda;

% Frequency for which the antenna is 3/2 the waveleght:
lamda = (2/3)*length;
f_2 = (c_0)/lamda;

% Frequency for which the antenna is twice the waveleght:
lamda = (1/2)*length;
f_3 = (c_0)/lamda;

% Take the fft of total current to go into frequency domain:
I_all_w = fft(transpose(I_all_z));

% Extract the frequencies of interest (frequencies above):
[~,idx] = min( abs(freq-f_1) );
I_f_1 = I_all_w(idx,:);

[~,idx] = min( abs(freq-f_2) ); 
I_f_2 = I_all_w(idx,:);

[~,idx] = min( abs(freq-f_3) ); 
I_f_3 = I_all_w(idx-2,:);


% Plot to visualize current distribution:

% Antenna distance array in mm:
dist = (-len:len)*dz*1000;
% Index of feed gap (to exclude from plot)
mid =ceil(size(I_f_1,2)/2);

% Plot current distributions:
figure;
plot(dist(1:mid-1),abs(I_f_1(1:mid-1)))
hold on 
plot(dist(mid+1:end),abs(I_f_1(mid+1:end)))
title('$ \bf{Current \: Distribution \: for \: l = \lambda/2 } $','Interpreter','latex')
xlabel('Lenght (mm)')
ylabel('Current (A)')
grid on 
grid minor 

figure;
plot(dist(1:mid-1),abs(I_f_2(1:mid-1)))
hold on 
plot(dist(mid+1:end),abs(I_f_2(mid+1:end)))
title('$ \bf{Current \: Distribution \: for \: l = 3\lambda/2 } $','Interpreter','latex')
xlabel('Lenght (mm)')
ylabel('Current (A)')
grid on 
grid minor 


figure;
plot(dist(1:mid-1),abs(I_f_3(1:mid-1)))
hold on 
plot(dist(mid+1:end),abs(I_f_3(mid+1:end)))
title('$ \bf{Current \: Distribution \: for \: l = 2\lambda } $','Interpreter','latex')
xlabel('Lenght (mm)')
ylabel('Current (A)')
grid on 
grid minor 


%% Post-processing: Return loss plot:
figure;
z0 = 50; % impedance of generator
gamma = (Z_w - z0)./(Z_w + z0);
plot(freq(2:end).*(length/c_0),20*log10(abs(gamma(2:end))))
title('Return Loss for Dipole l = 250 mm')
xlim([0 3])
ylabel('Return Loss (dB)')
xlabel('$\frac{l}{\lambda}$','Interpreter','latex')
grid on 
grid minor 

           