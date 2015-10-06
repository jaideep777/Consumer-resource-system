%close all; clear all;

% define these variables externally when running model iterator
EN__ = sprintf('Het%g-%g', 0.15,3);
KI__ = 'L';
Nc__ = 100;
hc__ = .5;
RTc__ = 15;
Kd_sd__ = 8;
cDisp__ = 0.1;
% bHarvest__ = .2000;
% rImit__ = 0.02;

b_imitate_hc__ = 1;
b_imitate_RT__ = 0;
b_imitate_Kd__ = 1;

b_graphics__ = 0;
b_video_out__ = 0;
% -------------------------------------------------------------------------

% sim params
L = 100;                        % domain size
nx = 200;   % grid size
dt = 0.1;   % time step
Tgen = 750000;

% resource 
r = 0.2 + 0.15*(2*perlin2D(nx,3)-1);       % growth rate (must be scalar or NxN matrix)
K = 50; %+ 25*(2*perlin2D(nx)-1);      % carrying capacity (must be scalar or NxN matrix)
D = 0.;                         % resource diffusion constant
R = ones(nx,nx).*K;             % resource grid
R(1,1) = 50;                    % delta function to test diffusion

% grid definition
dL = L/nx;                      % cell size
grid_centres = dL/2:dL:nx*dL;   % cell centre coordinates (same for X and Y)
[grid_x, grid_y] = meshgrid(grid_centres, grid_centres);

% consumers
Nc = Nc__;                        % number of consumers

xc = rand(Nc,1)*L;              % initial positions of consumers
yc = rand(Nc,1)*L;              % initial positions of consumers

hc = hc__*ones(Nc,1); % sort(randi([1,2], [Nc,1])*0.2);   % 0.25*ones(Nc,1); % 1+0.25*randn(Nc,1); % harvesting rates of consumers
RTc = RTc__*ones(Nc,1);            % randi([1,2], [Nc,1])*10; % 10+2*randn(Nc,1);  % threshold resource level for dispersal
Kd_sd = Kd_sd__*ones(Nc, 1);          % SD of dispersal kernel

Ke_sd = 4;                      % SD of exploitation kernel
Ke_cutoff = 12;                 % cutoff distance for exploitation kernel
Ke_nmax = floor(Ke_cutoff/dL);  % cutoff number of grid cells of exploitation kernel

Rc = zeros(Nc, 1);
nDc = zeros(Nc,1);

% Imitation
rImit = rImit__;                    % imitation rate (probability of immitation per unit time)

% payoff
Tw = 20;                        % Window size for moving average for payoff
Vc_window = zeros(Nc, Tw);      % payoff window (contents are not moved. It is updated in cyclic manner)
Vc = zeros(Nc,1);               % current averge payoff 
Vbase = 0;
CDisp = cDisp__;                      % cost of dispersal (per unit distance)
BRes = bHarvest__;                     % benefit from resource accummulation

% -------------------------------------------------------------------------

% Initialization
ixc = pos2id(xc, dL);           % initial grid indices of consumers
iyc = pos2id(yc, dL); 

% single exploitation kernel
[Ke_x, Ke_y] = meshgrid((-Ke_nmax:Ke_nmax)*dL, (-Ke_nmax:Ke_nmax)*dL);
Ke = exp(-(Ke_x.^2+Ke_y.^2)/Ke_sd^2);
% Ke = 1-(Ke_x.^2+Ke_y.^2 > 4*Ke_sd^2);
Ke = 1-(Ke_x.^2+Ke_y.^2)/Ke_sd^2;
Ke(Ke<0)=0;

% calculate resource consumption rate for entire grid
Ke_all = updateKernels2(ixc, iyc, hc, Ke, nx);

% -------------------------------------------------------------------------

% timeseries for storing data
Tskip = 100;
opsize = Tgen/Tskip;
h_vec = zeros(Nc, opsize);
RT_vec = zeros(Nc, opsize);
Kd_vec = zeros(Nc, opsize);
V_vec = zeros(Nc, opsize);
Rc_vec = zeros(Nc, opsize);
zD_vec = zeros(Nc, opsize);
nD_vec = zeros(Nc, opsize);
R_percap = zeros(opsize, 1);
dV_cumm = [];

len_d_cumm = zeros(Nc, 1);
n_d_cumm = zeros(Nc, 1);
% -------------------------------------------------------------------------

% figure; axis([0 L 0 L]); imagesc(grid_centres, grid_centres, R'); colormap(gray);
% aviobj = avifile('movie1.avi', 'fps', 20);
% Main run1
tic;
for t = 1:Tgen

    % resource consumed 
    Rc = calcResConsumed(R, Ke, ixc, iyc, hc)*dt;
    
    % resource left
    R = R + dt*(D*laplacian(R) + r.*R.*(1-R./K) - Ke_all.*R);
    R(R<0)=0;

    % check dispersal
    points_1d = sub2ind(size(R), ixc, iyc);
    p_disperse = 1./(1+exp(10*(R(points_1d) - RTc)));  % sign(R(points_1d) - RTc)
    l_disperse = rand(length(points_1d),1) < p_disperse;        % flag - disperse or not
%    l_disperse = sign(R(points_1d) - RTc) < 0;
    
    % disperse
    len_dispersal = abs(randn(Nc, 1)).*Kd_sd;           % distance to disperse (|gaussian with sd Kd_sd|) 
    theta_dispersal = rand(Nc, 1)*2*pi;                 % direction of dispersal    
    x_disp = len_dispersal.*cos(theta_dispersal);
    y_disp = len_dispersal.*sin(theta_dispersal);

    xc = xc + x_disp.*l_disperse;
    yc = yc + y_disp.*l_disperse;
    xc = make_periodic(xc, L);
    yc = make_periodic(yc, L);

    ixc = pos2id(xc, dL);
    iyc = pos2id(yc, dL);

    % update exploitation kernels
    Ke_all = updateKernels2(ixc, iyc, hc, Ke, nx);

    % calculate payoffs
    Vc_window(:, make_periodic_index(t,Tw)) = 1/dt*(Rc*BRes - len_dispersal.*l_disperse*CDisp - 0.08*hc.^2);
    Vc = mean(Vc_window, 2);
    Vc = Vc + Vbase;
    
    % Imitate
    n_imitation_events = rbinom(Nc, rImit*dt); % binornd(Nc, rImit*dt); % wrote my own because stats-toolbox licenses were over.
    for i=1:n_imitation_events
        % select individual who will copy
        id_who = randi(Nc);

        % LOCAL IMITATION
        % select closest individual to imitate
        dx_others1 = abs(ixc - ixc(id_who));
        dx_others2 = nx - dx_others1;
        dx_others = min(dx_others1, dx_others2);
        dy_others1 = abs(iyc - iyc(id_who));
        dy_others2 = nx - dy_others1;
        dy_others = min(dy_others1, dy_others2);
        r_others = sqrt(dx_others.^2+dy_others.^2);
        r_others = r_others([1:id_who-1, id_who+1:end]);
        [minr, idmr] = min(r_others);
        id_whom = idmr; 
        
        dV = Vc(id_whom) - Vc(id_who);          % payoff difference
%         dV_cumm = [dV_cumm, dV];
%         x1 = [-1e6,  -0.2,  0.2, 1e6]';
%         y1 = [   0,     0,    1,   1]';
%         imitation_prob = interp1q(x1,y1,dV); 
        imitation_prob =  0.5*(1+sign(dV));      % imitation probability (currently a step function)
        if (rand < imitation_prob)
            if b_imitate_hc__
                % imitate hc
                hc(id_who) = hc(id_whom) + 0.02*randn();
                if hc(id_who) < 0
                    hc(id_who) = 0;
                end
            end
            
            if b_imitate_Kd__
                % imitate Kd
                Kd_sd(id_who) = Kd_sd(id_whom) +  0.2*randn();
                if Kd_sd(id_who) < 0
                    Kd_sd(id_who) = 0;
                end
            end
            
            if b_imitate_RT__
                % imitate RT
                RTc(id_who) = RTc(id_whom) + 1*randn();
                if RTc(id_who) < 0
                    RTc(id_who) = 0;
                end
            end
                
        end
    end

    % update timeseries
    len_d_cumm = len_d_cumm + len_dispersal.*l_disperse;
    n_d_cumm = n_d_cumm + l_disperse;
    if mod(t,Tskip)==0
        it = ceil(t/Tskip);

        h_vec(:, it) = hc;
        RT_vec(:, it) = RTc;
        Kd_vec(:, it) = Kd_sd;

        Rc_vec(:, it) = Rc;
        V_vec(:, it) = Vc;
        R_percap(it) = sum(sum(R))/Nc;
        zD_vec(:, it) = len_d_cumm/Tskip;
        nD_vec(:, it) = n_d_cumm/Tskip;
        
        % reset cummulative variables
        len_d_cumm = zeros(Nc, 1);
        n_d_cumm = zeros(Nc, 1);
    end
    
    % visualize 
    if b_graphics__
        if mod(t,100)==0
            clf();
            axis([0 L 0 L]); imagesc(grid_centres, grid_centres, R'); colormap(gray);
            % figure; axis([0 L 0 L]); imagesc(grid_centres, grid_centres, Ke_all'); colormap(gray);
            hold on;
            % plot(xc,yc,'g*', 'markerSize', 5);
            plot((ixc-0.5)*dL,(iyc-0.5)*dL,'rs', 'markerSize', 5);
            hold off;
            pause(0.05)
%             frame = getframe(gca);
%             aviobj = addframe(aviobj,frame);
        end
    end
    
    % print progress
    if (mod(t,10000)==0)
        fprintf('t = %d\n', t)
    end

end
toc;

plot_stuff;
% aviobj = close(aviobj);