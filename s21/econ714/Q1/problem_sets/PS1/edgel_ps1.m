%{
    This file is used to conduct all necessary tasks for problem set 1 of
    the first quarter of Econ 714

    Date created:  27 Jan 2021
    Last modified: 28 Jan 2021
    Author: Danny Edgel
%}

% clean workspace
clc; clear

% set shooting method convergence variables
tol             = 1e-3;  % convergence tolerance
numcands        = 3e4;  % size of consumption grid increments
traj_periods    = 500;   % number of periods to test for convergence to SS


%% Define parameters
sigma   = 1;
alpha   = 1/3;
beta    = 0.99^(1/12);
delta   = 0.01;
T       = 12;
D       = 1;

%% Calculate model's steady-state

% capital, from Euler equation
k_ss = (((1/beta) - 1 + delta)/alpha)^(1/(alpha-1));

% consumption, from the resource constraint
c_ss = k_ss^alpha - delta*k_ss; % NOTE: assuming D=0 in the steady-state


%% Apply shooting method

% loop through consumption candidates, attempting convergence to the
% steady-state
mindist = 10;   % initialize minimum distance
j       = 0;    % convergence attempt count

while mindist > tol && j < 100
    fprintf('\nLoop %d; Current tolerance level: %d\n', j + 1, mindist)
    
    % after first loop, update consumption candidates
    if j > 0; cmin = .85*cmin + .15*c_cand(minind); else; cmin = 0.001; end
    
    % initialize vector of possible consumption levels
    inccon = (c_ss - cmin)/numcands;
    c_cand = cmin:inccon:c_ss;

    % initialize distance and trajectory vectors
    dist    = zeros(1, length(c_cand)); 
    c_traj  = zeros(length(c_cand), traj_periods);
    k_traj  = zeros(length(c_cand), traj_periods);
    
    % fill in every first period of the capital and consumption trajectories
    % with the consumption candidate and the capital steady state minus D,
    % respectively
    c_traj(:, 1) = c_cand';
    k_traj(:, 1) = k_ss;
    
    for i = 1:length(c_cand)
        fprintf('  Testing candidate %d of %d...\n',i,length(c_cand))

        % fill in trajectory of current consumption candidate
        for t = 2:traj_periods

            % extract 'current' capital and consumption from last trajectory
            % period
            c = c_traj(i, t-1);
            k = k_traj(i, t-1);

            % choose period j's capital and consumption using Euler equation
            % and resource constraint; if next period is T, subtract D from
            % knext
            
            knext = k^alpha + (1-delta)*k - c;
            if t == T; knext = knext - D; end
            
            cnext = (beta*(alpha*knext^(alpha-1)+1-delta))^(1/sigma)*c;
            
            if isreal(knext) && isreal(cnext) && ~isnan(knext) && ~isnan(cnext)
                k_traj(i,t) = knext;
                c_traj(i,t) = cnext;
            else
                k_traj(i, t) = k;
                c_traj(i, t) = c;
            end

        end % end j loop (trajectory periods)

        % calculate distance between final trajectory period and the
        % steady-state
        dist(i) = norm(abs([k_traj(i, end), c_traj(i, end)] - [k_ss, c_ss]));

    end % end i loop (consumption candidates)

    % extract consumption candidate that came closest to the steady-state
    [mindist, minind] = min(dist);
    
    j = j + 1;
    
end % end while loop (tolerance check)

%% Save consumption and capital trajectories on transition path

% save transition paths from chosen candidate
k_path = k_traj(minind, :);
c_path = c_traj(minind, :);

% add steady state to start of transition path
k_path = [k_ss, k_path(1:150)];
c_path = [c_ss, c_path(1:150)];

% save time vector
t_vec = 0:(length(k_path)-1);

%% Plot consumption and capital dynamics
figure(1)
    plot(t_vec,k_path,t_vec,k_ss*ones(151),'r--')
    text(t_vec(end)*.8,k_ss*1.0003,'Steady State','Color','red')
    xline(12,'--'); text(15,(max(k_path)+min(k_path))/2,'Event')
    title('Capital post-shock transition path')
    xlabel('t'); ylabel('K_t')
    saveas(gcf,'figure4a.png')
    
figure(2)
    plot(t_vec,c_path,t_vec,c_ss*ones(151),'r--')
    text(t_vec(end)*.8,c_ss*.9999,'Steady State','Color','red')
    xline(12,'--'); text(15,(max(c_path)+min(c_path))/2,'Event')
    title('Consumption post-shock transition path')
    xlabel('t'); ylabel('C_t')
    saveas(gcf,'figure4b.png')












