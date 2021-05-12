function [rmsError,rmsError_int1,rmsError_int2,time_elapsed,ESS,state_est,error] = ...
    pf_df(y, observation_mode,N,d,resamplingScheme,no_of_runs,state_true,A,sigma_u,sigma_v,ns,init_mean,init_std)
% Deterministic Fusion Method(df)
% Specify a model for each feature, combine results by averaging them (in another work, candidate models' weights are fixed and equal)
%**************************************************************************************
% SETUP BUFFERS TO STORE PERFORMANCE RESULTS
% ==========================================
rmsError = zeros(1,no_of_runs);
rmsError_int1 = zeros(1,no_of_runs);
rmsError_int2 = zeros(1,no_of_runs);
time_elapsed       = zeros(1,no_of_runs);
ESS           =zeros(no_of_runs,ns);
state_est    =ones(d,ns,no_of_runs);
error           =zeros(d,ns,no_of_runs);
%**************************************************************************************
M=2; % number of mixture components ; two features, namely bearing and range, of observations
if observation_mode==2
    idx_start_1=190;  % start time for missing sensor 1 observation
    idx_end_1=200;   % end time for missing sensor 1 observation
    idx_start_2=250;  % start time for missing sensor 2 observation
    idx_end_2=260;   % end time for missing sensor 2 observation
    for t=1:ns
        if t>=idx_start_1 && t<=idx_end_1
            y(1,t)=unifrnd(-pi,pi);
        end
        if t>=idx_start_2 && t<=idx_end_2
            y(2,t)=unifrnd(0,1e4);
        end
    end
end
for j=1:no_of_runs
    % INITIALISATION:
    % ==============
    xparticle_pf = repmat(init_mean,1,N)+repmat(init_std,1,N).*randn(d,N);        % These are the particles for the estimate
    xparticlePred_pf = ones(d,N);    % One-step-ahead predicted values of the states.
    yPred_pf = ones(2,N);            % One-step-ahead predicted values of y.
    w = ones(ns,N,M);                   % Importance weights.
    wm = ones(ns,N);                    % w after model averaging
    state_est(:,1,j)=mean(xparticle_pf,2);
    tic;
    for t=2:ns,
        fprintf('run = %i / %i :  PF-DF : t = %i / %i  \r',j,no_of_runs,t,ns);
        fprintf('\n')
        if observation_mode==2
            if t>=idx_start_1 && t<=idx_end_1
                y(1,t)=unifrnd(-pi,pi);
            end
            if t>=idx_start_2 && t<=idx_end_2
                y(2,t)=unifrnd(0,1e5);
            end
        end
        % PREDICTION STEP:
        % We use the transition prior as proposal.
        for i=1:N
            % xparticlePred_pf(:,i) = A*xparticle_pf(:,i)+B*[1.0*sigma_u(1)*randn;1.0*sigma_u(2)*randn];
            xparticlePred_pf(:,i) = A*xparticle_pf(:,i)+sigma_u'.*randn(d,1);
        end
        % EVALUATE IMPORTANCE WEIGHTS:
        for i=1:N
            yPred_pf(1,i) =solve_atan(xparticlePred_pf(3:4,i));
            yPred_pf(2,i) =norm(xparticlePred_pf(3:4,i));
        end
        for i=1:N
            loglik1= -1/2*(y(1,t)-yPred_pf(1,i))^2/(sigma_v(1)^2)-log(sigma_v(1)*sqrt(2*pi));
            loglik2= -1/2*(y(2,t)-yPred_pf(2,i))^2/(sigma_v(2)^2)-log(sigma_v(2)*sqrt(2*pi));
            w(t,i,1)=1/(sqrt(2*pi)*sigma_v(1))*exp(loglik1)+1e-323;%1e-99;  % 
            w(t,i,2)=1/(sqrt(2*pi)*sigma_v(2))*exp(loglik2)+1e-323;%1e-99;  % 
        end
        pai=[.5 .5];
        for m=1:M
            w(t,:,m) = w(t,:,m)/(sum(w(t,:,m)));  % Normalise the weights.
        end  
        wm(t,:) =  pai(1)*w(t,:,1)+pai(2)*w(t,:,2);  
        wm(t,:)=wm(t,:)/sum(wm(t,:));
        ESS(j,t)=1/sum(wm(t,:).^2)/N;
        % STATE ESTIMATED:
        state_est(:,t,j)=xparticlePred_pf*wm(t,:)';
        % RESAMPLING
        if resamplingScheme == 1
            outIndex = residualR(1:N,wm(t,:)');        % Residual resampling.
        elseif resamplingScheme == 2
            outIndex = systematicR(1:N,wm(t,:)');      % Systematic resampling.
        else
            outIndex = multinomialR(1:N,wm(t,:)');     % Multinomial resampling.
        end
        xparticle_pf = xparticlePred_pf(:,outIndex); % Keep particles with
    end
     %% Performance related statistics
    error(:,:,j)=state_est(:,:,j)-state_true;
    rmsError(j)  = sqrt(sum(sum(error(:,:,j).^2,1))/ns);
    rmsError_int1(j)= sqrt(sum(sum(error(:,1:100,j).^2,1))/100);  % first 100 time steps before algorithm convergence
    rmsError_int2(j)= sqrt(sum(sum(error(:,181:270,j).^2,1))/90); % interval with sensor failures
    time_elapsed(j)=toc;
end