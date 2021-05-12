function [rmsError,rmsError_int1,rmsError_int2,time_elapsed,ESS,state_est,error] = ...
    pf(y, observation_mode,N,d,resamplingScheme,no_of_runs,state_true,A,sigma_u,sigma_v,ns,init_mean,init_std)
% int represents interval 
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
if observation_mode==1
    for t=190:210  %  sensor 1 fails
        y(1,t)=unifrnd(-pi,pi);
    end
    for t=250:260  %  sensor 2 fails
        y(2,t)=unifrnd(0,1e4);
    end
end
if observation_mode==2
    idx_start_1=190;  % start time for missing sensor 1 observation
    idx_end_1=200;   % end time for missing sensor 1 observation
    idx_start_2=250;  % start time for missing sensor 2 observation
    idx_end_2=260;   % end time for missing sensor 2 observation
end
for j=1:no_of_runs
    % INITIALISATION:
    % ==============
    xparticle_pf = repmat(init_mean,1,N)+repmat(init_std,1,N).*randn(d,N);        % These are the particles for the estimate
    xparticlePred_pf = ones(d,N);    % One-step-ahead predicted values of the states.
    yPred_pf = ones(2,N);            % One-step-ahead predicted values of y.
    wm=zeros(ns,N);
    state_est(:,1,j)=mean(xparticle_pf,2);
    tic;
    for t=2:ns,
        fprintf('run = %i / %i : PF : t = %i / %i  \r',j,no_of_runs,t,ns);
        fprintf('\n')
        
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
            loglik1= -1/2*(y(1,t)-yPred_pf(1,i))^2/(sigma_v(1)^2); %-log(sigma_v(1));
            loglik2= -1/2*(y(2,t)-yPred_pf(2,i))^2/(sigma_v(2)^2); %-log(sigma_v(2));
            if observation_mode==2
                if t>=idx_start_1 && t<=idx_end_1
                      wm(t,i)=exp(loglik2)+1e-323; %1e-99;
                elseif t>=idx_start_2 && t<=idx_end_2
                      wm(t,i)=exp(loglik1)+1e-323; %1e-99;
                else
                      wm(t,i)=exp(loglik1+loglik2)+1e-323; %1e-99;
                end
            else
                wm(t,i)=exp(loglik1+loglik2)+1e-323; %1e-99;
            end
        end
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