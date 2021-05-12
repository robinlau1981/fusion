function [rmsError,rmsError_int1,rmsError_int2,time_elapsed,ESS,state_est,error,pai_ave] = ...
    dmmpf(y, observation_mode,N,d,resamplingScheme,no_of_runs,state_true,A,sigma_u,sigma_v,ns,init_mean,init_std)
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
M=4; % number of mixture components
% ff=1; % forgetting factor
pai_ave=zeros(ns,M); % averaged pai over Monte Carlo runs
thr_wm=1e-3;  % restrict the minimum model weight to be thr_wm
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
    pai=1/M*ones(1,M); %[.9 1/(M-1)*ones(1,M-1)];% weights of mixture components
    pai_ave(1,:)=pai_ave(1,:)+pai;
    mLik=ones(1,M); % marginal likelihood of mixture components
    
    % INITIALISATION:
    % ==============
    xparticle_pf =repmat(init_mean,1,N)+repmat(init_std,1,N).*randn(d,N);       % These are the particles for the estimate
    xparticlePred_pf = ones(d,N);    % One-step-ahead predicted values of the states.
    yPred_pf = ones(2,N);            % One-step-ahead predicted values of y.
    w = ones(ns,N,M);                   % Importance weights.
    wm = ones(ns,N);                    % w after model averaging
    state_est(:,1,j)=mean(xparticle_pf,2);
    tic;
    for t=2:ns,
        fprintf('run = %i / %i :  DMM PF : t = %i / %i  \r',j,no_of_runs,t,ns);
        fprintf('\n');
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
            loglik1= -1/2*(y(1,t)-yPred_pf(1,i))^2/(sigma_v(1)^2);%-log(sigma_v(1)*sqrt(2*pi));
            loglik2= -1/2*(y(2,t)-yPred_pf(2,i))^2/(sigma_v(2)^2);%-log(sigma_v(2)*sqrt(2*pi));
            w(t,i,1)=exp(loglik1+loglik2)/(2*pi*sigma_v(1)*sigma_v(2));%+1e-99;   % both sensors function well 
            w(t,i,2)=1/(sqrt(2*pi)*sigma_v(1))*exp(loglik1)*1e-4;%+1e-99;      % sensor 1 functions well & sensor 2 fails
            w(t,i,3)=1/(sqrt(2*pi)*sigma_v(2))*exp(loglik2)*1/(2*pi);%+1e-99;  % sensor 2 functions well & sensor 1 fails
            w(t,i,4)=1e-4*1/(2*pi);          % sensor 2 functions well & sensor 1 fails
        end
        for m=1:M
            mLik(m)=sum(w(t,:,m)); 
            if pai(m)<thr_wm
                pai(m)=thr_wm;
            end
        end
        pai=pai/sum(pai); 
        % use a forgetting mechanism to generate predictive distribution of
        % the candidate models 
        % pai_tmp=pai.^ff;
        % pai_pred=pai_tmp/sum(pai_tmp); 
        pai_pred=pai; % set the predictive distribution of the candidate models as the posterior at the previous time step 
        pai=pai_pred.*mLik/sum(pai_pred.*mLik);
        pai_ave(t,:)=pai_ave(t,:)+pai;
        
        for m=1:M
            w(t,:,m) = w(t,:,m)/(sum(w(t,:,m))+1e-99);  % Normalise the weights.
        end  
        wm(t,:) =  pai(1)*w(t,:,1)+pai(2)*w(t,:,2)+pai(3)*w(t,:,3)+pai(4)*w(t,:,4);  
        
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
pai_ave=pai_ave/no_of_runs;