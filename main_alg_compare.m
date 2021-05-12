clear;
close all;

load simu_data2;
y=z; % observations with Gaussian noise
no_of_runs = 1e2;     % number of experiments to generate statistical  averages
prob_fail=.8;
observation_mode=0; %  0: the normal case without Modality failure;
                    %  1: the case of Modality failures
                    %  2: the case of missing one Modality's observations
                    %  3: two modalitys fail simultaneously
if observation_mode==1
    for t=190:210  %  Modality 1 fails
        y(1,t)=unifrnd(-pi,pi);
    end
    
    for t=220:230  %  Modality 1 fails
        if rand<prob_fail
            y(1,t)=unifrnd(-pi,pi);
        end
    end
    
    for t=235:245  %  Modality 1 fails
        if rand<prob_fail
            y(2,t)=unifrnd(0,1e4);
        end
    end
    
    for t=250:260  %  Modality 2 fails
        y(2,t)=unifrnd(0,1e4);
    end
end
if observation_mode==2
    idx_start_1=190;  % start time for missing Modality 1 observation
    idx_end_1=200;    % end time for missing Modality 1 observation
    idx_start_2=250;  % start time for missing Modality 2 observation
    idx_end_2=260;    % end time for missing Modality 2 observation
end

if observation_mode==3
    for t=190:200  %  Modality 1 and 2 fail
        y(1,t)=unifrnd(-pi,pi);
        y(2,t)=unifrnd(0,1e4);
    end
    
    for t=210:240
        if rand<prob_fail
           y(1,t)=unifrnd(-pi,pi);
           y(2,t)=unifrnd(0,1e4); 
        end
    end
    
    for t=250:260  %  Modality 1 and 2 fail
        y(1,t)=unifrnd(-pi,pi);
        y(2,t)=unifrnd(0,1e4);
    end
end
N = 1e4;                    % Number of particles.
d=size(state_true,1);  % dimension of state vector
resamplingScheme = 1;

init_mean=state_true(:,1);
init_std=[1 1 10 10]';

[rmsError_pf,rmsError_int1_pf,rmsError_int2_pf,time_elapsed_pf,ESS_pf,state_est_pf,error_pf] = ...
    pf(y,observation_mode,N,d,resamplingScheme,no_of_runs,state_true,A,sigma_u,sigma_v,ns,init_mean,init_std);
[rmsError_pf_alpha,rmsError_int1_pf_alpha,rmsError_int2_pf_alpha,time_elapsed_pf_alpha,ESS_pf_alpha,state_est_pf_alpha,error_pf_alpha,alpha_pf_alpha] = ...
    pf_alpha(y,observation_mode,N,d,resamplingScheme,no_of_runs,state_true,A,sigma_u,sigma_v,ns,init_mean,init_std);
[rmsError_pf_df,rmsError_int1_pf_df,rmsError_int2_pf_df,time_elapsed_pf_df,ESS_pf_df,state_est_pf_df,error_pf_df] = ...
    pf_df(y,observation_mode,N,d,resamplingScheme,no_of_runs,state_true,A,sigma_u,sigma_v,ns,init_mean,init_std);
[rmsError_dmmpf,rmsError_int1_dmmpf,rmsError_int2_dmmpf,time_elapsed_dmmpf,ESS_dmmpf,state_est_dmmpf,error_dmmpf,pai_ave_dmmpf] = ...
    dmmpf(y,observation_mode,N,d,resamplingScheme,no_of_runs,state_true,A,sigma_u,sigma_v,ns,init_mean,init_std);
% calculate mean of RMSE errors
mean_RMSE_pf      = mean(rmsError_pf);
mean_RMSE_int1_pf      = mean(rmsError_int1_pf);
mean_RMSE_int2_pf      = mean(rmsError_int2_pf);
% calculate variance of RMSE errors
var_RMSE_pf      = var(rmsError_pf);
var_RMSE_int1_pf      = var(rmsError_int1_pf);
var_RMSE_int2_pf      = var(rmsError_int2_pf);
% calculate mean of execution time
mean_time_pf      = mean(time_elapsed_pf);
% display final results
disp('************* FINAL RESULTS *****************');
disp('RMSE of PF : mean and variance');
disp('---------');
disp(['PF           = ' num2str(mean_RMSE_pf) ' (' num2str(var_RMSE_pf) ')']);
disp(' ');
disp('RMSE of PF (interval 1 & 2) : mean and variance');
disp('---------');
disp(['PF           = ' num2str(mean_RMSE_int1_pf) ' (' num2str(var_RMSE_int1_pf) ')']);
disp('---------');
disp(['PF           = ' num2str(mean_RMSE_int2_pf) ' (' num2str(var_RMSE_int2_pf) ')']);
disp('Execution time of PF (seconds)');
disp('-------------------------');
disp(['PF           = ' num2str(mean_time_pf)]);


% calculate mean of RMSE errors
mean_RMSE_pf_alpha      = mean(rmsError_pf_alpha);
mean_RMSE_int1_pf_alpha      = mean(rmsError_int1_pf_alpha);
mean_RMSE_int2_pf_alpha      = mean(rmsError_int2_pf_alpha);
% calculate variance of RMSE errors
var_RMSE_pf_alpha      = var(rmsError_pf_alpha);
var_RMSE_int1_pf_alpha      = var(rmsError_int1_pf_alpha);
var_RMSE_int2_pf_alpha      = var(rmsError_int2_pf_alpha);
% calculate mean of execution time
mean_time_pf_alpha      = mean(time_elapsed_pf_alpha);
% display final results

disp('************* FINAL RESULTS *****************');
disp('RMSE of PF-alpha : mean and variance');
disp('---------');
disp(['PF           = ' num2str(mean_RMSE_pf_alpha) ' (' num2str(var_RMSE_pf_alpha) ')']);
disp('RMSE of PF-alpha (interval 1 & 2) : mean and variance');
disp('---------');
disp(['PF           = ' num2str(mean_RMSE_int1_pf_alpha) ' (' num2str(var_RMSE_int1_pf_alpha) ')']);
disp('---------');
disp(['PF           = ' num2str(mean_RMSE_int2_pf_alpha) ' (' num2str(var_RMSE_int2_pf_alpha) ')']);
disp('Execution time of PF-alpha (seconds)');
disp('-------------------------');
disp(['PF           = ' num2str(mean_time_pf_alpha)]);

% calculate mean of RMSE errors
mean_RMSE_pf_df      = mean(rmsError_pf_df);
mean_RMSE_int1_pf_df      = mean(rmsError_int1_pf_df);
mean_RMSE_int2_pf_df      = mean(rmsError_int2_pf_df);
% calculate variance of RMSE errors
var_RMSE_pf_df      = var(rmsError_pf_df);
var_RMSE_int1_pf_df      = var(rmsError_int1_pf_df);
var_RMSE_int2_pf_df      = var(rmsError_int2_pf_df);
% calculate mean of execution time
mean_time_pf_df      = mean(time_elapsed_pf_df);
% display final results
disp('************* FINAL RESULTS *****************');
disp('RMSE of PF-DF : mean and variance');
disp('---------');
disp(['PF-DF           = ' num2str(mean_RMSE_pf_df) ' (' num2str(var_RMSE_pf_df) ')']);
disp('RMSE of PF-DF (Interval 1 & 2): mean and variance');
disp('---------');
disp(['PF-DF           = ' num2str(mean_RMSE_int1_pf_df) ' (' num2str(var_RMSE_int1_pf_df) ')']);
disp(['PF-DF           = ' num2str(mean_RMSE_int2_pf_df) ' (' num2str(var_RMSE_int2_pf_df) ')']);
disp('Execution time of PF-DF (seconds)');
disp('-------------------------');
disp(['PF-DF           = ' num2str(mean_time_pf_df)]);

% calculate mean of RMSE errors
mean_RMSE_dmmpf      = mean(rmsError_dmmpf);
mean_RMSE_int1_dmmpf      = mean(rmsError_int1_dmmpf);
mean_RMSE_int2_dmmpf      = mean(rmsError_int2_dmmpf);
% calculate variance of RMSE errors
var_RMSE_dmmpf      = var(rmsError_dmmpf);
var_RMSE_int1_dmmpf      = var(rmsError_int1_dmmpf);
var_RMSE_int2_dmmpf      = var(rmsError_int2_dmmpf);
% calculate mean of execution time
mean_time_dmmpf      = mean(time_elapsed_dmmpf);
% display final results
disp('************* FINAL RESULTS *****************');
disp('RMSE of DMMPF : mean and variance');
disp('---------');
disp(['DMMPF           = ' num2str(mean_RMSE_dmmpf) ' (' num2str(var_RMSE_dmmpf) ')']);
disp('RMSE of DMMPF (Interval 1 & 2): mean and variance');
disp('---------');
disp(['DMMPF           = ' num2str(mean_RMSE_int1_dmmpf) ' (' num2str(var_RMSE_int1_dmmpf) ')']);
disp(['DMMPF           = ' num2str(mean_RMSE_int2_dmmpf) ' (' num2str(var_RMSE_int2_dmmpf) ')']);
disp('Execution time of DMMPF (seconds)');
disp('-------------------------');
disp(['DMMPF           = ' num2str(mean_time_dmmpf)]);

% figure,
% plot(2:ns, mean(error_pf(1,2:ns,:),3),'LineWidth',1,'Color',[.4 .4 .4]);
% xlabel('Time');ylabel('V_x error'); 
% grid on;
% hold on;
% plot(2:ns, mean(error_pf_alpha(1,2:ns,:),3),'LineWidth',1,'Color',[1 0 0]);
% hold on;
% plot(2:ns, mean(error_pf_df(1,2:ns,:),3),'LineWidth',1,'Color',[0 1 1]);
% hold on;
% plot(2:ns, mean(error_idmmpf(1,2:ns,:),3),'LineWidth',1,'Color',[0 0 .1]);
% hold on;
% plot(2:ns, mean(error_dmmpf(1,2:ns,:),3),'LineWidth',1,'Color',[0 1 0]);
% legend('PF','PF-alpha','PF-df','PF-iDMM','PF-DMM');

% figure,
% plot(2:ns, mean(error_pf(2,2:ns,:),3),'LineWidth',1,'Color',[.4 .4 .4]);
% xlabel('Time');ylabel('V_y error'); 
% grid on;
% hold on;
% plot(2:ns, mean(error_pf_alpha(2,2:ns,:),3),'LineWidth',1,'Color',[1 0 0]);
% hold on;
% plot(2:ns, mean(error_pf_df(2,2:ns,:),3),'LineWidth',1,'Color',[0 1 1]);
% hold on;
% plot(2:ns, mean(error_idmmpf(2,2:ns,:),3),'LineWidth',1,'Color',[0 0 1]);
% hold on;
% plot(2:ns, mean(error_dmmpf(2,2:ns,:),3),'LineWidth',1,'Color',[0 1 0]);
% legend('PF','PF-alpha','PF-df','PF-iDMM','PF-DMM');

% figure,
% plot(mean(state_est_pf(3,:,:),3), mean(state_est_pf(4,:,:),3),'LineWidth',1,'Color',[.4 .4 .4]);
% xlabel('X');ylabel('Y');
% grid on;
% hold on;
% plot(mean(state_est_pf_alpha(3,:,:),3), mean(state_est_pf_alpha(4,:,:),3),'LineWidth',1,'Color',[1 0 0]);
% hold on;
% plot(mean(state_est_pf_df(3,:,:),3), mean(state_est_pf_df(4,:,:),3),'LineWidth',1,'Color',[0 1 1]);
% hold on;
% plot(mean(state_est_idmmpf(3,:,:),3), mean(state_est_idmmpf(4,:,:),3),'LineWidth',1,'Color',[0 0 1]);
% hold on;
% plot(mean(state_est_dmmpf(3,:,:),3), mean(state_est_dmmpf(4,:,:),3),'LineWidth',1,'Color',[0 1 0]);
% hold on;
% plot(state_true(3,:), state_true(4,:),'k-');
% legend('PF','PF-alpha','PF-df','PF-iDMM','PF-DMM','true answer');
% 
% figure,
% plot(2:ns, mean(error_pf(3,2:ns,:),3),'LineWidth',1,'Color',[.4 .4 .4]);
% xlabel('Time');ylabel('R_x error'); 
% grid on;
% hold on;
% plot(2:ns, mean(error_pf_alpha(3,2:ns,:),3),'LineWidth',1,'Color',[1 0 0]);
% hold on;
% plot(2:ns, mean(error_pf_df(3,2:ns,:),3),'LineWidth',1,'Color',[0 1 1]);
% hold on;
% plot(2:ns, mean(error_idmmpf(3,2:ns,:),3),'LineWidth',1,'Color',[0 0 1]);
% hold on;
% plot(2:ns, mean(error_dmmpf(3,2:ns,:),3),'LineWidth',1,'Color',[0 1 0]);
% legend('PF','PF-alpha','PF-df','PF-iDMM','PF-DMM');
% 
% figure,
% plot(2:ns, mean(error_pf(4,2:ns,:),3),'LineWidth',1,'Color',[.4 .4 .4]);
% xlabel('Time');ylabel('R_y error'); 
% grid on;
% hold on;
% plot(2:ns, mean(error_pf_alpha(4,2:ns,:),3),'LineWidth',1,'Color',[1 0 0]);
% hold on;
% plot(2:ns, mean(error_pf_df(4,2:ns,:),3),'LineWidth',1,'Color',[0 1 1]);
% hold on;
% plot(2:ns, mean(error_idmmpf(4,2:ns,:),3),'LineWidth',1,'Color',[0 0 1]);
% hold on;
% plot(2:ns, mean(error_dmmpf(4,2:ns,:),3),'LineWidth',1,'Color',[0 1 0]);
% legend('PF','PF-alpha','PF-df','PF-iDMM','PF-DMM');
% 
% figure,
% plot(1:ns,pai_ave_idmmpf(:,1)',1:ns,pai_ave_idmmpf(:,2)');
% grid on;
% legend('Model 1','Model 2');
% xlabel('time');
% ylabel('Model Weight');
% title('Averaged model weights of PF-iDMM');
% 
% figure,
% plot(1:ns,mean(alpha_pf_alpha(1,:,:),3),1:ns,mean(alpha_pf_alpha(2,:,:),3));
% grid on;
% legend('Modality 1','Modality 2');
% xlabel('time');
% ylabel('Probability of Modality');
% title('Trustworthiness of Senor of PF-alpha');

figure,
plot(1:ns,pai_ave_dmmpf(:,1)',1:ns,pai_ave_dmmpf(:,2)',1:ns,pai_ave_dmmpf(:,3)',1:ns,pai_ave_dmmpf(:,4)');
grid on;
legend('Model 1','Model 2','Model 3','Model 4');
xlabel('time');
ylabel('Model Weight');
title('Averaged model weights of PF-DMM');

if observation_mode==0
    save res_alg_compare_observe_mode_0_2nd_setting ...
        observation_mode state_true error_pf state_est_pf  ...
        rmsError_pf time_elapsed_pf mean_RMSE_pf mean_RMSE_int1_pf mean_RMSE_int2_pf var_RMSE_pf var_RMSE_int1_pf var_RMSE_int2_pf...
        error_pf_alpha state_est_pf_alpha  ...
        rmsError_pf_alpha time_elapsed_pf_alpha mean_RMSE_pf_alpha mean_RMSE_int1_pf_alpha mean_RMSE_int2_pf_alpha var_RMSE_pf_alpha var_RMSE_int1_pf_alpha var_RMSE_int2_pf_alpha alpha_pf_alpha...
        error_dmmpf state_est_dmmpf rmsError_dmmpf time_elapsed_dmmpf ...
        mean_RMSE_dmmpf mean_RMSE_int1_dmmpf mean_RMSE_int2_dmmpf var_RMSE_dmmpf var_RMSE_int1_dmmpf var_RMSE_int2_dmmpf pai_ave_dmmpf...
        error_pf_df state_est_pf_df rmsError_pf_df time_elapsed_pf_df ...
        mean_RMSE_pf_df mean_RMSE_int1_pf_df mean_RMSE_int2_pf_df var_RMSE_pf_df var_RMSE_int1_pf_df var_RMSE_int2_pf_df;
end
if observation_mode==1
    save res_alg_compare_observe_mode_1_2nd_setting ...
        observation_mode state_true error_pf state_est_pf  ...
        rmsError_pf time_elapsed_pf mean_RMSE_pf mean_RMSE_int1_pf mean_RMSE_int2_pf var_RMSE_pf var_RMSE_int1_pf var_RMSE_int2_pf...
        error_pf_alpha state_est_pf_alpha  ...
        rmsError_pf_alpha time_elapsed_pf_alpha mean_RMSE_pf_alpha mean_RMSE_int1_pf_alpha mean_RMSE_int2_pf_alpha var_RMSE_pf_alpha var_RMSE_int1_pf_alpha var_RMSE_int2_pf_alpha alpha_pf_alpha...
        error_dmmpf state_est_dmmpf rmsError_dmmpf time_elapsed_dmmpf ...
        mean_RMSE_dmmpf mean_RMSE_int1_dmmpf mean_RMSE_int2_dmmpf var_RMSE_dmmpf var_RMSE_int1_dmmpf var_RMSE_int2_dmmpf pai_ave_dmmpf ...
        error_pf_df state_est_pf_df rmsError_pf_df time_elapsed_pf_df ...
        mean_RMSE_pf_df mean_RMSE_int1_pf_df mean_RMSE_int2_pf_df var_RMSE_pf_df var_RMSE_int1_pf_df var_RMSE_int2_pf_df;
end
if observation_mode==2
    save res_alg_compare_observe_mode_2_2nd_setting ...
        observation_mode state_true error_pf state_est_pf  ...
        rmsError_pf time_elapsed_pf mean_RMSE_pf mean_RMSE_int1_pf mean_RMSE_int2_pf var_RMSE_pf var_RMSE_int1_pf var_RMSE_int2_pf...
        error_pf_alpha state_est_pf_alpha  ...
        rmsError_pf_alpha time_elapsed_pf_alpha mean_RMSE_pf_alpha mean_RMSE_int1_pf_alpha mean_RMSE_int2_pf_alpha var_RMSE_pf_alpha var_RMSE_int1_pf_alpha var_RMSE_int2_pf_alpha alpha_pf_alpha...
        error_dmmpf state_est_dmmpf rmsError_dmmpf time_elapsed_dmmpf ...
        mean_RMSE_dmmpf mean_RMSE_int1_dmmpf mean_RMSE_int2_dmmpf var_RMSE_dmmpf var_RMSE_int1_dmmpf var_RMSE_int2_dmmpf pai_ave_dmmpf ...
        error_pf_df state_est_pf_df rmsError_pf_df time_elapsed_pf_df ...
        mean_RMSE_pf_df mean_RMSE_int1_pf_df mean_RMSE_int2_pf_df var_RMSE_pf_df var_RMSE_int1_pf_df var_RMSE_int2_pf_df;
end
if observation_mode==3
    save res_alg_compare_observe_mode_3_2nd_setting ...
        observation_mode state_true error_pf state_est_pf  ...
        rmsError_pf time_elapsed_pf mean_RMSE_pf mean_RMSE_int1_pf mean_RMSE_int2_pf var_RMSE_pf var_RMSE_int1_pf var_RMSE_int2_pf...
        error_pf_alpha state_est_pf_alpha  ...
        rmsError_pf_alpha time_elapsed_pf_alpha mean_RMSE_pf_alpha mean_RMSE_int1_pf_alpha mean_RMSE_int2_pf_alpha var_RMSE_pf_alpha var_RMSE_int1_pf_alpha var_RMSE_int2_pf_alpha alpha_pf_alpha...
        error_dmmpf state_est_dmmpf rmsError_dmmpf time_elapsed_dmmpf ...
        mean_RMSE_dmmpf mean_RMSE_int1_dmmpf mean_RMSE_int2_dmmpf var_RMSE_dmmpf var_RMSE_int1_dmmpf var_RMSE_int2_dmmpf pai_ave_dmmpf ...
        error_pf_df state_est_pf_df rmsError_pf_df time_elapsed_pf_df ...
        mean_RMSE_pf_df mean_RMSE_int1_pf_df mean_RMSE_int2_pf_df var_RMSE_pf_df var_RMSE_int1_pf_df var_RMSE_int2_pf_df;
end