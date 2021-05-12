
clear;
close all;

ts=1;                  % sampling period of observations, unit: second
tmax=5*60;       %  tracking lasts 15 minutes = 15*60 seconds 
ns=fix(tmax/ts);   % number of observations 

%%                     % initialization of stats of the target and the observer 
aox(1)=0;            % X-coordinate acceleration of observer 
aoy(1)=0;            % Y-coordinate acceleration of observer 
rox(1)=0;             % X-coordinate position of observer 
roy(1)=0;             % Y-coordinate position of observer 
vox(1)=10*sqrt(2); % X-coordinate velocity of observer 
voy(1)=1*sqrt(2); % Y-coordinate velocity of observer 
rtx(1)=1e3;         % X-coordinate position of target
rty(1)=5e3;              % Y-coordinate position of target
vtx(1)=15;            % X-coordinate velocity of target    %须保证 vtx0>vox(1),否则会出现角度从 pi到-pi突变的情况
vty(1)=-1;              % Y-coordinate velocity of target
atx(1)=0;              % X-coordinate acceleration of target 
aty(1)=0;              % Y-coordinate acceleration of target 
% sigma_u=[1e-2 1e-2];    % state noise std
sigma_u=[1 1 10 10];    % state noise std
sigma_v=[1e-2 1e1];    % standard error of observations  
A=[1 0 0 0;0 1 0 0;ts 0 1 0;0 ts 0 1];
% B=[ts 0;0 ts;.5*ts^2 0;0 .5*ts^2];
y(:,1)=[vtx(1)-vox(1) vty(1)-voy(1) rtx(1)-rox(1) rty(1)-roy(1)]';
for k=1:ns-1
%     aox(k+1)=0;
%     aoy(k+1)=0;
%     atx(k+1)=sigma_u(1)*randn;
%     aty(k+1)=sigma_u(2)*randn;
%    u=[atx(k+1)-aox(k+1) aty(k+1)-aoy(k+1)]';  %目标相对观测站的加速度
%    y(:,k+1)=A*y(:,k)+B*u;
     y(:,k+1)=A*y(:,k)+sigma_u'.*randn(4,1);
end
state_true=y;

z=zeros(2,ns); %存放理想测量值（无噪声）
for k=1:ns  
    z(1,k)=solve_atan(state_true(3:4,k)) + sigma_v(1)*randn;     
    z(2,k)=norm(state_true(3:4,k))  + sigma_v(2)*randn;    
end       

figure,
plot(ts*(1:ns),z(1,1:ns)/pi*180); grid on; title('Bearing observation');
figure,
plot(ts*(1:ns),z(2,1:ns)); grid on; title('Range observation');
figure,
plot(ts*(1:ns),state_true(1,1:ns)); grid on; title('V_x');
figure,
plot(ts*(1:ns),state_true(2,1:ns)); grid on; title('V_y');
figure,
plot(ts*(1:ns),state_true(3,1:ns)); grid on; title('R_x');
figure,
plot(ts*(1:ns),state_true(4,1:ns)); grid on; title('R_y');

%save simu_data A B ts ns state_true z sigma_u sigma_v;
save simu_data A ts ns state_true z sigma_u sigma_v;
