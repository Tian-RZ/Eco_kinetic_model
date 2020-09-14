function Eco_GlcNAc
% simulation of growth curve & GlcAc & HAc curve of Eco
close all
global par
def_pacf=[0.1:0.1:10];

% Random Sampling Loop
       for i = 1:100
          par.pacf=def_pacf(i);
          pacf(i)=par.pacf;
        % Km Values    
          par.Km = ones(1,6); %initial km
      
        % Vmax Values per OD(a.u.) per hour(a.u.) are fixed to 1
          par.vmax = 0.24*ones(1,6);

        % Influx per OD(a.u.) per hour(a.u.)
          par.vmax(2) = 1;
          
        % initial other parameters
          par.miu0 = 0.76;  % miu value without limitation
          par.n0 = 50;  % Toral nitrogen source
          par.Y_x = 0.7;  % 0.7 Unit cell per Unit nitrogen
          par.Y_GNA = 0.54;  %  0.54 Unit GlcNAc per Unit nitrogen
        
        % initial reaction parameters
          x0(1,9) = 0.1;  %initial biomass
          x0(1,7) = 0;
          [time,X] = ode23s(@pathw_fun,[0 10],x0);
          time_hat = [0:0.1:10];
          for k = 1:9
          X_hat(:,k) = interp1(time,X(:,k),time_hat);
    end

X_cella(i,:) = X_hat(:,9)';
X_Na(i,:) = X_hat(:,8)';
X_Acea(i,:) = X_hat(:,7)';
X_GlcNAca(i,:) = X_hat(:,6)';

% X_cell = mean(X_cella); X_cell_std = std(X_cella);
% X_N = mean(X_Na); X_N_std = std(X_Na);
% X_Ace = mean(X_Acea); X_HAc_std = std(X_Acea);
% X_GlcNAc = mean(X_GlcNAca); X_GlcNAc_std = std(X_GlcNAca);
 
X_cellb(i) = X_cella(end,end);
X_Nb(i) = X_Na(end,end);
X_Aceb(i) = X_Acea(end,end);
X_GlcNAcb(i) = X_GlcNAca(end,end);

Sim(i,1)=par.pacf;
 
figure(1)

subplot(2,2,1)
plot(time_hat,X_cella,'r-','LineWidth',2)
set(gca, 'FontSize', 14);
title('Cell growth','fontsize',12);
hold on

 end % of sampling loop
 Sim(:,2)= X_cellb';
 Sim(:,3)= X_Nb';
 Sim(:,4)= X_Aceb';
 Sim(:,5)= X_GlcNAcb';
 
subplot(2,2,2)
plot( pacf,X_cellb,'b-','LineWidth',2)
set(gca, 'xscale', 'log');
title('Maximum Biomass','fontsize',12);

subplot(2,2,3)
plot( pacf,X_Aceb,'g-','LineWidth',2)
set(gca, 'xscale', 'log');
title('Acetate','fontsize',12);

subplot(2,2,4)
plot( pacf,X_GlcNAcb,'r-','LineWidth',2)
set(gca, 'xscale', 'log');
title('GlcNAc','fontsize',12); 
     
end


%linear pathway ODEs - calculate rates and mass balances
function dx_dt = pathw_fun(~,x)

global par

% Glycolytic pathway parameters
I_glc = 0.152*log(par.pacf)+0.65;
if I_glc <0.6
    I_atp = 7/3*I_glc-0.4;
else
    I_atp = 1.5*I_glc+0.1;  
end

% Acetate biosynthesis
  % k2 = 10^(-0.06*x(7)-0.1192);
 k2=0.1+0.76*0.0534^x(7);
 I_ace = k2/par.miu0;
  % Ace biosynthesis rate/OD/hour
if I_glc <= 0.6
    v_ace = 0;
else
    v_ace = I_ace*(0.6*I_glc-0.36); 
end

dx_dt(7,1) = v_ace*x(9);

% GlcNAc biosynthesis
v_2 = par.vmax(2); % Glucose absorption rate/OD/hour
v_3 = par.vmax(3)*x(2)/(1*x(2)+par.Km(3)); % Reaction v_3
v_4 = par.vmax(4)*x(3)/(2*(x(3)+par.Km(4))); % Reaction v_4
v_5 = par.vmax(5)*x(4)/(2*(x(4)+par.Km(5))); % Reaction v_5 
v_6 = I_ace*par.vmax(6)*x(5)/(2*(x(5)+par.Km(6))); % Reaction v_6
v_7 = v_6+0.05; %Basic reaction rate/OD/hour

dx_dt(2,1) = x(9)*(v_2-v_3);
dx_dt(3,1) = x(9)*(v_3-v_4);
dx_dt(4,1) = x(9)*(v_4-v_5);
dx_dt(5,1) = x(9)*(v_5-v_7);
dx_dt(6,1) = x(9)*v_7;

% Nitrogen source
x(8) = par.n0-x(9)/par.Y_x-x(6)/par.Y_GNA;
% I_N = x(8)/par.n0;
I_N = x(8)/(x(8)+15); % KS=(5/10/15/20)

%Cell growth
miu = I_ace*I_atp*I_N*par.miu0;
% miu1 = (x(8)*par.miu0)/(x(8)+1);
dx_dt(9,1) = x(9)*miu;

if x(8) <= 0
    dx_dt(9,1) =0;
    
end

end
