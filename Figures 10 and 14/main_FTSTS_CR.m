clear -all
clearvars

% start timer
    tic
 
% -- Run Parameters --
    duration = 1000*5 % ms
    step = 0.1; %ms
    duration_step = duration/step;
    
% -- Neuron Parameters --
    mew_e = 20.8; 
    sigma_e = 1;
    mew_i = 18; 
    sigma_i = 3;
    mew_c = 0;
    N_E = 1600;
    N_I = 400;
    Ntot = N_E + N_I;
    Energy = zeros(duration_step,1);
    
% -- Synaptic Parameters --
    Weight_0 = 1;
    J_E = 100;
    J_EE = 50;    
    J_II = 50;
    N_i = 1;
    C_E = 0.3*Ntot; 
    C_EE = Ntot;
    C_I = 0.3*Ntot;
    C_II = Ntot;
    
    % TO MAKE FIGURE 10 UNCOMMENT
    tau_LTP = 20; 
    tau_LTD = 22; 
    J_I = 260
    
    % TO MAKE FIGURE 14 UNCOMMENT
%     tau_LTP = 20; 
%     tau_LTD = 20;  
%     J_I = 175
    
    
% -- Make Random Synaptic Conncetions ---
    epsilon_E = 0.1; 
    epsilon_I = 0.1;
    epsilon_EE = 0.1;
    epsilon_II = 0.1;
    S_key_IE = zeros(N_E,N_I);
    S_key_EI = zeros(N_I,N_E);
    S_key_EE = zeros(N_E,N_E);
    S_key_II = zeros(N_I,N_I);
    
    %--- I to E ---
    syn_count = 0;
    for pre_neuron = 1:N_I
        for post_neuron = 1:N_E
            x = rand;
            if x <= epsilon_I
                syn_count = syn_count + 1;
                S_key_EI(pre_neuron,post_neuron) = syn_count;
            else
                S_key_EI(pre_neuron,post_neuron) = 0;
            end
    
        end
    end
    
    % --- E to I ---
    syn_count = 0;
    for pre_neuron = 1:N_E
        for post_neuron = 1:N_I
            x = rand;
            if x <= epsilon_E
                syn_count = syn_count + 1;
                S_key_IE(pre_neuron,post_neuron) = syn_count;
            else
                S_key_IE(pre_neuron,post_neuron) = 0;
            end
    
        end
    end

    % --- E to E ---
    syn_count = 0;
    for pre_neuron = 1:N_E
        for post_neuron = 1:N_E
            x = rand;
            if x <= epsilon_EE && pre_neuron ~= post_neuron
                syn_count = syn_count + 1;
                S_key_EE(pre_neuron,post_neuron) = syn_count;
            else
                S_key_EE(pre_neuron,post_neuron) = 0;
            end
    
        end
    end
    
    % --- I to I ---
    syn_count = 0;
    for pre_neuron = 1:N_I
        for post_neuron = 1:N_I
            x = rand;
            if x <= epsilon_II && post_neuron ~= pre_neuron
                syn_count = syn_count + 1;
                S_key_II(pre_neuron,post_neuron) = syn_count;
            else
                S_key_II(pre_neuron,post_neuron) = 0;
            end
    
        end
    end
    
    num_synapses_IE = max(max(S_key_IE));
    num_synapses_EI = max(max(S_key_EI));
    W_IE = zeros(duration_step,1); 
    W_IE_std = zeros(duration_step,1); 
    
% -- Initial Conditions --
    vE0 = 14*ones(1,N_E);
    vI0 = 14*ones(1,N_I);
    S_EI0 = zeros(1,N_E);
    S_IE0 = zeros(1,N_I);
    S_EE0 = zeros(1,N_E);
    S_II0 = zeros(1,N_I);
    X_IE0 = zeros(1,N_I);
    X_EI0 = zeros(1,N_E);
    X_II0 = zeros(1,N_I);
    X_EE0 = zeros(1,N_E);
    spt_E0 = 0;
    phif = zeros(1,N_E);
    
    Apost0 = zeros(1,num_synapses_IE);
    Apre0 = zeros(1,num_synapses_IE);
    W_IE0 = Weight_0*ones(1,num_synapses_IE);   
    W_EI0 = Weight_0; 
    leftover_S_EI = zeros(5/step + 1 ,N_E);
    leftover_S_IE = zeros(5/step + 1 ,N_I);
    leftover_S_EE = zeros(5/step + 1 ,N_E);
    leftover_S_II = zeros(5/step + 1 ,N_I);
    ref_E = zeros(1,N_E);
    ref_I = zeros(1,N_I);
    tau_E_m = 10;
    tau_I_m = 10;
    
% -- Run --    
    sample_duration = 20;
    Synchrony = zeros(duration/sample_duration,1);
    time_syn = zeros(duration/sample_duration,1);

% Generate General Stimulation Pattern
% TO MAKE FIGURE 10 UNCOMMENT
V_stim = 1;   
T_stim = 1;     
x_neutral = 7; 
multiple = 1;   
t_pulse = T_stim*(x_neutral+multiple+1);
interval_time = 4*4*(x_neutral*T_stim+T_stim+multiple*T_stim);

% TO MAKE FIGURE 14 UNCOMMENT
% V_stim = 1;   
% T_stim = 0.5;     
% x_neutral = 3; 
% multiple = 1;   
% t_pulse = T_stim*(x_neutral+multiple+1);
% interval_time = 4*4*(x_neutral*T_stim+T_stim+multiple*T_stim);

[Ue Ui] = pulsatile_input_25(multiple,V_stim,T_stim,x_neutral,duration,step,interval_time);
[Ue2 Ui2] = pulsatile_input_25(multiple,V_stim,T_stim,x_neutral,duration,step,interval_time);
sum(Ue,2)
sum(Ui,2)
cross_100 = 1;
Energy0 = 0;

for i = 1:duration/sample_duration
    % run parameters
    comp_time = (i-1)*sample_duration;
    if mean(W_IE0)*J_I < 75
        cross_100 = 0;
    end
    ON = 1*(i*sample_duration >= 2000)*cross_100; % turn on control input
    plast_on = 1*(i*sample_duration >= 100); % past on/off
       
    % indexes
    a = 1 + (i>=2)*(i-1)*sample_duration/step;
    b = i*sample_duration/step;
    
    
    % TO MAKE FIGURE 10 W/ FTSTS-CR UNCOMMENT
    Vstim = 100; 
    ue = Vstim*Ue(:,a:b);
    ui = Vstim*Ui(:,a:b);
    percent_V_stim = 1;
    
    % TO MAKE FIGURE 10 W/ CR UNCOMMENT
%     Vstim = 100; 
%     ue = Vstim*Ue(:,a:b);
%     ui = Vstim*Ue2(:,a:b);
%     percent_V_stim = 1;

    % TO MAKE FIGURE 14 W/ FTSTS-CR UNCOMMENT
%     Vstim = 200; 
%     ue = Vstim*Ue(:,a:b);
%     ui = Vstim*Ui(:,a:b);
%     percent_V_stim = 1;
    
    % TO MAKE FIGURE 14 W/ CR UNCOMMENT
%     Vstim = 200; 
%     ue = Vstim*Ue(:,a:b);
%     ui = Vstim*Ue2(:,a:b);
%     percent_V_stim = 1;
    
    
   [timem v_Em v_Im S_EIm S_IEm S_EEm S_IIm X_EIm X_IEm X_EEm X_IIm Apostm Aprem W_IEm spike_Em spike_Im ref_Em ref_Im synchronym spt_Em phif] = ode_neuron_model(plast_on,ON,vE0,vI0,S_EI0,S_IE0,S_II0,S_EE0,X_EI0,X_IE0,X_EE0,X_II0,Apost0,Apre0,W_IE0,W_EI0,mew_e,sigma_e,ue,ui,mew_i,sigma_i,J_E,J_EE,J_I,J_II,C_E,C_I,C_EE,C_II,tau_LTP,tau_LTD,step,sample_duration,N_E,N_I,S_key_EI,S_key_IE,S_key_EE,S_key_II,leftover_S_EI,leftover_S_IE,leftover_S_EE,leftover_S_II,ref_E,ref_I,tau_E_m,tau_I_m, percent_V_stim,comp_time,spt_E0,phif);
    
    
    
    % recorded variables
    time(a:b,:) = timem(1:sample_duration/step,:) + (i-1)*sample_duration;
    W_IE(a:b,1) = mean(W_IEm(1:sample_duration/step,:),2);
%     spike_E(a:b,:) = spike_Em(1:sample_duration/step,:);
%     spike_I(a:b,:) = spike_Im(1:sample_duration/step,:);
    Synchrony(i) = synchronym;
    time_syn(i) = sample_duration*(i);
    spike_time_E(a:b,:) = spt_Em(1:sample_duration/step,:);
    
    
    % generate intial condition for next run
    sample_end = sample_duration/step ;
    vE0 = v_Em(sample_end,:);
    vI0 = v_Im(sample_end,:);
    S_EI0 = S_EIm(sample_end,:);
    S_IE0 = S_IEm(sample_end,:);
    S_II0 = S_IIm(sample_end,:);
    S_EE0 = S_EEm(sample_end,:);
    X_EI0 = X_EIm(sample_end,:);
    X_IE0 = X_IEm(sample_end,:);
    X_II0 = X_IIm(sample_end,:);
    X_EE0 = X_EEm(sample_end,:);
    Apost0 = Apostm(sample_end,:);
    Apre0 = Aprem(sample_end,:);
    W_IE0 = W_IEm(sample_end,:);
    W_EI0 = Weight_0; 
    left_sample_end = sample_end - 5/step;
    leftover_S_EI = S_EIm(left_sample_end:sample_end,:);
    leftover_S_IE = S_IEm(left_sample_end:sample_end,:);
    leftover_S_EE = S_EEm(left_sample_end:sample_end,:);
    leftover_S_II = S_IIm(left_sample_end:sample_end,:);
    spt_E0 = spt_Em(sample_end,:);
    
    
    
    if i == floor((5*duration/sample_duration)/100)
        disp('5% Complete')
    end    
    if i == floor((10*duration/sample_duration)/100)
        disp('10% Complete')
    end    
    if i == floor((15*duration/sample_duration)/100)
        disp('15% Complete')
    end    
    if i == floor((20*duration/sample_duration)/100)
        disp('20% Complete')
    end    
    if i == floor((duration/sample_duration)/4)
        disp('25% Complete')
    end    
    if i == floor(30*(duration/sample_duration)/100)
        disp('30% Complete')
    end    
    if i == floor(35*(duration/sample_duration)/100)
        disp('35% Complete')
    end    
    if i == floor(4*(duration/sample_duration)/10)
        disp('40% Complete')
    end    
    if i == floor(45*(duration/sample_duration)/100)
        disp('45% Complete')
    end    
    if i == floor((duration/sample_duration)/2)
        disp('50% Complete')
    end    
    if i == floor(55*(duration/sample_duration)/100)
        disp('55% Complete')
    end    
    if i == floor(6*(duration/sample_duration)/10)
        disp('60% Complete')
    end    
    if i == floor(65*(duration/sample_duration)/100)
        disp('65% Complete')
    end    
    if i == floor(70*(duration/sample_duration)/100)
        disp('70% Complete')
    end    
    if i == floor(75*(duration/sample_duration)/100)
        disp('75% Complete')
    end    
    if i == floor(8*(duration/sample_duration)/10)
        disp('80% Complete')
    end    
    if i == floor(8.5*(duration/sample_duration)/10)
        disp('85% Complete')
    end
    
    if i == floor(9*(duration/sample_duration)/10)
        disp('90% Complete')
    end    
    if i == floor(9.5*(duration/sample_duration)/10)
        disp('95% Complete')
    end
    
end
% -- run time --   
minute = toc/60     
 
figure(1)
plot(time/1000,J_I*W_IE,'k')
xlabel('Time (sec)')
ylabel('Average Synpatic Weight')
xlim([0 duration/1000])
ylim([0 300])


% calculated the Kuramoto Order Parameter
step = 0.1;
t = [0.1:step:duration];
[RE] = kuramoto_syn(spike_time_E,t,step,duration,N_E);

figure(2)
plot(t/1000,RE)
ylim([0 1.1])


function [RE] = kuramoto_syn(sptime,t,step,duration,N)
tspE = sptime;
phi = zeros(size(t,2),N);

% make vector of phases
for neuron = 1:N
        second_spike = 0;
    for i = 1:(duration)/step-1
        phi(i,neuron) = 2*pi*(t(i) - tspE(i,neuron));
        
        if tspE(i+1,neuron) ~= tspE(i,neuron)            
            if second_spike == 1
                delt = tspE(i+1,neuron) - tspE(i,neuron);
                a = floor(tspE(i,neuron)/step);
                b = floor(tspE(i+1,neuron)/step);
                unnorm = phi(a:b,neuron);
                phi(a:b,neuron) = unnorm./delt;
            end
            second_spike = 1;
        end
    end
end
RE = abs(mean((exp(j.*phi)),2));
end

