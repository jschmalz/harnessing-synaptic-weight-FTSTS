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
    J_E = 100; % NOTE: J_E = J_EI
    % TO MAKE FIGURE 12A uncomment
    J_I = 260; % NOTE: J_I = J_IE 
    
    % TO MAKE FIGURE 12B uncomment
%     J_I = 75; % NOTE: J_I = J_IE
    
    N_i = 1;
    C_E = 0.3*Ntot;
    C_I = 0.3*Ntot;
    tau_LTP = 20;
    tau_LTD = 22;   
    

    
% -- Make Random Synaptic Conncetions ---
    epsilon_E = 0.1; % connectivity
    epsilon_I = 0.1; % connectivity
    S_key_IE = zeros(N_E,N_I);
    S_key_EI = zeros(N_I,N_E);

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

    num_synapses_IE = max(max(S_key_IE));
    num_synapses_EI = max(max(S_key_EI));
    W_IE = zeros(duration_step,1); 
    W_IE_std = zeros(duration_step,1); 
    spE = zeros(duration_step,N_E);
    spI = zeros(duration_step,N_I);
    ISI_EI = zeros(duration_step, 1); 
    
% -- Initial Conditions --
    vE0 = 14*ones(1,N_E);
    vI0 = 14*ones(1,N_I);
    S_EI0 = zeros(1,N_E);
    S_IE0 = zeros(1,N_I);
    X_IE0 = zeros(1,N_I);
    X_EI0 = zeros(1,N_E);
    Apost0 = zeros(1,num_synapses_IE);
    Apre0 = zeros(1,num_synapses_IE);
    W_IE0 = Weight_0*ones(1,num_synapses_IE);   
    W_EI0 = Weight_0; 
    leftover_S_EI = zeros(5/step + 1 ,N_E);
    leftover_S_IE = zeros(5/step + 1 ,N_I);
    ref_E = zeros(1,N_E);
    ref_I = zeros(1,N_I);
    spike_E_time0 = 0;
    spike_I_time0 = 0;
    ISI_EI0 = 0;
    tau_E_m = 10; 
    tau_I_m = 10;
    
% -- Run --  
    sample_duration = 20;
    Synchrony = zeros(duration/sample_duration,1);
    time_syn = zeros(duration/sample_duration,1);

% Generate General Stimulation Patter
% TO MAKE FIGURE 12A uncomment
V_stim = 1;      
T_stim = 1;
x_neutral = 22;
multiple = 1;
deltaT_stim = 11;
t_pulse = T_stim*(x_neutral+multiple+1);

% TO MAKE FIGURE 12B uncomment
% V_stim = 1;      
% T_stim = 1;
% x_neutral = 10;
% multiple = 1;
% deltaT_stim = 5;
% t_pulse = T_stim*(x_neutral+multiple+1);

[Ue Ui] = pulsatile_input(multiple,V_stim,T_stim,x_neutral,duration,step,deltaT_stim);
sum(Ue)
sum(Ui)
cross_100 = 1;
Energy0 = 0;
 



for i = 1:duration/sample_duration
    % run parameters
    % TO MAKE FIGURE 12A uncomment
    if mean(W_IE0)*J_I < 75
        cross_100 = 0;
    end
    
    % TO MAKE FIGURE 12B uncomment
%     if mean(W_IE0)*J_I > 125
%         cross_100 = 0;
%     end

    ON = 1*(i*sample_duration >= 2000)*cross_100; % turn on control input
    plast_on = 1*(i*sample_duration >= 100); % past on/off
    comp_time = (i-1)*sample_duration;
    percent_V_stim = 1;
        
    % indexes
    a = 1 + (i>=2)*(i-1)*sample_duration/step;
    b = i*sample_duration/step;

    % TO MAKE FIGURE 12A uncomment
    Vstim = 200;
    
    
    % TO MAKE FIGURE 12B uncomment
%     Vstim = 100;
    
    ue = Vstim*Ue(1,a:b);
    ui = Vstim*Ui(1,a:b);
    
    
    [timem v_Em v_Im S_EIm S_IEm X_EIm X_IEm Apostm Aprem W_IEm spike_Em spike_Im ref_Em ref_Im synchronym spike_E_timem spike_I_timem ISI_EIm] = ode_neuron_model(plast_on,ON,vE0,vI0,S_EI0,S_IE0,X_EI0,X_IE0,Apost0,Apre0,W_IE0,W_EI0,mew_e,sigma_e,ue,ui,mew_i,sigma_i,J_E,J_I,C_E,C_I,tau_LTP,tau_LTD,step,sample_duration,N_E,N_I,S_key_EI,S_key_IE,leftover_S_EI,leftover_S_IE,ref_E,ref_I,Energy0,tau_E_m,tau_I_m, percent_V_stim,spike_E_time0,spike_I_time0,ISI_EI0,comp_time);
    
    
    
    % recorded variables
    time(a:b,:) = timem(1:sample_duration/step,:) + (i-1)*sample_duration;
    W_IE(a:b,1) = mean(W_IEm(1:sample_duration/step,:),2);
%     spike_E(a:b,:) = spike_Em(1:sample_duration/step,:);
%     spike_I(a:b,:) = spike_Im(1:sample_duration/step,:);
    Synchrony(i) = synchronym;
    time_syn(i) = sample_duration*(i);
    spE(a:b,:) = spike_E_timem(1:sample_duration/step,:);
    ISI_EI(a:b,1) = mean(ISI_EIm(1:sample_duration/step,:),2);
    
    % generate intial condition for next run
    sample_end = sample_duration/step ;
    vE0 = v_Em(sample_end,:);
    vI0 = v_Im(sample_end,:);
    S_EI0 = S_EIm(sample_end,:);
    S_IE0 = S_IEm(sample_end,:);
    X_EI0 = X_EIm(sample_end,:);
    X_IE0 = X_IEm(sample_end,:);
    Apost0 = Apostm(sample_end,:);
    Apre0 = Aprem(sample_end,:);
    W_IE0 = W_IEm(sample_end,:);
    W_EI0 = Weight_0; 
    left_sample_end = sample_end - 5/step;
    leftover_S_EI = S_EIm(left_sample_end:sample_end,:);
    leftover_S_IE = S_IEm(left_sample_end:sample_end,:);
    spike_E_time0 = spike_E_timem(sample_end,:);
    spike_I_time0 = spike_I_timem(sample_end,:);
    ISI_EI0 = ISI_EIm(sample_end,:);
    
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
    Kj = Energy0/1000
   

figure(1)
plot(time/1000,J_I*W_IE,'k')
ylim([0 300])
xlabel('Time (sec)')
ylabel('Average Synpatic Weight')
xlim([0 duration/1000])

% calculated the Kuramoto Order Parameter
step = 0.1;
t = [0.1:step:duration];
[RE] = kuramoto_syn(spE,t,step,duration,N_E);

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


