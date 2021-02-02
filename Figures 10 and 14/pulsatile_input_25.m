function [Ue Ui] = pulsatile_input(multi,V_stim,T_stim,x,duration,step,interval_time)
Ue = zeros(4,floor(duration/step) );
Ui = zeros(4,floor(duration/step));

% deals with any machine error in calcuation
epsilon = 0.001;
% biphasic pulse shape symetric = 1, asymetric > 1
pulse_shape = 1; 

% ue input
counter = 0;
counter2 = 0;
counter3 = 0;
counter4 = 0;
tk = 0;
A = 0;
for i = 1:duration/step
    
    if tk == 0
        rnm = randperm(4);
        stim_1 = rnm(1);
        stim_2 = rnm(2);
        stim_3 = rnm(3);
        stim_4 = rnm(4);
    end
     
tk = tk + step;
A = A +1;
    
% no. 1    
if tk > interval_time*(stim_1-1)/4 + epsilon  && tk < interval_time*(stim_1)/4 + step - epsilon 
    
    
    
    counter = counter + step;
    % anodic phase
    if counter >= step && counter < T_stim + step - epsilon
        Ue(1,i) = -V_stim/multi;
    end
    
    % cathodic phase
    if counter >= T_stim + step - epsilon && counter < 2*T_stim + step - epsilon
        Ue(1,i) = V_stim;
    end
    
    if counter >= 2*T_stim + step && counter < (2+x)*T_stim
        Ue(1,i) = 0;
    end
    
    if counter >= (2+x)*T_stim - epsilon && counter < (2+x+multi-1)*T_stim - epsilon
        Ue(1,i) = 0*-V_stim/multi;
    end
    
    if counter > (2+x+multi-1)*T_stim - epsilon
        counter = 0;
        Ue(1,i) = 0;
    end
end

     % no. 2
    
if tk > interval_time*(stim_2-1)/4 + epsilon   && tk < interval_time*(stim_2)/4 + step - epsilon   

    
    
    counter2 = counter2 + step;
    % anodic phase
    if counter2 >= step && counter2 < T_stim + step - epsilon
        Ue(2,i) = -V_stim/multi;
    end
    
    % cathodic phase
    if counter2 >= T_stim + step - epsilon && counter2 < 2*T_stim + step - epsilon
        Ue(2,i) = V_stim;
    end
    
    if counter2 >= 2*T_stim + step && counter2 < (2+x)*T_stim 
        Ue(2,i) = 0;
    end
    
     if counter2 >= (2+x)*T_stim - epsilon && counter2 < (2+x+multi-1)*T_stim - epsilon
        Ue(2,i) = 0*-V_stim/multi;
     end
    
    
    if counter2 > (2+x+multi-1)*T_stim - epsilon
        counter2 = 0;
        Ue(2,i) = 0;
    end
    
end   
     % no. 3
    
if tk > interval_time*(stim_3-1)/4 + epsilon && tk < interval_time*(stim_3)/4 + step - epsilon    
    

    
    counter3 = counter3 + step;
    % anodic phase
    if counter3 >= step && counter3 < T_stim + step - epsilon
        Ue(3,i) = -V_stim/multi;
    end
    
    % cathodic phase
    if counter3 >= T_stim + step - epsilon && counter3 < 2*T_stim + step - epsilon
        Ue(3,i) = V_stim;
    end
    
    if counter3 >= 2*T_stim + step && counter3 < (2+x)*T_stim + step
        Ue(3,i) = 0;
    end
    
     if counter3 >= (2+x)*T_stim - epsilon && counter3 < (2+x+multi-1)*T_stim - epsilon
        Ue(3,i) = 0*-V_stim/multi;
     end
    
    
    if counter3  > (2+x+multi-1)*T_stim - epsilon
        counter3 = 0;
        Ue(3,i) = 0;
    end
    
end
    
     % no.4

if tk > interval_time*(stim_4-1)/4 + epsilon && tk < interval_time*(stim_4)/4 + step - epsilon    
    

    
    
    counter4 = counter4 + step;
    % anodic phase
    if counter4 >= step && counter4 < T_stim + step - epsilon
        Ue(4,i) = -V_stim/multi;
    end
    
    % cathodic phase
    if counter4 >= T_stim + step - epsilon && counter4 < 2*T_stim + step - epsilon
        Ue(4,i) = V_stim;
    end
    
    if counter4 >= 2*T_stim + step && counter4 < (2+x)*T_stim + step
        Ue(4,i) = 0;
    end
    
     if counter4 >= (2+x)*T_stim - epsilon && counter4 < (2+x+multi-1)*T_stim - epsilon
        Ue(4,i) = 0*-V_stim/multi;
     end
    
    
    if counter4 > (2+x+multi-1)*T_stim - epsilon
        counter4 = 0;
        Ue(4,i) = 0*-V_stim/multi;
    end
     
   
end
  
    if tk > interval_time - epsilon
        tk = 0;
    end
 
    
end



%ui input
counter = 0;
counter1 = 0;
counter2 = 0;
counter3 = 0;
tk = 0;
for i = 1:duration/step
    
    if tk == 0
        rnm = randperm(4);
        stim_1 = rnm(1);
        stim_2 = rnm(2);
        stim_3 = rnm(3);
        stim_4 = rnm(4);
    end
  
 tk = tk + step;
    
% n0. 1
if tk > interval_time*(stim_1-1)/4 + epsilon && tk < interval_time*(stim_1)/4 + step - epsilon
    counter = counter + step;
    % cathodic phase
    if counter >= step && counter < T_stim + step - epsilon
        Ui(1,i) = V_stim;
    end
    
    % anodic phase
    if counter >= T_stim + step - epsilon && counter < (multi+1)*T_stim + step - epsilon
        Ui(1,i) = -V_stim/multi;
    end
    
    % neutral phase
    if counter >= (multi+1)*T_stim + step - epsilon  && counter < (multi+1+x)*T_stim - epsilon
        Ui(1,i) = 0;
    end
    
    if counter >= (multi+1+x)*T_stim - epsilon
       counter = 0;
       Ui(1,i) = 0*V_stim;
    end
end

% n0. 2

if tk > interval_time*(stim_2-1)/4 + epsilon && tk < interval_time*(stim_2)/4 + step - epsilon     
    counter1 = counter1 + step;
    % cathodic phase
    if counter1 >= step && counter1 < T_stim + step - epsilon
        Ui(2,i) = V_stim;
    end
    
    % anodic phase
    if counter1 >= T_stim + step - epsilon && counter1 < (multi+1)*T_stim + step - epsilon
        Ui(2,i) = -V_stim/multi;
    end
    
    % neutral phase
    if counter1 >= (multi+1)*T_stim + step - epsilon && counter1 < (multi+1+x)*T_stim - epsilon
        Ui(2,i) = 0;
    end
    
    if counter1 >= (multi+1+x)*T_stim - epsilon
       counter1 = 0;
       Ui(2,i) = 0*V_stim;
    end
end


% n0. 3

if tk > interval_time*(stim_3-1)/4 + epsilon && tk < interval_time*(stim_3)/4 + step - epsilon    
    counter2 = counter2 + step;
    % cathodic phase
    if counter2 >= step && counter2 < T_stim + step - epsilon
        Ui(3,i) = V_stim;
    end
    
    % anodic phase
    if counter2 >= T_stim + step - epsilon && counter2 < (multi+1)*T_stim + step - epsilon
        Ui(3,i) = -V_stim/multi;
    end
    
    % neutral phase
    if counter2 >= (multi+1)*T_stim + step - epsilon && counter2 < (multi+1+x)*T_stim - epsilon
        Ui(3,i) = 0;
    end
    
    if counter2 >= (multi+1+x)*T_stim - epsilon
       counter2 = 0;
       Ui(3,i) = 0*V_stim;
    end
end


% n0. 4

if tk > interval_time*(stim_4-1)/4 + epsilon && tk < interval_time*(stim_4)/4 + step - epsilon     
    counter3 = counter3 + step;
    % cathodic phase
    if counter3 >= step && counter3 < T_stim + step - epsilon
        Ui(4,i) = V_stim;
    end
    
    % anodic phase
    if counter3 >= T_stim + step - epsilon && counter3 < (multi+1)*T_stim + step - epsilon
        Ui(4,i) = -V_stim/multi;
    end
    
    % neutral phase
    if counter3 >= (multi+1)*T_stim + step - epsilon && counter3 < (multi+1+x)*T_stim - epsilon
        Ui(4,i) = 0;
    end
   
    if counter3 >= (multi+1+x)*T_stim - epsilon
       counter3 = 0;
       Ui(4,i) = 0*V_stim;
    end
end


 % reset overall counter
 if tk > interval_time - epsilon
     tk = 0;
 end
    
end
end