%2.c Dakota Ewing
clearvars
%                           Loop Parameters                              %
%------------------------------------------------------------------------%
Tf = 3;  % Tf is the simulation time in ms.
dt = .001;% dt is the time step.
% steps is the total number of timesteps.
Steps = Tf/dt + 1;
% T is the array in which we hold all the time values.
T = 0:dt:Tf;

%                          Variables & Arrays                            %
%------------------------------------------------------------------------%
%(microF/cm^2) The capacitance of the membrane.
Cm = 0.1; 
%(mV) The initial activation potential of the membrane.
Vi = -65; 
%(nA) The current applied to the neuron.
Iapp = 15; 
Ia = Iapp;

%(mV) The displacement from the equilibrium potential for K.
K_eq_dis = -77; 
%(mM/L) The concentration of K ions inside the mebrane.
K_conc_in = 150; 
%(mM/L) The concentration of K ions outside the mebrane.
K_conc_out = 5.5; 

% Ki is the array that holds all values of the potassium ions inside the 
% membrane at a given timestep.
Ki(1:Steps) = 0;
Ki(1) = K_conc_in;

% Ko is the array that holds all values of the potassium ions outside the 
% membrane at a given timestep.
Ko(1:Steps) = 0;
Ko(1) = K_conc_out;

%(mV) The displacement from the equilibrium potential for K.
Na_eq_dis = 50;
%(mM/L) The concentration of Na ions inside the mebrane.
Na_con_in = 15; 
%(mM/L) The concentration of Na ions outside the mebrane.
Na_conc_out = 150; 

% Nai is the array that holds all values of the sodium ions inside the 
% membrane at a given timestep.
Nai(1:Steps) = 0;
Nai(1) = K_conc_in;

% Nao is the array that holds all values of the sodium ions outside the 
% membrane at a given timestep.
Nao(1:Steps) = 0;
Nao(1) = K_conc_out;


%(mV) The displacement from the equilibrium potential for Leakage.
Leak_eq_dis = -54.4; 
%(mS/cm^2) The maximum K conductance.
K_m_cond = 36; 
%(mS/cm^2) The maximum Na conductance.
Na_m_cond = 120; 
%(mS/cm^2) The maximum Leakage conductance.
Leak_m_cond = 0.3; 


% (ms^-1)alpha n is the opening rate constant for n:
alpha_n = @(vt) (0.01*(vt + 55))/(1-exp(-1*((vt+55)/10)));
% (ms^-1)alpha m is the opening rate constant for m:
alpha_m = @(vt) (0.1*(vt + 40))/(1-exp(-1*((vt+40)/10)));
% (ms^-1)alpha h is the opening rate constant for h:
alpha_h = @(vt) 0.07*exp(-1*((vt+65)/20));

% (ms^-1)beta n is the closing rate constant for n:
beta_n = @(vt) 0.125*exp(-1*((vt+65)/80));
% (ms^-1)beta m is the closing rate constant for m:
beta_m = @(vt) 4*exp(-1*((vt+65)/18));
% (ms^-1)closing rate constant for h:
beta_h = @(vt) 1/(exp(-1*((vt+35)/10))+1);

%________________________________________________________________________%
%                      GATING VARIABLES                                  %

% n is the potassium activation gating variable, initial probability of K 
% gate being open:
n_initial = 0.317;
%(ms^-1)rate of change of n over time:
dNdt = @(vt,n) alpha_n(vt) * (1 - n) - beta_n(vt) * n; 
% N is the array that holds all values of the n at a given timestep.
N(1:Steps) = 0;
N(1) = n_initial;

% m is the sodium activation gating variable, initial probability of 
% Na gate being open:
m_initial = 0.05;
% (ms^-1)rate of change of m over time:
dMdt = @(vt,m) alpha_m(vt) * (1 - m) - beta_m(vt) * m;
% M is the array that holds all values of the m at a given timestep.
M(1:Steps) = 0;
M(1) = m_initial;

% h is the sodium inactivation gating variable, initial probability of Na 
% gate being deactivated:
h_initial = 0.6;                 
% (ms^-1)rate of change of h over time:
dHdt = @(vt,h) alpha_h(vt) * (1 - h) - beta_h(vt) * h;
% H is the array that holds all values of the h at a given timestep.
H(1:Steps) = 0;
H(1) = h_initial;




%The current from the potassium pump.
Ik = @(n,vt) K_m_cond * (n^4) * (vt - K_eq_dis);
% The current caused by leaks in the membrane.
Il = @(vt) Leak_m_cond * (vt - Leak_eq_dis);
% The current from the sodium pump.
INa = @(m,h,vt) Na_m_cond * (m^3) * h * (vt - Na_eq_dis);
% The current from the Na-K pump.
Ip = Il(Vi);

% The change in voltage over time for the membrane.
dVdt = @(n,m,h,vt,gk,gna,gnak,ca)(ca +gnak*Ip -gk*Ik(n,vt) -gna*INa(m,h,vt) -Il(vt))/Cm;

% The change in the sodium ion concentration inside the cell over time.
dNadt = @(vt, m, h) 3*Ip -INa(m,h,vt); 

% The change in the potassium ion concentration outside the cell over time.
dKdt = @(vt, n)  Ik(n,vt) - Il(vt) - 2*Ip;

% V is the array that holds all values of the voltage at a given time.
V(1:Steps) = 0;
V(1) = Vi;


%GATING VARIABLES:

% The gate for the potassium pump. when 0, the K pump is off, when 1, it is
% on.
Kg = 0;
% The gate for the sodium pump. when 0, the K pump is off, when 1, it is
% on.
Nag = 0;
% The gate for the Na-K pump. when 0, the Na-K pump is off, when 1, it is
% on.
NaKg = 1;
% The repolarization switch used to tell the simulation if the neuron is
% repolarizing. 
repolarize = 0;
% The sodium pump activation threshold.
Vna = -55;
% The potassium pump activation/sodium pump deactivation threshold.
Vk = 49.3;



%                             Simulation                                 %
%------------------------------------------------------------------------%

%loop:

%while the iterator is less than or equal to the total number of steps:
for s = 1:Steps
    
    % Applied current switch:
    if(T(s) >=0.5 && T(s) <=1)
        Ia = Iapp;
    else
        Ia = 0;
    end
    
    %RK4:
    n1 = dNdt(V(s), N(s))*dt;
    m1 = dMdt(V(s), M(s))*dt;
    h1 = dHdt(V(s), H(s))*dt;
    k1 = dKdt(V(s), N(s))*dt;
    na1= dNadt(V(s), M(s), H(s))*dt;
    v1 = dVdt(N(s), M(s), H(s), V(s), Kg, Nag, NaKg, Ia)*dt;
    
    
    n2 = dNdt(V(s)+v1/2, N(s)+n1/2)*dt;
    m2 = dMdt(V(s)+v1/2, M(s)+m1/2)*dt;
    h2 = dHdt(V(s)+v1/2, H(s)+h1/2)*dt;
    k2 = dKdt(V(s)+v1/2, N(s)+n1/2)*dt;
    na2= dNadt(V(s)+v1/2, M(s)+m1/2, H(s)+h1/2)*dt;
    v2 = dVdt(N(s)+n1/2, M(s)+m1/2, H(s)+h1/2, V(s)+v1/2, Kg, Nag, NaKg, Ia)*dt;
    
    n3 = dNdt(V(s)+v2/2, N(s)+n2/2)*dt;
    m3 = dMdt(V(s)+v2/2, M(s)+m2/2)*dt;
    h3 = dHdt(V(s)+v2/2, H(s)+h2/2)*dt;
    k3 = dKdt(V(s)+v2/2, N(s)+n2/2)*dt;
    na3= dNadt(V(s)+v2/2, M(s)+m2/2, H(s)+h2/2)*dt;
    v3 = dVdt(N(s)+n2/2, M(s)+m2/2, H(s)+h2/2, V(s)+v2/2, Kg, Nag, NaKg, Ia)*dt;
    
    n4 = dNdt(V(s)+v3, N(s)+n3)*dt;
    m4 = dMdt(V(s)+v3, M(s)+m3)*dt;
    h4 = dHdt(V(s)+v3, H(s)+h3)*dt;
    k4 = dKdt(V(s)+v3, N(s)+n3)*dt;
    na4= dNadt(V(s)+v3, M(s)+m3, H(s)+h3)*dt;
    v4 = dVdt(N(s)+n3, M(s)+m3, H(s)+h3, V(s)+v3, Kg, Nag, NaKg, Ia)*dt;
    
    N(s + 1) = N(s) +(n1 + 2*n2 + 2*n3 + n4)/6;
    M(s + 1) = M(s) +(m1 + 2*m2 + 2*m3 + m4)/6;
    H(s + 1) = H(s) +(h1 + 2*h2 + 2*h3 + h4)/6;
    Ko(s+ 1) = Ko(s)+(k1 + 2*k2 + 2*k3 + k4)/6;
    Nai(s+1) = Nai(s)+(na1 + 2*na2 + 2*na3 + na4)/6;
    V(s + 1) = V(s) +(v1 + 2*v2 + 2*v3 + v4)/6;
    
    
    % Gating channel switches:
    if(V(s)>= Vna && V(s)<=Vk && ~repolarize)
        Kg = 0;
        Nag = 1;
        repolarize = 0;
        G(s) = 10;
    elseif(V(s)>=Vk && ~repolarize)
        Kg = 1;
        Nag = 0;
        repolarize = 1;
        G(s) = 20;
    elseif((V(s+1)>V(s)) && repolarize)
        Kg = 0;
        Nag = 0;
        repolarize = 0;
        G(s) = 0;
    elseif(repolarize)
        Kg = 1;
        Nag = 0;
        repolarize = 1;
        G(s) = 20;
    end
    %Na-K Pump switch. If the sodium ion concentrain inside the membrane
    %and the potassium ion concentration outside the membrane are greater
    %than zero, the Na-K switch is on, else it is off.
    if  Nai(s)>0 && Ko(s)>0
        NaKg = 1;
    else
        NaKg = 0;
    end
   
end

%                               Plots                                    %
%------------------------------------------------------------------------%

plot(T(1:Steps), V(1:Steps))
title('2.c Graph of Membrane Potential Over Time:')
xlabel('Time(ms)') 
ylabel('Voltage(mV)') 
figure;


plot(T(1:Steps), M(1:Steps));
hold on;
plot(T(1:Steps), N(1:Steps));
hold on;
plot(T(1:Steps), H(1:Steps));
title('2.c Graph of n, m and h Over Time:')
xlabel('Time(ms)') 
ylabel('Probability(%)') 
legend('m','n','h','Location','northwest')

