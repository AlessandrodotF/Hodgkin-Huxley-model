%function project(model,T,tmax,xmax,stim,dur) 
%Hodgkin Huxley model
%"model" select the equation to be solved respect to the assigment, enter 1
%for solving the first two requests of the exercise or 2 for the
%propagating action potential.
%T: temperature (°C) of the experiment, it can be 18.5°C or 37°C 
%tmax: Simulation time (ms)
%xmax: Simulation space (cm), used only for the propagating action potential
%stim:Currentstimulus (uA/cm^2) 
%dur:Duration of the stimulus (ms)
model=2;
T=18.5;
tmax=5;
xmax=0.1;
stim=100;
dur=0.5;
%Takes into account temperature-induced changes
factor=3^((T-18.5)/10); 
% define constant parameters
Vrest = -65.0;  % Membrane voltage at rest (mV)
gbarNa = 120.0; % Max Na conductivity (mS/cm^2)
vNa = 50.0;     % Equilibrium voltage for Na (mV)
gbarK = 36.0;   % Max K conductivity (mS/cm^2)
vK = -77.0;     % Equilibrium voltage for K (mV)
gbarL = 0.3;    % Leakage conductivity (mS/cm^2)
vL = -54.4;     % Equilibrium voltage for L (mV)

if model==1
    Cm=1;           % uF/cm^2
    dt = 0.01;      %ms
    t = 0:dt:tmax;  %ms


    % Set up vectors for membrane voltage, gating variables, and stimulus
    [Vm, m, h,Istim, n, INa, IK,IL] = deal(zeros(1,length(t)));
    % Apply stimulus
    dur_var = 0:dt:dur;
    Istim(1:length(dur_var)) = stim;
    Istim(floor(length(t)/2):floor(length(t)/2)+length(dur_var)) = stim;

    
    % Initialize state variables
    Vm(1) = Vrest; % Membrane begins at resting potential
    m(1) = alpham(Vrest,factor)/(alpham(Vrest,factor)+betam(Vrest,factor));
    h(1) = alphah(Vrest,factor)/(alphah(Vrest,factor)+betah(Vrest,factor));
    n(1) = alphan(Vrest,factor)/(alphan(Vrest,factor)+betan(Vrest,factor));
    
    % time loop
    for i=1:length(t)-1   
        % compute ionic currents
        INa(i) = gbarNa*(m(i)^3)*h(i)*(Vm(i)-vNa);
        IK(i) = gbarK*(n(i)^4)*(Vm(i)-vK);
        IL(i) = gbarL*(Vm(i)-vL);
        
        % update state variables-Euler method
        Vm(i+1) = Vm(i) - (dt/Cm)*(INa(i)+IK(i)+IL(i)-Istim(i));
        m(i+1) = m(i) + dt*(alpham(Vm(i),factor)*(1-m(i))-betam(Vm(i),factor)*m(i));
        h(i+1) = h(i) + dt*(alphah(Vm(i),factor)*(1-h(i))-betah(Vm(i),factor)*h(i));
        n(i+1) = n(i) + dt*(alphan(Vm(i),factor)*(1-n(i))-betan(Vm(i),factor)*n(i));   
    end
    
    % Plot membrane voltage 
    figure(1)
    plot(t,Vm)
    title("Membrane potential for I_{stim} = " + stim + " \muA/cm^2 ")
    ylabel('Membrane voltage (mV)'), xlabel('Time (ms)')
    grid on
    
    % Plot IK and INa
    figure(2)
    plot(t,INa,t,IK)
    legend('Na','K')
    title("Sodium and potassium ionic currents for I_{stim} = " + stim + " \muA/cm^2")
    xlabel('Time (ms)'), ylabel('Current (mA/cm^2)')
    grid on
    
    figure(3)
    plot(t,m,t,n,t,h)
    legend('m','n','h')
    title('Gating variables (m, h, n) vs time')
    xlabel('Time (ms)')
    grid on

else

    dt = 0.001;     %ms
    t = 0:dt:tmax;  %ms
    dx=0.001;       %cm        
    x = 0:dx:xmax;  %cm     
    r=0.025;        %cm axon radius 
    Cm = (2*pi*r*1);            %uF/cm 
    deltax=2*dx;
    
    
    %Rc value according to the temperature
    if T==37
        Rc=(25/(pi*(r^2)));
    else
        Rc=(35.4/(pi*(r^2)));   %Ohm/cm
    end  

    % Initialize state variables
    [ V_final, m, h, n, Istim, INa, IK, IL] = deal(zeros(length(t),length(x)));
    V_final(:,:)=Vrest; %axon initial state
    
    %stimulus for the central portion of the axon for "dur" ms
    dur_var = 0:dt:dur;
    halfx=floor(length(x)/2);
    Istim(1:length(dur_var),halfx-5:halfx+5)=stim;
    
    m(1,:) = alpham(Vrest,factor)/(alpham(Vrest,factor)+betam(Vrest,factor));
    h(1,:) = alphah(Vrest,factor)/(alphah(Vrest,factor)+betah(Vrest,factor));
    n(1,:) = alphan(Vrest,factor)/(alphan(Vrest,factor)+betan(Vrest,factor));
    for l=1:length(t)-1
        for p=2:length(x)-1
            INa(l,p) = gbarNa*(m(l,p)^3)*h(l,p)*(vNa-V_final(l,p));
            IK(l,p)  = gbarK*(n(l,p)^4)*(vK-V_final(l,p));
            IL(l,p)  = gbarL*(vL-V_final(l,p)); 
    
            V_final(l+1,p)= V_final(l,p)+(dt/Cm)*( (1/((deltax*deltax*Rc)))*(V_final(l,p-1)-2*V_final(l,p) ...
               +V_final(l,p+1))+(INa(l,p)+IK(l,p)+IL(l,p)+Istim(l,p)));
            m(l+1,p) = m(l,p) + dt*(alpham(V_final(l,p),factor)*(1-m(l,p))-betam(V_final(l,p),factor)*m(l,p));
            h(l+1,p) = h(l,p) + dt*(alphah(V_final(l,p),factor)*(1-h(l,p))-betah(V_final(l,p),factor)*h(l,p));
            n(l+1,p) = n(l,p) + dt*(alphan(V_final(l,p),factor)*(1-n(l,p))-betan(V_final(l,p),factor)*n(l,p));                                                                          
        end
    end
    figure(1)
    %just for plotting procedure
    t1=((tmax/dt)*0.25)/tmax;
    t2=((tmax/dt)*0.5)/tmax;
    t3=((tmax/dt)*0.85)/tmax;
    t4=((tmax/dt)*1.6)/tmax;
    t5=((tmax/dt)*2.35)/tmax;
    
    figure(1)
    plot(x,V_final(t1,:),x,V_final(t2,:),x,V_final(t3,:),x,V_final(t4,:),x,V_final(t5,:),x,V_final(floor(3*length(t)/4),:),x,V_final(length(t),:),'LineWidth',1)
    legend("t="+dt*t1+"ms","t="+dt*t2+"ms","t="+dt*t3+"ms" ,"t="+dt*t4+"ms","t="+dt*t5+"ms","t="+dt*floor(3*length(t)/4)+"ms","t="+dt*length(t)+"ms",'Location','eastoutside','Orientation','vertical','Box','off','FontSize',14, "LineWidth",1)
    title("Membrane potential for I_{stim} = "+ stim +"  \muA/cm^2 ")
    ylabel('Membrane voltage (mV)'), xlabel('space (cm)')
    grid on
    figure(2)
    plot(t,V_final(:,450),'LineWidth',1)
    title("Membrane potential for I_{stim} = "+ stim +"  \muA/cm^2 ")
    ylabel('Membrane voltage (mV)'), xlabel('time (ms)')
    grid on


end


% Alpha and Beta functions for the gating variables 
% They include the temperature releted factor
function [aN] = alphan(v,factor)
aN = factor*((0.01*(v+55)) / (1-exp(-(v+55)/10)));
end

function [bN] = betan(v,factor)
bN =factor*( 0.055*exp(-v/80));
end

function [aM] = alpham(v,factor)
aM = factor*(0.1*(v+40)) / (1-exp(-(v+40)/10));
end

function [bM] = betam(v,factor)
bM =factor*(0.108*exp(-v/18));
end

function [aH] = alphah(v,factor)
aH = factor*(0.0027*exp(-v/20));
end

function [bH] = betah(v,factor)
bH = factor*(1/(1+exp(-(v+35)/10)));
end

