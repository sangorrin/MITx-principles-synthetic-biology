clear

% fetch the framework
urlwrite('http://web.mit.edu/20.305/www/part_composition_setup.m', ...
    'part_composition_setup.m');
rehash;
part_composition_setup('v5');

freqmultiplier = BioSystem();

% Species + initial conditions
A = freqmultiplier.AddCompositor('A', 0);
R = freqmultiplier.AddCompositor('R', 0); 
B = freqmultiplier.AddCompositor('B', 0); 
AB = freqmultiplier.AddCompositor('AB', 0); 
OUT = freqmultiplier.AddCompositor('OUT', 0); 

% Constants 
freqmultiplier.AddConstant('K_a', 1);
freqmultiplier.AddConstant('n_a', 4);
freqmultiplier.AddConstant('K_b', 1);
freqmultiplier.AddConstant('n_b', 4);
freqmultiplier.AddConstant('K_r', 1);
freqmultiplier.AddConstant('n_r', 4);
freqmultiplier.AddConstant('K_ab', 1);
freqmultiplier.AddConstant('n_ab', 4);
freqmultiplier.AddConstant('k_on', 0.2);
freqmultiplier.AddConstant('k_off', 0.8);
freqmultiplier.AddConstant('k_tln', 0.8);
freqmultiplier.AddConstant('k_trxo', 0.6);
freqmultiplier.AddConstant('k_trxb', 0.5);
freqmultiplier.AddConstant('k_trxr', 0.6);
freqmultiplier.AddConstant('k_mdeg', 0.5);
freqmultiplier.AddConstant('k_pdeg', 0.05);

% Parts
freqmultiplier.AddPart(Part('All interactions', [OUT A R B AB], ...
              [Rate('(k_trxo*k_tln/k_mdeg)*(K_ab^n_ab/(K_ab^n_ab+AB^n_ab))*((A^n_a/(A^n_a+K_a^n_a))+(B^n_b/(B^n_b+K_b^n_b)))-k_pdeg*OUT') ...  % dOUT/dt
               Rate('k_off*AB-k_on*A*B-k_pdeg*A') ... % dA/dt 
               Rate('(k_trxr*k_tln/k_mdeg)*(A^n_a/(A^n_a+K_a^n_a))-k_pdeg*R') ... % dR/dt 
               Rate('(k_trxb*k_tln/k_mdeg)*(K_r^n_r/(R^n_r+K_r^n_r))+k_off*AB-k_on*A*B-k_pdeg*B') ... % dB/dt 
               Rate('k_on*A*B-k_off*AB')]));  % dAB/dt 
           
% Iterations (comment out 1  group of iter_values and iter_constant)
% iter_values = 0.01:0.001:0.1;
% iter_constant = 'k_pdeg';
% iter_values = 0.4:0.1:0.9;
% iter_constant = 'k_trxo';
% iter_values = 0.01:0.01:0.9;
% iter_constant = 'k_trxb';
% iter_values = 0.1:0.01:0.9;
% iter_constant = 'k_trxr';
% iter_values = 2:1:6; 
% iter_constant = 'n_ab';
% iter_values = 2:1:6; 
% iter_constant = 'n_b';
% iter_values = 2:1:6; 
% iter_constant = 'n_a';
% iter_values = 2:1:6; 
% iter_constant = 'n_r';
% iter_values = 0.5:0.1:0.9; 
% iter_constant = 'k_off';
% iter_values = 0.0:0.1:0.5; 
% iter_constant = 'k_on';
% iter_values = 1:1:4; 
% iter_constant = 'K_a';
% iter_values = 1:1:4; 
% iter_constant = 'K_b';
% iter_values = 1:1:4; 
% iter_constant = 'K_ab';
% iter_values = 1:2:5; 
% iter_constant = 'K_r';

if not(exist('iter_values', 'var'))
    iter_values = -1;
    iter_constant = '';
    iterate = false;
else
    iterate = true;
end

avg_period = [];
std_deviation = [];

for val = iter_values
    if iterate
        freqmultiplier.ChangeConstantValue(iter_constant, val);
    end
    [t, y] = freqmultiplier.run_pulses([...
        Pulse(0, 'A', 10), ...
        Pulse(100, 'A', 0), ...
        Pulse(200, 'A', 10), ...
        Pulse(300, 'A', 0), ...
        Pulse(400, 'A', 10), ...
        Pulse(500, 'A', 0), ...
        Pulse(600, 'A', 10), ...
        Pulse(700, 'A', 0), ...
        Pulse(800, 'A', 10), ...
        Pulse(900, 'A', 0), ...
        Pulse(1000, 'A', 10), ...
        Pulse(1100, 'A', 0), ...
        Pulse(1200, 'A', 10), ...
        Pulse(1300, 'A', 0), ...
        Pulse(1400, 'A', 10), ...
        Pulse(1500, 'A', 0), ...
        Pulse(1600, 'A', 10), ...
        Pulse(1700, 'A', 0), ...
        Pulse(1800, 'A', 10), ...
        Pulse(1900, 'A', 0), ...
    ]);

    figure(1)
    hold on;
%     plot(t, y(:, freqmultiplier.CompositorIndex('B')));
%     plot(t, y(:, freqmultiplier.CompositorIndex('R')));
    plot(t, y(:, freqmultiplier.CompositorIndex('OUT')));
    
    out = y(:, freqmultiplier.CompositorIndex('OUT'));
    y2 = mean(out)*ones(1,length(t));
    [xout,yout] = intersections(t,out,t,y2,1);
    d = diff(xout(1:2:end));
    avg_period = [avg_period, mean(d)];
    std_deviation = [std_deviation, std(d)];
end

xlabel('Time (s)');
ylabel('Concentration (mM)');
if iterate
    title(sprintf('OUT signal %s', iter_constant), 'fontsize', 14)
    legend(cellstr(num2str(iter_values', 2)));
    plot(t, y(:, freqmultiplier.CompositorIndex('A')), '--');
else
    title(sprintf('A vs OUT signal'), 'fontsize', 14)
    plot(t, y(:, freqmultiplier.CompositorIndex('A')), '--');
    legend('OUT', 'A');
    avg_period
    std_deviation
end

if iterate
    figure(2)
    errorbar(iter_values,avg_period,std_deviation,'x-')
    xlabel(sprintf('value for %s', iter_constant));
    ylabel('mean period and standard deviation (s)');
    %legend(cellstr(num2str(iter_values', 2)));
end





