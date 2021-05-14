function plotJaynesCummings()

n = 3;

matrixSet = makeJCMatrices(n);

qubitFreq = 1;
coupling = .04;
detunings = linspace(-.3,.3,50);
harmonicFreqs = detunings * 0;
%detunings = 0;

figure('Name','Jaynes Cummings Hamiltonian');
for i = 1:length(detunings)  %length(H) is size of the basis we're using
    thisDetuning = detunings(i);
    harmonicFreq = qubitFreq - thisDetuning;
    harmonicFreqs(i) = harmonicFreq;
    H = JCHamiltonian(matrixSet,harmonicFreq,qubitFreq,coupling);
    allEigenvalues(i,:) = eig(H);
end

for i = 1:matrixSet.numStates
    plot(harmonicFreqs,allEigenvalues(:,i),'LineWidth',2);
    hold all;
end
title('Jaynes-Cummings Hamiltonian');
xlabel('Cavity Frequency');
text(0.8,3,sprintf('Qubit frequency %3.2f, coupling %3.2f',qubitFreq,coupling));
ylabel('Energy');
axis tight;
a = ylim;
ylim(a*1.1);
figure1 = gcf;
annotation(figure1,'arrow',[0.55 0.516071428571428],...
    [0.246619047619048 0.30952380952381]);

% Create arrow
annotation(figure1,'arrow',[0.544642857142857 0.516071428571428],...
    [0.444238095238095 0.35]);

% Create arrow
annotation(figure1,'arrow',[0.8375 0.8875],...
    [0.465666666666667 0.395238095238095]);

% Create arrow
annotation(figure1,'arrow',[0.205357142857143 0.160714285714286],...
    [0.194238095238095 0.254761904761905]);

% Create arrow
annotation(figure1,'arrow',[0.219642857142857 0.155357142857143],...
    [0.391857142857143 0.338095238095238]);

% Create arrow
annotation(figure1,'arrow',[0.846428571428571 0.883928571428571],...
    [0.706142857142857 0.65952380952381]);a
%labelstates(); only works for n = 3


function H = JCHamiltonian(matrixSet,harmonicFreq,qubitFreq,coupling)
    %see equation 1 on problem set
    H = (matrixSet.N * harmonicFreq) + (matrixSet.sigmaZ * qubitFreq * 0.5) + ...
        ((matrixSet.aDagger * matrixSet.sigmaMinus + matrixSet.a * matrixSet.sigmaPlus)  * coupling);
    
function matrixSet = makeJCMatrices(n)
    matrixSet.allstates = makeJCstates(n);
    numStates = length(matrixSet.allstates);
    aDagger = zeros(numStates);
    sigmaMinus = zeros(numStates);
    sigmaZ = zeros(numStates);
    %now constrict adagger1 and adagger2

    for (i = 1:numStates)
        for j = 1:numStates
            bra = matrixSet.allstates{i};
            ket = matrixSet.allstates{j};
            if (bra.n1 == ket.n1+1) && (bra.n2 == ket.n2)
                aDagger(i,j) = sqrt(ket.n1+1);
            end
            if (bra.n2 == ket.n2+1) && (bra.n1 == ket.n1)
                sigmaMinus(i,j) = sqrt(ket.n2+1);
            end
            if (bra.n2 == ket.n2) && (bra.n1 == ket.n1)
                sigmaZ(i,j) = -(2*ket.n2-1);
            end
        end
    end
    matrixSet.numStates = numStates;
    matrixSet.aDagger = aDagger;
    matrixSet.a = aDagger';
    
    matrixSet.sigmaMinus = sigmaMinus;
    matrixSet.sigmaPlus = sigmaMinus';
    
    matrixSet.sigmaZ = sigmaZ;
    matrixSet.N = matrixSet.aDagger * matrixSet.a;
    
    function allstates = makeJCstates(n)
    allstates = {};
    for (i = 1:n)
        for (j = 1:2)
            state.n1 = i-1;
            state.n2 = j-1;  %n2 is spin up or spin down
            if state.n2 == 2
                state.descriptor = sprintf('|%d,up>',state.n1);
            else
                state.descriptor = sprintf('|%d,down>',state.n1);
            end
            allstates{end+1} = state;
        end
    end
    
    
function labelstates()
            figure1 = gcf;
            ylim(axes1,[-0.55 3.85]);

annotation('textbox',dim,'String','|1,up>', 'FitBoxToText','on')

% Create textarrow
annotation(figure1,'textarrow',[0.255357142857143 0.196428571428571],...
    [0.201061855670103 0.136082474226804],'String',{'|0,\downarrow>'});

% Create textarrow
annotation(figure1,'textarrow',[0.214285714285714 0.160714285714286],...
    [0.442298969072165 0.360824742268041],'String',{'|1,\downarrow>'});

% Create textarrow
annotation(figure1,'textarrow',[0.782142857142857 0.841071428571429],...
    [0.508278350515464 0.445360824742268],'String',{'|1,\uparrow>'});