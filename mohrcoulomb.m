%% Automated script for evaluating triaxial compression tests
% Mohr-Coulomb Yield Criteria for brittle failure

% Version 0.6

%% Possible improvements:
% -> add GUI with GUIDE
% -> add legend to figure
% -> add option to load data via .csv or .txt
% -> atm test results for intermediate shear stresses are overrepresented
% in the calculation of the mean values for the cohesion and the angle of
% internal friction due to a very simplified approach at calculating the
% mean values.
% This leads to an increasing divergence of the tangents at the smallest
% and highest values of the shear stress, while the fit in the middle is
% relatively accurate.
% Introducing a loading in the calculation which takes into account that
% a) tangents for intermediate values are more abundant
% b) the mean error for the tangents calculated from  Mohr-circles close
%    to each other, that is circles representing intermediate values, is
%    expected to be bigger
% could take care of that divergence and improve the accuracy of the mean
% values a lot.
% -> after that has been taken care of, a statistical error analysis should
% be introduced, preferably by someone who knows his mathematics.

%% License and Copyright:
% (c) Matthias Doerfler, 2016
% By using this code you agree to the conditions stated here:
% For private, scientific or educational use only.
% A commercial use is not allowed.
% The code may be modified and distributed under the following conditions:
% ->    the license and copyright as it is here must be included in the
%       distributed or modified code
% ->    The code must be released under the same license
% ->    Changes made to the code must be indicated in "Changelog"
% ->    This code is provided without warranty. The author can't be held
%       accountable for any damage caused by the use of this code.
% University of Freiburg, Germany

%% Changelog
% Version 0.1:
%   first version
% Version 0.2:
%   - added option for conversion from kN to MPa
%   - Sigma 1 and sigma 3 swapped (now correct use)
% Version 0.3:
%   - added detection system for imaginary and negative results
%   - improved readability
% Version 0.4:
%   - added option for correction for affected piston surface
%   - clearing out redundant varibales
% Version 0.5:
%   - results are now saved in a .mat file
% Version 0.6:
%   - added option to predict ideal values of theta and compare with
%     measured values (basically a tool to help you determine if your
%     samples behave according to the Mohr-Coulomb theory or not)
%     This option is to be assumed to be still methodically flawed
%   - added comments for features added in 0.3
%   - improved readability
%   - added arrays to document interval size and circles for postprocessing

%% Options

clear
% convert sigma 1 (Piston Force) from kN to MPa ?
convert = 0;
% correct for piston surface affected by fluid pressure ?
correct = 0;
% compare your results for the angle theta with the values for theta
% predicted by the Mohr-Coulomb theory ?
predict = 1;


%% Specifing dimensions of sample cylinders and piston surface

diameter = 0.04;                      % Diameter of the sample cylinder [m]
pistondiameter = 0.1;                 % Diameter of the piston [m]

%% Data input (your test results go here)

% Number of measurements (at least 2)
N = 3;
% Angle tetha [°] (only necessary if you wish to use the tetha-comparision
% function)
theta = [65.3 62.0 59.7];
% Maximum confining pressure before failure [MPa]
s3 = [0 5 10.1];              % Create a vector with the maximum confining pressures
                        % from each measurement
% Maximum force of the piston before failure [kN]
s1 = [42.4 58.6 74.8];        % Create a vector with the maximum piston forces
                        % from each measurement
% NOTICE: the test results must be sorted from smallest to biggest shear
% stress or the script will break / produce nonsense results.
% CAUTION: make sure that corresponding s1 and s3 values are always higher
% than the preceding values. Otherwise one Mohr-Circle will end up being
% entirely surrounded by another one which leads to an error.

%% Correcting piston force

if correct == 1 && diameter ~= pistondiameter
    A = pi*(diameter/2)^2;              % Surface area of the cylinder [m^2]
    affA = pistonA-A;                   % Area of the piston affected by the fluid pressure
    pistonA = pi*(pistondiameter/2)^2;  % Surface area of the piston [m^2]
    if convert == 1
        counterForce = s3.*(1000*affA); % Force working on the piston
        trueForce = s1-counterForce;    % effective force of the piston
        s1 = trueForce;
        disp('Piston force corrected for fluid pressure affecting piston.')
        clear('counterForce','trueForce')
    else
        counterPressure = s3.*(affA/A); % Pressure working on the piston
        trueSigma1 = s1-counterPressure;% effective pressure of the piston
        s1 = trueSigma1;
        disp('Sigma 1 corrected for fluid pressure affecting piston.')
        clear('counterPressure','trueSigma1')
    end
    clear('A','affA','pistonA')
else
    disp('No correction for piston surface affected by fluid pressure needed.')
    clear('pistondiameter')
end
clear('correct')

%% kn to mpa

if convert==1
    A = pi*(diameter/2)^2;              % Surface area of the cylinder [m^2]
    s1 = ((s1 .* 1000)./A) ./ 1000000;
    disp('Sigma 1 converted from kN to MPa.')
    clear('A')
else
    disp('Sigma 1 is already in MPa. No conversion needed.')
    clear('diameter')
end
clear('convert')

%% Magic happens here

% Brief explanation on what is going on here:
% Say you have done 5 measurements so there will be 5 circles in the
% Mohr-Coulomb diagram.
% The algorithm will calculate all possible common tangents to those
% circles:
% interval 1: tangents to circles: 1+2, 2+3, 3+4, 4+5
% interval 2: tangents to circles: 1+3, 2+4, 3+5
% interval 3: tangents to circles: 1+4, 2+5
% interval 4: tangent  to circles: 1+5
% Then the mean values of the tangents' slope (angle of internal friction)
% and y-intercept (cohesion) will be calculated.

% determining no. of steps necessary for the calculation
Nd = N-1;                 % No. of differently sized intervals for tangents
Nxsteps = zeros(1,Nd);
for i = 1:Nd
    Nxsteps(i) = N-i;     % No of steps for each interval size
end
% creating vector necessary for calculation of the absolute step no.
NxstepsP = zeros(1,Nd);
for i = 1:(Nd-1)
    NxstepsP(i+1) = sum(Nxsteps(1:i));
end

% calculating N cirlce parameters and curves
% Initalizing variables
M = zeros(N,1);                   % array used to store the circle centres
R = zeros(N,1);                   % array used to store the circle radii
X = zeros(100,N);                 % array used to store the x-values of the curve
Y = zeros(100,N);                 % as above, y-values
t = linspace(0,pi,100);
for i = 1:N
    % circle parameters
    m = ((s1(i)-s3(i))/2)+s3(i);  % centre
    M(i) = m;
    r = (s1(i)-s3(i))/2;          % radius
    R(i) = r;
    % circle curve
    x = r*cos(t);
    x = x';
    X((1+100*(i-1)):(100*i)) = x;
    y = r*sin(t);
    y = y';
    Y((1+100*(i-1)):(100*i)) = y;
end
clear('x','y','r','m','i','t')

% evaluation
% precalculations
sumSteps = 0.5*(N^2-N);                % sumSteps=no. of individual steps, relationship to N
xbase = 0:1:(round(max(s1))+10);       % x-values for calculation of tangents from poly
xbaseN = round(max(s1)+11);            % size of arrays in x dimension
% Initalizing variables
CN = zeros(sumSteps,1);                % Storage for cohesion values
myN = zeros(sumSteps,1);               % Storage for friction coefficient values
tanregcurveN = zeros(xbaseN,sumSteps); % Storage for tangent curves
imaginaryIndex = zeros(sumSteps,1);    % Indication if results were imaginary
negativeMYindex = zeros(sumSteps,1);   % Indication if friction coefficient was negative
negativeCindex = zeros(sumSteps,1);    % Indication if cohesion was negative
discardedTrue = 0;                     % Indication if any step has been discarded
predictedThetas = zeros(sumSteps,1);   % Storage for predicted values of tetha
IntervalDoc = zeros(sumSteps,1);       % Documentation of interval size at individual step
CirclesDoc = zeros(sumSteps,2);        % Documentation of circles used for calculation

for j = 1:Nd              % there are Nd differently sized intervals
    for i = 1:Nxsteps(j)  % each interval size has Nxsteps(j) individual steps
        % determining no. of cumulative individual step
        pos = i + NxstepsP(j);
        % documenting interval size and circles
        IntervalDoc(pos) = j;
        CirclesDoc(pos,1) = i;
        CirclesDoc(pos,2) = i+j;
        % get circle parameters
        m1 = M(i); m2 = M(i+j);
        r1 = R(i); r2 = R(i+j);
        % calculating common tangent
        m21 = m2-m1;
        d2 = m21^2;
        r21 = (r2-r1) / d2;
        s21 = sqrt(d2 - (r2-r1)^2) / d2;
        u11 = [-m21*r21, m21*s21];
        L11 = [m1,0] + r1*u11; L21 = [m2,0] + r2*u11;
        tanx1 = [L11(1),L21(1)]; tany1 = [L11(2),L21(2)];
        tanreg = polyfit(tanx1,tany1,1);
        thetaP = 45 + (atand(tanreg(1)) / 2);
        if isreal(tanreg) == 1
            if tanreg(1) >= 0
                if tanreg(2) >= 0
                    % saving cohesion and aoif at cumulative step 'pos'
                    myN(pos) = tanreg(1);
                    CN(pos) = tanreg(2);
                    % calculating curve of tangent
                    tanregcurve = polyval(tanreg,xbase);
                    % saving curve of tangent at cumulative step 'pos'
                    tanregcurve = tanregcurve';
                    tanregcurveN((1 + xbaseN * (pos-1)):(xbaseN * pos)) = tanregcurve;
                    % saving predicted value of theta
                    predictedThetas(pos) = thetaP;
                else
                    % saving cohesion and aoif at cumulative step 'pos'
                    myN(pos) = tanreg(1);
                    CN(pos) = tanreg(2);
                    % calculating curve of tangent
                    tanregcurve = polyval(tanreg,xbase);
                    % saving curve of tangent at cumulative step 'pos'
                    tanregcurve = tanregcurve';
                    tanregcurveN((1 + xbaseN * (pos-1)):(xbaseN * pos)) = tanregcurve;
                    % saving predicted value of theta
                    predictedThetas(pos) = thetaP;
                    % Warning (negative cohesion)
                    W = ['#01: The circles ',num2str(i),' and ',num2str(i+j),' have a common tangent with a negative y-intercept.'];
                    warning(W)
                    disp('Negative cohesion does not occur in nature.')
                    disp('The tangent will be flagged yellow.')
                    I = ['The interval size was ',num2str(j),'.'];
                    disp(I)
                    negativeCindex(pos) = 1;
                    clear('I','W','tanreg')
                end
            else
                % Warning (negative friction coefficient)
                W = ['#02: The circles ',num2str(i),' and ',num2str(i+j),' have a common tangent with a negative slope.'];
                warning(W)
                disp('Negative friction coefficient does not occur in nature.')
                disp('Results have been discarded.')
                I = ['The interval size was ',num2str(j),'.'];
                disp(I)
                negativeMYindex(pos) = 1;
                discardedTrue = 1;
                clear('I','W','tanreg','thetaP')
            end
        else
            % Warning (no common tangent)
            W = ['#03: The circles ',num2str(i),' and ',num2str(i+j),' do not share a common tangent.'];
            warning(W)
            disp('Values for C and my are imaginary.')
            disp('Results have been discarded.')
            I = ['The interval size was ',num2str(j),'.'];
            disp(I)
            imaginaryIndex(pos) = 1;
            discardedTrue = 1;
            clear('tanreg','I','W','m21','d2','r21','s21')
            clear('u11','L11','L21','tanx1','thetaP')
        end
    end
end
% clearing out obsolete variables
clear('pos','i','j','Nd','L11','L21','d2','m1','m2','m21')
clear('Nxsteps','NxstepsP','r1','r2','r21','s21','tanreg')
clear('tanregcurve','tanx1','tany1','u11','R','thetaP')

%% Sorting out positions of discarded results

if discardedTrue == 1   % checking if any results were discarded
    % copying results
    myNd = myN;
    CNd = CN;
    tanregcurveNd = tanregcurveN;
    predictedThetasd = predictedThetas;
    % precalculations
    discarded = imaginaryIndex+negativeMYindex; % combined discarded values
    negativeCsum = sum(negativeCindex); % number of results with a negative cohesion
    % Information output about amount of negative values for C
    if negativeCsum > 0
        I = ['A total number of ',num2str(negativeCsum),' values for C are negative.'];
        disp(I)
        disp('The results are carried on for further calculations anyway, since discarding them would lead to an overestimation of the mean value for C.')
    else
    end
    clear('myN','CN','tanregcurveN')
    % Initalizing variables
    discardedsum = 0;
    discardedIndex = zeros(sumSteps,1);

    for i = 1:sumSteps
        if discarded(i)~=0
            % counting discarded results
            discardedsum = discardedsum+1;
            % positions of discarded results
            discardedIndex(i) = 1;
        else
        end
    end

    trueSum = sumSteps-discardedsum;    % total number of valid results
    percentDiscarded = round((discardedsum/sumSteps)*100);  % percentage of results which were discarded
    I = ['A total number of ', num2str(discardedsum), ' steps have been discarded. That is ',num2str(percentDiscarded),'% of all possible steps.'];
    disp(I)
    I = ['The new total number of steps is ',num2str(trueSum)];
    disp(I)

    % Initalizing variables
    myN = zeros(trueSum,1);
    CN = zeros(trueSum,1);
    tanregcurveN = zeros(xbaseN,trueSum);
    predictedThetas = zeros(trueSum,1);
    IntervalDocTrue = zeros(trueSum,1);
    CirclesDocTrue = zeros(trueSum,2);
    for i = 1:sumSteps
        if discardedIndex(i)==0 % checking if position wasn't discarded
            % saving valid results in new variable
            CN(i-sum(discardedIndex(1:i))) = CNd(i);
            myN(i-sum(discardedIndex(1:i))) = myNd(i);
            tanregcurveN(:,i-sum(discardedIndex(1:i))) = tanregcurveNd(:,i);
            predictedThetas(i-sum(discardedIndex(1:i))) = predictedThetasd(i);
            IntervalDocTrue(i-sum(discardedIndex(1:i))) = IntervalDoc(i);
            CirclesDocTrue(i-sum(discardedIndex(1:i)),:) = CirclesDoc(i,:);
        else
        end
    end
else
    trueSum = sumSteps;
    disp('No results had to be discarded.')
end
% clearing out obsolete variables
clear('discardedsum','I','negativeCsum','percentDiscarded','tanregcurveNd')
clear('CNd','myNd','predictedThetasd')

%% Searching for interval leaps

IntSize = 1;
IntExtend = zeros((2*(N-1)),1);
IntExtend(1) = 1;
IntExtend(end) = sumSteps;
for i = 1:sumSteps
    if IntSize == IntervalDoc(i)
    else
        IntExtend((IntervalDoc(i)*2)-2) = i-1;
        IntExtend((IntervalDoc(i)*2)-1) = i;
        IntSize = IntervalDoc(i);
    end
end

if discardedTrue == 1
    IntSize = 1;
    IntExtendd = zeros((2*(N-1)),1);
    IntExtendd(1) = 1;
    IntExtendd(end) = trueSum;
    for i = 1:trueSum
        if IntSize == IntervalDocTrue(i)
        else
            IntExtendd((IntervalDocTrue(i)*2)-2) = i-1;
            IntExtendd((IntervalDocTrue(i)*2)-1) = i;
            IntSize = IntervalDocTrue(i);
        end
    end
else
end
clear('IntSize','i')


%% Mean values for theta
if predict == 1
    ThetaNsum = zeros(N,1);
    CircleCalls = zeros(N,1);
    if discardedTrue == 1
        for i = 1:trueSum
            ThetaNsum(CirclesDocTrue(i,1)) = ThetaNsum(CirclesDocTrue(i,1)) + predictedThetas(i);
            ThetaNsum(CirclesDocTrue(i,2)) = ThetaNsum(CirclesDocTrue(i,2)) + predictedThetas(i);
            CircleCalls(CirclesDocTrue(i,1)) = CircleCalls(CirclesDocTrue(i,1)) +1;
            CircleCalls(CirclesDocTrue(i,2)) = CircleCalls(CirclesDocTrue(i,2)) +1;
        end
        ThetaNmean = ThetaNsum ./ CircleCalls;
    else
        for i = 1:sumSteps
            ThetaNsum(CirclesDoc(i,1)) = ThetaNsum(CirclesDoc(i,1)) + predictedThetas(i);
            ThetaNsum(CirclesDoc(i,2)) = ThetaNsum(CirclesDoc(i,2)) + predictedThetas(i);
            CircleCalls(CirclesDoc(i,1)) = CircleCalls(CirclesDoc(i,1)) +1;
            CircleCalls(CirclesDoc(i,2)) = CircleCalls(CirclesDoc(i,2)) +1;
        end
        ThetaNmean = ThetaNsum ./ CircleCalls;
    end
else
end

%% Mean value for coefficient of internal friction and cohesion

% calculating mean values
myMean = mean(myN);
CMean = mean(CN);
disp('Mean friction coefficient (my):')
disp(myMean)
PhiMean = atand(myMean);
disp('Mean angle of internal friction (Phi) [�]:')
disp(PhiMean)
disp('Mean cohesion (C) [MPa]:')
disp(CMean)
% tangent calculated from mean values
Pmean = [myMean,CMean];
tanregcurveMean = polyval(Pmean,xbase);
yMax = max(tanregcurveMean);
clear('Pmean')

%% Plotting results

figure(1)
hold
axis equal
xlim([0 xbaseN])
ylim([0 yMax])
grid on
title('Mohr circles and failure envelope')
xlabel('\sigma_{n} [MPa]')
ylabel('\sigma_{s} [MPa]')

% tangents
for i = 1:trueSum
    tanregcurve = tanregcurveN(:,i);
    xbase;  % get
    figure(1)
    if CN(i) >= 0
        plot(xbase,tanregcurve,'green')
    else
        plot(xbase,tanregcurve,'yellow')
    end
end

% circles
for i = 1:N
    x = X((1+100*(i-1)):(100*i))+M(i);
    y = Y((1+100*(i-1)):(100*i));
    figure(1)
    plot(x,y,'blue')
end

% Mean tangent
figure(1)
plot(xbase,tanregcurveMean,'red')

clear('x','y','tanregcurve','i')

%% Theta graphs

n = 1:1:N;
figure(2)
ylim([0 90])
xlim([0 N+1])
plot(n,theta,'+',n,ThetaNmean,'.')
legend('measured','predicted')
title('Measured and predicted values for theta in comparision')
ylabel('\theta [°]')
xlabel('sample No.')
grid on

%% Saving

C = clock;
datename = [num2str(C(3)),'-',num2str(C(2)),'-',num2str(C(1)),'_',num2str(C(4)),'-',num2str(C(5))];
savename = ['mohrcoulomb_results',datename,'.mat'];
save(savename,'N','s1','s3','CMean','myMean')
