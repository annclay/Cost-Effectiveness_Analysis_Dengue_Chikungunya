function [t, popChik, popDen, chikTable, denTable, incidenceChik, incidenceDen, cumulativeChik, cumulativeDen, population] = IV_Chik_Dengue_SIR_CombinedModel_VacStateInfs_March30_20v3(r,TranChik,TranDen,recChik,recDen,nu,mu,X0Chik, X0Den, Y0Chik, Y0Den, Z0Chik, StartTime, MaxTime, pSeqChik, pSeqRecovery, symptChik, symptDen0, symptDen1, HRChik, HRDen,limRepeatInfDen,colombiaPop, vaxDenCoverage, vaxChikCoverage, vaxDenEfficacy, vaxChikEfficacy, denSens, denSpec, vaxDenInitialVax)

%Chikungunya Dengue Transmission models
% written by Anneke Claypool
%August 30, 2021

% This version includes the serotype test for dengue as a pre-rec for
% receiving the dengue vaccine. It also includes the sensitivity and
% specificity for that test. 

% This version of the code was updated to only vaccinate people who have
% already been tested. 

% This model is adapted from the MATLAB version 4.2 from page 123 of 
% "Modeling Infectious Disease in humans and animals" 
% by Keeling & Rohani.
% 
%It has been combined to make it a dual Chikungunya Dengue model

%% Combined Model
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
nsteps= MaxTime;
tspan = linspace(1,MaxTime, nsteps); %create timespan with results for each time step
options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);%options = odeset('RelTol', 1e-1, 'AbsTol', 1e-3);%'RelTol', 1e-1, 'AbsTol', 1e-3);%, 'RelTol', 1, 'AbsTol', 1e-2)
%options = odeset('RelTol', 1e-4);

VaxBounceBack = 1-0.6;
%lastwarn(''); % Clear last warning message
%Prepare the population dataset
%where population = [S_HchikS_Hden S_Mchik I_H0chik I_H1chik I_Mchik R0chik R1chik cumulativeSymptInfectionschik totalDeaths chikDeaths S_H0den S_H1den S_Mden I_H0den I_H1den I_Mden cumulativeSymptInfectionsden denDeaths]

%% Initialize States

totsPop = X0Chik(1);
nVaccinatedDen = vaxDenCoverage*totsPop; %round(
nVaccinatedChik = round(vaxChikCoverage*totsPop,0);
initialTests = (vaxDenInitialVax*totsPop); 
        VaxDenPropS0ini = vaxDenInitialVax*(1-denSpec);
        VaxDenPropS1ini = vaxDenInitialVax*(denSens);
        
% Susceptible Humans 0 Dengue
%pop(1) = S_HchikS_H0den
S_HchikS_H0den0 = ((1-vaxChikCoverage)*(X0Chik(1)-X0Den(2)-Y0Den(1)-Y0Den(2))-Y0Chik(1)-Y0Chik(2)-Z0Chik(1)-Z0Chik(2))*(1-VaxDenPropS0ini);%-nVaccinatedChik;
%S_HchikS_H0den0 = round((X0Chik(1)-Y0Chik(1)-Y0Chik(2)-Z0Chik(1)-Z0Chik(2)-X0Den(2)-Y0Den(1)-Y0Den(2)-nVaccinatedChik),0);%*(1-vaxDenCoverage);


%pop(2) = I_H0chikS_H0den
I_H0chikS_H0den0 = Y0Chik(1)*(1-VaxDenPropS0ini);
%pop(3) = I_H1chikS_H0den
I_H1chikS_H0den0 = Y0Chik(2)*(1-VaxDenPropS0ini);
%pop(4) = R_H0chikS_H0den
R_H0chikS_H0den0 = Z0Chik(1)*(1-VaxDenPropS0ini);
%pop(5) = R_H1chikS_H0den
R_H1chikS_H0den0 = Z0Chik(2)*(1-VaxDenPropS0ini);

% Susceptible Humans 1 Dengue
%pop(6) = S_HchikS_H1den
S_HchikS_H1den0 = (1-vaxChikCoverage)*(X0Den(2))*(1-VaxDenPropS1ini);
%pop(7) = I_H0chikS_H1den
I_H0chikS_H1den0 = 0;
%pop(8) = I_H1chikS_H1den
I_H1chikS_H1den0 = 0;
%pop(9) = R_H0chikS_H1den
R_H0chikS_H1den0 = 0;
%pop(10) = R_H1chikS_H1den
R_H1chikS_H1den0 = 0;

% Infected Humans 0 Dengue
%pop(11) = S_HchikI_H0den
S_HchikI_H0den0 = (1-vaxChikCoverage)*Y0Den(1);
%pop(12) = I_H0chikI_H0den
I_H0chikI_H0den0 = 0;
%pop(13) = I_H1chikI_H0den
I_H1chikI_H0den0 = 0;
%pop(14) = R_H0chikI_H0den
R_H0chikI_H0den0=0;
%pop(15) = R_H1chikI_H0den
R_H1chikI_H0den0=0;

% Infected Humans 1 Dengue
%pop(16) = S_HchikI_H1den
S_HchikI_H1den0 = (1-vaxChikCoverage)*Y0Den(2);
%pop(17) = I_H0chikI_H1den
I_H0chikI_H1den0=0;
%pop(18) = I_H1chikI_H1den
I_H1chikI_H1den0=0;
%pop(19) = R_H0chikI_H1den
R_H0chikI_H1den0=0;
%pop(20) = R_H1chikI_H1den
R_H1chikI_H1den0=0;

%Mosquitos
%pop(21)=S_MchikS_Mden
S_MchikS_Mden0=X0Chik(2);
%pop(22)=I_MchikS_Mden
I_MchikS_Mden0 = Y0Chik(3);
%pop(23)=S_MchikI_Mden
S_MchikI_Mden0 = Y0Den(3);
%pop(24)=I_MchikI_Mden
I_MchikI_Mden0=0;

%Vaccine Initial States
%Dengue Vaccine   
% S_HchikV_Hden = pop(30);
S_HchikV_Hden0 = ((1-vaxChikCoverage)*(X0Chik(1)-X0Den(2)-Y0Den(1)-Y0Den(2))-Y0Chik(1)-Y0Chik(2)-Z0Chik(1)-Z0Chik(2))*(VaxDenPropS0ini)+ (1-vaxChikCoverage)*(X0Den(2))*(VaxDenPropS1ini);
%         I_H0chikV_Hden = pop(31); 
I_H0chikV_Hden0 = Y0Chik(1)*VaxDenPropS0ini;
%         I_H1chikV_Hden = pop(32);
I_H1chikV_Hden0 = Y0Chik(2)*VaxDenPropS0ini;
%         R_H0chikV_Hden = pop (33);
R_H0chikV_Hden0 = Z0Chik(1)*VaxDenPropS0ini;
%         R_H1chikV_Hden = pop(34); 
R_H1chikV_Hden0 = Z0Chik(2)*VaxDenPropS0ini; 
%    
%         %Chikungunya Vaccine
%         V_HchikS_H0den = pop(35);
V_HchikS_H0den0 = vaxChikCoverage*(X0Chik(1)-X0Den(2)-Y0Den(1)-Y0Den(2)); %-Y0Chik(1)-Y0Chik(2)    (1-vaxDenCoverage)*
%V_HchikS_H0den0 = nVaccinatedChik- (vaxChikCoverage*nVaccinatedDen); %(1-vaxDenCoverage)*
%         V_HchikS_H1den = pop(36);
V_HchikS_H1den0 =vaxChikCoverage*(X0Den(2))*(1-VaxDenPropS1ini);
%         V_HchikI_H0den = pop(37);
V_HchikI_H0den0 = Y0Den(1)*vaxChikCoverage; %*
%         V_HchikI_H1den = pop(38);
V_HchikI_H1den0= Y0Den(2)*vaxChikCoverage;
%         V_HchikV_Hden = pop(39);
V_HchikV_Hden0 = (vaxChikCoverage*(X0Chik(1)-X0Den(2)-Y0Den(1)-Y0Den(2))-Y0Chik(1)-Y0Chik(2)-Z0Chik(1)-Z0Chik(2))*(VaxDenPropS0ini) + (vaxChikCoverage)*(X0Den(2))*(VaxDenPropS1ini);%vaxDenCoverage*(X0Chik(1)-Y0Chik(1)-Y0Chik(2)-Z0Chik(1)-Z0Chik(2)-X0Den(2)-Y0Den(1)-Y0Den(2));
% 

    
% Extra states
%pop(25) = deathsTot
deathsTot0 = 0;
%pop(26) = deathsChik
deathsChik0 = 0;
%pop(27) = deathsDen
deathsDen0 = 0;
%pop(28) = incChik
incChik0= round(Y0Chik(2),10);
%pop(29) = incDen
incDen0= round(Y0Den(2),10);
%pop(40) = totalTests
totalTests0 = initialTests;
%pop(41) = totalVax;
totalVax0 = ((X0Chik(1)-X0Den(2)-Y0Den(1)-Y0Den(2))-Y0Chik(1)-Y0Chik(2)-Z0Chik(1)-Z0Chik(2))*(VaxDenPropS0ini)+ X0Den(2)*(VaxDenPropS1ini) +Y0Chik(1)*VaxDenPropS0ini +Y0Chik(2)*VaxDenPropS0ini +Z0Chik(1)*VaxDenPropS0ini+Y0Den(2)*vaxChikCoverage;

%Organize initial states
%Assume that all inital infections are not overlapped between Chikungunya
%and Dengue


initialStates = [S_HchikS_H0den0 I_H0chikS_H0den0 I_H1chikS_H0den0 R_H0chikS_H0den0 R_H1chikS_H0den0 S_HchikS_H1den0 I_H0chikS_H1den0 I_H1chikS_H1den0 R_H0chikS_H1den0 R_H1chikS_H1den0 S_HchikI_H0den0 I_H0chikI_H0den0 I_H1chikI_H0den0 R_H0chikI_H0den0 R_H1chikI_H0den0 S_HchikI_H1den0 I_H0chikI_H1den0 I_H1chikI_H1den0 R_H0chikI_H1den0 R_H1chikI_H1den0 S_MchikS_Mden0 I_MchikS_Mden0 S_MchikI_Mden0 I_MchikI_Mden0 deathsTot0 deathsChik0 deathsDen0 incChik0 incDen0 S_HchikV_Hden0 I_H0chikV_Hden0 I_H1chikV_Hden0 R_H0chikV_Hden0 R_H1chikV_Hden0 V_HchikS_H0den0 V_HchikS_H1den0 V_HchikI_H0den0 V_HchikI_H1den0 V_HchikV_Hden0 totalTests0 totalVax0];

% if double(sum(sum(initialStates(1:20))+sum(initialStates(30:39))))~= X0Chik(1)
%     difference = X0Chik(1) - double(sum(sum(initialStates(1:20))+sum(initialStates(30:39))));
%     initialStates(1) = initialStates(1) + difference;
% end

[t, population]=ode23t(@Diff,tspan,initialStates,options); % Use ode23t, The main ODE, Originally ode23tb


%% Organize Resulting Data

        %Chik initial variables
        S_HchikResults = population(:,1) + population(:,6) + population(:,11) + population(:,16) + population(:,30);
        S_MchikResults = population(:,21) + population(:,23);
        I_H0chikResults = population(:,2) +  population(:,7) + population(:,12)+ population(:,17) + population(:,31);
        I_H1chikResults = population(:,3) + population(:,8) + population(:,13) +population(:,18)+ population(:,32); 
        I_MchikResults= population(:,22) + population(:,24);
        R_H0chikResults = population(:,4) + population(:,9) +population(:,14) + population(:,19) + population(:,33); 
        R_H1chikResults = population(:,5) + population(:,10) + population(:,15) + population(:,20) + population(:,34);
        chikCasesCumulative = population(:,28);
        deathsCumulative = population(:,25);
        deathsChikCumulative = population(:,26);
        V_HchikResults = sum(population(:,35:39),2);
        
        %Den initial variables
        S_H0denResults = sum(population(:,1:5),2) + population(:,35);
        S_H1denResults = sum(population(:,6:10),2)+ population(:,36);
        S_MdenResults = sum(population(:,21:22),2);
        I_H0denResults = sum(population(:,11:15),2)+ population(:,37);
        I_H1denResults = sum(population(:,16:20),2)+ population(:,38); 
        I_MdenResults= sum(population(:,23:24),2);
        denCasesCumulative = population(:,29);
        deathsDenCumulative = population(:,27);
        V_HdenResults = sum(population(:,30:34),2)+ population(:,39);

%Separate Chik and Den populations
popChik = [S_HchikResults S_MchikResults I_H0chikResults I_H1chikResults I_MchikResults R_H0chikResults R_H1chikResults chikCasesCumulative deathsCumulative deathsChikCumulative V_HchikResults];
popDen = [S_H0denResults S_H1denResults S_MdenResults I_H0denResults I_H1denResults I_MdenResults denCasesCumulative deathsCumulative deathsDenCumulative V_HdenResults] ; %Add back in total deaths

%Calculate the resulting Chikungunya incidence
incidenceChik = IncidenceChik(population, colombiaPop);
cumulativeChik = cumsum(incidenceChik);

%Create a table of the results of the SIR models with differeing symptoms:Chik
colNamesChik = {'Susceptible_Humans','Susceptible_Mosquitos','Infected_Asympt_Humans', 'Infected_Sympt_Humans', 'Infected_Mosquitos', 'Recovered_Humans_No_Seq', 'Recovered_Humans_Seq' 'Cumulative_Cases', 'Cumulative_Deaths', 'Cumulative_Chik_Deaths', 'Vax_Humans_Chik'};
chikTable = array2table(popChik,'VariableNames',colNamesChik);

%Create a table of the results of the SIR models with differeing symptoms:Den
colNamesDen = {'Susceptible_Humans_noPrevInf','Susceptible_Humans_PrevInf','Susceptible_Mosquitos','Infected_Asympt_Humans', 'Infected_Sympt_Humans', 'Infected_Mosquitos','Cumulative_Cases', 'Cumulative_Deaths', 'Cumulative_Den_Deaths', 'Vax_Humans_Den'};
denTable = array2table(popDen,'VariableNames',colNamesDen);

%Calculates the incidence rate over time for Dengue
incidenceDen = IncidenceRatesDen(population,colombiaPop);
cumulativeDen = cumsum(incidenceDen);

% [warnMsg, warnId] = lastwarn;
%                 if ~isempty(warnMsg)
%                     error(warnMsg);
%                 end


%% Functions

% Calculates the differential rates used in the integration for Chikungunya.
    function dPop=Diff(t, pop)    

        
 % Susceptible Humans 0 Dengue
        S_HchikS_H0den = pop(1);
        I_H0chikS_H0den = pop(2);
        I_H1chikS_H0den = pop(3);
        R_H0chikS_H0den = pop(4);
        R_H1chikS_H0den = pop(5);

% Susceptible Humans 1 Dengue
        S_HchikS_H1den= pop(6);
        I_H0chikS_H1den= pop(7);
        I_H1chikS_H1den= pop(8);
        R_H0chikS_H1den = pop(9);
        R_H1chikS_H1den = pop(10);

% Infected Humans 0 Dengue
        S_HchikI_H0den = pop(11);
        I_H0chikI_H0den = pop(12);
        I_H1chikI_H0den = pop(13);
        R_H0chikI_H0den = pop(14);
        R_H1chikI_H0den = pop(15);

% Infected Humans 1 Dengue
        S_HchikI_H1den = pop(16);
        I_H0chikI_H1den = pop(17);
        I_H1chikI_H1den = pop(18);
        R_H0chikI_H1den = pop(19);
        R_H1chikI_H1den = pop(20);

%Total Human Population
        TotalHumanPop = sum(pop(1:20))+ sum(pop(30:39));

%Mosquitos
        S_MchikS_Mden = pop(21);
        I_MchikS_Mden = pop(22);
        S_MchikI_Mden = pop(23);
        I_MchikI_Mden = pop(24);
        
%Total Mosquito Population
        TotalMosPop = sum(pop(21:24));
        
% Extra states
        deathsTot = pop(25);
        deathsChik = pop(26);
        deathsDen = pop(27);
        incChik = pop(28);
        incDen = pop(29);
        
%Vaccine States
        %Dengue Vaccine
         
        S_HchikV_Hden = pop(30);
        I_H0chikV_Hden = pop(31);
        I_H1chikV_Hden = pop(32);
        R_H0chikV_Hden = pop (33);
        R_H1chikV_Hden = pop(34);        
        
        %Chikungunya Vaccine
        V_HchikS_H0den = pop(35);
        V_HchikS_H1den = pop(36);
        V_HchikI_H0den = pop(37);
        V_HchikI_H1den = pop(38);
        %Chikungunya and Dengue vaccine
        V_HchikV_Hden = pop(39);
        

%Total Human Population
        TotalHumanPop = sum(pop(1:20))+ sum(pop(30:39));

        %Chik initial variables
        S_Hchik = S_HchikS_H0den + S_HchikS_H1den + S_HchikI_H0den + S_HchikI_H1den + S_HchikV_Hden;
        S_Mchik = S_MchikS_Mden + S_MchikI_Mden;
        I_H0chik = I_H0chikS_H0den +  I_H0chikS_H1den + I_H0chikI_H0den+ I_H0chikI_H1den + I_H0chikV_Hden;
        I_H1chik = I_H1chikS_H0den + I_H1chikS_H1den + I_H1chikI_H0den +I_H1chikI_H1den + I_H1chikV_Hden; 
        I_Mchik= I_MchikS_Mden + I_MchikI_Mden;
        R_H0chik = R_H0chikS_H0den + R_H0chikS_H1den +R_H0chikI_H0den + R_H0chikI_H1den + R_H0chikV_Hden; 
        R_H1chik = R_H1chikS_H0den + R_H1chikS_H1den + R_H1chikI_H0den + R_H1chikI_H1den + R_H1chikV_Hden;
        
        %Den initial variables
        S_H0den = S_HchikS_H0den + I_H0chikS_H0den + I_H1chikS_H0den + R_H0chikS_H0den + R_H1chikS_H0den + V_HchikS_H0den;
        S_H1den = S_HchikS_H1den + I_H0chikS_H1den + I_H1chikS_H1den + R_H0chikS_H1den + R_H1chikS_H1den + V_HchikS_H1den; 
        S_Mden = S_MchikS_Mden + I_MchikS_Mden;
        I_H0den = S_HchikI_H0den + I_H0chikI_H0den + I_H1chikI_H0den +R_H0chikI_H0den + R_H1chikI_H0den + V_HchikI_H0den;
        I_H1den = S_HchikI_H1den + I_H0chikI_H1den + I_H1chikI_H1den + R_H0chikI_H1den +R_H1chikI_H1den + V_HchikI_H1den; 
        I_Mden= S_MchikI_Mden + I_MchikI_Mden;
        
        V_Hden = S_HchikV_Hden + I_H0chikV_Hden + I_H1chikV_Hden + R_H0chikV_Hden + R_H1chikV_Hden;
        V_Hchik = V_HchikS_H0den + V_HchikS_H1den + V_HchikI_H0den + V_HchikI_H1den;

        V_Hboth = V_HchikV_Hden;
        dPop=zeros(39,1);
        
        
        %Total Population
        Htot = sum(pop(1:20))+ sum(pop(30:39));%S_Hchik + I_H0chik + I_H1chik + R_H0chik + R_H1chik;
        %Mtotchik = S_Mchik + I_Mchik;
        %Htotden = S_H0den + S_H1den + I_H0den + I_H1den;
        Mtot = sum(pop(21:24));

        %Transmission variables
        Tchik= TranChik;
        gammachik = recChik;
        Tden= TranDen;
        gammaden = recDen;
        HRden = HRDen;
        lden=limRepeatInfDen;

        %% Combined Model
        
        %Humans
        %S_Hchik S_Hden
        infRateChik = r*(Tchik(1,2)*I_Mchik)/Htot;
        infRateDen = r*(Tden(1,2)*I_Mden)/Htot;
        infRateDen2 = r*(Tden(1,2)*I_Mden)*lden/Htot;
        recoverRateChik = gammachik(1);
        recoverRateDen = gammaden(1);
        
        %Vaccination Variables
        vaxRateWeek = vaxDenCoverage/52;
        sensitivity = denSens;
        specificity = denSpec;
        VaxProp = (vaxRateWeek*Htot)/(S_H0den + S_H1den);
        VaxPropS0 = VaxProp*(1-specificity);
        VaxPropS1 = VaxProp*(sensitivity);
        
        %% Susceptible Humans 0 Dengue
        %dS_Hchik S_H0den/dt 
        dPop(1) = nu(1)*TotalHumanPop - (infRateChik + infRateDen - infRateChik*infRateDen)*S_HchikS_H0den- VaxPropS0*S_HchikS_H0den - mu(1)*S_HchikS_H0den;  %checked
        %I_H0chik S_H0den
        dPop(2) = infRateChik*(1-symptChik)*(1-infRateDen)*(S_HchikS_H0den + (1-vaxChikEfficacy)*V_HchikS_H0den) - (infRateDen + recoverRateChik - infRateDen*recoverRateChik)*I_H0chikS_H0den - VaxPropS0*I_H0chikS_H0den- mu(1)*I_H0chikS_H0den; 
        %dI_H1chikS_H0den/dt
        dPop(3)= infRateChik*symptChik*(1-infRateDen)*(S_HchikS_H0den + (1-vaxChikEfficacy)*V_HchikS_H0den)- (infRateDen + recoverRateChik - infRateDen*recoverRateChik)*I_H1chikS_H0den - VaxPropS0*I_H1chikS_H0den- mu(1)*HRChik*I_H1chikS_H0den;  
        %dR_H0chikS_H0den/dt
        dPop(4)= recoverRateChik*(I_H0chikS_H0den+(1-pSeqChik)*I_H1chikS_H0den)*(1-infRateDen)+ pSeqRecovery*R_H1chikS_H0den*(1-infRateDen)- infRateDen*R_H0chikS_H0den - VaxPropS0*R_H0chikS_H0den-mu(1)*R_H0chikS_H0den;
        %dR_H1chikS_H0den/dt
        dPop(5)= recoverRateChik*pSeqChik*I_H1chikS_H0den*(1-infRateDen)-(pSeqRecovery + infRateDen - pSeqRecovery*infRateDen)*R_H1chikS_H0den- VaxPropS0*R_H1chikS_H0den-mu(1)*R_H1chikS_H0den;
        
        %% Susceptible Humans 1 Dengue
        %dS_HchikS_H1den/dt
        dPop(6) =recoverRateDen*(S_HchikI_H0den + S_HchikI_H1den)*(1-infRateChik) - (infRateChik + infRateDen2 - infRateChik*infRateDen2)*S_HchikS_H1den - VaxPropS1*S_HchikS_H1den - mu(1)*S_HchikS_H1den;  
        %dI_H0chikS_H1den/dt
        dPop(7) =recoverRateDen*(1-recoverRateChik)*(I_H0chikI_H0den + I_H0chikI_H1den) + infRateChik*(1-symptChik)*recoverRateDen*(S_HchikI_H0den + S_HchikI_H1den + (1-vaxChikEfficacy)*(V_HchikI_H0den+V_HchikI_H1den)) + infRateChik*(1-symptChik)*(S_HchikS_H1den+ (1-vaxChikEfficacy)*V_HchikS_H1den)*(1-infRateDen2) - (infRateDen2 +recoverRateChik - infRateDen2*recoverRateChik)*I_H0chikS_H1den - VaxPropS1*I_H0chikS_H1den - mu(1)*I_H0chikS_H1den; 
        %dI_H1chikS_H1den/dt
        dPop(8)= recoverRateDen*(1-recoverRateChik)*(I_H1chikI_H0den + I_H1chikI_H1den) + infRateChik*symptChik*recoverRateDen*(S_HchikI_H0den + S_HchikI_H1den+ (1-vaxChikEfficacy)*(V_HchikI_H0den+V_HchikI_H1den)) + infRateChik*symptChik*(S_HchikS_H1den+ (1-vaxChikEfficacy)*V_HchikS_H1den)*(1-infRateDen2) - (infRateDen2 + recoverRateChik -infRateDen2*recoverRateChik)*I_H1chikS_H1den - VaxPropS1*I_H1chikS_H1den - mu(1)*HRChik*I_H1chikS_H1den;  
        %dR_H0chikS_H1den/dt
        dPop(9)= recoverRateDen*(R_H0chikI_H0den + R_H0chikI_H1den) + recoverRateChik*(1-infRateDen2)*(I_H0chikS_H1den+(1-pSeqChik)*I_H1chikS_H1den)+ pSeqRecovery*R_H1chikS_H1den*(1-infRateDen2) + pSeqRecovery*recoverRateDen*(R_H1chikI_H0den + R_H1chikI_H1den) + recoverRateChik*recoverRateDen*(I_H0chikI_H0den + I_H0chikI_H1den + (1-pSeqChik)*(I_H1chikI_H0den+ I_H1chikI_H1den)) - infRateDen2*R_H0chikS_H1den - VaxPropS1*R_H0chikS_H1den -mu(1)*R_H0chikS_H1den;
        %dR_H1chikS_H1den/dt
        dPop(10)= recoverRateDen*(1-pSeqRecovery)*(R_H1chikI_H0den + R_H1chikI_H1den) + recoverRateChik*(1-infRateDen2)*pSeqChik*I_H1chikS_H1den + recoverRateChik*pSeqChik*recoverRateDen*(I_H1chikI_H0den + I_H1chikI_H1den) -(pSeqRecovery + infRateDen2 - pSeqRecovery*infRateDen2)*R_H1chikS_H1den- VaxPropS1*R_H1chikS_H1den -mu(1)*R_H1chikS_H1den;
       
        
        %% Infected Humans 0 Dengue
        %r*(T(1,2)*I_M + T(1,1)*I_H0)*(S_H0*(1-symptDen0)+ S_H1*(1-symptDen1)*l)
        
        %dS_HchikI_H0den/dt
        dPop(11) =  (1-infRateChik)*infRateDen*(1-symptDen0)*S_HchikS_H0den+ (1-infRateChik)*infRateDen2*S_HchikS_H1den*(1-symptDen1) + (VaxBounceBack)*(1-infRateChik)*infRateDen2*(1-symptDen1)*(1-vaxDenEfficacy)*S_HchikV_Hden -(recoverRateDen + infRateChik - recoverRateDen*infRateChik)*S_HchikI_H0den - mu(1)*S_HchikI_H0den;  
        %dI_H0chikI_H0den/dt
        dPop(12) = infRateChik*(1-symptChik)*(S_HchikI_H0den + (1-vaxChikEfficacy)*V_HchikI_H0den)*(1-recoverRateDen) + (1-recoverRateChik)*infRateDen*I_H0chikS_H0den*(1-symptDen0)+ (1-recoverRateChik)*infRateDen2*I_H0chikS_H1den*(1-symptDen1)+ (1-recoverRateChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*I_H0chikV_Hden*(1-symptDen1) + infRateChik*(1-symptChik)*infRateDen*(1-symptDen0)*(S_HchikS_H0den + (1-vaxChikEfficacy)*V_HchikS_H0den) +infRateChik*(1-symptChik)*(VaxBounceBack)*infRateDen2*(1-symptDen1)*((1-vaxDenEfficacy)*(1-vaxChikEfficacy)*V_HchikV_Hden+ (1-vaxDenEfficacy)*S_HchikV_Hden) + infRateChik*(1-symptChik)*infRateDen2*(1-symptDen1)*(S_HchikS_H1den+(1-vaxChikEfficacy)*V_HchikS_H1den) -(recoverRateDen + recoverRateChik - recoverRateDen*recoverRateChik)*I_H0chikI_H0den  - mu(1)*I_H0chikI_H0den; 
        %dI_H1chikI_H0den/dt
        dPop(13)=  infRateChik*symptChik*(S_HchikI_H0den+ (1-vaxChikEfficacy)*V_HchikI_H0den)*(1-recoverRateDen) + infRateChik*symptChik*infRateDen*(S_HchikS_H0den+ (1-vaxChikEfficacy)*V_HchikS_H0den)*(1-symptDen0)+ infRateChik*symptChik*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*S_HchikV_Hden*(1-symptDen1)+ infRateChik*symptChik*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*(1-vaxChikEfficacy)*V_HchikV_Hden*(1-symptDen1)+ infRateChik*symptChik*infRateDen2*(S_HchikS_H1den+ (1-vaxChikEfficacy)*V_HchikS_H1den)*(1-symptDen1) +(1-recoverRateChik)*infRateDen*I_H1chikS_H0den*(1-symptDen0)+(1-recoverRateChik)*(infRateDen2*(1-vaxDenEfficacy)*VaxBounceBack*I_H1chikV_Hden*(1-symptDen1) + infRateDen2*I_H1chikS_H1den*(1-symptDen1)) -(recoverRateDen + recoverRateChik - recoverRateDen*recoverRateChik)*I_H1chikI_H0den - mu(1)*HRChik*I_H1chikI_H0den;  
        %dR_H0chikI_H0den/dt
        dPop(14)=  recoverRateChik*(1-recoverRateDen)*(I_H0chikI_H0den+(1-pSeqChik)*I_H1chikI_H0den) + recoverRateChik*infRateDen*(1-symptDen0)*(I_H0chikS_H0den +(1-pSeqChik)*I_H1chikS_H0den)+ recoverRateChik*(VaxBounceBack)*infRateDen2*(1-symptDen1)*((1-vaxDenEfficacy)*I_H0chikV_Hden +(1-pSeqChik)*(1-vaxDenEfficacy)*I_H1chikV_Hden) + recoverRateChik*infRateDen2*(1-symptDen1)*(I_H0chikS_H1den+(1-pSeqChik)*I_H1chikS_H1den)+ pSeqRecovery*R_H1chikI_H0den*(1-recoverRateDen)+ infRateDen*R_H0chikS_H0den*(1-symptDen0)+ infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*R_H0chikV_Hden*(1-symptDen1)+ infRateDen2*R_H0chikS_H1den*(1-symptDen1)+ infRateDen*(1-symptDen0)*R_H1chikS_H0den*pSeqRecovery+ (VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*R_H1chikV_Hden*(1-symptDen1)*pSeqRecovery+ infRateDen2*R_H1chikS_H1den*pSeqRecovery*(1-symptDen1) -recoverRateDen*(R_H0chikI_H0den) -mu(1)*R_H0chikI_H0den;
        %dR_H1chikI_H0den/dt
        dPop(15)=  recoverRateChik*(pSeqChik*I_H1chikI_H0den)*(1-recoverRateDen) +(1-pSeqRecovery)*(infRateDen*R_H1chikS_H0den*(1-symptDen0)+ infRateDen2*R_H1chikS_H1den*(1-symptDen1)) +(1-pSeqRecovery)*(infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*R_H1chikV_Hden*(1-symptDen1)) + recoverRateChik*pSeqChik*(infRateDen*(1-symptDen0)*I_H1chikS_H0den + infRateDen2*(1-symptDen1)*I_H1chikS_H1den) + recoverRateChik*pSeqChik*infRateDen2*(VaxBounceBack)*(1-symptDen1)*(1-vaxDenEfficacy)*I_H1chikV_Hden-(recoverRateDen + pSeqRecovery - recoverRateDen*pSeqRecovery)*R_H1chikI_H0den - mu(1)*R_H1chikI_H0den;
        
        
        %% Infected Humans 1 Dengue
        %r*(T(1,2)*I_M + T(1,1)*I_H0)*(S_H0*symptDen0 + S_H1*symptDen1*lden)/Htot
        %dS_HchikI_H1den/dt
        dPop(16) =  (1-infRateChik)*infRateDen*S_HchikS_H0den*symptDen0 + (1-infRateChik)*infRateDen2*S_HchikS_H1den*symptDen1 + (1-infRateChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*S_HchikV_Hden*symptDen1 -(recoverRateDen + infRateChik -recoverRateDen*infRateChik)*S_HchikI_H1den  - mu(1)*HRDen*S_HchikI_H1den;  
        %dI_H0chikI_H1den/dt
        dPop(17) =  infRateChik*(1-symptChik)*(S_HchikI_H1den+ (1-vaxChikEfficacy)*V_HchikI_H1den)*(1-recoverRateDen) + infRateChik*(1-symptChik)*infRateDen*(S_HchikS_H0den + (1-vaxChikEfficacy)*V_HchikS_H0den)*symptDen0 + infRateChik*(1-symptChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*S_HchikV_Hden*symptDen1 + infRateChik*(1-symptChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*(1-vaxChikEfficacy)*V_HchikV_Hden*symptDen1 + infRateChik*(1-symptChik)*infRateDen2*(S_HchikS_H1den + (1-vaxChikEfficacy)*V_HchikS_H1den)*symptDen1 + (1- recoverRateChik)*infRateDen*I_H0chikS_H0den*symptDen0 + (1- recoverRateChik)*(infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*I_H0chikV_Hden*(symptDen1) + infRateDen2*I_H0chikS_H1den*symptDen1) -(recoverRateDen + recoverRateChik - recoverRateDen*recoverRateChik)*I_H0chikI_H1den - mu(1)*HRDen*I_H0chikI_H1den; 
        %dI_H1chikI_H1den/dt
        dPop(18)=   infRateChik*symptChik*S_HchikI_H1den*(1-recoverRateDen) + infRateChik*symptChik*(1-vaxChikEfficacy)*V_HchikI_H1den*(1-recoverRateDen) + infRateChik*(symptChik)*(VaxBounceBack)*infRateDen2*(1-vaxDenEfficacy)*S_HchikV_Hden*symptDen1 + infRateChik*symptChik*(infRateDen*S_HchikS_H0den*symptDen0  + infRateDen*(1-vaxChikEfficacy)*V_HchikS_H0den*symptDen0+ infRateDen2*S_HchikS_H1den*symptDen1 + infRateDen2*(1-vaxChikEfficacy)*V_HchikS_H1den*symptDen1 + infRateDen2*(VaxBounceBack)*(1-vaxChikEfficacy)*(1-vaxDenEfficacy)*V_HchikV_Hden*symptDen1) + (1- recoverRateChik)*(infRateDen*I_H1chikS_H0den*symptDen0 + infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*I_H1chikV_Hden*symptDen1 + infRateDen2*I_H1chikS_H1den*symptDen1)  -(recoverRateDen + recoverRateChik- recoverRateDen*recoverRateChik)*I_H1chikI_H1den - mu(1)*HRDen*HRChik*I_H1chikI_H1den;  
        %dR_H0chikI_H1den/dt
        dPop(19)=   recoverRateChik*(1-recoverRateDen)*(I_H0chikI_H1den+(1-pSeqChik)*I_H1chikI_H1den) + infRateDen*R_H0chikS_H0den*symptDen0 + infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*R_H0chikV_Hden*symptDen1 + infRateDen2*R_H0chikS_H1den*symptDen1 + pSeqRecovery*R_H1chikI_H1den*(1-recoverRateDen) + pSeqRecovery*(infRateDen*symptDen0*R_H1chikS_H0den + (VaxBounceBack)*infRateDen2*symptDen1*(1-vaxDenEfficacy)*R_H1chikV_Hden +infRateDen2*symptDen1*R_H1chikS_H1den)   + recoverRateChik*(symptDen0*infRateDen*(I_H0chikS_H0den+(1-pSeqChik)*I_H1chikS_H0den)+ symptDen1*(VaxBounceBack)*infRateDen2*((1-vaxDenEfficacy)*I_H0chikV_Hden+(1-pSeqChik)*(1-vaxDenEfficacy)*I_H1chikV_Hden)  + infRateDen2*symptDen1*(I_H0chikS_H1den+(1-pSeqChik)*I_H1chikS_H1den))   -recoverRateDen*(R_H0chikI_H1den) -mu(1)*HRDen*R_H0chikI_H1den;
        %dR_H1chikI_H1den/dt
        dPop(20)=   recoverRateChik*pSeqChik*I_H1chikI_H1den*(1-recoverRateDen) + (1-pSeqRecovery)*(infRateDen*R_H1chikS_H0den*symptDen0 + infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*R_H1chikV_Hden*symptDen1 + infRateDen2*R_H1chikS_H1den*symptDen1) + recoverRateChik*pSeqChik*(infRateDen*I_H1chikS_H0den*symptDen0 + infRateDen2*(VaxBounceBack)*(1-vaxDenEfficacy)*I_H1chikV_Hden*symptDen1 + infRateDen2*I_H1chikS_H1den*symptDen1) -(recoverRateDen + pSeqRecovery - recoverRateDen*pSeqRecovery)*R_H1chikI_H1den -mu(1)*HRDen*R_H1chikI_H1den;
        
        %% Vaccinated Humans Dengue
        %dS_HchikV_Hden/dt
        dPop(30) = VaxPropS0*S_HchikS_H0den + VaxPropS1*S_HchikS_H1den -(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + infRateChik - VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*infRateChik)*S_HchikV_Hden  - mu(1)*S_HchikV_Hden;  
        %dI_H0chikV_Hden/dt
        dPop(31) = VaxPropS0*I_H0chikS_H0den + VaxPropS1*I_H0chikS_H1den + (1-vaxChikEfficacy)*infRateChik*(1-symptChik)*(1-(infRateDen2*(1-vaxDenEfficacy)))*V_HchikV_Hden + infRateChik*(1-symptChik)*(1-(infRateDen2*(1-vaxDenEfficacy)))*S_HchikV_Hden-(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + recoverRateChik - VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*recoverRateChik)*I_H0chikV_Hden - mu(1)*I_H0chikV_Hden; 
        %%dI_H1chikV_Hden/dt
        dPop(32) = VaxPropS0*I_H1chikS_H0den + VaxPropS1*I_H1chikS_H1den + (1-vaxChikEfficacy)*infRateChik*(symptChik)*(1-(infRateDen2*(1-vaxDenEfficacy)))*V_HchikV_Hden +infRateChik*symptChik*(1-(infRateDen2*(1-vaxDenEfficacy)))*S_HchikV_Hden-(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + recoverRateChik - VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*recoverRateChik)*I_H1chikV_Hden- mu(1)*HRChik*I_H1chikV_Hden;  
        %dR_H0chikV_Hden/dt
        dPop(33) = VaxPropS0*R_H0chikS_H0den + VaxPropS1*R_H0chikS_H1den + recoverRateChik*(1-(infRateDen2*(1-vaxDenEfficacy)))*I_H0chikV_Hden + recoverRateChik*(1-pSeqChik)*(1-(infRateDen2*(1-vaxDenEfficacy)))*I_H1chikV_Hden + pSeqRecovery*(1-(infRateDen2*(1-vaxDenEfficacy)))*R_H1chikV_Hden -VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*R_H0chikV_Hden -mu(1)*R_H0chikV_Hden;
        %dR_H1chikV_Hden/dt
        dPop(34) = VaxPropS0*R_H1chikS_H0den + VaxPropS1*R_H1chikS_H1den + recoverRateChik*pSeqChik*(1-(infRateDen2*(1-vaxDenEfficacy)))*I_H1chikV_Hden-(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + pSeqRecovery - VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*pSeqRecovery)*R_H1chikV_Hden -mu(1)*R_H1chikV_Hden;
        
        %% Vaccinated Humans Chikungunya
        %dV_Hchik S_H0den/dt 
        dPop(35) =  -((1-vaxChikEfficacy)*infRateChik + infRateDen - (1-vaxChikEfficacy)*infRateChik*infRateDen)*V_HchikS_H0den -VaxPropS0*V_HchikS_H0den - mu(1)*V_HchikS_H0den;  
        %dV_HchikS_H1den/dt
        dPop(36) = (1-(infRateChik*(1-vaxChikEfficacy)))*recoverRateDen*V_HchikI_H0den + (1-(infRateChik*(1-vaxChikEfficacy)))*recoverRateDen*V_HchikI_H1den- ((1-vaxChikEfficacy)*infRateChik + infRateDen2 - ((1-vaxChikEfficacy)*infRateChik*infRateDen2))*V_HchikS_H1den -VaxPropS1*V_HchikS_H1den - mu(1)*V_HchikS_H1den;  
        %dV_HchikI_H0den/dt
        dPop(37) =  (1-(infRateChik*(1-vaxChikEfficacy)))*VaxBounceBack*infRateDen2*(1-symptDen1)*(1-vaxDenEfficacy)*V_HchikV_Hden + (1-(infRateChik*(1-vaxChikEfficacy)))*infRateDen*(1-symptDen0)*V_HchikS_H0den + (1-(infRateChik*(1-vaxChikEfficacy)))*infRateDen2*(1-symptDen1)*V_HchikS_H1den-(recoverRateDen + (1-vaxChikEfficacy)*infRateChik - (recoverRateDen*(1-vaxChikEfficacy)*infRateChik))*V_HchikI_H0den - mu(1)*V_HchikI_H0den;  
        %dV_HchikI_H1den/dt
        dPop(38) =  (1-(infRateChik*(1-vaxChikEfficacy)))*VaxBounceBack*infRateDen2*(symptDen1)*(1-vaxDenEfficacy)*V_HchikV_Hden + (1-(infRateChik*(1-vaxChikEfficacy)))*infRateDen*symptDen0*V_HchikS_H0den + (1-(infRateChik*(1-vaxChikEfficacy)))*infRateDen2*symptDen1*V_HchikS_H1den -(recoverRateDen + (1-vaxChikEfficacy)*infRateChik -(recoverRateDen*(1-vaxChikEfficacy)*infRateChik))*V_HchikI_H1den  - mu(1)*HRDen*V_HchikI_H1den;  
        %dV_HchikV_Hden/dt
        dPop(39) = VaxPropS0*V_HchikS_H0den + VaxPropS1*V_HchikS_H1den -(VaxBounceBack*infRateDen2*(1-vaxDenEfficacy) + (1-vaxChikEfficacy)*infRateChik -VaxBounceBack*infRateDen2*(1-vaxDenEfficacy)*(1-vaxChikEfficacy)*infRateChik)*V_HchikV_Hden  - mu(1)*V_HchikV_Hden;  

        
        %% Mosquitos
        infChikMos = r*(Tchik(2,1)*(I_H0chik +I_H1chik))/Mtot;
        infDenMos = r*(Tden(2,2)*I_Mden + Tden(2,1)*(I_H0den +I_H1den))/Mtot;
        
        %dS_MchikS_Mden/dt
        dPop(21)=nu(2)*TotalMosPop - (infChikMos + infDenMos -infChikMos*infDenMos)*S_MchikS_Mden - mu(2)*S_MchikS_Mden;
        %dI_MchikS_Mden/dt
        dPop(22)= infChikMos*(1-infDenMos)*S_MchikS_Mden-infDenMos*I_MchikS_Mden- mu(2)*I_MchikS_Mden;
        %dS_MchikI_Mden/dt
        dPop(23)= infDenMos*(1-infChikMos)*S_MchikS_Mden-infChikMos*S_MchikI_Mden- mu(2)*S_MchikI_Mden;
        %dI_MchikI_Mden/dt
        dPop(24)= infDenMos*I_MchikS_Mden + infChikMos*S_MchikI_Mden + infChikMos*infDenMos*S_MchikS_Mden- mu(2)*I_MchikI_Mden;
        
        
        %% Extra data
        %deathsTot
         dPop(25) = mu(1)*(S_HchikS_H0den + I_H0chikS_H0den  + I_H1chikS_H0den*HRChik  + R_H0chikS_H0den  + R_H1chikS_H0den  + S_HchikS_H1den + I_H0chikS_H1den + I_H1chikS_H1den*HRChik + R_H0chikS_H1den  + R_H1chikS_H1den  + S_HchikI_H0den  + I_H0chikI_H0den  + I_H1chikI_H0den*HRChik  + R_H0chikI_H0den + R_H1chikI_H0den  + HRDen*(S_HchikI_H1den  + I_H0chikI_H1den  + I_H1chikI_H1den*HRChik  + R_H0chikI_H1den  + R_H1chikI_H1den)+S_HchikV_Hden + I_H0chikV_Hden + HRChik*I_H1chikV_Hden + R_H0chikV_Hden + R_H1chikV_Hden +V_HchikS_H0den +  V_HchikS_H1den + V_HchikI_H0den + HRDen*V_HchikI_H1den + V_HchikV_Hden);
        % deathsChik
         dPop(26) = mu(1)*(HRChik-1)*(I_H1chikS_H0den  + I_H1chikS_H1den + I_H1chikI_H0den + I_H1chikI_H1den + I_H1chikV_Hden);
        % deathsDen
         dPop(27) = mu(1)*(HRDen-1)*(S_HchikI_H1den + I_H0chikI_H1den  + I_H1chikI_H1den  + R_H0chikI_H1den  + R_H1chikI_H1den + V_HchikI_H1den);
        %cumulativeInfectionChik
         dPop(28) = (infRateChik*S_Hchik*symptChik)+ V_Hchik*infRateChik*symptChik*(1-vaxChikEfficacy) + V_Hboth*infRateChik*symptChik*(1-vaxChikEfficacy);
        %cumulativeInfectionDen
         dPop(29) = infRateDen*S_H0den*symptDen0+infRateDen2*S_H1den*symptDen1+ V_Hden*infRateDen2*symptDen1*(1-vaxDenEfficacy)+ V_Hboth*infRateDen2*symptDen1*(1-vaxDenEfficacy);
      
        % Number of people tested
        totalTests = (vaxRateWeek*Htot);
        dPop(40)=totalTests;
        % Number of people vaccinated
        totalVax = VaxPropS0*S_HchikS_H0den + VaxPropS1*S_HchikS_H1den +VaxPropS0*I_H0chikS_H0den + VaxPropS1*I_H0chikS_H1den + VaxPropS0*I_H1chikS_H0den + VaxPropS1*I_H1chikS_H1den +VaxPropS0*R_H0chikS_H0den + VaxPropS1*R_H0chikS_H1den +VaxPropS0*R_H1chikS_H0den + VaxPropS1*R_H1chikS_H1den +  VaxPropS0*V_HchikS_H0den + VaxPropS1*V_HchikS_H1den;
        dPop(41)=totalVax;
    
    end


%Calculates the weekly incidence rates rate per 100,000 of population from
%cumulative incidence rate
    function incidence = IncidenceChik(pop,colombiaPop)
        incidence = NaN(MaxTime ,1);
        initialPop = sum(pop(1,1:20))+sum(pop(1,30:39));%+1;
        incidence(1) = (pop(1,28))*100000/initialPop;%colombiaPop(1);%(pop(1,3)+pop(1,8)+pop(1,13)+pop(1,18))*100000/initialPop;
        for i = 2:(MaxTime)
            totalPop = sum(pop(i,1:20))+sum(pop(i,30:39));%+1;
            incidence(i) = (pop(i,28)-pop(i-1,28))*100000/totalPop;%colombiaPop(i);%((pop(i,3)+pop(i,8)+pop(i,13)+pop(i,18))-(pop(i-1,3)+pop(i-1,8)+pop(i-1,13)+pop(i-1,18)))*100000/totalPop;
        end
    end

%Calculates the weekly incidence rates rate per 100,000 of population for
%Den
    function incidence = IncidenceRatesDen(pop,colombiaPop)
        incidence = NaN(MaxTime ,1);
        initialPop = sum(pop(1,1:20))+sum(pop(1,30:39));%+1; %Population of modelled area used
        incidence(1) = (pop(1,29)*100000)/initialPop;%colombiaPop(1);%sum(pop(1,16:20))*100000/initialPop;
        for i = 2:(MaxTime)
            totalPop = sum(pop(i,1:20))+sum(pop(i,30:39));%+1;
            incidence(i) = ((pop(i,29)-pop(i-1,29))*100000)/totalPop;%colombiaPop(i);%(sum(pop(i,16:20))-sum(pop(i-1,16:20)))*100000/totalPop;
        end
    end



end
