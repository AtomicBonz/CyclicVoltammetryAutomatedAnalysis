% v.4 - single Scan - 
%% README
% in order to analyze Cyclic Voltammetry data with this code, the
% VersaStudio ".par" file must first be converted directly into a ".txt" format
% this can be done by opening ".par" using Notepad and then "save as: .txt"
% for batches of ".par" files, it can be done by:
%i) open CMD
%ii) navigate to parent directory of where your ".par" files are stored, where you intend for recursive operation to happen 
%iib) make sure a backup copy of ".par" files are stored somewhere else
%iii) >> for /r %g in (*.par) do ren "%g" "%~ng.txt"
%iv) your ".par" files will turn into ".txt" files which are ready for
%inputinto CVAA4

%this code has 2 outputs:
%1) a data ".txt" file with cyclic voltammetry parameters of interest: peak
%current value, apparent reduction potential, no. of electrons exchanged,
%analyte concentration, analyte diffusion coefficient, charging current
%slope value, R-squared (Rsq) value of the Berzins-Delahay fit. 
%2) a plot of the fit, saved as a ".fig" file or ".png" file. 

% the following distinctions are made:
% "experimenter" is the person who physically performed the CV, so
% parameters like Vmax_set, Vmin_set and v, A and T were determined by them. 
% "analyst" is the person who's running the CVAA 4 code, so they can set parameters like
% Vmax_clip, Vmin_clip, but they are responsible for inputting the correct,
% v, A and T set by experimenter into this code. Analyst is also
% responsible for choosing the correct n, alpha and initialized C_0 OR Do. 

%line 563, 564: switch on/off which param C_0 or Do you desire to search
%for, where the other is defined in line 94,95
%%

fclose all;

clear all
close all

CCstartWidScale = 1;
CCendWidScale = 1;

codeVers = 'CVAA4_redOnly--'; %version of the code, appended to output name, "redOnly" means only reduction side is analyzed
folder = 'C:\Users\bonit\MATLAB Drive\CV_Analysis\2020_11_10_miniProbeCorrExpt\CVAAout\'; %folder where your output files are saved into

%%input data
tic
fid = fopen('BatchConcFix.txt'); % text file that contains a list of all the files in the input folder you wish to analyze. this line is legacy. this means not all files in the folder will get analyzed, only the ones in this list. filenames in this list must include input folder directory.
while ~feof(fid) % feof(fid) is true when the file ends 
%% open each file and pull in the data 
filename = fgetl(fid); % file from list to be analyzed
datain = importdata(filename);%importing data from txt into matrix
%% defining the variables in which to pull data. The data in these variables will not be processed. 
% these variables are retained for convenient debugging (compare
% pre-processed data with post-processed data)
filename_indiv = extractAfter(filename,'CVAAin\');%extracts filename from list without input folder directory
dataname = strtok(filename_indiv,'.'); %creating a string excluding '.txt' to name the fig
  
%% ===========EXTRACTED DATA FROM EXPERIMENT OUTPUT.txt FILE specific to VS.txt format
WE1PotentialV = datain.data(:,5);%col of data for working electrode potential in [V]
WE1PotentialV(1:5) = [];% chopping off front matter
WE1CurrentA = datain.data(:,6);%col of data for working electrode current in [A]
WE1CurrentA(1:5) = [];% chopping off front matter
Vmax_set = datain.data(1,1);% upper limit of potential window in [V], set by experimenter
Vmin_set = datain.data(1,2);% lower limit of potential window in[V], set by experimenter
% OCPvalue = datain.data(1,); %VS data has no OCP value output. this is legacy from Gamry.  
StartScanPotential = datain.data(1,10); %VS data has no start scan potential value. thisis legacy from Gamry. 
V_sweepRate = datain.data(1,3); %potential scan rate, set by experimenter.
NoOfStopCross = datain.data(1,4);%number of times the potential scan switches direction  
t = datain.data(:,7); %timestamp in [s] for each datapoint of current & potential 
t(1:5) = [];% chopping off front matter
ScanIndex = datain.data(:,11);%VS calls this "Segment" but we have to +1 to make the index 1-based not 0-based 
ScanIndex(1:5) = [];% chopping off front matter
DataIndex = datain.data(:,12);% they call this "point" but we have to +1 to make the index 1-based not 0-based 
DataIndex(1:5) = [];% chopping off front matter
%% ========== HARD CODED PARAMS
%function tolerances
tol_dV = 0.001 ;
ScanDir = 0;% direction of potential scan (-1 for cathodic, +1 for anodic, 0 for agnostic upon initialization)

% numerical (aka brute force) params
int_limit = 21; %legacy, for earlier version of the code using least-squares minimization 
% physical fixed constants
F = 96500;%Faraday's constant in [C/mol]
R = 8.31;%ideal gas const in [J/molK]
v = V_sweepRate; %CV sweep rate [V/s], this is extracted from your data file
T = 273+500; % absolute temperature [K], input by analyst					   
A = 0.47885; %_cmsq electrode surface area [cm^2], 
WEType = 'W'; % corrosion time in seconds;
Vmin_clip = -1.5; %region of interest to analyze and plot. 
Vmax_clip = 0.35; %region of interest to analyze and plot.

%other experimental params
n = 2; %number of electrons in reaction couple of interest
alpha = 0.5; % electron transfer symmetry coefficient
%C_0 = 3.8E-5;% analyte concentration in [mol/cm^3]
Do = 1.13E-4; % analyte diffusion coefficient in [cm^2/s]
%the analysis will solve for either C_0 or Do, but the other must be known.
titl = 't = 4.3 h'; %plot title text. 

%after specifying up to this point, the analyst finally need only specify in line
%XXX which of C_0 or Do they would like to solve for, then just run the code.  
% "time jump": some choices of Vmin_clip & Vmax_clip can cause there to be discontinuities in the Time vector
% which will cause an analysis error and plotting error.
% if you get a plotting error where "Etot and Itot are not the
% same length", check the values of the Time vector and ensure that the
% timestamps run consecutively in intervals of dt[s] (usually about 0.1s)
% with no large gaps in them aka "time jump" e.g. if the timestamp skips
% from 17s immediately to 62s. if this happens you have to manually chop
% off the jump. "Ctr+G" to line 247 & 447 to manually input the value of
% the vector element to truncate at. Many errors can be traced back to an
% inappropriate choice of Vmin_clip/Vmax_clip.
%% initializing TotalData matrix where all data will be stored before/after operation

TotalData = [DataIndex ScanIndex t WE1PotentialV WE1CurrentA];
TotalData(5964:end,:) = []; % cut off NaN values. 
[TotalDataPtCount,~] = size(TotalData);

%% OTHER VARIABLES these are the variables that will be processed
V = smooth(TotalData(:,4));%put them into variables, don't act on them directly.
I = smooth(TotalData(:,5));
t = TotalData(:,3) - TotalData(1,3);
time = t;
ScanNo = TotalData(:,2);
dat_index = TotalData(:,1);

%% grad module: calculates local gradients dV, dt etc. 
Grad = zeros(TotalDataPtCount,4);%the reason i have them all in one matrix is so that it's easier to visually compare
	for i = 1:TotalDataPtCount-1
	Grad(i,1) = I(i+1)-I(i);
	Grad(i,2) = V(i+1)- V(i);
	Grad(i,3) = (I(i+1) - I(i))/(V(i+1) - V(i));% dI/dE, to convert to scan rate
	Grad(i,4) = (V(i+1)- V(i))/(t(i+1) - t(i)); %dV/dt, to verify that scan rate is the same as what the expeimenter input
	dt(i) = (time(i+1) - time(i));
	end
dt = mean(min(dt));
%% adding Grad data to Total Data
TotalData = [TotalData Grad];
% TotalData = [(1)DataIndex (2)ScanIndex (3)time (4)WE1PotentialV (5)WE1CurrentA (6)dI (7)dV (8)dI/dE (9)dV/dt];
%% Sawtooth fit
Sweep = ones(TotalDataPtCount,1); %the sweep vector is the length of all the data pts 
	if (TotalData(1,7)> 0) && (TotalData(2,7)> 0) && (TotalData(3,7)> 0) && (TotalData(4,7)> 0) && (TotalData(5,7)> 0) 
	%text(0,0,'scan started in positive direction ')%TotalData(x,7) = (7)dV
	ScanDir = 1;
	InitSweepDirxn = 1;
	elseif (TotalData(1,7)< 0) && (TotalData(2,7)< 0) && (TotalData(3,7)< 0) && (TotalData(4,7)< 0) && (TotalData(5,7)< 0) 
	%text(0,0,' scan started in negative direction ')
	ScanDir = -1;
	InitSweepDirxn = -1;
	else
	fprintf('unstable initialization of scan. please clip TotalData before proceeding ')
	end

%% Exclusion based on fit to Sawtooth
sawtooth_amplitude = 0.5*(Vmax_set-Vmin_set);
sawtooth_period = 2*(Vmax_set-Vmin_set)/V_sweepRate;

	switch ScanDir
	case 1
	sawtooth_shift = abs(V(1)-Vmin_set)/V_sweepRate;
	case -1
	sawtooth_shift = abs(V(1)-Vmax_set)/V_sweepRate;
	end
sawModel = 0.5*(Vmax_set + Vmin_set ) + (ScanDir)*sawtooth_amplitude*sawtooth(2*pi*(t+sawtooth_shift)/sawtooth_period,0.5);

%plot(sawModel,'r.')
devFrmSweep = sawModel - V;
TotalData = [TotalData devFrmSweep];
% TotalData = [(1)DataIndex (2)ScanIndex (3)time (4)WE1PotentialV (5)WE1CurrentA (6)dI (7)dV (8)dI/dE (9)dV/dt (10)devFrmSweep];

%% Assigning Sweep coordinate
[VmaxTurningPtValues, VmaxTurningPtValue_indices] = findpeaks(sawModel,'MinPeakHeight',Vmax_set - V_sweepRate);%,'MinPeakDistance',);% what is the instrument LOQ? 
[VminTurningPtValues, VminTurningPtValue_indices] = findpeaks(-sawModel,'MinPeakHeight',abs(Vmin_set + V_sweepRate));%,'MinPeakDistance',);% what is the instrument LOQ? 
VminTurningPtValues = -VminTurningPtValues;%the negative sign has to be put back

[noOfMaxTurningPts,~] = size(VmaxTurningPtValues);
%noOfMaxTurningPts = noOfMaxTurningPts +1; %fix for depper22
[noOfMinTurningPts,~] = size(VminTurningPtValues);

AllTurningPts = zeros(1,2);%preallocated for speed
	switch ScanDir
	case 1
		for p = 1:noOfMaxTurningPts %this came from "findpeaks" which you don't have anymore 
		AllTurningPts = [AllTurningPts; VmaxTurningPtValue_indices(p) VmaxTurningPtValues(p)];
		AllTurningPts = [AllTurningPts; VminTurningPtValue_indices(p) VminTurningPtValues(p)];
		end
	case  -1
		for p = 1:noOfMinTurningPts
		AllTurningPts = [AllTurningPts; VminTurningPtValue_indices(p) VminTurningPtValues(p)];
            %AllTurningPts = [AllTurningPts; VmaxTurningPtValue_indices(p) VmaxTurningPtValues(p)]
            try AllTurningPts = [AllTurningPts; VmaxTurningPtValue_indices(p) VmaxTurningPtValues(p)];
            catch AllTurningPts = [AllTurningPts; 5962 TotalData(5962,4)];
            end %this try catch sequence is for stuff before depper22
		end
	end
%AllTurningPts(AllTurningPts(:,1)== 0,:) = []; 

%scatter(AllTurningPts(:,1),AllTurningPts(:,2),'ms')

switch InitSweepDirxn
  case -1
	for q = 2:NoOfStopCross
		for p = AllTurningPts(q,1):AllTurningPts(q+1,1)
            p;
            Sweep(p) = ScanDir;
		end
	ScanDir = -ScanDir;
    end
    Sweep = Sweep*(-ScanDir);%i don't know why we have to put another (-1) factor in here. i don't know when the direction got reversed 
    case 1
        for q = 1:NoOfStopCross
            for p = AllTurningPts(q,1):AllTurningPts(q+1,1)
                p = p+1;
                Sweep(p) = ScanDir;
            end
            ScanDir = -ScanDir;
        end
    
    end 
        
%plot(Vmax_set*Sweep,'g.')    
TotalData = [TotalData Sweep];% this works
% TotalData = [(1)DataIndex (2)ScanIndex (3)time (4)WE1PotentialV (5)WE1CurrentA (6)dI (7)dV (8)dI/dE (9)dV/dt (10)devFrmSweep (11)Sweep];
%% Sawtooth Exclusion condition
TotalData_V = TotalData; %save the old data for later plotting
TotalData_VI = TotalData;
%% Noise reduction & potential overload exclusion step
TotalData_VI(abs(TotalData_VI(:,10))>0.05*sawtooth_amplitude,:) = []; %exclusion condition: if devSweep is some absolute value larger than modeled sawtooth amplitude, exclude

noOfScans = max(ScanNo);

				%% Plot setup
				for x = 1:noOfScans
				%x = 2;
				%fprintf(strcat('\n this is scan number ',num2str(x)))
				t = 1;
				Ecoupl = [];
                tiler = figure;%tiledlayout(1,noOfScans);
				title(titl) 
				%nexttile 
				xlabel('Potential [V]');
				ylabel('Current [A]');
				%% Analyzing both REDOX scans to find Ered 
					for y = -1:2:1 %we need to check both scans to find Eguess accurately 
					%fprintf(strcat('\n this is scan direction ',num2str(y)))
					% 
					hold on 
					%     
					initialScan = TotalData;
					initialScan(initialScan(:,2)~=x,:) = []; 
					initialScan(initialScan(:,11)~=y,:) = [];
					
					%% if you want to analyze a sub-window of your potential scan window
					initialScan(initialScan(:,4)< Vmin_clip,:) = [];
					initialScan(initialScan(:,4)> Vmax_clip,:) = [];
                    %% get rid of time jump
                    initialScan(564:end,:) = [];
                    %%
					
					origData = [initialScan(:,4) initialScan(:,5)];
					Time = initialScan(:,3); %this is necessary for BD functions, not so for MO
					Current= initialScan(:,5);
					Potentials = initialScan(:,4);                                         

					[ptsPerPass,~] = size(Time);
					h1 = plot(Potentials,Current,'--b','LineWidth',2);
                    [peaks_expData, loc_expData, width_expData, prom_expData] = findpeaks(smooth(y*Current),'MinPeakProminence',0.0001,'MaxPeakWidth',300,'Sortstr','descend','Npeaks',1,'Annotate','extents')%play with params here, no specifications work decently but could be optimized
                        try isTherePeak = peaks_expData(1);
						% do this if there are peaks found
						isCleanSalt = 0;
						% =========end pass selection
							switch y
							case -1
							isCleanSalt_red = isCleanSalt;
							peaks_red = peaks_expData;
							peaklocs_red = loc_expData;
							[NoOfPeaks,~] = size(peaks_red);			 																			 
							peakRiseBegin = round(loc_expData-(0.5+0.5)*width_expData)
							peakFallEnd = round(loc_expData+(1+0.5)*width_expData);
							risingLegStart_red = Potentials(peakRiseBegin);
								if peakFallEnd > ptsPerPass
								peakFallEnd = ptsPerPass;
								end

							case 1
							isCleanSalt_ox = isCleanSalt;
							peaks_ox = peaks_expData;
							peaklocs_ox = loc_expData;
							%writecell()
							[NoOfPeaks,~] = size(peaks_ox);			 																													  																	  
							peakRiseBegin = round(loc_expData-(1+0.6)*width_expData);
							peakFallEnd = round(loc_expData+width_expData);
							risingLegStart_ox = Potentials(peakRiseBegin);
								if peakFallEnd > ptsPerPass
								peakFallEnd = ptsPerPass;
								end	
							end

						plot(Potentials(loc_expData(1)),y*peaks_expData(1),'rx','MarkerSize',16)%

						%begin charge current extrapolation																											  
                        switch y
                            case -1
                            [value,location] = min(abs(Potentials-0.2)); % we assume  Eu2+/Eu3+ hump has subsided by 0.2V. 
                            chargeCurrentStart = location;
                            case 1
                            chargeCurrentStart = 1;
                         end
						chargeCurrentEnd = peakRiseBegin
						chargeCurrent_E = Potentials(chargeCurrentStart:chargeCurrentEnd);
						chargeCurrent_i = Current(chargeCurrentStart:chargeCurrentEnd);

							switch y 
							case -1 
							chargeCurrent_i_red = chargeCurrent_i;
							chargeCurrent_E_red = chargeCurrent_E;
							%integralUnderPeak_red = integralUnderPeak;
							peakRiseBegin_red = peakRiseBegin;
							peakFallEnd_red = peakFallEnd;
							peakSelection_Potentials_red_exp = Potentials(peakRiseBegin:peakFallEnd);
							peakSelection_Current_red_exp = Current(peakRiseBegin:peakFallEnd);
                            e_pChargeCurrent18 = polyfit(chargeCurrent_E,chargeCurrent_i,1);
                            
                            chargeCurrent_i_fit_red = e_pChargeCurrent18(1)*Potentials + e_pChargeCurrent18(2); %extrapolate																				
                            chargeCurrentSlope_red = e_pChargeCurrent18(1);
                            h2 = plot(Potentials,chargeCurrent_i_fit_red,':r','LineWidth',1)
                            %h3 = plot(chargeCurrent_E, chargeCurrent_i,'gd','MarkerSize',10);    
                            chargeCurrentRegion = -1.3*peaks_red*ones(size(chargeCurrent_i_red));
                            h3 = area(chargeCurrent_E,chargeCurrentRegion,'EdgeColor','none','FaceColor','black','FaceAlpha',0.45,'ShowBaseLine','off')
                            %h4 = plot(Potentials(peakRiseBegin:peakFallEnd), Current(peakRiseBegin:peakFallEnd),'m+','MarkerSize',7);
                            peakRegion = -1.3*peaks_red*ones(size(Current(peakRiseBegin:peakFallEnd)));
                            h4 = area(Potentials(peakRiseBegin:peakFallEnd),peakRegion,'EdgeColor','none','FaceColor','black','FaceAlpha',0.15,'ShowBaseLine','off');
							case 1
							chargeCurrent_i_ox = chargeCurrent_i;
							chargeCurrent_E_ox = chargeCurrent_E;
							%integralUnderPeak_ox = integralUnderPeak;
							peakRiseBegin_ox = peakRiseBegin;
							peakFallEnd_ox = peakFallEnd;
							peakSelection_Potentials_red_ox = Potentials(peakRiseBegin:peakFallEnd);
							peakSelection_Current_red_ox = Current(peakRiseBegin:peakFallEnd);
                            e_pChargeCurrent18 = polyfit(chargeCurrent_E,chargeCurrent_i,1);
						    chargeCurrent_i_fit_ox = e_pChargeCurrent18(1)*Potentials + e_pChargeCurrent18(2); %extrapolate																				
                            chargeCurrentSlope_ox = e_pChargeCurrent18(1);
                            %plot(Potentials,chargeCurrent_i_fit_ox)
							end
                                                        
                                                                                             
						%plot(Potentials,chargeCurrent_i_fit)
						%text(min(Potentials),y*peaks_expData,strcat('Integral under pink - bgChargeCurrent = ', num2str(integralUnderPeak)))
						[noOfPeaksDetectedPerScan, ~]= size(peaks_expData);
						%peaks = y*peaks
						%Ecoupl is our storage matrix.
						Ecoupl(t) = initialScan(loc_expData(1),4);%"loc(1) hard coded because we're only interested in biggest peak for now, ie "most prominent peak". 
						widCoupl(t) = width_expData(1); %width(1) hard coded because we're only interested in biggest peak for now, ie "most prominent peak"
						% eventually fix this hard coded bit, "loc(1)". 
                        t = t+1;

						% do this if there are no peaks found 
						catch 
						isCleanSalt = 1;
						%% write a bit of code to calculate the clean salt charging current. 
						chargeCurrentStart =round(ptsPerPass*0.5)-30;
						chargeCurrentEnd = round(ptsPerPass*0.5)+30;
						chargeCurrent_E = Potentials(chargeCurrentStart:chargeCurrentEnd);
						chargeCurrent_i = Current(chargeCurrentStart:chargeCurrentEnd);
						scatter(chargeCurrent_E,chargeCurrent_i,'g.')
						e_pChargeCurrent18 = polyfit(chargeCurrent_E,chargeCurrent_i,1);
						chargeCurrent_i_fit = e_pChargeCurrent18(1)*Potentials + e_pChargeCurrent18(2); %extrapolate
						plot(Potentials,chargeCurrent_i_fit)
							switch y
							case -1
							isCleanSalt_red = isCleanSalt;
							ccStart_red = chargeCurrentStart;
							ccEnd_red = chargeCurrentEnd;
							chargeCurrentSlope_red = e_pChargeCurrent18(1);
							text(-0.9,y*2E-4,strcat('Charge Current slope, red =', num2str(e_pChargeCurrent18(1))))
							case 1
							isCleanSalt_ox = isCleanSalt;
							ccStart_ox = chargeCurrentStart;
							ccEnd_ox = chargeCurrentEnd;
							chargeCurrentSlope_ox = e_pChargeCurrent18(1);
							text(-0.9,y*2E-4,strcat('Charge Current slope, ox =', num2str(e_pChargeCurrent18(1))))
							end

						continue 
						end 							


                    end
                    
                    if (isCleanSalt_red==0 && isCleanSalt_ox==0)

                        electronConvConst = 2.3*R*T/F;
                        e_Ered = mean(Ecoupl);
                    
                    end
                    %% ==== Nernstian Shift
                    
                    %% ====End Nernstian Shift                           
                    


				dV= mean(abs(Grad(:,2)));% voltage step

				%% ==============START FITTING MODULE HERE through to line 439==============
				%% ==========peak segmentation module
				% =======pass selection 
				%tic

				r = 1;


					if isCleanSalt_red == 1 && isCleanSalt_ox ==1

					text(-0.2,0,'Clean Salt')

					dataout = [folder strcat('DATOUT-CLEANSALT',codeVers,dataname,'.txt')];
					dattabl = table(T,isCleanSalt_red,isCleanSalt_ox,ccStart_red,ccEnd_red,chargeCurrentSlope_red, ccStart_ox,ccEnd_ox, chargeCurrentSlope_ox);
					writetable(dattabl,dataout,'Delimiter','comma','WriteRowNames',true)
					saveas(tiler, [folder strcat('out-CLEANSALT',codeVers,dataname,'.png')]);%the name of the image file saved as is the dataname concatenated with the .png extension
					hold off
					close 
                    else
                    e2_alpha = (peaks_red/peaks_ox)/(1+peaks_red/peaks_ox);
                    extractedPeakAmpl_red = peaks_red;
                    ccControlledExtPeakAmpl_red = peaks_red - abs(chargeCurrent_i_fit_red(ptsPerPass - peaklocs_red));
                    extractedPeakAmpl_ox = peaks_ox;
                    ccControlledExtPeakAmpl_ox = peaks_red - abs(chargeCurrent_i_fit_red(peaklocs_ox));


                    e_noOfeExchanged = electronConvConst/abs(diff(Ecoupl)); % estimation of the number of electrons exchanged in this redox couple. 
                    noOfSpGuess = 1;
                    n = 2 ;%assumed %e_noOfeExchanged;%
		
					%% begin analysis
                    y = -1;
                 iBDlossFnValue_min = 1E8; %just some arbitrarily large no
                       % for y = -1:2:1% ends (786)

						initialScan = TotalData_VI;
						initialScan(initialScan(:,2)~=x,:) = []; 
						initialScan(initialScan(:,11)~=y,:) = []; 
						
						%% for latest batch of test data, snip off V < -2.5 & V> -0.4
						initialScan(initialScan(:,4)< Vmin_clip,:) = [];
						initialScan(initialScan(:,4)> Vmax_clip,:) = [];
						%
                        %% another way to get rid of time jump
                        initialScan(564:end,:) = [];
                        %%
						
						Time = initialScan(:,3); %this is necessary for BD functions, not so for MO
						Current= y*initialScan(:,5); %************************************this is differeing from your original Current = initialscan()
						Potentials = initialScan(:,4);
                        

						[ptsPerPass,~] = size(Time);
						[peaks_expData, loc_expData, width_expData] = findpeaks(smooth(Current),'MinPeakProminence',0.0001,'MaxPeakWidth',300,'Sortstr','descend','Npeaks',1,'Annotate','extents');%play with params here, no specifications work decently but could be optimized

							try peakfound = peaks_expData(1);
							catch 
							continue 
							end 
						%chargeCurrentStart = 1; %round(loc_expData-CCstartWidScale*width_expData);%
						%chargeCurrentEnd = round(loc_expData-CCendWidScale*width_expData);
						risingLegStart = round(loc_expData-width_expData);
						risingLegEnd = loc_expData;					
                        for b = 1:2

						switch y
						case -1

						% diagnostic
						redScan_Time = Time;
						redScan_Pot = Potentials;
						redScan_Current = Current;
						%

						
                            
                        switch b
                            case 1
                                t0_red = (Vmax_clip - e_Ered)/v; % this uses Ered, so it must come after the Ered definition statement
                            case 2
                                t0_red = t0_red+arrayShift*dt;
                        end
                        
                        Time2 = Time - Time(1);
						[loc_Ered_rows, loc_Ered_cols] = find(abs(Time2 - t0_red )< dt);% what is this value % 0.0305? it is the value of the time step
							try t0_red = Time2(loc_Ered_rows(1),loc_Ered_cols(1));% we hope there aren't more than 2 values sometimes there's one
							catch 
							continue
							end 


						tred1 = Time2(1:loc_Ered_rows(1));
						tred2 = Time2(loc_Ered_rows(1)+1:end);
						[tred1Len,~] = size(tred1);
						[tred2Len,~] = size(tred2);
						tmax = (Vmax_set-Vmin_set)/v;

						%chargeCurrent_i_fit_seg = chargeCurrent_i_fit_red(ptsPerPass - (tred2Len-1):end);

						Ered1 = Vmax_clip - v*tred1;
						Ered2 = Vmax_clip - v*tred2;

						case 1

						% diagnostic
						oxScan_Time = Time;
						oxScan_Pot = Potentials;
						oxScan_Current = Current;
						%
                       
                            switch b
                                case 1
                                t0_red = (e_Ered - Vmin_clip)/v; %Vmin_set
                                case 2
                                t0_red = t0_red+arrayShift*dt;
                            end
						Time2 = Time - Time(1); % for the 2nd scan you have to frame shift the time, regardless of if it's negative or positive. 
						[loc_Ered_rows, loc_Ered_cols] = find(abs(Time2 - t0_red )< dt);

							try t0_red = Time2(loc_Ered_rows(1),loc_Ered_cols(1));% we hope there aren't more than 2 values sometimes there's one
							catch
							continue
							end 

						tred1 = Time2(1:loc_Ered_rows(1));
						tred2 = Time2(loc_Ered_rows(1)+1:end);

						[tred1Len,~] = size(tred1);
						[tred2Len,~] = size(tred2);
						tmax = (Vmax_set-Vmin_set)/v;

						%chargeCurrent_i_fit_seg = chargeCurrent_i_fit_ox(ptsPerPass - (tred2Len-1):end);

						Ered1 = Vmin_clip + v*tred1; %Vmin_set
						Ered2 = Vmin_clip + v*tred2; %Vmin_set

						end 

							if tred2(1) < t0_red %a control in case the midpoint definition is off or something 
							tred2(1) = t0_red; % could this happen in points beyond tred2(1)? like tred(1:n) or something? 
							end


						%% ===============2020-04-08 this works, and this seems to work for the BD functions. now you just have to get the fitting right. 
						%toc
						t = 1;
						k = 1;
						
                      %% Gradient search for Concentration or Diffusion coeff
						% for k = 1:noOfSpGuess 
  
                       %% I_bd component functions =====finding the prefactor
                       THET = n*F*v/(R*T);
                       fun = @(x) exp(x.^2);
                       fun2 = @(x2) integral(fun,0,x2);
                       fun3 = @(x3) exp(-x3.^2).*fun2(x3);
                       for i = 1:tred2Len
                           i_bd_fn(i) = 2/sqrt(pi)*fun3(sqrt(THET)*sqrt(tred2(i)-tred2(1)));
                       end
                       bd_ampl_factor = max(i_bd_fn)%this bd_ampl_factor = function_peak
                       % ==========
                       ccCorrectedPeaks_red = peaks_red - abs(chargeCurrent_i_fit_red(loc_expData));
                       %Do = (ccCorrectedPeaks_red/(alpha*bd_ampl_factor*A*C_0*F*n*sqrt(F*n*v/(R*T))))^2
                       C_0 = ccCorrectedPeaks_red/(alpha*bd_ampl_factor*A*F*n*sqrt(F*n*Do*v/(R*T)))
                       ampl_SCALE = alpha*A*C_0*F*n*sqrt(F*n*v*Do/(R*T));
                       
                       %%
                       switch y 
                       case -1
                       %model cathodic scan by Berzins-Delahay function
                       for i = 1:loc_Ered_rows(1)
                           itot(i) = chargeCurrent_i_fit_red(i);
                       end
                       for i = loc_Ered_rows(1)+1:ptsPerPass
                           itot(i) = chargeCurrent_i_fit_red(i) + y*ampl_SCALE*i_bd_fn(i-loc_Ered_rows(1));
                       end
                       
                       chargeCurrent_i_fit_seg = chargeCurrent_i_fit_red(ptsPerPass - (tred2Len-1):end);
                       case 1
                       % anodic scan not modeled by Berzins-Delahay function
                       end
    
                        [size_itot,~] = size(itot);
						Etot = [Ered1; Ered2];% Eox1 ECx2];
					
						%% =====================================


                    
						[peaks_mod, loc_mod] = findpeaks(smooth(abs(itot)),'MaxPeakWidth',300,'Sortstr','descend','Npeaks',1,'Annotate','extents');%play with params here, no specifications work decently but could be optimized
						if b ==2
                        h5 = plot(Potentials(loc_mod),itot(loc_mod),'rx','LineWidth',2,'MarkerSize',24);%
                        end
						arrayShift = loc_expData - loc_mod;
							switch y
							case -1
							peakRiseBegin = peakRiseBegin_red;
							peakFallEnd = peakFallEnd_red;
							arrayShift_red = arrayShift;
							case 1
							peakRiseBegin = peakRiseBegin_ox;
							peakFallEnd = peakFallEnd_ox;
							arrayShift_ox = arrayShift;
							end
                        end

						SSres = 0;
						SStot = 0;
						% calculating the R^2 value
						i_datamean = mean(initialScan(:,5)); %this is why we need the same number of values in exp data & fit ==> to calculate R^2
							for i = peakRiseBegin:peakFallEnd %% the problem here is that the peak for (raw data) y*Current and (modeled) itot do not occur at the same *vector location*
							%you need a fix for this. 
							SStot = SStot + (abs(Current(i))-abs(i_datamean))^2; % this number is small in general 
							%plot(Potentials(i),y*Current(i),'r*')
								try SSres = SSres + (abs(itot(i))-abs(Current(i)))^2; %explained sum of sq || assuming loc_mod > loc_expData, peak_mod is forwards. 
								catch 
								arrayShift = 0;
								SSres = SSres + (abs(itot(i))-abs(Current(i)))^2; %explained sum of sq || assuming loc_mod > loc_expData, peak_mod is forwards. 
								end
							%plot(Etot(i-arrayShift),itot(i-arrayShift),'k.')
                            
							end

						SSres;
						SStot;
						Rsq = 1 - SSres/SStot;
						%text(-0.2,0.3*y*extractedPeakAmpl_ox,strcat('R^2 = ', num2str(Rsq)))
							switch y
							case -1
							Rsq_red = Rsq;
							case 1
							Rsq_ox = Rsq;
                            end
                        h6 = plot(Potentials,itot,'k','LineWidth',1);
						r = r+1;
                        % end %starts (791)

					legend([h1(1) h2(1) h3(1) h4(1) h5(1) h6(1)],'experimental data','charging current extrapolation','charging current segment','peak ROI','peak current (cathodic)','modeled CV (cathodic)','Location','northwest','FontSize',9)
                    ylim([-1.7*peaks_red,0.12])
                    xlim([-1.5,1.2*max(Potentials)])
                    dataout =[folder strcat('DATOUT-',codeVers,dataname,'x=',num2str(x),'alpha=',num2str(alpha),'Do=',num2str(Do),'C_0=',num2str(C_0),'n=',num2str(n),'.txt')];
					fileCVdata = convertCharsToStrings(dataname);
					WEMat = convertCharsToStrings(WEType);
					dattabl = table(fileCVdata,T,WEMat,alpha,Do,C_0,n,isCleanSalt_red,extractedPeakAmpl_red,ccControlledExtPeakAmpl_red,e_Ered,e_noOfeExchanged,chargeCurrentSlope_red,Rsq_red);
					writetable(dattabl,dataout,'Delimiter','comma','WriteRowNames',true)
					% end % ends iteration over all species k
					% 

					% %======================END FITTING MODULE FOR SINGLE PASS

					

					%% save data as matrix, save plots as file.png

 					saveas(tiler, [folder strcat('plot',codeVers,dataname,'x=',num2str(x),'.fig')]);%
                    saveas(tiler, [folder strcat(titl,' scan#',num2str(x),'.png')]);%
                    end % For x = ... ends iteration over number of scans.

					end %if isCleanSalt_red/ox ==1

end %this is the end of the while fid open loop
toc
fclose(fid); % close the file

fprintf(strcat(codeVers,' bien.'));