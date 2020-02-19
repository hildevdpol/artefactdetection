%% Alarmering & Artefact detectie
% Artefact detectie: Beslisboom
% Hilde van der Pol

clear all; close all; clc;

%% Parameters and reading excel documents

path= '/Users/hildevanderpol/Desktop/artefactdetectie/KT3401_dataWeek2_2019'; 
folder= 'Delft04'; 
filename= 'D04Cal (2).xls'; 
fs= 100; 

[t, ABP, Paw, CVP] = readArtefacts(path, folder, filename, fs); 

%% Insparitory hold Paw (artefact 1) 
abs_diff_Paw = abs(diff(Paw)); 
fs = 100; 
t_absinsp = (1/fs:1/fs:length(ABP)/fs)';
heightinsp = mean(abs_diff_Paw) + 1.5*std(abs_diff_Paw); 
[pks,locs] = findpeaks(abs(diff(Paw)), fs,'MinPeakHeight',heightinsp,'MinPeakDistance',1);
findLoc_Paw = find(diff(locs) >5);

%start and endtime artefact
inspwaarde = (t>locs(findLoc_Paw) & t<locs(findLoc_Paw+1));
movsuminsp=movsum(inspwaarde,2)>=1;
tfindloc_inspstart1= t(diff(movsuminsp)==1);
tfindloc_inspend1= t(diff(movsuminsp)==-1);

if length(tfindloc_inspstart1)>length(tfindloc_inspend1); 
    tfindloc_inspend1(length(tfindloc_inspend1)+1)= t(length(t)); 
end
if length(tfindloc_inspstart1)<length(tfindloc_inspend1); 
    tfindloc_inspstart1= [t(1) tfindloc_inspstart1]; 
end

% eliminating short artefacts 
tfindloc_inspstart= tfindloc_inspstart1((tfindloc_inspend1- tfindloc_inspstart1)>3);
tfindloc_inspend= tfindloc_inspend1((tfindloc_inspend1- tfindloc_inspstart1)>3);

%putting the results in a celarray
cellarray_insp={};
if ~isempty((tfindloc_inspstart)) && ~isempty((tfindloc_inspend))
    for i = 1:1:length(tfindloc_inspstart)
    cellarray_insp(i,:) = {tfindloc_inspstart(i), tfindloc_inspend(i), 'Calibratie', 'ABP&CVP'};
    end
end

%% Calibratie ABP & CVP (artefact 2) 
M_ABP_Cal= movmean(ABP,100);
M_CVP_Cal= movmean(CVP,100);

%start and endtime artefact
calwaarde = (M_CVP_Cal > -5 & M_CVP_Cal < 5 & M_ABP_Cal<5 & M_ABP_Cal> -5);
movsumcal=movsum(calwaarde,2)>=1;
tfindloc_calstart1= t(diff(movsumcal)==1);
tfindloc_calend1= t(diff(movsumcal)==-1);


if length(tfindloc_calstart1)>length(tfindloc_calend1); 
    tfindloc_calend1(length(tfindloc_calend1)+1)= t(length(t)); 
end
if length(tfindloc_calstart1)<length(tfindloc_calend1); 
    tfindloc_calstart1= [t(1) tfindloc_calstart1]; 
end


tfindloc_calstart= tfindloc_calstart1((tfindloc_calend1- tfindloc_calstart1)>3);
tfindloc_calend= tfindloc_calend1((tfindloc_calend1- tfindloc_calstart1)>3);


%celarray
cellarray_cal={};
if ~isempty((tfindloc_calstart)) && ~isempty((tfindloc_calend))
    for i = 1:1:length(tfindloc_calstart)
    cellarray_cal(i,:) = {tfindloc_calstart(i), tfindloc_calend(i), 'Calibratie', 'ABP&CVP'};
    end
end

%% Infuus CVP (artefact 3) 

infuus_threshold = mean(CVP) + 0.5*std(CVP);
infuusx = CVP > infuus_threshold;
CVP_M = movmean(infuusx ,35);
tinfuus = CVP_M==1;
movmeaninfuus= movmean(tinfuus,400); 

Minfuus = movmean(CVP,600);
infuusbovengrens = 3*mean(Minfuus)+ 1*std(Minfuus); 

%start and endtime artefact
infuuswaarde = (movmeaninfuus >= 0.03 & Minfuus < infuusbovengrens);
movsuminfuus=movsum(infuuswaarde,2)>=1;
tfindloc_infuusstart1= t(diff(movsuminfuus)==1);
tfindloc_infuusend1= t(diff(movsuminfuus)==-1);



if length(tfindloc_infuusstart1)>length(tfindloc_infuusend1); 
    tfindloc_infuusend1(length(tfindloc_infuusend1)+1)= t(length(t)); 
end
if length(tfindloc_infuusstart1)<length(tfindloc_infuusend1); 
    tfindloc_infuusstart1= [t(1) tfindloc_infuusstart1];
end


tfindloc_infuusstart= tfindloc_infuusstart1((tfindloc_infuusend1- tfindloc_infuusstart1)>4);
tfindloc_infuusend= tfindloc_infuusend1((tfindloc_infuusend1- tfindloc_infuusstart1)>4);

%celarray
cellarray_infuus={};
if ~isempty((tfindloc_infuusstart)) && ~isempty((tfindloc_infuusend))
    for i = 1:1:length(tfindloc_infuusstart)
    cellarray_infuus(i,:) = {tfindloc_infuusstart(i), tfindloc_infuusend(i), 'Infuus', 'CVP'};
    if  tfindloc_infuusstart > tfindloc_infuusend
        cellarray_infuus={}
    end
    end
end

%% Transducer hoog (artefact 4) 

%ABP
M_ABP_transd= movmean(ABP,800);
threshold_ABP_transd = mean(M_ABP_transd)-0.1*std(ABP); 
artLoc_ABP_transd = M_ABP_transd<threshold_ABP_transd ;

% CVP
M_CVP_transd= movmean(CVP,700);
threshold_CVP_transd = mean(M_CVP_transd)-0.2*std(CVP); 
artLoc_CVP_transd = M_CVP_transd<threshold_CVP_transd ;

%start and endtime artefact
transdwaarde = (M_ABP_transd < threshold_ABP_transd  & M_CVP_transd < threshold_CVP_transd & M_ABP_transd > threshold_CVP_transd);
movsumtransd=movsum(transdwaarde,2)>=1;
tfindloc_transdstart1= t(diff(movsumtransd)==1);
tfindloc_transdend1= t(diff(movsumtransd)==-1);


if length(tfindloc_transdstart1)>length(tfindloc_transdend1); 
    tfindloc_transdend1(length(tfindloc_transdend1)+1)= t(length(t)); 
end
if length(tfindloc_transdstart1)<length(tfindloc_transdend1);  
    tfindloc_transdstart1= [t(1) tfindloc_transdstart1];
end


tfindloc_transdstart= tfindloc_transdstart1((tfindloc_transdend1- tfindloc_transdstart1)>3);
tfindloc_transdend= tfindloc_transdend1((tfindloc_transdend1- tfindloc_transdstart1)>3);

%celarray
cellarray_transd={};
if ~isempty((tfindloc_transdstart)) && ~isempty((tfindloc_transdend))
    for i = 1:1:length(tfindloc_transdstart)
    cellarray_transd(i,:) = {tfindloc_transdstart(i), tfindloc_transdend(i), 'Transducer', 'ABP&CVP'};
    end
end


%% Slinger ABP (artefact 5) 
resolution = 0.5; 
[s, frx, tm] = spectrogram(ABP, resolution*fs,0,[0.5:0.1:50],fs);   
frange_slinger = [15 40]; 
f_amp_slinger = mean(abs(s(frx >= frange_slinger(1) & frx <= frange_slinger(2),:)));
M_f_amp_slinger = movmean(f_amp_slinger,5)';
threshold_slinger = mean(M_f_amp_slinger) + 0.2*std(f_amp_slinger);

%start and endtime artefact
slingerwaarde = (M_f_amp_slinger > threshold_slinger);
movsumslinger=movsum(slingerwaarde,2)>=1;
tfindloc_slingerstart1= tm(diff(movsumslinger)==1);
tfindloc_slingerend1= tm(diff(movsumslinger)==-1);

if length(tfindloc_slingerstart1)>length(tfindloc_slingerend1); 
    tfindloc_slingerend1(length(tfindloc_slingerend1)+1)= t(length(t)); 
end
if length(tfindloc_slingerstart1)<length(tfindloc_slingerend1); 
    tfindloc_slingerstart1= [t(1) tfindloc_slingerstart1]; 
end


tfindloc_slingerstart= tfindloc_slingerstart1((tfindloc_slingerend1- tfindloc_slingerstart1)>4);
tfindloc_slingerend= tfindloc_slingerend1((tfindloc_slingerend1- tfindloc_slingerstart1)>4);

%celarray
cellarray_slinger={};
if ~isempty((tfindloc_slingerstart1)) && ~isempty((tfindloc_slingerend))
    for i = 1:1:length(tfindloc_slingerstart)
    cellarray_slinger(i,:) = {tfindloc_slingerstart(i), tfindloc_slingerend(i), 'Slinger', 'ABP'};
    end
end
%% Flush (artefact 6)
%ABP
M_ABPflush= movmean(ABP,200);
threshold_ABPflush = mean(M_ABPflush)+0.9*std(ABP); 

%CVP
M_CVPflush= movmean(CVP,200);
threshold_CVPflush = 2.5*mean(M_CVPflush)+ 1.5*std(CVP); 

%ABP
flushwaardeABP = (M_ABPflush > threshold_ABPflush);
movsumflushABP = movsum(flushwaardeABP,2)>=1;
tfindloc_flushstartABP1= t(diff(movsumflushABP)==1);
tfindloc_flushendABP1= t(diff(movsumflushABP)==-1);

%CVP
flushwaardeCVP = (M_CVPflush > threshold_CVPflush);
movsumflushCVP = movsum(flushwaardeCVP,2)>=1;
tfindloc_flushstartCVP1= t(diff(movsumflushCVP)==1);
tfindloc_flushendCVP1= t(diff(movsumflushCVP)==-1);

if length(tfindloc_flushstartCVP1)>length(tfindloc_flushendCVP1); 
    tfindloc_flushendCVP1(length(tfindloc_flushendCVP1)+1)= t(length(t)); 
end
if length(tfindloc_flushstartCVP1)<length(tfindloc_flushendCVP1); 
    tfindloc_flushstartCVP1= [t(1) tfindloc_flushstartCVP1]; 
end


if length(tfindloc_flushstartABP1)>length(tfindloc_flushendABP1); 
    tfindloc_flushendABP1(length(tfindloc_flushendABP1)+1)= t(length(t)); 
end
if length(tfindloc_flushstartABP1)<length(tfindloc_flushendABP1); 
    tfindloc_flushstartABP1= [t(1) tfindloc_flushstartABP1]; 
end

tfindloc_flushstartABP= tfindloc_flushstartABP1((tfindloc_flushendABP1- tfindloc_flushstartABP1)>3);
tfindloc_flushendABP= tfindloc_flushendABP1((tfindloc_flushendABP1- tfindloc_flushstartABP1)>3);


tfindloc_flushstartCVP= tfindloc_flushstartCVP1((tfindloc_flushendCVP1- tfindloc_flushstartCVP1)>3);
tfindloc_flushendCVP= tfindloc_flushendCVP1((tfindloc_flushendCVP1- tfindloc_flushstartCVP1)>3);


cellarray_flushABP={};
if ~isempty((tfindloc_flushstartABP)) && ~isempty((tfindloc_flushendABP))
    for i = 1:1:length(tfindloc_flushstartABP)
    cellarray_flushABP(i,:) = {tfindloc_flushstartABP(i), tfindloc_flushendABP(i), 'Slinger', 'ABP'};
    end
end

cellarray_flushCVP={};
if ~isempty((tfindloc_flushstartCVP)) && ~isempty((tfindloc_flushendCVP))
    for i = 1:1:length(tfindloc_flushstartCVP)
    cellarray_flushCVP(i,:) = {tfindloc_flushstartCVP(i), tfindloc_flushendCVP(i), 'Slinger', 'ABP'};
    end
end

%% Gasbel (gass artefact) (artefact 7)
Movgasbel = movmean(ABP, 80);
meangasbel= mean(ABP); 
thresholdgasbel = meangasbel +std(ABP);

vargasbel = 20*movvar(ABP, 80); 
vargasbel1= movvar(ABP,80);
gasbelgrensvar= 1*mean(vargasbel1);

%start and endtime artefact
gasbelwaarde = (Movgasbel >= 20 & Movgasbel <= thresholdgasbel & vargasbel <= gasbelgrensvar);
movsumgasbel=movsum(gasbelwaarde,200)>=15;
tfindloc_gasbelstart1= t(diff(movsumgasbel)==1);
tfindloc_gasbelend1= t(diff(movsumgasbel)==-1);

if length(tfindloc_gasbelstart1)>length(tfindloc_gasbelend1); 
    tfindloc_gasbelend1(length(tfindloc_gasbelend1)+1)= t(length(t)); 
end
if length(tfindloc_gasbelstart1)<length(tfindloc_gasbelend1); 
    tfindloc_gasbelstart1 =[t(1) tfindloc_gasbelstart1];
end


tfindloc_gasbelstart= tfindloc_gasbelstart1((tfindloc_gasbelend1- tfindloc_gasbelstart1)>3);
tfindloc_gasbelend= tfindloc_gasbelend1((tfindloc_gasbelend1- tfindloc_gasbelstart1)>3);

%celarray
cellarray_gasbel={};
if ~isempty((tfindloc_gasbelstart)) && ~isempty((tfindloc_gasbelend))
    for i = 1:1:length(tfindloc_gasbelstart)
    cellarray_gasbel(i,:) = {tfindloc_gasbelstart(i), tfindloc_gasbelend(i), 'Gasbel', 'ABP'};
       end
end

%% making a central cell array to combine all 7 artefacts 
Cellarraytotaal = [cellarray_insp' cellarray_cal' cellarray_infuus' cellarray_transd' cellarray_slinger' cellarray_flushABP' cellarray_flushCVP' cellarray_gasbel']
