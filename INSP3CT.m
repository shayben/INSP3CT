%{
    Copyright 2012 Shay Ben Elazar ©
    This file is part of INSP3CT.

    INSP3CT is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    INSP3CT is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with INSP3CT.  If not, see <http://www.gnu.org/licenses/>.
%}
%% INSP3CT pipeline
%Runtime parameters. Make sure to set these according to the README before running the code.
global PARAMETERS
PARAMETERS.ChromosomeSizes = [];
PARAMETERS.interfrequencies = '';
PARAMETERS.intrafrequencies = '';
PARAMETERS.cores = 1;
PARAMETERS.binningresolution = 10000;
PARAMETERS.analysis_shuffle = false;
PARAMETERS.analysis_cyclic = false;

%example parameters for YeastDatasetExample.zip input.
%PARAMETERS.ChromosomeSizes =  [230218, 813184, 316620, 1531933, 576874,...
%270161, 1090940, 562643, 439888, 745751, 666816, 1078177, 924431, 784333,...
%1091291, 948066]; %Yeast genome example
%PARAMETERS.interfrequencies = 'interactions_HindIII_MspI_beforeFDR_inter.txt';
%PARAMETERS.intrafrequencies = 'interactions_HindIII_MspI_beforeFDR_intra_all.txt';

% Assert input and parameters
assert(~isempty(PARAMETERS.interfrequencies) || ~isempty(PARAMETERS.intrafrequencies),...
    ~isempty(PARAMETERS.ChromosomeSizes) || ~isempty(PARAMETERS.binningresolution),...
    'ERROR: Critical parameters are missing. Follow the README for instructions.');

if (PARAMETERS.cores > 1)
    try
        matlabpool ('open', PARAMETERS.cores)
    catch e
        display (e)
    end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Load and parse datasets
obj=DataClass();
obj=obj.M_ParseIntraRaw();
obj=obj.M_ParseInterRaw();

%% Generate interpolated and binned frequencies
obj=obj.M_GenomeTriScatteredInterp();

%% Analysis on annotation
%Example of annotation and loci datasets:
%Loci=randi(sum(cellfun(@length,obj.Interpolated(1,:))),1,100); %100 Randomly sampled loci.
%AnnotationData=randi(ceil(length(Loci)/PARAMETERS.binningresolution)+1,length(Loci),10)-1; %10 Random annotations.
assert(exist('AnnotationData','var') || exist('Loci','var'),...
'ERROR: You need to define an annotation matrix and a Loci vector to continue. Check README.');

[allpvals3d, allpvals2d, allbestspheresizes, allbestintsizes]=...
obj.AnalyzeCoLocalization(AnnotationData, Loci);
[corrected3d,corrected1d]=obj.CorrectResults(AnnotationData, allpvals3d, allpvals2d);

%% Show a result figure for each annotation
assert(exist('AnnotationData','var') || exist('Loci','var'),...
'ERROR: You need to define an annotation matrix and a Loci vector to continue. Check README.');
assert(exist('corrected3d','var') || exist('corrected1d','var'),...
'ERROR: You need to first run the co-localization analysis before attempting to plot the results.');

fig=figure('Color','w');
set(fig, 'defaulttextfontsize',18,'defaultlinelinewidth',.1,'defaultaxesfontsize',18);
for currannotation=1:size(AnnotationData,2)
    obj.ShowFigure(Loci, corrected3d(currannotation,:), corrected1d(currannotation,:));
    title(horzcat('Annotation #',num2str(currannotation)));
    input('Press any key to continue');
    cla;
end