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
classdef DataClass < dynamicprops
    properties
        Interpolated; %Resulting dataset in a cell of chrom. pair matrices
        Inter;
        Intra;
    end
    methods
        %constructor
        function obj=DataClass(~)
        end
        
        %% Parse raw data
        %Updates the main-diagonal of obj.P_M to hold the raw data
        function obj=M_ParseIntraRaw(obj)
            display('Parsing intrachromosomal interactions.')
            global PARAMETERS
            obj.Intra=dlmread(PARAMETERS.intrafrequencies,'\t',1,0);
            display('Done.');
        end
        
        %Updates the off-diagonal of obj.P_M to hold the raw data
        function obj=M_ParseInterRaw(obj)
            display('Parsing interchromosomal interactions.')
            global PARAMETERS
            obj.Inter=dlmread(PARAMETERS.interfrequencies,'\t',1,0);
            display('Done.')
        end
        
        %% Preprocessing
        %Scattered interpolation function for converting fragments
        %to 1kb bins in the frequency matrix "obj.P_Gdist".
        function obj=M_GenomeTriScatteredInterp(obj)
            display('Interpolating chromosomal interactions.')
            global PARAMETERS
            chromsizes=PARAMETERS.ChromosomeSizes;
            binsize=PARAMETERS.binningresolution;
            tmax=max([obj.Inter(:,5);obj.Intra(:,5)]);
            
            M=[obj.Inter;obj.Intra];
            [C,~,ic]=unique(M(:,[1, 3]),'rows');
            C=num2cell(C,2);
            cM=cell(size(C,1),1);
            for i=1:size(C,1)
                cM{i}=M(ic==i,:);
            end
            supmat=cell(length(cM),1);
            parfor chrom=1:length(cM) %multi-core
                chromx=C{chrom}(2); chromy=C{chrom}(1);
                if (chromx>=chromy) %multi-core
                    lenx=chromsizes(chromx); leny=chromsizes(chromy);
                    ti=floor(binsize/2):binsize:lenx; tj=floor(binsize/2):binsize:leny;
                    [qx,qy]=meshgrid(ti,tj);
                    display(horzcat('Chromosome ',num2str(chromx),' vs ',num2str(chromy)))
                    x=cM{chrom}(:,2);
                    y=cM{chrom}(:,4);
                    z=cM{chrom}(:,5);
                    hasvals=and(~isnan(z),z>0);
                    tx=x(hasvals); ty=y(hasvals); tz=z(hasvals);
                    ftx=[]; fty=[]; ftz=[];
                    tlenx=ceil(lenx/10); tleny=ceil(leny/10);
                    if (chromx==chromy) %ALG: force diagonal
                        ftx=([(floor(binsize/2):binsize:lenx)']);
                        fty=([(floor(binsize/2):binsize:leny)']);
                        ftz=[repmat(tmax,1,ceil((leny-floor(binsize/2))/binsize))'];
                    end
                    %ALG: insert shadow points for extrapolation
                    ftx=([ftx; (-tlenx:binsize:(lenx+tlenx))';...
                        repmat(-tlenx,length(-tleny:binsize:(leny+tleny)),1);...
                        repmat((lenx+tlenx),length(-tleny:binsize:(leny+tleny)),1);...
                        (-tlenx:binsize:(lenx+tlenx))']);
                    fty=([fty; repmat(-tleny,length(-tlenx:binsize:(lenx+tlenx)),1);...
                        (-tleny:binsize:(leny+tleny))';...
                        (-tleny:binsize:(leny+tleny))';...
                        repmat((leny+tleny),length(-tlenx:binsize:(lenx+tlenx)),1)]);
                    ftz=[ftz; zeros(2*length(-tlenx:binsize:(lenx+tlenx))+2*length(-tleny:binsize:(leny+tleny)),1)];
                    
                    tx=[tx; ftx]; ty=[ty; fty]; tz=[tz; ftz];
                    uvals=unique([tx ty tz],'rows');
                    %ALG: actual interpolation part:
                    F=TriScatteredInterp(uvals(:,1),uvals(:,2),uvals(:,3),'natural');
                    
                    qz=F(qx,qy);
                    if (chromx==chromy)
                        qz=(qz+qz')./2; %fixes nonlinearity invariance of interpolation to generate a valid dissimilarity matrix.
                        qz(logical(eye(size(qz))))=tmax;
                        
                    end
                    supmat{chrom}=qz'; %multi-core
                end
            end
            finmat=cell(length(chromsizes));
            for i=1:length(cM)
                chromx=C{i}(2); chromy=C{i}(1);
                finmat{chromx,chromy}=supmat{i};
                finmat{chromy,chromx}=supmat{i}';
            end
            obj.Interpolated=finmat;
            display('Done.');
        end
        
        %% Analysis
        %calculates Spatial and Genomic co-localization enrichments.
        function [allpvals3d, allpvals2d, allbestspheresizes, allbestintsizes]=...
                AnalyzeCoLocalization(obj, annotationmat, loci)
            display('Starting analysis on co-localization of markers');
            global PARAMETERS
            tM=cell2mat(obj.Interpolated);
            M=mat2cell(tM,ones(1,length(tM)),length(tM));
            cumbins=cumsum(unique(cellfun(@length,(cellfun(@(x) x(1,:), obj.Interpolated,'uniformoutput',0))),'rows'));
            tcumbins=[1 cumbins];
            tmat=zeros(length(tM),size(annotationmat,2));
            tmat(loci,:)=annotationmat;

            if PARAMETERS.analysis_shuffle
                shufl=randperm(length(unique(loci)));
                tmat(unique(loci),:)=tmat(shufl,:);
            end
            
             if PARAMETERS.analysis_cyclic
                for i=1:length(PARAMETERS.ChromosomeSizes)
                    ShiftSize=randi([ceil(.1*cumbins(i)),.9*cumbins(i)]);
                    tmat(tcumbins(i):tcumbins(i+1),:)=circshift(tmat(tcumbins(i):tcumbins(i+1),:),[0 ShiftSize]);
                end
            end
            tmat=num2cell(tmat,1);
            numterms=length(tmat);
            binloci=false(length(M),1); binloci(loci)=1;
            tmpM=M(loci);
            
            % Check gene 3d enrichment vs 1d
            allpvals3d=ones(numterms,length(loci)); allbestspheresizes=zeros(numterms,length(loci));
            allpvals2d=ones(numterms,length(loci)); allbestintsizes=zeros(numterms,length(loci));
            
            parfor o=1:length(loci)
                %3d
                currbin=loci(o);
                [~,idx_i]=sort(tmpM{o},'descend');
                idx_i=idx_i(binloci(idx_i)); %bins w/genes ordered by 3d dist
                tsphere=zeros(1,numterms); p=ones(1,numterms);
                interval=zeros(1,numterms); p2=ones(1,numterms);
                %interval
                Chrom=find(arrayfun(@(x) currbin<=x, cumbins),1); ChromVec=[tcumbins(Chrom):tcumbins(Chrom+1)];
                [~,idx_int]=sort(abs(ChromVec-currbin),'ascend');
                ChromVec=ChromVec(idx_int);
                randrest=setdiff(1:tcumbins(end),ChromVec);
                intvec=[ChromVec randrest(randperm(length(randrest)))];
                updintvec=intvec(binloci(intvec));
                tChromVec=ChromVec(binloci(ChromVec));
                for i=1:numterms
                    [tsphere(i),p(i)]=mHGthresh(full(tmat{i}(idx_i)),length(tChromVec));
                    [interval(i),p2(i)]=mHGthresh(full(tmat{i}(updintvec)),length(tChromVec));
                end
                allpvals3d(:,o)=p; allbestspheresizes(:,o)=tsphere;
                allpvals2d(:,o)=p2; allbestintsizes(:,o)=interval;
                
                display(horzcat('Finished locus #',num2str(o)));
            end
            display('Done.');
        end
        %corrects for multiple hypothesis testing.
        function [corrected3d,corrected1d]=CorrectResults(obj, AnnotationData, allpvals3d, allpvals2d)
            lipsoncorrection=repmat(cellfun(@(x) sum(full(x)),num2cell(AnnotationData,1))',1,size(AnnotationData,1));
            lipsoncorrection(lipsoncorrection==0)=nan;
            corrected3d=(allpvals3d.*lipsoncorrection);
            corrected1d=(allpvals2d.*lipsoncorrection);
        end
        %shows the resulting genomic landscape.
        function ShowFigure(obj, Loci, scores3d, scores2d)
            cumbins=cumsum(unique(cellfun(@length,(cellfun(@(x) x(1,:), obj.Interpolated,'uniformoutput',0))),'rows'));
            tcumbins=[1 cumbins];
            
            h=zeros(1,2);
            h_1=scatter(Loci,-log10(scores3d),50,'r','filled');
            hold on;
            h_2=scatter(Loci,-log10(scores2d),70,'b');
            ylabel('P-value(-log_1_0, corrected)'); xlabel('Bin number across the genome');
            
            ptch=patch([-100 -100 cumbins(end)+100 cumbins(end)+100]',[min(-log10([scores3d scores2d]))-.5 -log10(.05) -log10(.05) min(-log10([scores3d scores2d]))-.5]',1);
            
            alpha(ptch,.2); set(ptch,'facecolor',[.7 .7 .7]);
            wrn=warning('off');
            [pks3,locs3]=findpeaks(-log10(scores3d),'minpeakdistance',10,'minpeakheight',-log10(.05));
            [pks2,locs2]=findpeaks(-log10(scores2d),'minpeakdistance',10,'minpeakheight',-log10(.05));
            warning(wrn);
            arrayfun(@(x) arrow([Loci(x) max(pks3(locs3==x),pks2(locs2==x))+.45],[Loci(x) max(pks3(locs3==x),pks2(locs2==x))+.1],12,45,'facecolor','k'),intersect(locs3,locs2));
            [tlocs,tlocsid]=setdiff(locs3,locs2);
            arrayfun(@(x) arrow([Loci(x) pks3(locs3==x)+.45],[Loci(x) pks3(locs3==x)+.1],12,45,'facecolor','r'),tlocs);
            
            [tlocs,tlocsid]=setdiff(locs2,locs3);
            arrayfun(@(x) arrow([Loci(x) pks2(locs2==x)+.45],[Loci(x) pks2(locs2==x)+.1],10,45,'facecolor','b'),tlocs)
            axis tight;
            legend('3D enrichment','1D enrichment','P-value > .05 (corrected)','location','best');
        end
    end
    % helper functions
end