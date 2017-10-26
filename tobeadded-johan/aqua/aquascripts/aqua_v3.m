function QA=aqua(flags)
%
%AQUA: Automated Quality Assurance protocol for fMRI data.
%      --        --      -
%
%Provides several attributes about statistical properties 
%and noise-level in an EPI time-series.
%
%Usage in four different modes:
%
% >> aqua             : INTERACTIVE MODE: prompts for fMRI data and other info,
%                       computes QA structure and saves it to local directory,
%                       (QA.mat), displays results and saves plots in qa.ps
%
% >> QA=aqua(flags)   : BATCH MODE: computes QA structure.
%             
% >> aqua(QA)         : GRAPHIC ONLY: display/save plots of existing QA struct.
%8
% >> flags=aqua       : FLAG OUTPUT: just returns the default flag structure.
%
% The routine is able to treat several subjects and runs in one call.
%
%
% Input Structure 'flags'
% =======================
%
%  flags.data         : N times M cell array, where N is the number of subjects
%                       and M is the number of runs. Each cell contains a
%                       string array of filenames.
%  flags.threshold    : Only voxels with "intensity/sigma > threshold" are
%                       taken into account. default = 20
%  flags.eye_remove   : If true, tries to find a single connected region in the
%                       mask, which should exclude the eyes.
%                       default = 1
%  flags.n_discard    : Number of pre-scans to discard.
%                       default = 0
%  flags.realign      : If true, do realign images.
%                       (Writes affine-transformation to mat-file.)
%                       default = 0
%  flags.reslice      : If true, do reslice images.
%                       (Only internally, not saved to disk.)
%                       default = 1
%  flags.scans        : An index-vector specifying the scans for analysis.
%                       Omitted, if empty. Default = []
%  flags.save_result  : If true, results will be stored to QA.mat in the pwd.
%                       default = 1
%  flags.plot_result  : If true, graphical output of results are provided.
%                       default = 1
%  flags.save_plot    : If true, graph will be printed to qa.ps in the pwd.
%                       default = 1
%
%
% Output Structure 'QA'
% =====================
%
% QA is an array, if analysis is done for several subjects in one call.
%
% QA.name : name of subject directory (wihtout full path)
% QA.stat : structure of statistical results
%
% QA.stat is an array, if analysis is done for several runs.
%
% QA.stat.run         : name of subdirectory (run)
% QA.stat.sigma_bg    : noise level computed by 'background method'
% QA.stat.mean_img    : the mean image
% QA.stat.threshold   : threshold as defined in the input flags
% QA.stat.mean_intensity : along all voxels in the masked mean image
% QA.stat.mask        : The mask. It is only stored, if flags.eye_remove=1,
%                       since otherwise the mask can be reobtained.
% QA.stat.D           : Structure of some distribution distance measures
%                       .D_ks = Kolmogorov-Smirnov
%                       .D_ku = Kuiper Statistic
%                       .D_ad = Anderson-Darling
%                       .D_ts = my own measure ...
% QA.stat.r_qq        : Corellation coefficient of q-q plot
% QA.stat.sigma       : Noise level computed by slope of q-q plot
% QA.stat.psc         : Global percentage signal change
% QA.stat.qq_plot     : Samples of the q-q plot (only 500 points)
% QA.stat.temp        : Structure of temporal variations
%                       .r_qq   = cor.-coeff. for each timepoint (scan)
%                       .psc    = percentage signal change for each scan
%                       .pmv    = percentage mean variation for each scan
%                       .psc_sl = psc for each scan and slice
%                       .rp     = SPM2 realignment paramters
%
 
% Important remarks:
%
% - The routine needs SPM2 functions. SPM5 works as well, if the marsbar spm_get.m
%   routine is additionally somewhere in the searchpath.
%   see: http://marsbar.sourceforge.net/doc-devel/latest/marsbar/spm5/spm_get.html 
%
% - The routine should be applied to realigned data. If the images are resliced
%   already, the realign- and reslice-flag should be zero. If the data is
%   realigned but not resliced, only the realign-flag should be zero.

%
% Further information:
%
% Reference:
% Stoecker T. et al 2005. Automated Quality Assurance Routines for fMRI Data
% Applied to a Multicenter Study. Human Brain Mapping 25(2):237-246. 
%
% Copyright (C) 2005 Tony Stoecker
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the condition that the above
% copyright notice and this permission notice shall be included in all copies
% or substantial portions of the Software.
%
% Further, any publication partly arising from the use of the Software
% should reference to the article given above.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED,  INCLUDING BUT NOT  LIMITED TO THE WARRANTIES  OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR  PURPOSE AND  NONINFRINGEMENT. IN NO  EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN  AN ACTION OF CONTRACT,  TORT OR OTHERWISE,  ARISING
% FROM, OUT  OF OR  IN CONNECTION  WITH  THE SOFTWARE  OR  THE USE  OR OTHER
% DEALINGS IN THE SOFTWARE.
%

global defaults

%return default-flags and quit
if(nargin<1 && nargout==1) 
 flags.data        = []; 
 flags.threshold   = 20;
 flags.eye_remove  = 1;
 flags.n_discard   = 0;
 flags.scans       = [];
 flags.realign     = 0; 
 flags.reslice     = 0; %defalt=1
 flags.save_result = 1;
 flags.plot_result = 1;
 flags.save_plot   = 1;
 QA=flags;
 return
end 

%fill flags-structure in interactive mode
if(nargin<1 && nargout<1) 
 flags=qa_input;
end

%loop over subjects and runs, if input is the 'flags' structure
if isfield(flags,'data')

 fprintf('\n')
 [N,M]=size(flags.data);
 for i=1:N 
  for j=1:M
  if (~isempty(flags.data{i,j}))

   %get names for subjects and runs
   %simply from the directory structure
   %of the first image
   p=fileparts(flags.data{i,j}(1,:));
   [p,run]=fileparts(p);
   if j==1
    [tmp,~]=fileparts(p);
    [~,subject]=fileparts(tmp); % quick fix for me
   end
   QA(i).name = subject;
 
   fprintf('\n\nQA for subject %i (%s) , run %i (%s)\n',i,subject,j,run)
   %images of currrent run
   P = flags.data{i,j}([1+flags.n_discard:end],:); 
   %take only defined scans into account (e.g. block design)
   if (~isempty(flags.scans))
     P = P(flags.scans,:);
   end
   %calculate statistics 
      QA(i).stat(j) = qa_stats( P , flags.realign    , flags.reslice   ,...
                                    flags.threshold  , flags.eye_remove,...
                                    flags.plot_result, flags.save_plot , subject , run );
  end
  end
 end
 fprintf('\n')
 
 %% VB: change here to save differently
 %save results (if QA.mat exists, append date&time to filename)
 if flags.save_result
  QAMATFILE='QA';
  if exist([pwd,'/QA.mat'])==2
   QAMATFILE=[QAMATFILE,'_',strrep(strrep(strrep(datestr(now),' ','_'),':','_'),'-','_')];
  end
  try
      eval(['save ',QAMATFILE,' QA flags']);
  catch
      keyboard;
  end
 end
 
 % edit VB 22.08.2014
 disp('r')
 r_qq = QA.stat.r_qq
 disp('PSC')
 PSC  = QA.stat.psc
 
%display and save plot results, if input is QA structure
elseif isfield(flags,'stat')
  qa_plot(flags,1)
end
 
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% implementation of local functions %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------
%fill flags structure in interactive mode
%----------------------------------------
function flags=qa_input

 %get image data info
 nsub = spm_input('number of subjects',1,'i','1',1);
 nrun = spm_input('number of runs per subject',1,'i','1',1);
 bsame=1;
 if (nsub>1 || nrun>1),bsame = spm_input('Same num. of scans for all subjects/runs?',1,'b',['Yes';'No '],[1 0],1);end
 if (bsame), N = spm_input('number of scans per run',1,'i','1',1); end
 for i=1:nsub
  for j=1:nrun
   if (~bsame), N = spm_input(sprintf('subject %d, run %d: number of scans',i,j),1,'i','1',1); end
   % flags.data{i,j}=spm_get(N,'*.img',sprintf('select images for subject %d, run %d',i,j));
   flags.data{i,j}=spm_select([1 Inf],'image',sprintf('select images for subject %d, run %d',i,j),{},pwd,'.*',1);
  end
 end

 %get preprocessing info
 flags.n_discard = spm_input('number of Prescans to omit',1,'i','0',1);
 N = N-flags.n_discard;
 flags.realign = ~spm_input('Is data realigned?',1,'b',['Yes';'No '],[1 0],1);
 if ~flags.realign
   flags.reslice = ~spm_input('Is the realigned data resliced?',1,'b',['Yes';'No '],[1 0],0);
 else
    flags.reslice=1;
 end
 
 % keyboard;

 %get default params for qa analysis
 flags.eye_remove = spm_input('Apply eye removal?',1,'b',['Yes';'No '],[1 0],1);
 flags.threshold = spm_input('threshold for image masking',1,'r','20',1,[0 inf]);
 flags.scans=[];%default can not be overwritten, if number of scans vary with subjects/runs
 if (bsame),flags.scans = find(spm_input('vector specifying scans for QA analysis',1,'i',['ones(',num2str(N),',1)'],N));end

 %set other flags to default values
 flags.save_result = 1;
 flags.plot_result = 1;
 flags.save_plot   = 1;

return


%------------------------------------------------
%calculation of qa analysis for a single fMRI run
%------------------------------------------------

function stat=qa_stats(P,REALIGN,RESLICE,THRESHOLD,EYE_REMOVE,PLOT,PRINTPLOT,NAME,RUN)
%
%This function does all calculation of the quality assurance routine.
%More or less no documentation at all. Sorry. (Write some??)
%

%TS 09/2004

 global defaults
 stat.run=RUN;
 N=size(P,1); %number of scans

 %
 % realign images, if flag is set
 %
 if REALIGN
   fprintf('\n   -> SPM realignment')
   %graphics-output
    [Finter,Fgraph,CmdLine] = spm('FnUIsetup','Realign');
    spm('FigName',['realigning run ' deblank(RUN)],Finter,CmdLine);
   V=spm_vol(P);
   realign_flags = struct('quality',defaults.realign.estimate.quality,'fwhm',5,'rtm',0);
   V             = spm_realign(V,realign_flags);
   for i=1:N
    spm_get_space(V(i).fname,V(i).mat);
   end
 end 


 %
 % get realignment parametedfrs for comparison
 % BUT: if imgs are already resliced, we can only hope to find
 %      the params in the standard SPM text file 
 %
 
 % for our data - get RP from rp file ANYWAY. 
 if ~RESLICE
  try
   %     keyboard;
   RP=spm_get('Files',fileparts(P(1,:)),'rp_anomoco*.txt');
   fprintf('found rp file : %s\n',regexprep(RP(end,:),'\s*$',''));
   RP=load(regexprep(RP(end,:),'\s*$','')); % I wanted to nomoco RP files for this... but then it'll reslice anyway... hmm.   
   if size(RP,1)>N, RP=RP(end-N+1:end,:);end
  catch
   warning('could not determine realignment parameters')
   RP=zeros(N,6);
  end
 else %images are not resliced, but realigned, so get params from mat files
  fprintf('\n   -> retrieve realignment parameters ')
  for i=1:N
   V=spm_vol(P(i,:));
   M=V.mat;
   if ( spm_flip_analyze_images )
    M=diag([-1 1 1 1])*M;
   end
   if (i==1),M1=M; end
   rp(i,:)=spm_imatrix(M1/M);
   fprintf('%03d\b\b\b',N-i)
  end
  RP=[rp(:,1) -rp(:,2) -rp(:,3) -rp(:,4) rp(:,5) rp(:,6)];
 end
  

 %
 % read images and reslice internally, if flag is set
 %
 fprintf('\n   -> read images ')
 for j=1:N
  V=spm_vol(P(j,:));
  if (j==1)
   nx=V.dim(1); ny=V.dim(2); nz=V.dim(3);
   X=zeros(nx*ny*nz,N); F=ones(nx,ny,nz);
   if RESLICE
       disp('(internally reslicing...)');
     V1=V; d=[4 1; 4 1; 4 0];
     [x1,x2] = ndgrid(1:V1.dim(1),1:V1.dim(2));
   end
  end
  if ( RESLICE && j>1)
      % disp('(internally reslicing, II...)');
   C = spm_bsplinc(V, d);
   x=zeros(nx,ny,nz);
   for x3 = 1:nz
    [tmp,y1,y2,y3] = getmask(inv(V1.mat\V.mat),x1,x2,x3,V.dim(1:3),[0 0 0]);
    F(:,:,x3)=F(:,:,x3).*tmp;
    x(:,:,x3)=spm_bsplins(C, y1,y2,y3, d);
   end
  else
   x=spm_read_vols(V);
  end
  x(isnan(x))=0; %we don't like NaNs
  X(:,j)=x(:);
  fprintf('%03d\b\b\b',N-j)
 end

 for j=1:N; X(:,j)=X(:,j).*F(:); end %null all voxels where reslicing had problems

 %
 %compute background noise, thrsholds, and image masking 
 %
 fprintf('\n   -> background noise estimate & image masking ')
 %mean image
  MX=mean(X,2);   
  stat.mean_img=reshape(MX,nx,ny,nz);
 %track some minimal-signal corner-voxels through time-series
  dummy=zeros(nx,ny,nz);
  dummy([1:4,nx-3:nx],[1:4,ny-3:ny],[1:nz])=1;
  Ic=intersect(find(dummy),find(MX>1)); %index of corner-voxels (if mean <= 1 , these are 'strange' voxels)
  Y=[];
  for j=1:N
    [dummy,J]=sort(X(Ic,j)); % 1. sort gray values in corner voxels for the current scan
    try
    y=X(Ic(J(1:nz)),:);      % 2. nz voxels whith smallest signal for current scan => background !
    Y=[Y;y(:)];              % 3. track these voxels through all scans
    catch
        keyboard;
    end
  end
  Y(Y<0)=[]; %reslicing might produce artifical negative signal
 %noise estimate from all these background-voxels in all scans
  stat.sigma_bg=std(Y)/0.655; %Rayleigh correction
 %compute mask: it should cover at least 15 % of the nonzero elements in the mean image
  J=[];
  while (length(J)/length(find(MX)) < 0.15)  
   J=find(MX>THRESHOLD*stat.sigma_bg(1));
   if (length(J)/(nx*ny*nz) < 0.15)
    warning('can not find enough voxels above threshold ... reducing threshold!')
    THRESHOLD=THRESHOLD-5;
   end
  end
  MASK=zeros(nx,ny,nz);
  MASK(J)=1;
  stat.threshold=THRESHOLD;
  %remove eyes from mask
   if EYE_REMOVE
    try
     MASK=fcm(MASK);
    catch
     disp(' >>> could not remove eyes - check output !!!')
    end
    J=find(MASK);
   end
  stat.mask=MASK;
 %mean signal intensity
  stat.mean_intensity=mean(MX(J));
 %remove mean image
  for j=1:size(X,2)
     X(J,j)=X(J,j)-MX(J);
  end
 %reduce to mask
  X=X(J,:);


 %
 %similarity measures to standard normal distribution
 %
 fprintf('\n   -> distribution properties and similarity to gaussian')
 %q-q plot analysis
  [stat.r_qq, stat.sigma, dummy, stat.qq_plot]=qqa(X(:));
  stat.psc=100*stat.sigma/stat.mean_intensity;
 %compute some distance measures on z-scores
  stat.D=fit_normality_stats( (X(:)-mean(X(:))) / std(X(:)) ); 

 
 %
 %analysis of temporal variations
 %
 fprintf('\n   -> temporal variations')
 %qq-analyis per scan
 [stat.temp.r_qq , SIGMA, MEAN]=qqa(X);
 stat.temp.psc   = 100*SIGMA/stat.mean_intensity;
 stat.temp.pmv   = 100*MEAN/stat.mean_intensity;
 %time series of PSC per slice
 PSC_SL=ones(nz,N)*stat.psc; %set background to mean PSC for a meaningful colormap
 for k=1:nz
  Jk=find(J>(k-1)*nx*ny & J<k*nx*ny);  %the indices in J belonging to the k-th slice!
  if (length(Jk)>400) %if less than 20x20 voxels per slice, analysis is useless
   [RQQ, SIGMA, MEAN]=qqa(X(Jk,:));
   PSC_SL(k,:)=100*SIGMA/stat.mean_intensity;
  end
 end
 stat.temp.psc_sl = PSC_SL;
 %put realignment params to stat-struct
 stat.temp.rp=RP;

 %
 %plot result
 %
 if (PLOT)
  qa.name=NAME;
  qa.stat=stat;
  qa_plot(qa,PRINTPLOT);
 end

return


%---------------------------------------------------------------------
%compute distance measures of z-scores to standard normal distribution 
%---------------------------------------------------------------------
function D=fit_normality_stats(zs)
% D=fit_normality_stats(ZS);
%Some measures for testing a sample of z-scores ZS
%against standard normal distribution (mean zero, std one).
%Output D is a struct, having
%	D.N    = sample-size (length(ZS)),
%	D.D_ts = my own dist. measure (see m-file),
%	D.D_ks = Kolmogorov-Smirnov distance,
%	D.D_ku = Kuiper distance,
%	D.D_ad = Anderson-Darling distance.
 n=length(zs);
 %TS distance (my own distance: just area distance between prob.-distr.)
   Nbin=round(10*n^0.2);
   [G,X]=hist(zs,Nbin);
   g=exp(-X.^2/2);
   g=g/norm(g);
   G=G/norm(G);
   D_ts=trapz(X,abs(g-G))/trapz(X,g);
 %CDFs for other dist.-measures
   x=sort(zs);
   Sn=[1:n]/n;
   F=spm_Ncdf(x);
 % Kolmogorov-Smirnov distance
   D_ks=max(abs(Sn-F'));
   % Kuiper distance
   D_ku=max(Sn-F')+max(F'-Sn);
%Anderson-Darling distance
   J=find(F>0 & F<1);
   D_ad=length(J)^(-1)*max(abs(Sn(J)-F(J)')./sqrt(F(J)'.*(1-F(J)')));
 %all results in one structure
 D=struct('N',n,'D_ts',D_ts,'D_ks',D_ks,'D_ku',D_ku,'D_ad',D_ad);
return


%-------------------
%compute qq analysis
%-------------------
function [R,S,M,Zplot]=qqa(Z)
 %[R,S,M,Zplot]=qqa(Z)
 %Computes q-q analysis to standard normal distrbution for data vector Z.
 %Output are correlation coefficient R, standard deviation S, mean M
 %and 500 sample points Zplot for plotting. If Z is a matrix, analysis
 %is performed per column and the results are vectors.
 [n,m]=size(Z);
 N=ceil(n/500);
 q=([1:n]'-0.375)/(n+0.25);
 Zn=spm_invNcdf(q,0,1);
 for i=1:m
  z=sort(Z(:,i));
  r=corrcoef([z Zn]);
  p=polyfit(Zn,z,1);
  R(i)=r(1,2); %corelation coefficient
  S(i)=p(1);   %standard deviation (slope of q-q plot)
  M(i)=p(2);   %mean (y-axis crossing of q-q plot)
  if (nargout==4)
   Zplot(:,:,i)=[Zn(1:N:end) z(1:N:end)];
  end
 end
return



%-------------------------------------------------------------
%helper function for reslicing (copy-paste from spm_reslice.m)
%-------------------------------------------------------------
function [Mask,y1,y2,y3] = getmask(M,x1,x2,x3,dim,wrp)
tiny = 5e-2; % From spm_vol_utils.c
y1   = M(1,1)*x1+M(1,2)*x2+(M(1,3)*x3+M(1,4));
y2   = M(2,1)*x1+M(2,2)*x2+(M(2,3)*x3+M(2,4));
y3   = M(3,1)*x1+M(3,2)*x2+(M(3,3)*x3+M(3,4));
Mask = logical(ones(size(y1)));
if ~wrp(1), Mask = Mask & (y1 >= (1-tiny) & y1 <= (dim(1)+tiny)); end;
if ~wrp(2), Mask = Mask & (y2 >= (1-tiny) & y2 <= (dim(2)+tiny)); end;
if ~wrp(3), Mask = Mask & (y3 >= (1-tiny) & y3 <= (dim(3)+tiny)); end;
return;



%---------------------------------------
%eye removal (find conected mask region)
%---------------------------------------
function Mc=fcm(M_in)
%
%input: 3D-mask volume
%output: Mask containing a single connected region.

 %smooth mask
 [Nx,Ny,Nz]=size(M_in);
 %cluster mask indices
 I=find(M_in);
 [Ix,Iy,Iz]=ind2sub(size(M_in),I);
 J=spm_clusters([Ix Iy Iz]');
 N=hist(J,[1:max(I)]);        % count entries per cluster
 [dummy,j]=max(N);            % find index of maximum
 J=find(J==j);                % indices belonging to the largest cluster
 %mask of largest cluster (hopefully the brain with eyes removed)
 Mc=zeros(size(M_in));
 Mc(I(J))=1;

return


%------------------------------
%plot quality assurance results
%------------------------------
function qa_plot(QA,PRINTFLAG)
%
%qa_plot(QA,PRINTFLAG)
%
%Plot quality assurance results and write results to qa.ps.
%
%input: QA structure
%

 global defaults;
 %make spm-graph
 Fgraph   = spm_figure('FindWin','Graphics');
 if (isempty(Fgraph))
  Fgraph   = spm_figure('GetWin','Graphics');
 end

%display results 
 for ii=1:length(QA)
  IQA=QA(ii);
  for i=1:length(IQA.stat)
   figure(Fgraph)
   spm_clf;
   colormap([jet(100);gray(100)])
   h=subplot(5,1,1);
    p=get(h,'position');
    % p=p+[-0.025 0 0 0]; % move a bit to the RIGHT, plz.
    A=IQA.stat(i).temp.psc_sl;
    A(isnan(IQA.stat(i).temp.psc_sl(:)))=0;
    imagesc(A);
    CL=get(gca,'clim');
    set(gca,'clim',newclim(1,100,CL(1),CL(2),200))
    set(gca,'xticklabel',[])
    ylabel('slice','FontSize',10,'FontWeight','Bold');
    title('PSC  per Slice','FontSize',10,'FontWeight','Bold');
    Xp=get(gca,'xlim'); Xp=Xp(1)-.1*(Xp(2)-Xp(1));
    Yp=get(gca,'ylim'); Yp=Yp(1)-.1*(Yp(2)-Yp(1));
    text(Xp,Yp,[IQA.name ': run ' IQA.stat(i).run],'FontSize',12,'FontWeight','Bold');
    hc=colorbar;
    set(hc,'ylim',CL,'position',get(hc,'position')+[0.1 0 -1*.025 0]);
    axes(hc),title('%','FontSize',10,'FontWeight','Bold')
    set(h,'position',p) % restore axes - because it'll be shifted. good.
   subplot(5,1,2)
    N=length(IQA.stat(i).temp.psc);
    t=[1:N];
    plot(t,IQA.stat(i).temp.psc,'k','linewidth',1)
    set(gca,'xticklabel',[],'xlim',[1 N])
    ylabel('%','FontSize',10,'FontWeight','Bold')
    title(sprintf('PSC (mean=%4.2f%%)',IQA.stat(i).psc),'FontSize',10,'FontWeight','Bold');
   h3=subplot(5,1,3)
    [h,h1,h2]=plotyy(t,IQA.stat(i).temp.rp(:,1:3),t,IQA.stat(i).temp.rp(:,4:6));
    for j=1:3
     set(h1(j),'color',[0 0 1],'linewidth',1);
     set(h2(j),'color',[1 0 0],'linewidth',1);
    end
    set(h1(2),'linestyle','--');
    set(h2(2),'linestyle','--');
    set(h1(3),'linestyle',':');
    set(h2(3),'linestyle',':');
    axes(h(1))
     ylabel('mm','FontSize',10,'FontWeight','Bold');
     set(gca,'xlim',[1 N])
     legend('x','y','z',3)
    axes(h(2))
     xlabel('scans','FontSize',10,'FontWeight','Bold');
     ylabel('degree','FontSize',10,'FontWeight','Bold');
     set(gca,'xlim',[1 N])
     legend('p','r','y',4)
     title('SPM Realignment Parameters','FontSize',10,'FontWeight','Bold');
   subplot(5,2,7)
    plot(IQA.stat(i).qq_plot(:,1),IQA.stat(i).qq_plot(:,2),'.r',[-4 4],[-4*IQA.stat(i).sigma 4*IQA.stat(i).sigma],'--k')
    axis([-4 4 -4*IQA.stat(i).sigma 4*IQA.stat(i).sigma])
    xlabel('z-scores','FontSize',10,'FontWeight','Bold');
    ylabel('data quantiles','FontSize',10,'FontWeight','Bold');
    title(sprintf('qq-plot r_qq = %.2f',IQA.stat(i).r_qq),'FontSize',10,'FontWeight','Bold')
    set(gca,'position',get(gca,'position')+[0 -.03 0 0])
   subplot(5,2,8)
    text(.6,.9,sprintf('= %5.4f',IQA.stat(i).r_qq),'FontSize',10,'FontWeight','Bold')
    text(.6,.7,sprintf('= %4.2f (gray values)',IQA.stat(i).sigma),'FontSize',10,'FontWeight','Bold')
    text(.6,.5,sprintf('= %4.2f (gray values)',IQA.stat(i).sigma_bg),'FontSize',10,'FontWeight','Bold')
    text(.6,.3,sprintf('= %5.4f ',IQA.stat(i).D.D_ks),'FontSize',10,'FontWeight','Bold')
    text(.6,.1,sprintf('= %5.4f ',IQA.stat(i).D.D_ad),'FontSize',10,'FontWeight','Bold')
    text(-.2,.9,'CorrCoeff','FontSize',10,'FontWeight','Bold')
    text(-.2,.7,'StdDev','FontSize',10,'FontWeight','Bold')
    text(-.2,.5,'StdDev (background)','FontSize',10,'FontWeight','Bold')
    text(-.2,.3,'Kolmogorov-Smirnov D','FontSize',10,'FontWeight','Bold')
    text(-.2,.1,'Anderson-Darling D','FontSize',10,'FontWeight','Bold')
   set(gca,'position',get(gca,'position')+[0 -.03 0 0],'visible','off')
   subplot(5,2,9)
   % keyboard;
    vdami(IQA.stat(i).mean_img);
    CL=get(gca,'clim');
    set(gca,'clim',newclim(101,200,CL(1),CL(2),200))
    ylabel('mean image','FontSize',10,'FontWeight','Bold')
    set(gca,'position',get(gca,'position')+[-.05 -.08 .05 .05])
   subplot(5,2,10)
    vdami(IQA.stat(i).mask);
    CL=get(gca,'clim');
    set(gca,'clim',newclim(101,200,CL(1),CL(2),200))
    ylabel('brain mask','FontSize',10,'FontWeight','Bold')
    set(gca,'position',get(gca,'position')+[-.05 -.08 .05 .05])
   colormap([jet(100);gray(100)])

  %append footnote
  FNote = sprintf('%s, AQUA %s: %s',spm('ver'),spm('GetUser',' (%s)'),spm('time'));
  delete(findobj(gcf,'Tag','SPMprintFootnote'));
  axes('Position',[0.005,0.005,0.1,0.1],'Visible','off','Tag','SPMprintFootnote')
  text(0,0,FNote,'FontSize',6);

  %print figure to qa.ps
   if PRINTFLAG
    print -dpsc2 -painters -append -noui qa.ps;
   end

   pause(1) %time to refresh display 
  end %over runs
 end %over subjects in QA structure

return


%volume display as mosaic image
function [Mi,Mj]=vdami(X,range,s)

% leave out several slices... like factor 4...
% X=X(:,:,(1:4:end));
% remove 31st slice... for plotting purposes, only!
X(:,:,31)=[]; % otherwise you won't have a nice 5x6 fit/

[Nx,Ny,Nz]=size(X);
% keyboard;
I=[1:10];
J=Nz./I;
K=find(J==fix(J));
[dummy,i]=min(abs(J(K)-I(K)));

Mi=min([I(K(i)) J(K(i))]);
Mj=max([I(K(i)) J(K(i))]);

A=zeros(Mi*Nx,Mj*Ny);

for i=1:Mi
 for j=1:Mj
  k=i+(j-1)*Mi;
  A([1:Nx]+(i-1)*Nx,[1:Ny]+(j-1)*Ny)=X(:,:,k);
 end
end

imagesc(A);
set(gca,'xtick',[],'ytick',[])
axis image

return

%---------------------------------------------------
%helper function to compute clims for mutliple
%colormaps in one figure (see matlab documentation)
%---------------------------------------------------
function CLim = newclim(BeginSlot,EndSlot,CDmin,CDmax,CmLength)
 %                Convert slot number and range
 %                to percent of colormap
 PBeginSlot    = (BeginSlot - 1) / (CmLength - 1);
 PEndSlot      = (EndSlot - 1) / (CmLength - 1);
 PCmRange      = PEndSlot - PBeginSlot;
 %                Determine range and min and max
 %                of new CLim values
 DataRange     = CDmax - CDmin;
 ClimRange     = DataRange / PCmRange;
 NewCmin       = CDmin - (PBeginSlot * ClimRange);
 NewCmax       = CDmax + (1 - PEndSlot) * ClimRange;
 CLim          = [NewCmin,NewCmax];
return


