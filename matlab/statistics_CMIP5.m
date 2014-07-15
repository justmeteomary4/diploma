clear; clc
indir = 'files';
goutdir = 'graphs/graphs_matlab/graphs_matlab_pr-evspsbl';
nyr = 1000;
format = '%f';
model = {'BCC','CCSM4','CSIRO','GISS_121','GISS_122','GISS_123','GISS_124','GISS_125',...
              'GISS_126','GISS_127','GISS_128','GISS_1221','HadCM3','MIROC','MPI','MRI'};
%model = {'BCC','GISS_121','GISS_122','GISS_123','GISS_124','GISS_125','GISS_126','GISS_127',...
%              'GISS_128','GISS_1221','MIROC','MPI','MRI'};
maparea = '1'; % 1 - Volgabasin, 2 - CaspianSea
gname = 'mrros_volume'; % pr-evspsbl

vname1 = 'mrros_volume';
for nm = 1:size(model(:))
 if strcmp(vname1,'mrros_accano') == 1
  fname{nm} = strcat(indir,'/','f_',vname1,'_',maparea,'_',model{nm},'.dat');
  outdir = 'graphs/graphs_matlab/graphs_matlab_mrros_accano';
 end
 if strcmp(vname1,'mrros_volume') == 1
  fname{nm} = strcat(indir,'/','f_',vname1,'_',maparea,'_',model{nm},'.dat');
  outdir = 'graphs/graphs_matlab/graphs_matlab_mrros_volume';
 end
 if strcmp(vname1,'evspsbl') == 1
  fname{nm} = strcat(indir,'/','f_',vname1,'_2_',model{nm},'.dat');
  outdir = 'graphs/graphs_matlab/graphs_matlab_evspsbl';
 end
 if strcmp(vname1,'pr') == 1
  fname{nm} = strcat(indir,'/','f_',vname1,'_2_',model{nm},'.dat');
  outdir = 'graphs/graphs_matlab/graphs_matlab_pr';
 end
  fileexist = exist(strcat('./',fname{nm}),'file');
  if (fileexist == 2)
   fid(nm) = fopen(fname{nm},'r');
   arr1(nm) = textscan(fid(nm),format);
   fclose(fid(nm));
    for yr = 1:nyr
     arrcell1{nm,yr} = arr1{nm}(yr);
    end
  else
    for yr = 1:nyr
     arrcell1{nm,yr} = NaN;
    end
  end
end
arrmat1 = cell2mat(arrcell1);

%vname2 = 'evspsbl';
% for nm = 1:size(model(:))
%  if strcmp(vname2,'mrros_accano') == 1
%   fname{nm} = strcat(indir,'/','f_',vname2,'_',maparea,'_',model{nm},'.dat');
%   outdir = 'graphs/graphs_matlab/graphs_matlab_mrros_accano';
%  end
%  if strcmp(vname2,'mrros_volume') == 1
%   fname{nm} = strcat(indir,'/','f_',vname2,'_',maparea,'_',model{nm},'.dat');
%   outdir = 'graphs/graphs_matlab/graphs_matlab_mrros_volume';
%  end
%  if strcmp(vname2,'evspsbl') == 1
%   fname{nm} = strcat(indir,'/','f_',vname2,'_2_',model{nm},'.dat');
%   outdir = 'graphs/graphs_matlab/graphs_matlab_evspsbl';
%  end
%  if strcmp(vname2,'pr') == 1
%   fname{nm} = strcat(indir,'/','f_',vname2,'_2_',model{nm},'.dat');
%   outdir = 'graphs/graphs_matlab/graphs_matlab_pr';
%  end
%   fileexist = exist(strcat('./',fname{nm}),'file');
%   if (fileexist == 2)
%    fid(nm) = fopen(fname{nm},'r');
%    arr2(nm) = textscan(fid(nm),format);
%    fclose(fid(nm));
%     for yr = 1:nyr
%      arrcell2{nm,yr} = arr2{nm}(yr);
%     end
%   else
%     for yr = 1:nyr
%      arrcell2{nm,yr} = NaN;
%     end
%   end
% end
% arrmat2 = cell2mat(arrcell2);

for nm = 1:size(model(:))
garrmat(nm,:) = arrmat1(nm,:); % arrmat1(nm,:) - arrmat2(nm,:);
end

for nm = 1:size(model(:))
 plot(garrmat(nm,:));
 imgname = strcat(outdir,'/','timeseries_',gname,'_',model{nm},'.png');
 title(strcat(model{nm},' time series'),'Interpreter','none');
 xlabel('Time, years');
 ylabel('km3');
 ymin = min(garrmat(:));
 ymax = max(garrmat(:));
 ylim([ymin ymax]);
 set(gcf,'visible','off')
 %print(gcf,'-dpng','-r300',imgname);
 
 %General estimates
 vmin(nm) = min(garrmat(nm,:));
 vmax(nm) = max(garrmat(nm,:));
 vdelta(nm) = vmax(nm) - vmin(nm);
 vmean(nm) = mean(garrmat(nm,:),2); %среднее (мат. ожидание)
 vmedian(nm) = median(garrmat(nm,:),2); %медиана
 vvar(nm) = var(garrmat(nm,:),0); %дисперсия (std)2
 vstd(nm) = std(garrmat(nm,:),0); %стандартное (среднеквадратическое) отклонение
 vcov(nm) = vvar(nm)/vmean(nm); %коэффициент вариации
 
 %ECDF
 [f,x] = ecdf(garrmat(nm,:));
 plot(x,f)
 imgname = strcat(outdir,'/','ecdf_',gname,'_',model{nm},'.png');
 title(strcat(model{nm},' empirical cumulative distribution function'),'Interpreter','none');
 ylabel('ECDF');
 set(gcf,'visible','off')
 %print(gcf,'-dpng','-r300',imgname);
 
 %Histfit with frequency of occurrence as Yaxis
 h = histfit(garrmat(nm,:),[],'normal');
 imgname = strcat(outdir,'/','histfit_',gname,'_',model{nm},'.png');
 title(strcat(model{nm},' histogram with a distribution fit'),'Interpreter','none');
 ylabel('PDF');
 set(h(1),'FaceColor',[.8 .8 1]);
 set(gcf,'visible','off')
 %print(gcf,'-dpng','-r300',imgname);
 
 %Normplot: check whether the data could come from a normal distribution
 n = normplot(garrmat(nm,:));
 imgname = strcat(outdir,'/','normplot_',gname,'_',model{nm},'.png');
 title(strcat(model{nm},' normal probability plot'),'Interpreter','none');
 set(gcf,'visible','off')
 %print(gcf,'-dpng','-r300',imgname);
 
 %Tests
 testchi2(nm) = chi2gof(garrmat(nm,:),'Alpha',0.05);
 testks(nm) = kstest(garrmat(nm,:));
 testlillie(nm) = lillietest(garrmat(nm,:)); %stationarity
 testlmctest(nm) = lmctest(garrmat(nm,:));

% ACF
 acf(nm,:) = autocorr(garrmat(nm,:),250,0,2);
 xx = 0:250;
 plot(xx,acf(nm,:),'r','LineWidth',2); grid on
 %autocorr(garrmat(nm,:),nyr/4,0,2);
 imgname = strcat(outdir,'/','acf_',gname,'_',model{nm},'.png');
 title(strcat(model{nm},' autocorrelation function for runoff volume'),'Interpreter','none');
%runoff volume,accumulated anomalies
 ylabel('ACF');
 ylim([-1 1])
 hlineh = hline(0,'k');
 set(gcf,'visible','off')
 %print(gcf,'-dpng','-r300',imgname);
 
%Periodogram and PSD
   acfs(nm,:) = smooth(acf(nm,:),2,'lowess');
   [pxx1,f1] = periodogram(garrmat(nm,:),[],[],1);
   [pxx2,f2] = periodogram(acfs(nm,:),[],[],1);
   plot(f1,10*log10(pxx1),'k'); hold on;
   plot(f2,10*log10(pxx2),'r'); grid on; hold off;
   imgname = strcat(outdir,'/','periodogram&psd_',gname,'_',model{nm},'.png');
   title(strcat(model{nm},' periodogram and power spectral density for runoff volume'),'Interpreter','none');
   xlabel('1/year'); ylabel('dB');
   set(gcf,'visible','off')
 % print(gcf,'-dpng','-r300',imgname);
end