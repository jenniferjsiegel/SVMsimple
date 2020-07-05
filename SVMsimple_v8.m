
%%  Simple SVM to categorize neural patterns across time bins within a repeating interval
%   SET USER DEFINED FLAGS AT START OF CODE, BELOW
%       Options include multi-class and/or single SVMs, bootstrapping of
%       single SVMs by resampling train and test sets from data matrix,
%       randomization of data or labels as control
%   NOTE: Structured data ('dataX' and 'labelsY') is provided in 'sampledata' folder, 
%   if loaded will skip the first section of code. If not loaded, error will result because data file
%   structure is not present to pull data for structuring. Code is included regardless here as sample of data pull.

%   First steps: 
%   Pull interval data from each neuron's file and normalize within each interval 
%   Add to data matrix such that each column is activity from one neuron as repeating intervals (8 repetitions x 50 trials = 400 samples/neuron)
%   Each row is a time bin within a given interval, labels are repeating time bins. 

%   Second step:
%   Parse data into train and test sets (user defined)
%   Train and test SVMs using oneVone strategy, with bootstrapped resampling of train and test sets if invoked
%   Run built-in Multi-Class SVM
%   Plot results (bootstrapped single SVMs) or print to command window (Multi-Class)

warning('off'); % to disable SVM warning regarding perfectly separated classes

%% Usr define option flags: 
%   Capitalized must match data set!

pulldata=~exist('dataX'); %#ok<EXIST>
prop_trained=0.75;  % proportion of data set to use for training SVM
interval_clock=1;   % '1' for analyzing data within repeating stimulus intervals,'0' for trial ramp not yet invoked
randomize=0;        % '1' to randomize labels, '2' to randomize data bins but SVM will have a hard time
labelswap=0;        % test for inverse labeling
if pulldata==0
    numbins=size(labelY,1)/400;
elseif pulldata==1  % if pulling data from files, 
    numbins=6;      % USER DEFINE NUMBER OF BINS FOR INTERVAL
end
stim_isi=300;       % INTERVAL LENGTH IN MSEC
singleSVM=1;        % choose single class SVM, set to '1'
if singleSVM==1
    boot_iterate=1; % invoke bootstrap
    numi=100;      % num iterations for bootstrap
    CI=0.995;         % confidence level for boot, assume 2-tailed
    plotgraph=numbins;         % equals number of interval bins for correct graphing   
    % OR
    choose_bins="bin1 bin2";    % ignored if boot_iterate ==1
end
multiSVMecoc=1;     % error-correcting output, multi classifier
                    % bootstrap not invoked
MYmultiSVM=0;       % my invocation of binary class voting, currently null

%%  Pull data
%   Import data list from excel spreadsheet as 'folderpath' structure, and
%   label as 'mouseID', 'sessionID', and 'cellID'
%   This section for pulling interval data from original data files
%   Section will be skipped if data matrix is already loaded

runvars=who;

if pulldata==1
    task=0;
    numcycles=8;        % number of cycles per trial
    binsize=stim_isi/numbins;
    smoothspikes=0;     % '1' for gaus smooth, '0' for no smoothing
    normalize=1;        % normalizes each interval after data is compiled, "0" for off
    num_trials=50;
    spiketime_align=0;  % '0' for first CS-aligned, '1' for US-aligned
    % create initial interval labels based on user defined number of interval bins
    labels=strings(numbins,1);
    for labelcount=1:numbins
        makelabel=sprintf("labels(%d,1)='bin%d';",labelcount,labelcount);
        eval(makelabel);
    end
    % for normalization range
    a=0;
    b=1;

    % Define Parameters
    windowlength=stim_isi*numcycles;          
    temp_numbins=windowlength/binsize;
    if interval_clock==1                    % define window length depending on interval or trial ramp
        win_shift=stim_isi/binsize;        
    else
        win_shift=windowlength/binsize;
    end
    if spiketime_align==1
        laststim=-(235+stim_isi);            % for US aligned data set only
        firststim=-(abs(laststim)+(stim_isi*(numcycles-1)));  % for US aligned
    elseif spiketime_align==0
        laststim=stim_isi*(numcycles+2);    % to skip first 2 cycles
        firststim=(stim_isi*2);             % to skip first 2 cycles
    end

    num_cells=size(folderpath,1);

    runvars=who;

    for cell_count=1:num_cells      % for each neuron in list 
        % find and load previously processed trial data
        mouse=folderpath.mouseID(cell_count);
        session=folderpath.sessionID(cell_count);
        cellfile=sprintf('%svariables.mat',folderpath.cellID(cell_count));
        if task==1
            fileload=fullfile('/Users/jennifersiegel/aJenniWork/Data/aOddball_PFCunits',mouse,session,'SpikeResponseProfile',cellfile);
        else
            fileload=fullfile('/Users/jennifersiegel/aJenniWork/Data/aFreqExploration',mouse,session,'SpikeResponseProfile',cellfile);
        end
        load(fileload,'TrialSpikes_USon');
        load(fileload,'TrialSpikes_CSon');

        % initialize and pull spike intervals from relevant window and concatenate to session matrix
        if spiketime_align==0  % for CS-aligned, uses interval #3 - ? numcycles
            trialnames=fieldnames(TrialSpikes_CSon);
            concatspikes=zeros(num_trials,temp_numbins);
            baserate=zeros(num_trials,1);
            for trial_count=1:num_trials
               trialname=string(trialnames(trial_count));
               calltrialspikes=sprintf('trialspikes=TrialSpikes_CSon.%s*1000;',trialname);  % convert spiketimes to sec
               eval(calltrialspikes); 
               startwin=find(trialspikes>=firststim,1);
               do_it=~isempty(startwin);      % are there are spikes within trial window
               if do_it==1                    % if so, then bin using histogram count function and concatenate trial
                   stopwin=find(trialspikes<=laststim,1,'last');
                   temp=round(trialspikes(startwin:stopwin));
                   temphist=histcounts(temp,'BinLimits',[firststim,laststim],'BinMethod','integers','BinWidth',binsize);
                   if smoothspikes==1
                       [temphist,window]=smoothdata(temphist,'gaussian',100/binsize);    
                   end
                   concatspikes(trial_count:trial_count,:)=temphist; 
               end 
            end
        elseif spiketime_align==1      % align to US, uses last numcycles
            % initialize and pull spike intervals from relevant window and concatenate to session matrix
            trialnames=fieldnames(TrialSpikes_USon);
            concatspikes=zeros(num_trials,temp_numbins);
            baserate=zeros(num_trials,1);
            for trial_count=1:num_trials
               trialname=string(trialnames(trial_count));
               calltrialspikes=sprintf('trialspikes=TrialSpikes_USon.%s*1000;',trialname);
               eval(calltrialspikes); 
               trialspikes=trialspikes-firststim;
               startwin=find(trialspikes>=0,1);
               do_it=~isempty(startwin);            % are there are spikes within trial window
               if do_it==1                          % if so, then concatenate trial
                   stopwin=find(trialspikes<=windowlength,1,'last');
                   temp=round(trialspikes(startwin:stopwin));
                   temphist=histcounts(temp,'BinLimits',[1,windowlength],'BinMethod','integers','BinWidth',binsize);
                   concatspikes(trial_count:trial_count,:)=temphist; 
               end 
            end
        end

        % flip concatenated intervals to column and make single vector array
        concatspikes=transpose(concatspikes);
        concatspikes=reshape(concatspikes,[],1);
        
        % normalize each interval using min-max (0-1)
        if normalize==1
            concatsize=size(concatspikes,1);
            startindex=1;  
            while (startindex+win_shift-1)<=concatsize
                endindex=startindex+win_shift-1;
                % only min-max transform if there are neuron spikes in interval
                if sum(concatspikes(startindex:endindex,1))~=0
                    concatspikes(startindex:endindex,1)=(b-a)*((concatspikes(startindex:endindex,1)-min(concatspikes(startindex:endindex,1)))/(max(concatspikes(startindex:endindex,1))-min(concatspikes(startindex:endindex,1))))+a;
                end
                % catch if transformation results in NaN
                TF=isnan(concatspikes(startindex:startindex,1));
                if TF==1
                    concatspikes(startindex:endindex,1)=0;
                end
                % build an interval average as data check (will reflect gaussian function if data pulled correctly) 
                iavgtemp=transpose(concatspikes(startindex:endindex,1));
                if startindex>1
                    iavg_eachcell(cell_count:cell_count,:)=iavg_eachcell(cell_count:cell_count,:)+iavgtemp;
                    icount=icount+1;
                elseif startindex==1
                    iavg_eachcell(cell_count:cell_count,:)=iavgtemp;
                    icount=1;
                end
                startindex=startindex+(win_shift);
            end
            iavg_eachcell(cell_count:cell_count,:)=iavg_eachcell(cell_count:cell_count,:)/icount;
            iavg_eachcell(cell_count:cell_count,:)=(b-a)*((iavg_eachcell(cell_count:cell_count,:)-min(iavg_eachcell(cell_count:cell_count,:)))/(max(iavg_eachcell(cell_count:cell_count,:))-min(iavg_eachcell(cell_count:cell_count,:))))+a;
        elseif normalize==0
            % build an interval average as data check (will reflect gaussian function if data pulled correctly) 
            concatsize=size(concatspikes,1);
            startindex=1;  
            while (startindex+win_shift-1)<=concatsize
                endindex=startindex+win_shift-1;
                iavgtemp=transpose(concatspikes(startindex:endindex,1));
                if startindex>1
                    iavg_eachcell(cell_count:cell_count,:)=iavg_eachcell(cell_count:cell_count,:)+iavgtemp;
                    icount=icount+1;
                elseif startindex==1
                    iavg_eachcell(cell_count:cell_count,:)=iavgtemp;
                    icount=1;
                end
                startindex=startindex+(win_shift);
            end
            iavg_eachcell(cell_count:cell_count,:)=iavg_eachcell(cell_count:cell_count,:)/icount;
            iavg_eachcell(cell_count:cell_count,:)=(b-a)*((iavg_eachcell(cell_count:cell_count,:)-min(iavg_eachcell(cell_count:cell_count,:)))/(max(iavg_eachcell(cell_count:cell_count,:))-min(iavg_eachcell(cell_count:cell_count,:))))+a;
        end
        
        dataX(:,cell_count:cell_count)=concatspikes;
        clearvars('-except',runvars{:},'dataX','iavg_eachcell');
        runvars=who;
    end

    % visualize interval average data as check for correct data pull
    figure;
    imagesc(iavg_eachcell); 
    colorbar;
    colormap('jet');
    axis([1 (stim_isi/binsize) 1 num_cells])
    set(gca,'TickDir','out')
    set(gcf,'Position',[1   527   297   428])
    
    % Make time bin labels for dataX matrix for SVM
    labelY=strings(size(dataX,1),1);
    for i=1:win_shift:size(dataX,1)
        labelY(i:i-1+size(labels,1),1)=labels(:,1);
    end

    clearvars('-except',runvars{:},'labelY');
    runvars=who;
end

%% Parse data based on bin pair to test, and divide into train vs test sets
%   Note that bootstrapping will resample train and test sets to avoid
%   possible biases from data train/test split

if singleSVM==1
    if boot_iterate==0
        % Usr defined:
        choose_bins=split(choose_bins);
        binA=choose_bins(1);
        binB=choose_bins(2);

        % get correct data for single SVM
        keepbinA=strcmp(labelY,binA);
        XA=dataX(keepbinA,:);
        YA=labelY(keepbinA,:);
        keepbinB=strcmp(labelY,binB);
        XB=dataX(keepbinB,:);
        YB=labelY(keepbinB,:);
            % parse into train and test sets
            trainbin=randi(size(XA,1),round(size(XA,1)*prop_trained),1);
            train=zeros(size(XA,1),1);
            train(trainbin,1)=1;
            testbin=find(train==0);
            XAtrain=XA(trainbin,:);
            XAtest=XA(testbin,:);
            XBtrain=XB(trainbin,:);
            XBtest=XB(testbin,:);
            YAtrain=YA(trainbin,1);
            YAtest=YA(testbin,1);
            YBtrain=YB(trainbin,1);
            YBtest=YB(testbin,1);
        Xtrain=cat(1,XAtrain,XBtrain);
        Xtest=cat(1,XAtest,XBtest);
        Ytrain=cat(1,YAtrain,YBtrain);
        Ytest=cat(1,YAtest,YBtest);
           % for randomization test
           if randomize==1
             scramble=randi(size(Ytrain,1),size(Ytrain,1),1);
             Ytrain=Ytrain(scramble);
           end
           if randomize==2
             scramble=randi(size(Xtrain,1)*size(Xtrain,2),size(Xtrain,1),size(Xtrain,2));
             Xtrain=Xtrain(scramble);
           end

        % train single SVM
            SVMModel=fitcsvm(Xtrain,Ytrain);
            probSVMModel=fitPosterior(SVMModel,Xtrain,Ytrain);
            [labelsvm,PostProb]=predict(probSVMModel,Xtest);
            
            % tabulate guessing and remove svm guess strategy (split goes to Bin A)
            CertaintyScore=round(PostProb(:,2),1);
            guessing=find(CertaintyScore==0.5);
            guessrate=size(guessing,1)/size(CertaintyScore,1);
            labelsvm(guessing,1)={''};
            
            % creat output table
            T=table(Ytest,labelsvm,CertaintyScore);
            
            % tabulate actual performance
            hardvote=strcmp(T.Ytest,T.labelsvm);
            
            binA_accuracy=sum(hardvote(1:size(XAtest,1)))/size(XAtest,1);
            binB_accuracy=sum(hardvote(size(XAtest,1)+1:size(Xtest,1)))/size(XBtest,1);
            binA_avgCertainty=1-mean(T.CertaintyScore(1:size(XAtest,1)));
            binB_avgCertainty=mean(T.CertaintyScore(size(XAtest,1)+1:size(Xtest,1)));
            svmaccur=sum(hardvote(1:size(Xtest,1)))/size(Xtest,1);
            
            fprintf('%s accuracy = %.2f\n',binA,binA_accuracy)
            fprintf('%s accuracy = %.2f\n',binB,binB_accuracy)
            fprintf('%s certainty = %.2f\n',binA,binA_avgCertainty)
            fprintf('%s certainty = %.2f\n',binB,binB_avgCertainty)
            fprintf('guess rate = %.2f\n',guessrate)
            fprintf('svm accuracy = %.2f\n',svmaccur)
            
    elseif boot_iterate==1 
        % re-randomly assign train and test data to create distribution of results not dependent on train-test split
        if pulldata==0
            labels=strings(numbins,1);
            for labelcount=1:numbins
                makelabel=sprintf("labels(%d,1)='bin%d';",labelcount,labelcount);
                eval(makelabel);
            end
        end
        labels=split(labels);
        numbins=size(labels,1);
        
        % create output vectors to compile each iteration of bootstrap
        medianobs=NaN(numbins,numbins);
        CIupper=NaN(numbins,numbins);
        CIlower=NaN(numbins,numbins);
        medianguess=NaN(numbins,numbins);
        CIupper_guess=NaN(numbins,numbins);
        CIlower_guess=NaN(numbins,numbins);
        mediancertainty=NaN(numbins,numbins);
        CIupper_certainty=NaN(numbins,numbins);
        CIlower_certainty=NaN(numbins,numbins);
        mediansvmacc=NaN(numbins,numbins);
        CIupper_svmacc=NaN(numbins,numbins);
        CIlower_svmacc=NaN(numbins,numbins);
        
        % for each possible bin x bin (OvO) comparison...
        binAloopstart=2;
        for binBcount=1:numbins-1
            for binAcount=binAloopstart:numbins
                % creat output arrays to compile iterative bootstrap results
                fprintf('bin%d   bin%d\n',binBcount,binAcount)
                i_outputB=zeros(numi,1);
                i_outputA=zeros(numi,1);
                i_guessrate=zeros(numi,1);
                i_Acertainty=zeros(numi,1);
                i_Bcertainty=zeros(numi,1);
                i_svmaccur=zeros(numi,1);
                
                for i=1:numi    % for each bootstrap iteration
                    %fprintf('i %d\n',i)
                    binA=sprintf('bin%d',binAcount);    % context is neg example for SVM
                    binB=sprintf('bin%d',binBcount);    % context is pos example for SVM

                    % get correct data for single SVM based on bins to be compared
                    keepbinA=strcmp(labelY,binA);
                    XA=dataX(keepbinA,:);
                    YA=labelY(keepbinA,:);
                    keepbinB=strcmp(labelY,binB);
                    XB=dataX(keepbinB,:);
                    YB=labelY(keepbinB,:);
                        % parse into train and test sets
                        trainbin=randi(size(XA,1),round(size(XA,1)*prop_trained),1);
                        train=zeros(size(XA,1),1);
                        train(trainbin,1)=1;
                        testbin=find(train==0);
                        XAtrain=XA(trainbin,:);
                        XAtest=XA(testbin,:);
                        XBtrain=XB(trainbin,:);
                        XBtest=XB(testbin,:);
                        YAtrain=YA(trainbin,1);
                        YAtest=YA(testbin,1);
                        YBtrain=YB(trainbin,1);
                        YBtest=YB(testbin,1);
                    Xtrain=cat(1,XAtrain,XBtrain);
                    Xtest=cat(1,XAtest,XBtest);
                    Ytrain=cat(1,YAtrain,YBtrain);
                    Ytest=cat(1,YAtest,YBtest);
                          % for randomization test
                          if randomize==1
                            scramble=randi(size(Ytrain,1),size(Ytrain,1),1);
                            Ytrain=Ytrain(scramble);
                          end
                          if randomize==2
                            scramble=randi(size(Xtrain,1)*size(Xtrain,2),size(Xtrain,1),size(Xtrain,2));
                            Xtrain=Xtrain(scramble);
                          end
                          if labelswap==1
                              Ytest=cat(1,YBtest,YAtest);
                          end
                    % train single SVM
                    SVMModel=fitcsvm(Xtrain,Ytrain,'ClassNames',{binA,binB});
                    probSVMModel=fitPosterior(SVMModel,Xtrain,Ytrain);
                    [labelsvm,PostProb]=predict(probSVMModel,Xtest);
                    
                    % tabulate guessing and remove svm guess strategy
                    CertaintyScore=round(PostProb(:,2),1);
                    guessing=find(CertaintyScore==0.5);
                    guessrate=size(guessing,1)/size(CertaintyScore,1);
                    labelsvm(guessing,1)={''}; 

                    % creat output table
                    T=table(Ytest,labelsvm,CertaintyScore);

                    % tabulate actual performance
                    hardvote=strcmp(T.Ytest,T.labelsvm);

                    binA_accuracy=sum(hardvote(1:size(XAtest,1)))/size(XAtest,1);
                    binB_accuracy=sum(hardvote(size(XAtest,1)+1:size(Xtest,1)))/size(XBtest,1);
                    binA_avgCertainty=1-mean(T.CertaintyScore(1:size(XAtest,1)));
                    binB_avgCertainty=mean(T.CertaintyScore(size(XAtest,1)+1:size(Xtest,1)));
                    svmaccur=sum(hardvote(1:size(Xtest,1)))/size(Xtest,1);
                    
                    i_outputB(i,1)=binB_accuracy;
                    i_outputA(i,1)=binA_accuracy;
                    i_guessrate(i,1)=guessrate;
                    i_Acertainty(i,1)=binA_avgCertainty;
                    i_Bcertainty(i,1)=binB_avgCertainty;
                    i_svmaccur(i,1)=svmaccur;
                      
                end
                
                % determine confidence intervals
                CIlo=round(numi*(1-CI));
                CIup=round(numi*CI);
                
                % calculate median and CI limits for bin B
                i_outputB=sort(i_outputB);
                medianobs(binBcount,binAcount)=median(i_outputB,'omitnan');
                CIupper(binBcount,binAcount)=i_outputB(CIup,1);
                CIlower(binBcount,binAcount)=i_outputB(CIlo,1);
                
                % calculate median and CI limits for bin A
                i_outputA=sort(i_outputA);
                medianobs(binAcount,binBcount)=median(i_outputA,'omitnan');
                CIupper(binAcount,binBcount)=i_outputA(CIup,1);
                CIlower(binAcount,binBcount)=i_outputA(CIlo,1);
                
                % for SVM guess-rate
                i_guessrate=sort(i_guessrate);
                medianguess(binAcount,binBcount)=median(i_guessrate,'omitnan');
                CIupper_guess(binAcount,binBcount)=i_guessrate(CIup,1);
                CIlower_guess(binAcount,binBcount)=i_guessrate(CIlo,1);
                
                % for certainty measures
                i_Acertainty=sort(i_Acertainty);
                mediancertainty(binAcount,binBcount)=median(i_Acertainty,'omitnan');
                CIupper_certainty(binAcount,binBcount)=i_Acertainty(CIup,1);
                CIlower_certainty(binAcount,binBcount)=i_Acertainty(CIlo,1);
                
                i_Bcertainty=sort(i_Bcertainty);
                mediancertainty(binBcount,binAcount)=median(i_Bcertainty,'omitnan');
                CIupper_certainty(binBcount,binAcount)=i_Bcertainty(CIup,1);
                CIlower_certainty(binBcount,binAcount)=i_Bcertainty(CIlo,1);
                
                % for overall accuracy of SVM independent of train/test split
                i_svmaccur=sort(i_svmaccur);
                mediansvmacc(binAcount,binBcount)=median(i_svmaccur,'omitnan');
                CIupper_svmacc(binAcount,binBcount)=i_svmaccur(CIup,1);
                CIlower_svmacc(binAcount,binBcount)=i_svmaccur(CIlo,1);
                %return
            end
            binAloopstart=binAloopstart+1;
        end
        
        % Plot bootstrapped results depending on number of interval bins
        if plotgraph==3
            xbin=1:3;
            zero=0:plotgraph+1;
            zero(:)=0.5;

            figure; plot(0:plotgraph+1,zero,'--k');
            hold on;
            p1=plot(xbin(1),mediansvmacc(2,1),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin(1),CIupper_svmacc(2,1),'+k','LineWidth',1);
            plot(xbin(1),CIlower_svmacc(2,1),'+k','LineWidth',1);
            p2=plot(xbin(2),mediansvmacc(3,1),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin(2),CIupper_svmacc(3,1),'+k','LineWidth',1);
            plot(xbin(2),CIlower_svmacc(3,1),'+k','LineWidth',1);
            p3=plot(xbin(3),mediansvmacc(3,2),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin(3),CIupper_svmacc(3,2),'+k','LineWidth',1);
            plot(xbin(3),CIlower_svmacc(3,2),'+k','LineWidth',1);

            axis([0 plotgraph+1 0 1.1]);
            xticks(0:plotgraph+1);
            xticklabels({'' 'Bin1vs2' 'Bin1vs3' 'Bin2vs3' ''});
            xlabel('Bin Comparison');
            ylabel('Accuracy');
            set(gca,'TickDir','out');
            %legend([p1 p2 p3],{'vs Bin1','vs Bin2','vs Bin3'});
            title(['Single SVM OvO comparisons, ',num2str(stim_isi),' ISI, binsize = ',num2str(stim_isi/numbins),' msec'])
            hold off;
        end
        
        if plotgraph==4
            xbin=1:6;
            zero=0:size(xbin,2)+1;
            zero(:)=0.5;

            figure; plot(1:size(zero,2),zero,'--k');
            hold on;
            plot(xbin(1),mediansvmacc(2,1),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin(1),CIupper_svmacc(2,1),'+k','LineWidth',1);
            plot(xbin(1),CIlower_svmacc(2,1),'+k','LineWidth',1);
            plot(xbin(2),mediansvmacc(3,1),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin(2),CIupper_svmacc(3,1),'+k','LineWidth',1);
            plot(xbin(2),CIlower_svmacc(3,1),'+k','LineWidth',1);
            plot(xbin(3),mediansvmacc(4,1),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin(3),CIupper_svmacc(4,1),'+k','LineWidth',1);
            plot(xbin(3),CIlower_svmacc(4,1),'+k','LineWidth',1);
            plot(xbin(4),mediansvmacc(3,2),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin(4),CIupper_svmacc(3,2),'+k','LineWidth',1);
            plot(xbin(4),CIlower_svmacc(3,2),'+k','LineWidth',1);
            plot(xbin(5),mediansvmacc(4,2),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin(5),CIupper_svmacc(4,2),'+k','LineWidth',1);
            plot(xbin(5),CIlower_svmacc(4,2),'+k','LineWidth',1);
            plot(xbin(6),mediansvmacc(4,3),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin(6),CIupper_svmacc(4,3),'+k','LineWidth',1);
            plot(xbin(6),CIlower_svmacc(4,3),'+k','LineWidth',1);

            axis([0 plotgraph+1 0 1.1]);
            set(gca,'TickDir','out');
            xticklabels({'' 'Bin1 vs' 'Bin2 vs' 'Bin3 vs' 'Bin4 vs' ''});
            xlabel('Bin Label');
            ylabel('Accuracy');
            legend([p1 p2 p3 p4],{'vs Bin1','vs Bin2','vs Bin3','vs Bin4'})
            title(['Single SVM OvO comparisons, ',num2str(stim_isi),' ISI, binsize = ',num2str(stim_isi/numbins),' msec'])
            hold off;
        end


        if plotgraph==6
            xbin1=[nan 2 3 4 5 6];
            xbin2=[1 nan 3 4 5 6];
            xbin3=[1 2 nan 4 5 6];
            xbin4=[1 2 3 nan 5 6];
            xbin5=[1 2 3 4 nan 6];
            xbin=0:plotgraph+1;
            zero=xbin;
            zero(:)=0.5;

            figure; plot(xbin,zero,'--k');
            hold on;
            p1=plot(xbin1,medianobs(1:1,:),'or','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r');
            plot(xbin1,CIupper(1:1,:),'+r','LineWidth',1);
            plot(xbin1,CIlower(1:1,:),'+r','LineWidth',1);
            p2=plot(xbin2,medianobs(2:2,:),'ok','Color',[1 0.65 0],'LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[1 0.65 0]);
            plot(xbin2,CIupper(2:2,:),'+k','Color',[1 0.65 0],'LineWidth',1);
            plot(xbin2,CIlower(2:2,:),'+k','Color',[1 0.65 0],'LineWidth',1);
            p3=plot(xbin3,medianobs(3:3,:),'og','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','g');
            plot(xbin3,CIupper(3:3,:),'+g','LineWidth',1);
            plot(xbin3,CIlower(3:3,:),'+g','LineWidth',1);
            p4=plot(xbin4,medianobs(4:4,:),'ob','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','b');
            plot(xbin4,CIupper(4:4,:),'+b','LineWidth',1);
            plot(xbin4,CIlower(4:4,:),'+b','LineWidth',1);
            p5=plot(xbin5,medianobs(5:5,:),'ok','Color',[0.5 0.25 0.5],'LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[0.5 0.25 0.5]);
            plot(xbin5,CIupper(5:5,:),'+k','Color',[0.5 0.25 0.5],'LineWidth',1);
            plot(xbin5,CIlower(5:5,:),'+k','Color',[0.5 0.25 0.5],'LineWidth',1);

            axis([0 plotgraph+1 0 1.1]);
            set(gca,'TickDir','out');
            xticklabels({'' 'Bin1 vs' 'Bin2 vs' 'Bin3 vs' 'Bin4 vs' 'Bin5 vs' 'Bin6 vs' ''});
            xlabel('Bin Label');
            ylabel('Accuracy');
            legend([p1 p2 p3 p4 p5],{'vs Bin1','vs Bin2','vs Bin3','vs Bin4','vs Bin5'})
            title(['Single SVM OvO comparisons, ',num2str(stim_isi),' ISI, binsize = ',num2str(stim_isi/numbins),' msec'])
            hold off;
        end
        
        if plotgraph==8
            xbin1=[nan 2 3 4 5 6 7 8];
            xbin2=[1 nan 3 4 5 6 7 8];
            xbin3=[1 2 nan 4 5 6 7 8];
            xbin4=[1 2 3 nan 5 6 7 8];
            xbin5=[1 2 3 4 nan 6 7 8];
            xbin6=[1 2 3 4 5 nan 7 8];
            xbin7=[1 2 3 4 5 6 nan 8];
            xbin=0:9;
            zero=xbin;
            zero(:)=0.5;

            figure; plot(xbin,zero,'--k');
            hold on;
            plot(xbin1,medianobs(1:1,:),'or','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','r');
            plot(xbin1,CIupper(1:1,:),'+r','LineWidth',1);
            plot(xbin1,CIlower(1:1,:),'+r','LineWidth',1);
            plot(xbin2,medianobs(2:2,:),'ok','Color',[1 0.65 0],'LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[1 0.65 0]);
            plot(xbin2,CIupper(2:2,:),'+k','Color',[1 0.65 0],'LineWidth',1);
            plot(xbin2,CIlower(2:2,:),'+k','Color',[1 0.65 0],'LineWidth',1);
            plot(xbin3,medianobs(3:3,:),'og','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','g');
            plot(xbin3,CIupper(3:3,:),'+g','LineWidth',1);
            plot(xbin3,CIlower(3:3,:),'+g','LineWidth',1);
            plot(xbin4,medianobs(4:4,:),'ob','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','b');
            plot(xbin4,CIupper(4:4,:),'+b','LineWidth',1);
            plot(xbin4,CIlower(4:4,:),'+b','LineWidth',1);
            plot(xbin5,medianobs(5:5,:),'ok','Color',[0.5 0.25 0.5],'LineWidth',1,'MarkerSize',10,'MarkerFaceColor',[0.5 0.25 0.5]);
            plot(xbin5,CIupper(5:5,:),'+k','Color',[0.5 0.25 0.5],'LineWidth',1);
            plot(xbin5,CIlower(5:5,:),'+k','Color',[0.5 0.25 0.5],'LineWidth',1);
            plot(xbin6,medianobs(6:6,:),'om','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','m');
            plot(xbin6,CIupper(6:6,:),'+m','LineWidth',1);
            plot(xbin6,CIlower(6:6,:),'+m','LineWidth',1);
            plot(xbin7,medianobs(7:7,:),'ok','LineWidth',1,'MarkerSize',10,'MarkerFaceColor','k');
            plot(xbin7,CIupper(7:7,:),'+k','LineWidth',1);
            plot(xbin7,CIlower(7:7,:),'+k','LineWidth',1);

            axis([0 plotgraph+1 0 1.1]);
            set(gca,'TickDir','out');
            xticklabels({'' 'Bin1 vs' 'Bin2 vs' 'Bin3 vs' 'Bin4 vs' 'Bin5 vs' 'Bin6 vs' 'Bin7 vs' 'Bin8 vs' ''});
            xlabel('Bin Label');
            ylabel('Accuracy');
            legend([p1 p2 p3 p4 p5 p6 p7],{'vs Bin1','vs Bin2','vs Bin3','vs Bin4','vs Bin5','vs Bin6','vs Bin7'})
            title(['Single SVM OvO comparisons, ',num2str(stim_isi),' ISI, binsize = ',num2str(stim_isi/numbins),' msec'])
            hold off;
        end
    end
end

% train multi-class one-vs-one, bootstrap not invoked here
if multiSVMecoc==1
    % get data and parse into train and test sets
    X=dataX;
    Y=labelY;
        trainbin=randi(size(X,1),round(size(X,1)*prop_trained),1);
        train=zeros(size(X,1),1);
        train(trainbin,1)=1;
        testbin=find(train==0);
        Xtrain=X(trainbin,:);
        Xtest=X(testbin,:);
        Ytrain=Y(trainbin,1);
        Ytest=Y(testbin,1);
    % train multi-class model, print results to command window
    Mdl = fitcecoc(Xtrain,Ytrain);
    error = resubLoss(Mdl)          %#ok
    CVMdl = crossval(Mdl);          
    genError = kfoldLoss(CVMdl)     %#ok
    
end

if plotgraph==3
    f=msgbox({'The multi-class SVM suggests good performance with a low error rate and little information loss. (See Cmd Window result)';'The single SVM analysis (see Plot) shows that the neural data can be used to decode interval time for all interval time bins with a high level of accuracy.';'For different results, try running the 6-bin data!'},'RESULTS');
elseif plotgraph==6
    f=msgbox({'The multi-class SVM indicates moderate performance and a notable level of information loss. (See Cmd Window result)';'The single SVM analysis (see Plot) reveals that as interval time proceeds, later time bins are not discriminated locally above chance levels, while early interval time bins are still discriminated well. The moderate SVM performance is specifically due to a decreased ability to discrimate momentary time late in the interval, and not moderate performance over all time bins.'});
end
    
warning('on'); % to re-enable warnings 
    
