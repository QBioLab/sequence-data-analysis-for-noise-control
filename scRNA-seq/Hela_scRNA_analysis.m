clear all
tic
% 1) load ChIP-atlas gene lists for Cpp/P300 and NfkB/p65
[chip]=load_chip_atlas_human();
[0 toc/1e3]
% 2) load and convert three scRNA data
[sc]=list_scRNA2(1);
for i1=1:length(sc)
    % 3) load scRNMA-seq UMI count Matrix in csv
    [a,gene]=convert_scRNA_to_mat(sc(i1));
    [1 i1 toc/1e3]
    % 4) separate human and mouse cells
    [hg,geneh,ah,mg,genem,am,jh,jm]=isolate_mixed_sample(a,gene);
        HeLa(i1).ah=ah;
        HeLa(i1).gene=geneh;    
    [2 i1 toc/1e3]
end
% 5) normalized UMI counts to library size and calculate mean and CV
[HeLa]=Calcu_Hela_A2_A3_C2(HeLa);
[3 toc/1e3]
% 6) find genes with CVs significantly above the basal level 
[HeLa,Fig6C,index_outlier,cv_threashold,CPM_t]=Hela_find_CV_outlier(HeLa);
[4 toc/1e3]
% 7) statistic analysis of effects of A485+LMK on CV with outlier genes with large CVs 
[HeLa,bin, Fig6D,Fig6E,Fig6F]=Hela_bin_CVs_chip(HeLa,index_outlier,chip);
% 8) calculate histograms for UMI counts and for "SAVER-denoised' UMI counts
    list(1).name='HMGCS1'; list(2).name='HMOX1';list(3).name='mRuby-2A-bla';
[Fig6G,list]=histograms_hela_UMI(HeLa,list);
save('HeLa_A2_A3_C2_final.mat','HeLa','chip','geneh','index_outlier','cv_threashold','CPM_t','Fig6C','Fig6D','Fig6E','Fig6F','Fig6G','-v7.3');
% end of main program

function [sc]=list_scRNA2(j1)
    % sample infomations 
    sc(j1).folder='.\';
    sc(j1).source='halpos+huada'; % sequencing platforms, data conbined
    sc(j1).name='A1A2_withIndex';sc(j1).species='human_mouse';
    j1=j1+1;sc(j1).name='B1A3_withIndex';sc(j1).folder=sc(1).folder;sc(j1).source=sc(1).source;sc(j1).species='human_mouse';
    j1=j1+1;sc(j1).name='C1C2_withIndex';sc(j1).folder=sc(1).folder;sc(j1).source=sc(1).source;sc(j1).species='human_mouse';
end

function [chip]=load_chip_atlas_human()
    % downloaded from ChIP-atlas
    fileP300='EP300-human.tsv';
    fileCbp='CREBBP-human.tsv';
    fileNfkb='RELA-human.tsv';

    temp=importdata(fileP300);
    ng=length(temp.data(:,1));
    for i1=1:ng
        P300.gene{i1}=temp.textdata{i1+1,1};
    end

    temp=importdata(fileCbp);
    ng=length(temp.data(:,1));
    for i1=1:ng
        Cbp.gene{i1}=temp.textdata{i1+1,1};
    end

    temp=importdata(fileNfkb);
    ng=length(temp.data(:,1));
    for i1=1:ng
        Nfkb.gene{i1}=temp.textdata{i1+1,1};
    end

    % merge ChIP list and scores for P300 and Cbp
    g1=string(P300.gene);
    g2=string(Cbp.gene);
    g12=unique([g1 g2]);
    CbpP300.gene=cellstr(g12);

    chip(1).name='Cbp/P300';
    chip(1).gene=CbpP300.gene;
    chip(2).name='Rela (NfkB P65)';
    chip(2).gene=Nfkb.gene;
end

function [a,gene]=convert_scRNA_to_mat(sc)
    filecsv=sprintf('%s%s.csv',sc.folder,sc.name)
     % load count matrix csv files
     temp=importdata(filecsv);
     clear gene
     % generate gene lists
    for i1=1:length(temp.textdata)-1
        temp0=char(temp.textdata{1,i1+1});
        if temp0(1)=='"'
            gene{i1}=temp0(2:end-1);
        else
            gene{i1}=temp0;
        end
    end
    a=temp.data;
end

function [hg,geneh,ah,mg,genem,am,jh,jm]=isolate_mixed_sample(a,gene)
    jh=0;
    jm=0;
    for i1=1:length(gene)
        temp=char(gene{i1});
        if temp(1:2)=='hg' % human genes
            jh=jh+1;
            hg(jh)=i1; % indices for human genes
            geneh{jh}=temp(6:end);
        elseif temp(1:2)=='mm' % mouse genes
            jm=jm+1;
            mg(jm)=i1; % indices for mouse genes
            genem{jm}=temp(6:end);
        elseif temp(1:5)=='mRuby'  % mRuby-blas transcripts in HeLa cell
            jh=jh+1;
            hg(jh)=i1; % indices for mouse genes
            geneh{jh}=temp;
        end
    end
    for i1=1:length(mg)
        temp=char(gene{mg(i1)});
        genem{i1}=temp(6:end);
    end
    clear ah am
    ah=a(:,hg);
    am=a(:,mg);
    sh=sum(ah');
    sm=sum(am');
    clear ah am
    jh=find(sm<0.01*sh);
    jm=find(sh<0.018*sm);
    ah=a(jh,hg);
    am=a(jm,mg);
end

function [HeLa]=Calcu_Hela_A2_A3_C2(HeLa)
    b1h=HeLa(1).ah;
    b2h=HeLa(2).ah;
    b3h=HeLa(3).ah;
    geneh=HeLa(1).gene;
    ngh=length(geneh);

    TR1=sum(b1h')'; % total UMI counts for each cells
    TR2=sum(b2h')'; % total UMI counts for each cells
    TR3=sum(b3h')'; % total UMI counts for each cells

    % counts per million (CPM) normalized by the library size
    b1hn=1e6*b1h./repmat(TR1,1,ngh); 
    b2hn=1e6*b2h./repmat(TR2,1,ngh); 
    b3hn=1e6*b3h./repmat(TR3,1,ngh); 
    % means of normalized UMI
    m10=mean(b1hn);
    m20=mean(b2hn);
    m30=mean(b3hn);
    % CVs for normalized UMI
    cv10=std(b1hn)./m10;
    cv20=std(b2hn)./m20;
    cv30=std(b3hn)./m30;

    HeLa(1).name='A2';
    HeLa(1).ah=b1h;
    HeLa(1).ahn=b1hn;
    HeLa(1).m=m10;
    HeLa(1).cv=cv10;

    HeLa(2).name='A3';
    HeLa(2).ah=b2h;
    HeLa(2).ahn=b2hn;
    HeLa(2).m=m20;
    HeLa(2).cv=cv20;

    HeLa(3).name='C2';
    HeLa(3).ah=b3h;
    HeLa(3).ahn=b3hn;
    HeLa(3).m=m30;
    HeLa(3).cv=cv30;
end

function [HeLa,Fig6C,index_outlier,cv_threashold,CPM_t]=Hela_find_CV_outlier(HeLa)
    geneh=HeLa(1).gene;
    ngh=length(geneh);
    b1hn=HeLa(1).ahn; 
    m10=HeLa(1).m;
    cv10=HeLa(1).cv;
    b2hn=HeLa(2).ahn;
    m20=HeLa(2).m;
    cv20=HeLa(2).cv;
    b3hn=HeLa(3).ahn;
    m30=HeLa(3).m;
    cv30=HeLa(3).cv;

% % % merge two control samples A2, A3
    b12hn=[b1hn;b2hn];
    nc12=length(b12hn(:,1));
    m12=mean(b12hn);
    cv12=std(b12hn)./m12;
    HeLa(4).name='A2 & A3';
    HeLa(4).ahn=b12hn;
    HeLa(4).m=m12;
    HeLa(4).cv=cv12;
% % % filter out the one cell dramatically influnces CV for each gene
    [cv11,c_1]=filter_onecell(b1hn,cv10,ngh);
    [cv21,c_2]=filter_onecell(b2hn,cv20,ngh);
    [cv31,c_3]=filter_onecell(b3hn,cv30,ngh);
    [cv121,c_12]=filter_onecell(b12hn,cv12,ngh);
    
    HeLa(1).cv1=cv11;HeLa(2).cv1=cv21;HeLa(3).cv1=cv31;HeLa(4).cv1=cv121;
    
% filter out low expression genes 
% find the power function for baseline for CV vs mean 
% find the power function for boundery for abnormal CV vs mean 
    j1=find(m10>1 & m20>1 & m30>1); % genes with mean > 1CPM in all of the three samples

% merge three data sets
    mm=[m10(j1) m20(j1) m30(j1)];
    ccv=[cv11(j1) cv21(j1) cv31(j1)];
    
% rank the means of the genes, and group the genes into 200 per groups
% and calculate average mean and average cv
    n_group_average=200;
    [temp,j2]=sort(mm);
    ngtotal=length(j2);
        j1x=1:n_group_average:ngtotal;
        j1y=j1x+n_group_average-1; j1y(end)=length(j2);
    npart=length(j1x);
    for i1=1:npart
        mean_mm(i1)=mean(mm(j2(j1x(i1):j1y(i1)))); %  means of the means
        mean_ccv(i1)=mean(ccv(j2(j1x(i1):j1y(i1)))); %  means of the CVs
        sd_ccv(i1)=std(ccv(j2(j1x(i1):j1y(i1)))); %  SDs of the CVs
    end
    
% calculate and plot the cutoff for genes with large CVs

    % power function and starting parameters
    cv_mean_model=@(a,b,x) (a*x.^b); % fitting function
    startPoints=[10 -0.5];
    % fit mean CCs to means
    ft1=fit(mean_mm',mean_ccv',cv_mean_model,'Start', startPoints);
    % fit mean CV + 3*sd CVs to means
    ft2=fit(mean_mm',mean_ccv'+3*sd_ccv',cv_mean_model,'Start', startPoints);
    
    %
    x0=(1:1:1e5)';
    cvbase=cv_mean_model(ft1.a,ft1.b,x0);
    cvbound=cv_mean_model(ft2.a,ft2.b,x0);

    figure(1); hold off
    plot(log10(mm),ccv,'b.','MarkerSize',8);hold on
    plot(log10(x0),cvbase,'m-','lineWidth',2);
    plot(log10(x0),cvbound,'r--','lineWidth',2);
    hold off
    legend('C.V.s','baseline','cutoff')
    xlabel({'single cell expressions','log10(CPM)'})
    ylabel('C.V.');
    Fig6C.xdata=log10(mm); Fig6C.ydata=ccv;
    Fig6C.xfit=log10(x0); Fig6C.yfit=cvbase;
    Fig6C.xcutoff=log10(x0); Fig6C.ycutoff=cvbound;
    Fig6C.fit1=ft1;Fig6C.fit2=ft2;
% finding outlier genes (abover the mean_CV+3sd vs mean curve)
    % finding abover the curve, other criteria 
    % 1) minimal mean expression: 5CPM
    % 2) minimal CV: 1.0
    cv_threashold=1;
    CPM_t=5;
    % calculated CV cutoff for each genes
    cv121_t=cv_mean_model(ft2.a,ft2.b,m12(j1));
    cv31_t=cv_mean_model(ft2.a,ft2.b,m30(j1));

    % find genes with CVs above the cutoff in any of the samples
    j12t=find(cv121(j1)>cv_threashold & cv121(j1)>=cv121_t & m12(j1)>CPM_t);
    j3t=find(cv31(j1)>cv_threashold & cv31(j1)>=cv31_t & m30(j1)>CPM_t);
    jj=sort([j12t j3t]);
    % removing redundant gene list
    temp=find(diff(jj)==0);
    jj(temp)=[];    
    index_outlier=j1(jj);
end

function [HeLa,bin, Fig6D,Fig6E,Fig6F]=Hela_bin_CVs_chip(HeLa,index_outlier,chip) 
    cv11=HeLa(1).cv1;
    cv21=HeLa(2).cv1;
    cv31=HeLa(3).cv1;
    cv121=HeLa(4).cv1;
    jj0=index_outlier;
    geneh=HeLa(1).gene;

    % CV threshold for outlier genes with low and high CV
        cvmn=[0  2  ];
        cvmx=[2  100];
    % find outliear genes with low and high CVS
    for i1=1:2
        k1=find(cv121(jj0)>cvmn(i1) & cv121(jj0)<cvmx(i1));
        bin(i1).cv12=cv121(jj0(k1));
        bin(i1).cv1=cv11(jj0(k1));
        bin(i1).cv2=cv21(jj0(k1));
        bin(i1).cv3=cv31(jj0(k1));
        bin(i1).jj0=jj0(k1);
        nbin(i1)=length(k1);
    end
    % find outliear genes with high CVS in ChiP target gene lists 
    % of Cbp/P300 or P65 
    for i1=1:2
       [jj1]=find_genes_in_chip_0(geneh,bin(2).jj0,chip,i1);
        bin(i1+2).jj0=jj1;
        bin(i1+2).cv12=cv121(jj1);
        bin(i1+2).cv1=cv11(jj1);
        bin(i1+2).cv2=cv21(jj1);
        bin(i1+2).cv3=cv31(jj1);
        nbin(i1+2)=length(jj1);
    end

    figure(3);hold off;
    % t-test and plot between merged control (A2+A3) and A+L (C2)
    for i1=1:4
        % paired t-test
        [h,p(i1)]=ttest(bin(i1).cv3,bin(i1).cv12);
        mb12(i1)=mean(bin(i1).cv12);
        mb3(i1)=mean(bin(i1).cv3);
        [i1 mb12(i1) mb3(i1) p(i1)*100]
        if i1<=2
            if i1==1
                plot(bin(i1).cv12,bin(i1).cv3,'mp','MarkerSize',8);hold on
            else
                plot(bin(i1).cv12,bin(i1).cv3,'r^','MarkerSize',8);hold on
            end
            xlabel('C.V. ctrl.');
            ylabel('C.V. A+L');
        end
    end
    % t-test and plot between 2 controls A2 vs A3
    for i1=1:4
        % paired t-test
        [h,p12(i1)]=ttest(bin(i1).cv1,bin(i1).cv2);
        mb1(i1)=mean(bin(i1).cv1);
        mb2(i1)=mean(bin(i1).cv2);
        [i1 mb1(i1) mb2(i1) p12(i1)*100]
         figure(3);%subplot 122
         if i1==1
            plot(bin(i1).cv1,bin(i1).cv2,'b+','MarkerSize',8);hold on
         else 
            plot(bin(i1).cv1,bin(i1).cv2,'ko','MarkerSize',8);hold on
         end
    end
    figure(3);%subplot 121
    plot([-1 100],[-1 100],'k--','LineWidth',1);hold off
    axis([0 12 0 12]);
    title('HeLa');
    legend('low CV, ctrl. vs A+L','high CV, ctrl. vs A+L','low CV, ctrl.-1 vs -2','high CV, ctrl.-1 vs -2');
        Fig6D.x1=bin(1).cv12;Fig6D.y1=bin(1).cv3;
        Fig6D.x2=bin(2).cv12;Fig6D.y2=bin(2).cv3;
        Fig6D.x3=bin(1).cv1;Fig6D.y3=bin(1).cv2;
        Fig6D.x4=bin(2).cv1;Fig6D.y4=bin(2).cv2;

    figure(5);subplot 111
    clear x g
    x=[bin(1).cv12'; bin(1).cv3' ;bin(2).cv12'; bin(2).cv3'];
    x=[x;bin(3).cv12';bin(3).cv3';bin(4).cv12';bin(4).cv3'];
    g11=repmat({'low ctrl.'},nbin(1),1);
    g12=repmat({'low A+L'},nbin(1),1);
    g21=repmat({'high ctrl.'},nbin(2),1);
    g22=repmat({'high A+L'},nbin(2),1);
    g31=repmat({'high ctrl. CBP/P300'},nbin(3),1);
    g32=repmat({'high A+L CBP/P300'},nbin(3),1);
    g41=repmat({'high ctrl. P65'},nbin(4),1);
    g42=repmat({'high A+L P65'},nbin(4),1);
    g=[g11;g12;g21;g22;g31;g32;g41;g42];

    boxplot(x,g,'PlotStyle','traditional')
    ylabel('C.V');
    title('HeLa');
    ax=gca;
    set(ax,'XTickLabel',...
        {'low ctrl.','low A+L','high ctrl.','high A+L','high ctrl. CBP/P300','high A+L CBP/P300','high ctrl. P65','high A+L P65'},...
        'XTickLabelRotation',45);
 
    hold off
        Fig6E.g1='low ctrl.';   Fig6E.x1=bin(1).cv12;
        Fig6E.g2='low A+L';     Fig6E.x2=bin(1).cv3;
        Fig6E.g3='high ctrl.';  Fig6E.x3=bin(2).cv12;
        Fig6E.g4='high A+L';    Fig6E.x4=bin(2).cv3;
        Fig6E.g5='high ctrl. CBP/P300'; Fig6E.x5=bin(3).cv12;
        Fig6E.g6='high A+L CBP/P300';   Fig6E.x6=bin(3).cv3;
        Fig6E.g7='high ctrl. P65';      Fig6E.x7=bin(4).cv12;
        Fig6E.g8='high A+L P65';        Fig6E.x8=bin(4).cv3;

        Fig6E.p=p;% paired t-test
        Fig6E.mb12=mb12; 
        Fig6E.mb3=mb3;    

    figure(6);subplot 111
    clear x12 gm
    x12=[bin(1).cv1'; bin(1).cv2' ;bin(2).cv1'; bin(2).cv2'];
    x12=[x12;bin(3).cv1';bin(3).cv2';bin(4).cv1';bin(4).cv2'];
    gm11=repmat({'low ctrl.-1'},nbin(1),1);
    gm12=repmat({'low ctrl.-2'},nbin(1),1);
    gm21=repmat({'high ctrl.-1'},nbin(2),1);
    gm22=repmat({'high ctrl.-2'},nbin(2),1);
    gm31=repmat({'high ctrl.-1 CBP/P300'},nbin(3),1);
    gm32=repmat({'high ctrl.-2 CBP/P300'},nbin(3),1);
    gm41=repmat({'high ctrl.-1 P65'},nbin(4),1);
    gm42=repmat({'high ctrl.-2 P65'},nbin(4),1);
    gm=[gm11;gm12;gm21;gm22;gm31;gm32;gm41;gm42];

    boxplot(x12,gm,'PlotStyle','traditional')
    ylabel('C.V');
    title('HeLa');
    ax=gca;
    set(ax,'XTickLabel',...
     {'low ctrl.-1','low ctrl.-2','high ctrl.-1','high ctrl.-2','high ctrl.-1 CBP/P300','high ctrl.-2 CBP/P300','high ctrl.-1 P65','high ctrl.-2 P65'},...
        'XTickLabelRotation',45);
    hold off   
    Fig6F.g1='low ctrl.-1.';    Fig6F.x1=bin(1).cv1;
    Fig6F.g2='low ctrl.-2';     Fig6F.x2=bin(1).cv2;
    Fig6F.g3='high ctrl.-1';    Fig6F.x3=bin(2).cv1;
    Fig6F.g4='high ctrl.-2';    Fig6F.x4=bin(2).cv2;
    Fig6F.g5='high ctrl.-1 CBP/P300';   Fig6F.x5=bin(3).cv1;
    Fig6F.g6='high ctrl.-2 CBP/P300';   Fig6F.x6=bin(3).cv2;
    Fig6F.g7='high ctrl.-1 P65';        Fig6F.x7=bin(4).cv1;
    Fig6F.g8='high ctrl.-2 P65';        Fig6F.x8=bin(4).cv2;
    
    Fig6F.p=p12;% paired t-test
    Fig6F.mb1=mb1; 
    Fig6F.mb2=mb2; 
end
    
function [Fig6G,list]=histograms_hela_UMI(HeLa,list)
    geneh=string(HeLa(1).gene);
    b3hn=HeLa(3).ahn;
    b12hn=HeLa(4).ahn;
    b3h=HeLa(3).ah;
    b12h=[HeLa(1).ah;HeLa(2).ah];
    x=-0.1:0.1:4.7;
    for i1=1:length(list)
        j1=find(geneh==list(i1).name);
        list(i1).UMI_control=b12hn(:,j1)+1;
        list(i1).UMI_AL=b3hn(:,j1)+1;
        csv_name=sprintf('HeLa-A2-%s-saver.csv',list(i1).name)
        b1=importdata(csv_name);
        csv_name=sprintf('HeLa-A3-%s-saver.csv',list(i1).name)
        b2=importdata(csv_name);
        b12=[b1;b2];
        csv_name=sprintf('HeLa-C2-%s-saver.csv',list(i1).name)
        b3=importdata(csv_name);
        list(i1).sUMI_control=b12*1e6./(sum(b12h')');
        list(i1).sUMI_AL=b3*1e6./(sum(b3h')');
        list(i1).x=x;
        list(i1).y1=log10(list(i1).UMI_control);
        list(i1).h1=hist(list(i1).y1,x)/length(list(i1).y1);
        list(i1).y2=log10(list(i1).UMI_AL);
        list(i1).h2=hist(list(i1).y2,x)/length(list(i1).y2);

        list(i1).sy1=log10(list(i1).sUMI_control);
        list(i1).sh1=hist(list(i1).sy1,x)/length(list(i1).sy1);
        list(i1).sy2=log10(list(i1).sUMI_AL);
        list(i1).sh2=hist(list(i1).sy2,x)/length(list(i1).sy2);
        Fig6G(i1).x=x;
        Fig6G(i1).h1=list(i1).h1;
        Fig6G(i1).h2=list(i1).h2;
        Fig6G(i1).sh1=list(i1).sh1;
        Fig6G(i1).sh2=list(i1).sh2;
    end
end 

function [cv11,c_1]=filter_onecell(b1hn,cv10,ngh)
    cv11=cv10;
    c_1=zeros(1,ngh);% the outlier cell for each gene
    for i1=1:ngh
        [temp,k1]=sort(b1hn(:,i1));
        cv_0=cv10(i1);
        cv_1=std(temp(1:end-1))/mean(temp(1:end-1));
        cv_2=std(temp(1:end-2))/mean(temp(1:end-2));
        % if removing the highest expressed cell, 
        % drastic reduce CV twices as much as removing the second cells
        if (cv_0-cv_1)>2*abs(cv_1-cv_2)
            cv11(i1)=cv_1;
            c_1(i1)=k1(end);
        end
    end
end

function [jj1]=find_genes_in_chip_0(geneh,jj0,chip,i0)
    gene_chip=chip(i0).gene;
    nj0=length(jj0);
    jj1=[];
    j1=0;
    for i1=1:nj0
        temp1=char(geneh{jj0(i1)});
        for i2=1:length(gene_chip)
            temp2=char(gene_chip{i2});
            if length(temp1)==length(temp2)
                if temp1==temp2
                    j1=j1+1;
                    jj1(j1)=jj0(i1);
                end
            end
        end
    end
end