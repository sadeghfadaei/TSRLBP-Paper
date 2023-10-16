clc
clear all

Noise = 'noo';
mappMode = 'riu2';
histMode = 'nh';
block_size = 3;

for indDataset = 1:6
    if indDataset == 1
        database = 'Brodatz';
    elseif indDataset == 2
        database = 'VisTex';
    elseif indDataset == 3
        database = 'STex-512-Splitted';
    elseif indDataset == 4
        database = 'UIUC';
    elseif indDataset == 5
        database = 'KTH';
    else
        database = 'DTD';
    end
    
    ind1 = 0;
    
    if strcmp(database,'GTface')
        subdatabase = 'GT1';
        ImageFormat = '.jpg';
        TT = [1:1:11];
        TVisComp = 10;
    end
    if strcmp(database,'KDEFface')
        subdatabase = 'KDEF1';
        ImageFormat = '.jpg';
        TT = [1:1:11];
        TVisComp = 10;
    end
    if strcmp(database,'Yaleface')
        subdatabase = 'Yale1';
        ImageFormat = '.gif';
        TT = [1:1:11];
        TVisComp = 10;
    end
    if strcmp(database,'ORLface')
        subdatabase = 'ORL1';
        ImageFormat = '.jpg';
        TT = [1:1:11];
        TVisComp = 10;
    end
    if strcmp(database,'Brodatz')
        TT = [5:5:40];
        Subd = 'NB';
        subdatabase = 'Brodatz4480';
        ImageFormat = '.png';
        TVisComp = 10;
        Ind1 = 6; Ind2 = 5;
    end
    if strcmp(database,'VisTex')
        TT = [5:5:40];
        Subd = 'VT';
        subdatabase = 'VisTex2160';
        ImageFormat = '.png';
        TVisComp = 10;
        Ind1 = 12; Ind2 = 16;
    end
    if strcmp(database,'Corel-1k')
        subdatabase = 'Corel-1k-N';
        ImageFormat = '.jpg';
        TT = [10:10:100];
        TVisComp = 10;
        Ind1 = 10; Ind2 = 29;    
    end
    if strcmp(database,'Corel-10k')
        subdatabase = 'Corel-10k-1';
        ImageFormat = '.jpg';
        TT = [10:10:100];
        TVisComp = 10;
        Ind1 = 2; Ind2 = 8;
    end
    if strcmp(database,'STex-512-Splitted')
        subdatabase = 'STex-512-Splitted-1';
        ImageFormat = '.jpeg';
        TT = [2:2:16];
        TVisComp = 10;
        Ind1 = 6; Ind2 = 3;
    end
    if strcmp(database,'UIUC')
        TT = [3:3:24];
        Subd = 'UIUC';
        subdatabase = 'UIUC_natural';
        ImageFormat = '.jpg';
        TVisComp = 10;
        Ind1 = 6; Ind2 = 3;
    end
    if strcmp(database,'KTH')
        TT = [10:10:80];
        Subd = 'KTH';
        subdatabase = 'KTH_TIPS_Ver1';
        ImageFormat = '.png';
        TVisComp = 10;
        Ind1 = 6; Ind2 = 3;
    end
    if strcmp(database,'DTD')
        TT = [15:15:120];
        Subd = 'DTD';
        subdatabase = 'DTD_Ver1';
        ImageFormat = '.jpg';
        TVisComp = 10;
        Ind1 = 6; Ind2 = 3;
    end

    AD = ['G:\Dropbox\Data\databases\',database,'\',subdatabase];
    Nd1 = dir(AD);
    NumCat = nnz(~ismember({Nd1.name},{'.','..'})&[Nd1.isdir]);
    IMperCAT = 0;
    for i = 1:NumCat
        AD = ['G:\Dropbox\Data\databases\',database,'\',subdatabase,'\',int2str(i),'\*',ImageFormat];
        files = dir(AD);
        IMperCAT(i) = numel(files);
    end
    NumIm = sum(IMperCAT);

    for radius = 1:3

        neighbors = radius*8;
        MaxComb = neighbors/2;
        samples = neighbors;
        [table,newMax] = getmapping(samples,mappMode);
       
        for combMode = 1:MaxComb
            PP = zeros(64,length(TT));
            RR = zeros(64,length(TT));
            for Coding = 1:31
                % Feature Extraction
                Features = cell(NumIm,1);
                for i1 = 1:NumCat
                    AD = ['G:\Dropbox\Data\databases\',database,'\',subdatabase,'\',int2str(i1),'\*',ImageFormat];
                    files = dir(AD);
                    for i2 = 1:numel(files)
                        iii = sum(IMperCAT(1:(i1-1)))+i2;
                        fprintf('\n Image number = %d',iii);
                        filepath = fullfile(files(i2).folder, files(i2).name);
                        im = imread(filepath);
                        [M,N,C] = size(im);
                        if C == 1
                            im = cat(3,im,im,im);
                            C = 3;
                        end
                        if Noise == 'yes'
                            im = rgb2gray(im);
                            im = double(im)/255;
                            varIm = var(im(:));
                            im = imnoise(im, 'gaussian', 0, varIm / (10^(NdB/10)));
                        end
                        Features{iii,1} = TSRLBP(im,radius,neighbors,mappMode,histMode,Coding,combMode,table,newMax,block_size);
                    end
                end

                filename1 = ['Featurs','_',database,'_',mappMode,'_R',int2str(radius),'_','Comb',int2str(combMode),'_','Coding',int2str(Coding)];
                save(filename1,'Features')

                Precision = zeros(NumCat,max(IMperCAT),length(TT));
                Recall = zeros(NumCat,max(IMperCAT),length(TT));
                Similarity1 = zeros(1,NumIm);
                Similarityall = zeros(NumIm);

                % Enter the Image Query number between 1:NumIm
                for i1 = 1:NumCat
                    t=0;
                    for i2 = 1:IMperCAT(i1)
                        QuaryNum = sum(IMperCAT(1:(i1-1)))+i2;
                        t = t+1;
                        %Load Features of query image and all other images
                        for i = QuaryNum:NumIm
                            FirstImage = Features{QuaryNum,1};   
                            SecondImage = Features{i,1};
                            Similarity1(1,i) = SimilarityMeasure(FirstImage,SecondImage); 
                        end
                        Similarityall(QuaryNum,QuaryNum:NumIm) = Similarity1(QuaryNum:NumIm);
                        Similarityall(QuaryNum:NumIm,QuaryNum) = Similarity1(QuaryNum:NumIm)';
                        Similarity = Similarityall(QuaryNum,:);

                        %Retrive first T images
                        [MaxSimilarities,MaxSimilaritiesIndex] = sort(Similarity,'ascend');
                        %ReturnedImages = MaxSimilaritiesIndex(1:T);

                        %Compute Precision and Recall      
                        for ggggg = 1:length(TT)
                            T = TT(ggggg);
                            ReturnedImages = MaxSimilaritiesIndex(1:T);
                            QueryClass = ceil(QuaryNum/IMperCAT(i1));
                            ReturnedCorrect = 0;
                            for k = 1:length(ReturnedImages)
                                if (ReturnedImages(k) > sum(IMperCAT(1:(i1-1)))) & (ReturnedImages(k) <= sum(IMperCAT(1:(i1-1)))+IMperCAT(i1))
                                    ReturnedCorrect = ReturnedCorrect+1; 
                                end
                            end
                            Precision(i1,i2,ggggg) = ReturnedCorrect/T;
                            Recall(i1,i2,ggggg) = ReturnedCorrect/IMperCAT(i1);
                        end
                    end
                end
                for i = 1:length(TT)
                    PP(Coding+1,i) = 100*sum(sum(Precision(:,:,i)))/NumIm;
                    RR(Coding+1,i) = 100*sum(sum(Recall(:,:,i)))/NumIm;
                    PP(65-(Coding+1),i) = PP(Coding+1,i);
                    RR(65-(Coding+1),i) = RR(Coding+1,i);
                end

                fprintf('\n Dataset = %s, mappMode = %s, Comb = %d, Coding = %2d,  P = %4.2f%%,  R = %4.2f%%',database,mappMode,combMode,Coding,PP(Coding+1,4),RR(Coding+1,4));
                if Noise == 'yes'
                    filename = ['PreRec','_',database,'_',mappMode,'_R',int2str(radius),'_','Comb',int2str(combMode),'_NdB_',int2str(NdB)];
                else
                    filename = ['PreRec','_',database,'_',mappMode,'_R',int2str(radius),'_','Comb',int2str(combMode)];
                end
                save(filename,'PP','RR')
            end
        end

    end
 end