

% Code developed by Andrej Bicanski
% andrej.bicanski@gmail.com
% www.andrejbicanski.com
%
% model published in Current Biology
%
% Bicanski A, Burgess N. - A computational model of recognition 
% memory via grid cells. Current Biology, 2019, 29, 1–12. 
% DOI: 10.1016/j.cub.2019.01.077
%
% This is the main srcipt of the model with enough data to get you started.
% See readme file for instructions

% Note, with the limited number of stimuli in the repo (10), there is less
% ambiguity across stimuli. You have to scale up the number to see drops in
% performance. In the paper we used 99 stimuli. Also note, a number substaintially 
% higher than 99 stimuli will liekly require a more sophisticated visual 
% system or other changes to the model, since the interference from similar local
% features will increase.


% preparation
setupDCs = 0; % learn distance cell connections, do only once
setupGCs = 0; % make grid cells, do only once
getIMpts = 0; % once for each stimulus select targets for foveation, manually ;-(, I know
addIMpts = 0; % if you want to add more points later
makeWTS  = 0; % make weigths for first stimulus encounter (learning phase), can take long with many stimuli.

% config to simulate, set all to 0 if you just doing preparatory stuff above
cleanIMs = 0;   smallNfeat = 0;   % unperturnbed stimuli 
occluIMs = 0;   % noise occlusions
realoIMs = 0;   % realistic occlusions (other stimuli)
smallIMs = 1;   % for downscaled stimuli

ATTonly  = 0;   % no grid cell computations, select feature to attend randomly, for prosopagnosia-like symptoms, combine with cleanIMs=1
addedPTs = 0;   % include added distractors, i.e. targets for foveation, not encoded in the learning phase

occfix1  = 0;   % allow only one cosecutive fixation in the occulsion, as if inferring only one location, combine with occluIMs=1 or realoIMs=1, or famocIMs=1

savedata = 0;   % set to 1 to save data
PLT      = 1;   % set to 0 to suppress plots even for Ninst=1;

br    = 5;   % amount of blurr, better odd number because then centering blurred image on canvas works better
brs   = 3;   % approximaely the same for small IMs, 

small = 0.5;
Cfac  = 1/2;
l     = 440;
ls    = l*Cfac;

N_IM = 10;

N_PTSperIM = 9;    % start with 9 pts per image, can use fewer, or in future make it variable
if smallNfeat
    N_PTSperIM = 4;
end

if setupGCs
    DOE_prep_makeGCs %#ok<UNRCH>
else
    load DOE_GC_FRmaps GC_FRmaps
end

N_scales  = length(squeeze(GC_FRmaps(1,1,1,:)));
N_offsets = length(squeeze(GC_FRmaps(1,1,:,1)));
res       = length(squeeze(GC_FRmaps(:,1,1,1)));   % 441
N_GCs     = N_scales*N_offsets;
GCs       = reshape(GC_FRmaps,[res,res,N_offsets*N_scales]);
GCavg     = mean(mean(sum(sum(GC_FRmaps,4),3)));



%%%%%
%%%%%
%%%%%



if setupDCs
    
    N_DCs = res;                  %#ok<UNRCH> % doesn't make more sense to have higher resolution if the smallest unit is a pixel
    
    y_curr_dist = zeros(res,1);
    y_goal_dist = zeros(res,1);
    x_curr_dist = zeros(res,1);   % non-co-linear to y, here and 60 deg offset from y but we call it x for notational simplicity
    x_goal_dist = zeros(res,1);
    
    GC2DyA_wts    = zeros(N_DCs,N_GCs);   % two distance arrays A,B per axis,
    
    DyA2rdout_wts   = zeros(1,N_DCs);     % weights to readout cells
    DyB2rdout_wts   = zeros(1,N_DCs);
    Dy60A2rdout_wts = zeros(1,N_DCs);
    Dy60B2rdout_wts = zeros(1,N_DCs);
    
    % this actually assumes we already know the mapping of GCs to phases
    for i = 1:N_DCs                              % GCs with common phase project to DCs
        GC_rates        = squeeze(GCs(i,1,:));   % squeeze(GC_FRmaps(i,1,:)); % j=1 because we only do this for one column, i.e. for one axis
        GC2DyA_wts(i,:) = GC_rates;              % GC2Dy_wts = GC2Dy_wts + GC_rates*DC_rate';
    end
    
    % same as above for all positions
    AllPosGC2Dy_wts = zeros(N_DCs,N_GCs);
    AllPosGC2Dx_wts = zeros(N_DCs,N_GCs);
    for i = 1:N_DCs
        AllPosGC2Dy_wts(i,:) = AllPosGC2Dy_wts(i,:) + squeeze(sum(GCs(i,:,:),2))';     %
        AllPosGC2Dx_wts(i,:) = AllPosGC2Dx_wts(i,:) + squeeze(sum(GCs(:,i,:),1))';     %
    end

    AllPosGC2Dy_wts./repmat(max(AllPosGC2Dy_wts,[],2),[1,N_GCs]);
    AllPosGC2Dx_wts./repmat(max(AllPosGC2Dx_wts,[],2),[1,N_GCs]);

    GC2D_wts   = GC2DyA_wts;              % because of symmetry along 60 deg axes we only need one
    
    LINwts = linspace(0,res-1,res);         % template for weigths to readout cells
    
    yDcurr2yUPrdout_wts = flip(LINwts);   % weights to readout cells
    yDgoal2yUPrdout_wts = LINwts;
    yDcurr2yDOrdout_wts = LINwts;
    yDgoal2yDOrdout_wts = flip(LINwts);
    
    xDcurr2yUPrdout_wts = flip(LINwts);   % for x, "UP" means going right
    xDgoal2yUPrdout_wts = LINwts;
    xDcurr2yDOrdout_wts = LINwts;
    xDgoal2yDOrdout_wts = flip(LINwts);
    
    save DOE_DCs GC2D_wts AllPosGC2Dy_wts AllPosGC2Dx_wts yDcurr2yUPrdout_wts ...
        yDgoal2yUPrdout_wts yDcurr2yDOrdout_wts yDgoal2yDOrdout_wts xDcurr2yUPrdout_wts ...
        xDgoal2yUPrdout_wts xDcurr2yDOrdout_wts xDgoal2yDOrdout_wts
    
else
    load DOE_DCs
end



%%%%%
%%%%%
%%%%%



if getIMpts
    figure(34668)
    for i = 1:N_IM %#ok<UNRCH>
        eval(['IM = imread(''stimlib/testimage' int2str(i) '.png'');']);
        if length(size(IM)) == 3 % means it's RGB, has a 3rd dim
            IM = rgb2gray(IM);
        end
        IM = imadjust(IM);
        IM = double(IM);
        imagesc(IM);colormap(gray);axis square
        %load clown
        %image(X);colormap(gray);
        eval(['IM' int2str(i) '_PTS = nan(N_PTSperIM,2);']);
        prompt = 'export points of interest via datatips';  % this is a bit clunky, sorry
        dummy = input(prompt);                                  
        % promt user for input, user must export datatips as p structs.
        % place cursor on image, rightclick and 'export cursor to workspace' as p1, p2 etc.   
        % do all N_PTSperIM before hitting the promt to continue
        for j = 1:N_PTSperIM
            eval(['IM' int2str(i) '_PTS(j,:) = p' int2str(j) '.Position;']);
            eval(['clear p' int2str(j) '']);
        end
        eval(['save IMcoords/IM' int2str(i) '_PTS IM' int2str(i) '_PTS']);
    end
    
end
if addIMpts
    toadd = 5;
    figure(34668)
    for i = 1:N_IM %#ok<UNRCH>
        eval(['load IMcoords/IM' int2str(i) '_PTS IM' int2str(i) '_PTS']);   % load targets
        eval(['IMpts = IM' int2str(i) '_PTS;']);
        eval(['IM = imread(''stimlib/testimage' int2str(i) '.png'');']);
        if length(size(IM)) == 3 % means it's RGB, has a 3rd dim
            IM = rgb2gray(IM);
        end
        IM = imadjust(IM);
        IM = double(IM);
        imagesc(IM);colormap(gray);axis square; hold on        
        plot(IMpts(:,1),IMpts(:,2),'co','LineWidth',2,'MarkerSize',15);
        prompt = 'export points of interest via datatips';
        dummy = input(prompt);    
        % promt user for input, user must export datatips as p structs. 
        % As above, use p10, p11 etc, if you roiginally had 9
        for j = 10:10+toadd-1
            eval(['IM' int2str(i) '_PTS(j,:) = p' int2str(j) '.Position;']);
            eval(['clear p' int2str(j) '']);
        end
        eval(['save IMcoords/IM' int2str(i) '_PTS IM' int2str(i) '_PTS']);
    end
    
end



%%%%%
%%%%% 
%%%%%



lfov    = 61;
N_fov   = 61*61;   % make quadratic fovea for simplicity
N_fov_s = 31*31;   % for small stimuli, assume attention can modulate surround 
df      = 30;      % ...
df_s    = 15;

if makeWTS % bottom-up phase/learning, no vector calculations here, these weights are specific to you stimulus library

    % sensory cells, a poor man's visual system
    SC_2_PLC_wts   = zeros(N_PTSperIM*N_IM,N_fov*256);   % 256 possible gray scale values, could make this more coarse
    SC_2_PLC_wts_s = zeros(N_PTSperIM*N_IM,N_fov_s*256);
    
    PRC_2_PLC_wts = zeros(N_PTSperIM*N_IM,N_IM);
    FPC_2_PLC_wts = zeros(N_PTSperIM*N_IM,N_fov);
    GC_2_PLC_wts  = zeros(N_PTSperIM*N_IM,N_GCs);
    
    FeatDetRefVals   = zeros(N_fov,N_IM*N_PTSperIM);
    FeatDetRefVals_s = zeros(N_fov_s,N_IM*N_PTSperIM);   % for smaller stimuli
    % feature detector reference values, this encapsulates the assumption,
    % that there are neurons in the visual systme with a preferred
    % gray-scale value, i.e. a preferred stimulus. Explicit simulation
    % of these neurons could be replaced by a Gaussian kernel: exp(-(x-ref).^2)
    FeatDetRefVals_br   = zeros(N_fov,N_IM*N_PTSperIM);
    FeatDetRefVals_brs  = zeros(N_fov_s,N_IM*N_PTSperIM);

    canvasbr  = zeros(440,440);
    canvasbrs = zeros(330,330); % this is just to place the downscaled image, it will later be repositioned
    
    for i=1:N_IM
        
        eval(['load IMcoords/IM' int2str(i) '_PTS IM' int2str(i) '_PTS']);   % load targets
        eval(['IMpts = IM' int2str(i) '_PTS;']);
        
        PRCs    = zeros(N_IM,1);              % number of perirhinal/concept cells
        PRCs(i) = 1;   % corresonding perirhinal cell active while an object is sampled
        % in the BUphase this is simply stipulated, i.e. a cell is assigned
        
        eval(['IM = imread(''stimlib/testimage' int2str(i) '.png'');']);   % load image
        if length(size(IM)) == 3 % means it's RGB, has a 3rd dim
            IM = rgb2gray(IM);
        end
        IM = imadjust(IM);
        IM = double(IM);
        IMs = imresize(IM, small);      % downscaled images for size-inv
        
        intIM = integralImage(IM);
        IMbr  = integralFilter(intIM, integralKernel([1 1 br br], 1/br^2)); % blurr images, basically deterministic noise
        IMbr  = uint8(IMbr);
        IMbr  = double(IMbr);
        
        IMbr_off = ceil((length(IM(:,1))-length(IMbr(:,1)))/2);
        tmp_br   = canvasbr;
        tmp_br(IMbr_off+1:IMbr_off+length(IMbr(:,1)),IMbr_off+1:IMbr_off+length(IMbr(:,1))) = IMbr;
        IMbr = tmp_br;
        
        intIMs = integralImage(IMs);
        IMbrs  = integralFilter(intIMs, integralKernel([1 1 brs brs], 1/brs^2));
        IMbrs  = uint8(IMbrs);
        IMbrs  = double(IMbrs);
        
        IMbrs_off = ceil((length(IMs(:,1))-length(IMbrs(:,1)))/2);
        tmp_brs   = canvasbrs;
        tmp_brs(IMbrs_off+1:IMbrs_off+length(IMbrs(:,1)),IMbrs_off+1:IMbrs_off+length(IMbrs(:,1))) = IMbrs;
        IMbrs = tmp_brs;
        
        order = randperm(N_PTSperIM);
        
        for j = order                           % pick targets at random, assumed to given by Vsys
            
            FPCs  = zeros(N_fov,1);             % Foveal patch cells, re-init for each pt
            PLCs  = zeros(N_IM*N_PTSperIM,1);   % feature/patch label cells, re-init for each pt
            cGCs  = zeros(N_GCs,1);             % current location grid cells, re-init for each pt
            tGCs  = zeros(N_GCs,1);             % target location grid cells, re-init for each pt
            
            SCs   = zeros(N_fov*256,1);         % sensory cells
            
            % each position in field of view can be target or goal
            x = IMpts(j,1);
            y = IMpts(j,2);
            
            xs = round(x*Cfac);  % cut in half for half size image
            ys = round(y*Cfac);
            
            cGCs    = reshape(squeeze(GC_FRmaps(y,x,:,:)),[N_offsets*N_scales,1]);
            tGCs    = cGCs;
            PRCs(i) = 1;
            PLCs((i-1)*N_PTSperIM+j) = 1;
            
            FPCs   = reshape(IM(y-df:y+df,x-df:x+df),[N_fov,1]); % should check for max/min of IM but our pts are not that close to the borders
            FPCs_s = reshape(IMs(ys-df_s:ys+df_s,xs-df_s:xs+df_s),[N_fov_s,1]);
            
            FPCs_br  = reshape(IMbr(y-df:y+df,x-df:x+df),[N_fov,1]); % should check for max/min of IM but our pts are not that close to the borders
            FPCs_brs = reshape(IMbrs(ys-df_s:ys+df_s,xs-df_s:xs+df_s),[N_fov_s,1]);
            
            
            SCs = exp( -( ( reshape(repmat(FPCs_br,[1,256])',[256*N_fov,1])-repmat((0:255)',[N_fov,1]) )/10 ).^2 );  % response of all sensory neurons to pixel k 
            SCs(SCs<0.01) = 0;   % for sparseness
            SC_2_PLC_wts = SC_2_PLC_wts + PLCs*SCs';
            
            SCs_s = exp( -( ( reshape(repmat(FPCs_brs,[1,256])',[256*N_fov_s,1])-repmat((0:255)',[N_fov_s,1]) )/10 ).^2 );  % response of all sensory neurons to pixel k 
            SCs_s(SCs_s<0.01) = 0;   % for sparseness
            SC_2_PLC_wts_s = SC_2_PLC_wts_s + PLCs*SCs_s';
            
            
            FeatDetRefVals(:,(i-1)*N_PTSperIM+j) = FPCs;%/max(max(IM));  %max(FPCs);
            FeatDetRefVals_s(:,(i-1)*N_PTSperIM+j) = FPCs_s;%/max(max(IM));  %max(FPCs);
        
            FeatDetRefVals_br(:,(i-1)*N_PTSperIM+j) = FPCs_br;%/max(max(IM));  %max(FPCs);
            FeatDetRefVals_brs(:,(i-1)*N_PTSperIM+j) = FPCs_brs;%/max(max(IM));  %max(FPCs);
            
            % memorize/familiarize
            PRC_2_PLC_wts = PRC_2_PLC_wts + (PLCs*PRCs'*(1-rand/10));   % small random fluctuation for random selection of next target, see below
            FPC_2_PLC_wts = FPC_2_PLC_wts + PLCs*FPCs';
            GC_2_PLC_wts  = GC_2_PLC_wts + PLCs*cGCs';    %same for current and target
            
        end
        foreignPLCs = ones(N_IM*N_PTSperIM,1);
        foreignPLCs((i-1)*N_PTSperIM+1:i*N_PTSperIM) = 0;
        foreignPLCs = logical(foreignPLCs);
        PRC_2_PLC_wts(foreignPLCs,i) = 0;  
        % this last bit was to test the effect of inhibiting all labels not belonging to a given stim, 
        % but was not used in the paper
    end
    
    % normalize to 1
    norm = repmat(sum(SC_2_PLC_wts,2),[1,256*N_fov]);
    SC_2_PLC_wts = SC_2_PLC_wts./norm;
    
    norm_s = repmat(sum(SC_2_PLC_wts_s,2),[1,256*N_fov_s]);
    SC_2_PLC_wts_s = SC_2_PLC_wts_s./norm_s;
    
    SC_2_PLC_wts   =  sparse(SC_2_PLC_wts);
    SC_2_PLC_wts_s =  sparse(SC_2_PLC_wts_s);
    
    PLC_2_PRC_wts = PRC_2_PLC_wts';
    PLC_2_FPC_wts = FPC_2_PLC_wts';
    PLC_2_GC_wts  = GC_2_PLC_wts';
    PLC_2_PLC_wts = (diag(ones(N_IM*N_PTSperIM,1))-1)*0;%.1;   % *1 *100 *750 lateral inhib, but we might just hardcode winner take all
    
    % stronger lateral inhib is justified because if stimuli are locally very
    % similar we want to stick to the current hypothesis and see how it pans
    % out
    
    PRC_2_PRC_wts = (diag(ones(N_IM,1))-1)*0;%0.1; % also no inibition
    
    % as feature detectors
    FD_FPC_2_PLC_wts = ones(N_PTSperIM*N_IM,N_fov);
    
    if N_PTSperIM == 9
        save DOE_WTS_blurr_FDwide10_001 PLC_2_PRC_wts PLC_2_FPC_wts PLC_2_GC_wts PLC_2_PLC_wts PRC_2_PRC_wts ...
            PRC_2_PLC_wts FD_FPC_2_PLC_wts FD_FPC_2_PLC_wts FPC_2_PLC_wts GC_2_PLC_wts SC_2_PLC_wts SC_2_PLC_wts_s...
            FeatDetRefVals FeatDetRefVals_s FeatDetRefVals_br FeatDetRefVals_brs '-v7.3'
    else
        save DOE_WTS_blurr_FDwide10_4PTSperIM_001 PLC_2_PRC_wts PLC_2_FPC_wts PLC_2_GC_wts PLC_2_PLC_wts PRC_2_PRC_wts ...
            PRC_2_PLC_wts FD_FPC_2_PLC_wts FD_FPC_2_PLC_wts FPC_2_PLC_wts GC_2_PLC_wts SC_2_PLC_wts SC_2_PLC_wts_s...
            FeatDetRefVals FeatDetRefVals_s FeatDetRefVals_br FeatDetRefVals_brs '-v7.3'
    end
    
    
else
    
    if N_PTSperIM == 9
        load DOE_WTS_blurr_FDwide10_001
    else
        load DOE_WTS_blurr_FDwide10_4PTSperIM_001
    end
    
end



%%%%%
%%%%%
%%%%%



% now test the model
tic;
Ninst = 1;   % number of random instances to run, did a few thousand to test the failure rate reported in the paper is representative

% NOTE ON TERMINOLOGY: 
% Stimulus Identity (ID) neurons in the paper are here called PRCs
% since I was first considering they might all reside in perirhinal cortex while I
% was developing the model. Feature label cells are called PLCs (Patch label cells)

IDvTvP     = 0;   %IDvTvP  = zeros(N_IM,200,N_IM,Ninst); 
RECvP      = zeros(N_IM,Ninst);  
OtakevP    = zeros(N_IM,Ninst);
vectors    = [];
vect2ID    = [];   % will be 2D, (startPT, endPT, resetNo, ID, perm)
maxN_reset = 10;   % arbitary choice, I think humans make far fewer attempt to try to recognize someone
resetsP    = zeros(N_IM,maxN_reset+1,Ninst);

trackocc   = [];   % track occlusions
trackoccX  = [];
trackoccY  = [];
trackIMs   = [];

% key parameters
reset_thresh = 5;      % left over from an alternative reset criterion, how much ID firing is allowed in total
Nstd         = 0.8;    % 2.8 in paper - number of standard deviations for initial pruning
% to see the effects of interference from other stimuli with a smaller amount of stimuli, this number can be
% reduced. We can imagine it scaling up with the number of stimuli. I
% tested 0.8 with 10 stimuli. 

if smallNfeat
    FLC2ID   = 0.5;    % 0.4 - strength of feature label cell to ID cell connections, higher with few features (not in paper)
else
    FLC2ID   = 0.4; 
end

% Note, "Nstd" has a strong influence. A high treshold for feature label cell
% firing means fewer co-active feature label cells at a given time, which
% has the consequence that the subsequent soft-max operation assigns the
% leading feature label cell a higher fraction of overall firing. As a
% consequence fewer saccades are needed on average. We want to preserve
% some alternative hypothesis, so shouldn't be too high.

% Finally, scaling down the contribution from feature label cells to ID neurons 
% (FLC2ID) should increase the median number of saccades with the exception
% of 1-saccade recogntion. If only one FLC fires, the softmax gives us
% instant recognition.

for q = 1:Ninst % if we want to run more than one iteration
 
    if smallIMs
        N_fov = N_fov_s;
        df    = df_s;
    end

    IDvT    = zeros(N_IM,500,N_IM);    % length 30, in case we do reset search later
    REC     = zeros(N_IM,1);           % record recogntion performance
    Otake   = zeros(N_IM,1);           % see how often a correct hypothesis overtakes an intially leading false one 
     
    dec_thresh   = 0.9;                % arbitrary recognition threshold, note, this does not mean recognition occurs then, it more like a confidence level
    resets       = zeros(N_IM,maxN_reset+1);   % 2nd dim is number of reset, entry how many saccades till reset, sum across 2nd dim, total saccade number
    RESET        = 0;                          % flag

    tol = 4;   % tolerance for vector error in number of pixels,
    % supposed to be adjusted by microsaccades (4 pixels, less than 1 %)

    
    % LOOP ACROSS ALL STIMULI
    for i = 1:N_IM      
        
        no_reset  = 0;    % number of resets 
        prevSTpts = [];   % starting points
        
        eval(['load IMcoords/IM' int2str(i) '_PTS IM' int2str(i) '_PTS']);   % load targets
        eval(['IMpts = IM' int2str(i) '_PTS;']);
        IMpts_bkp = IMpts;
        
        FPCs  = zeros(N_fov,1);             % Foveal patch cells, re-init for each pt, analogous to PW
        PLCs  = zeros(N_IM*N_PTSperIM,1);   % feature/patch label cells, re-init for each pt
        cGCs  = zeros(N_GCs,1);             % current location grid cells, re-init for each pt
        tGCs  = zeros(N_GCs,1);             % target location grid cells, re-init for each pt
        PRCs  = zeros(N_IM,1);              % number of perirhinal/concept cells/Identity neurons        
        
        eval(['IM = imread(''stimlib/testimage' int2str(i) '.png'');']);
        if length(size(IM))==3; IM=rgb2gray(IM); end   % if 3rd dim exists, is RGB
        IM = imadjust(IM);
        IM = double(IM);
        
        % with occlusions have half the features visible in each picture for good comparison across stimuli
        % note, this is not always half the image, can be more, can be less
        xotmp  = sort(IMpts(1:9,1)); 
        xocc = xotmp(5)+30;              % add quarter of foveal extent
        yotmp  = sort(IMpts(1:9,2)); 
        yocc = yotmp(5)+30;
        Nunocc = 5;                      % unoclcuded 
        
        if occluIMs || realoIMs
            startPT = ceil(rand*9);
            if occluIMs
                if rand<0.5   % add occlusion
                    occ = 1;
                    IM(:,xocc:end) = round(255*rand(440,440-xocc+1));
                else
                    occ = 2;
                    IM(yocc:end,:) = round(255*rand(440-yocc+1,440));
                end
                eval(['tmp_IMpts = IM' int2str(i) '_PTS;']);
                % SUBROUTINE CALL: select a starting pt, if we have reset
                % (in loop below) we don't want the same starting pt again
                [startPT,prevSTpts] = DOE_subr_SelectStartPT(occ,startPT,tmp_IMpts,xocc,yocc,prevSTpts);
            end
            if realoIMs
                occlNO = 300+ceil(rand*9);
                eval(['occlIM = imread(''stimlib/selected_occlusions/testimage' int2str(occlNO) '.png'');']);
                if length(size(occlIM))==3; occlIM=rgb2gray(occlIM); end   % if 3rd dim exists, is RGB
                occlIM = imadjust(occlIM);
                occlIM = double(occlIM);
                if rand<0.5   % add occlusion
                    occ = 1;
                    IM(:,xocc:end) = occlIM(:,xocc:end);
                else
                    occ = 2;
                    IM(yocc:end,:) = occlIM(yocc:end,:);
                end
                eval(['tmp_IMpts = IM' int2str(i) '_PTS;']);
                % SUBROUTINE CALL: select a starting pt, if we have reset
                % (in loop below) we don't want the same starting pt again
                [startPT,prevSTpts] = DOE_subr_SelectStartPT(occ,startPT,tmp_IMpts,xocc,yocc,prevSTpts);
            end
        else
            occ = 0;
        end
        if smallIMs
            IMs     = imresize(IM, small);
            canvas  = zeros(size(IM));
            startPT = ceil(rand*9);
        end
        if cleanIMs
            startPT = ceil(rand*9);
        end
        if ATTonly && ~addedPTs
            startPT = ceil(rand*9);
        end
        if ATTonly && addedPTs
            startPT = ceil(rand*14);
        end


        noInfoC   = 0;
        count     = 1;
        ATTcount  = 1; 
        ATTperm   = randperm(9);   % always 9 points
        if addedPTs
            ATTperm = randperm(14);   % use added attention distractors
        end
        visited   = [];
        pcount    = 1;   % for plotting

        PredResCount    = 0;
        flatprior       = ones(size(PLCs));  % for upweighting of predicted sensory outcome 
        prediction_bias = flatprior;
        
        conseqocc = 0;
        
        h1 = [];   % for fancy plotting
        h2 = [];

        
        % WITHIN STIMULUS LOOP
        while ~(any(PRCs>dec_thresh))   % compute vectors until one IDn/PRC reaches decision/confidence threshold

            
            if RESET 
                disp('RESET');
                no_reset = no_reset + 1;
                if no_reset>maxN_reset
                    str = (['NO MATCH AFTER ' int2str(maxN_reset) ' RESETS.']);
                    disp(str);
                    %pause(0.25);
                    REC(i) = -1;
                    trackocc  = [trackocc occ];
                    trackoccX = [trackoccX xocc];
                    trackoccY = [trackoccY yocc];
                    break
                end
                eval(['tmp_IMpts = IM' int2str(i) '_PTS;']);
                if (length(prevSTpts)==9 && (cleanIMs || smallIMs)) || (length(prevSTpts)==5 && ( occluIMs || occfix1 || realoIMs))
                    prevSTpts = [];
                end
                % SUBROUTINE CALL: select a starting pt, if we have reset
                % we don't want the same starting pt again
                [startPT,prevSTpts] = DOE_subr_SelectStartPT(occ,startPT,tmp_IMpts,xocc,yocc,prevSTpts);
                ATTperm  = randperm(9);         % always 9 points
                if addedPTs
                    ATTperm   = randperm(14);   % use added attention distractors
                end
                count    = 1;
                ATTcount = 1;
                visited  = [];
                pcount   = pcount+1;
                resets(i,no_reset) = pcount-1;
                PRCs     = PRCs*0;             % reset identity neurons with each reset
                RESET    = 0;
                PredResCount    = 0;
                prediction_bias = flatprior;
            end

            if count == 1           % start somewhere
                xIM = IMpts(startPT,1);   % to get size invariance, the motor output must be scaled, separate GC distance from image distance
                yIM = IMpts(startPT,2);
                x   = IMpts(startPT,1);   %
                y   = IMpts(startPT,2);
                IMpts_B = IMpts;
                if smallIMs   % align grids to stimulus, i.e. global remapping
                    % placing the downscaled image correctly (anchor to 
                    % foveated feature is equivalent to remapping the grids, 
                    % this way we don't need to actually shift all grids
                    x2      = round(xIM*Cfac);
                    y2      = round(yIM*Cfac);
                    lys_up  = length(1:y2-1);
                    lys_do  = min(length(y2:ls),ls-(lys_up+1));
                    lxs_le  = length(1:x2-1);
                    lxs_ri  = min(length(x2:ls),ls-(lxs_le+1));
                    canvas(yIM-lys_up:yIM+lys_do,xIM-lxs_le:xIM+lxs_ri) = IMs;
                    IM        = canvas;
                    transPTS  = round([IMpts(:,1)-xIM IMpts(:,2)-yIM]*Cfac); % translate reference values for 1 pixel correction
                    IMpts_bkp = [transPTS(:,1)+xIM transPTS(:,2)+yIM];
                    trackIMs  = [trackIMs; lys_up lys_do lxs_le lxs_ri xIM yIM, no_reset, i];
                end
            end
 
            
            % STEP ONE: computations at current position
            
            cGCs = reshape(squeeze(GC_FRmaps(y,x,:,:)),[N_offsets*N_scales,1]);   
            % current position GCs given by foveal position (self-motion update)
            % that is, we assume correct updating of the grid cell attractor
            % based on eye-velocity, then look up GC population vector
            FPCs = reshape(IM(yIM-df:yIM+df,xIM-df:xIM+df),[N_fov,1]);  
            
            if smallIMs % sensory discrimination with gaussian tuning
                SCs      = exp( -( ( reshape(repmat(FPCs,[1,256])',[256*N_fov,1])-repmat((0:255)',[N_fov,1]) )/10 ).^2 );
                FeatComp = SC_2_PLC_wts_s*SCs;           
            else
                SCs      = exp( -( ( reshape(repmat(FPCs,[1,256])',[256*N_fov,1])-repmat((0:255)',[N_fov,1]) )/10 ).^2 );
                FeatComp = SC_2_PLC_wts*SCs;
            end

            if max(FeatComp)>(mean(FeatComp) + Nstd*std(FeatComp)) % pruning, only update hypothesis in this case, but continue saccades otherwise
                noInfoC = 0;
                % map sensory patch to a given (feature) label
                % no input from PRC to PLC at this stage, only sensory evidence should select the next PRC
                PLCs = FeatComp/max(FeatComp);  % could maybe replace this by a softmax too
                if count>1   % only after the first feature there's a prediction
                    if find(PLCs==max(PLCs))~=find(prediction_bias==2)  % && length(find(PLCs==max(PLCs)))<length(PLCs) % second conditon does not increment reset counter when no info at all gathered
                        PredResCount = PredResCount + 1;
                    end
                end
                if PredResCount>2 % reset if many false predctions. this is very effective and allows for early reset 
                    RESET = 1;
                end
                if ~RESET
                    % integrate prediction before or after pruning? before, to double weak evidence
                    PLCs = PLCs.*prediction_bias;                 % 2 is the value of the prediction bias for the single predicted PLC
                    PLCs(PLCs<(mean(PLCs)+Nstd*std(PLCs))) = 0;   % no winner take all, but some thresholding
                    PLCs(find(PLCs)) = softmax(PLCs(PLCs>0));     % fnally softmax
                    PRCs_prev = PRCs;
                    incr = FLC2ID*(PLC_2_PRC_wts*PLCs);           % increment for PRCs
                    PRCs = PRCs + incr;                           % Label cell drives perithinal/stimulus cell = hypothesis generated
                    PRCs = max(PRCs,0);                           % map negative values to zero
                    PRCs = min(PRCs,1);                           % saturation
                    if Ninst == 1 && PLT                          % look at evolving actiivty
                        figure(54); cla;
                        figure(54);bar(PRCs,'r');title('ID neurons');hold on
                        figure(54);bar(PRCs_prev,'k');
                        figure(55);bar(PLCs,'k');title('Feature Label Cells');
                    end
                    if sum(PRCs)>reset_thresh                     % optional, other reset criterion acts much earlier
                        RESET = 1;
                    end
                end
            else
                disp('not enough evidence to update');
                noInfoC = noInfoC + 1;
                if count == 1
                    RESET = 1;
                    disp('First feature not informative enough');   % this is very rare
                end
            end

            
            % STEP TWO: computations for target (* OR ** depending on type of simulation)
            
            
            %  *  STEP TWO in the attention only case, default below
            if ATTonly && ~RESET   
                if length(visited)>=N_PTSperIM-1
                    visited = [visited(end); find(PLCs==max(PLCs))];
                else
                    visited = [visited; find(PLCs==max(PLCs))];
                end
                PLCs          = PLCs*0;                           % reset and inhibit current patch cell
                WeightNoise   = 1-rand(size(PRC_2_PLC_wts));
                PLCs          = WeightNoise.*PRC_2_PLC_wts*PRCs; 
                PLCs(visited) = 0;                                % without this if no PLC is selected
                PLCs(PLCs<max(PLCs)) = 0;                         % winner take all, just for target selection
                if max(PLCs)>0
                    PLCs = PLCs/max(PLCs);
                end
                prediction_bias = flatprior+PLCs;           
                ATTcount = ATTcount + 1;            % attention counter
                nextPT   = ATTperm(ATTcount);       % randomly selected via bottom-up attention
                next_xIM = IMpts_bkp(nextPT,1);
                next_yIM = IMpts_bkp(nextPT,2);
                x_diff = next_xIM-x;
                y_diff = next_yIM-y;
                xIM_diff = x_diff;  % separate motor output from difference on grid
                yIM_diff = y_diff;
                if Ninst == 1 && PLT
                    figure((i)*1000+no_reset);   % plotting 1
                    if count == 1
                        imagesc(IM);colormap(gray);axis square; hold on
                        plot(IMpts_bkp(1:9,1),IMpts_bkp(1:9,2),'co','LineWidth',2,'MarkerSize',15); hold on
                        if addedPTs
                            plot(IMpts_bkp(10:14,1),IMpts_bkp(10:14,2),'yo','LineWidth',2,'MarkerSize',15); hold on
                        end
                    end
                    if ~isempty(h1) set(h1,'Visible','off'); end
                    if ~isempty(h2) set(h2,'Visible','off'); end
                    ofStim     = ceil(find(prediction_bias==2)/N_PTSperIM);
                    intendedPT = find(prediction_bias==2)-(ofStim-1)*N_PTSperIM;
                    if ofStim~=i
                        eval(['load ''IMcoords/IM' int2str(ofStim) '_PTS'' IM' int2str(ofStim) '_PTS']);
                        eval(['TMPpts = IM' int2str(ofStim) '_PTS;']);
                        h1 = plot(TMPpts(intendedPT,1),TMPpts(intendedPT,2),'go','LineWidth',3,'MarkerSize',15); hold on
                    else
                        h2 = plot(IMpts_bkp(intendedPT,1),IMpts_bkp(intendedPT,2),'ro','LineWidth',3,'MarkerSize',10); hold on
                    end
                    hold on;
                    len = sqrt(x_diff^2+y_diff^2);
                    quiver(xIM,yIM,xIM_diff,yIM_diff,'Color','r','LineWidth',2,'MaxHeadSize',40.0/len,'AutoScale','off'); hold on
                end
                vect2ID = [vect2ID; xIM yIM xIM_diff yIM_diff no_reset i q];   % everthing we need to redraw the saccades
                vectors = cat(3, vectors, [x x+x_diff; y y+y_diff]);
                yIM = yIM+yIM_diff;  % update postion on image (scaled version of above)
                xIM = xIM+xIM_diff;
                y   = y+y_diff;      % update postion on grid cell sheet
                x   = x+x_diff;
            end
            
            
            % STEP TWO:  **  DEFAULT
            if ~ATTonly && ~RESET && ~(any(PRCs>dec_thresh))   % first part of if is nonsense any(IDvT(:,pcount+1,i)<dec_thresh) && 
                x_diff = 0; 
                y_diff = 0;
                while abs(x_diff)<=tol && abs(y_diff)<=tol % if same GCs active, repeat selection process
                    if length(visited)>=N_PTSperIM-1
                        visited = [visited(end); find(PLCs==max(PLCs))];
                    else
                        visited = [visited; find(PLCs==max(PLCs))];
                    end 
                    PLCs          = PLCs*0;             % reset and inhibit current patch cell
                    WeightNoise   = 1-rand(size(PRC_2_PLC_wts));
                    PLCs          = WeightNoise.*PRC_2_PLC_wts*PRCs;   
                    PLCs(visited) = 0;                  % without this if no PLC is selected
                    PLCs(PLCs<max(PLCs)) = 0;           % winner take all, just for target selection
                    if occfix1 && conseqocc==0 && (occluIMs || realoIMs)
                        PredPLC = mod(find(PLCs),N_PTSperIM);   if PredPLC==0; PredPLC=9; end
                        LastPLC = mod(visited(end),N_PTSperIM);   if LastPLC==0; LastPLC=9; end
                        if occ == 1
                            % if previous fixation was within occlusion and current one is as well
                            if tmp_IMpts(LastPLC,1)>xocc && tmp_IMpts(PredPLC,1)>xocc 
                                conseqocc = 1;
                                while conseqocc==1
                                    PLCs          = PLCs*0;             % reset and inhibit current patch cell
                                    WeightNoise   = 1-rand(size(PRC_2_PLC_wts));
                                    PLCs          = WeightNoise.*PRC_2_PLC_wts*PRCs;
                                    PLCs(visited) = 0;                  % without this if no PLC is selected
                                    PLCs(PLCs<max(PLCs)) = 0;           % winner take all, just for target selection
                                    PredPLC  = mod(find(PLCs),N_PTSperIM);   if PredPLC==0; PredPLC=9; end
                                    if tmp_IMpts(PredPLC,1)<xocc 
                                        conseqocc = 0;
                                    end
                                end 
                            end
                        end
                        if occ == 2
                            if tmp_IMpts(LastPLC,2)>yocc && tmp_IMpts(PredPLC,2)>yocc %&& CurrLEAD== 
                                conseqocc = 1;
                                while conseqocc==1
                                    PLCs          = PLCs*0;             % reset and inhibit current patch cell
                                    WeightNoise   = 1-rand(size(PRC_2_PLC_wts));
                                    PLCs          = WeightNoise.*PRC_2_PLC_wts*PRCs;
                                    PLCs(visited) = 0;                  % without this if no PLC is selected
                                    PLCs(PLCs<max(PLCs)) = 0;           % winner take all, just for target selection
                                    PredPLC  = mod(find(PLCs),N_PTSperIM);   if PredPLC==0; PredPLC=9; end
                                    if tmp_IMpts(PredPLC,2)<yocc %&& CurrLEAD==
                                        conseqocc = 0;
                                    end
                                end
                            end
                        end 
                    end
                    if max(PLCs)>0
                        PLCs = PLCs/max(PLCs);
                    end
                    prediction_bias = flatprior+PLCs;
                    tGCs = max(PLC_2_GC_wts*PLCs,0);
                    % SUBROUTINE CALL: calculate vector components with
                    % distance model
                    [x_diff,y_diff] = DOE_subr_GCs2Dist(cGCs,tGCs,AllPosGC2Dy_wts,AllPosGC2Dx_wts, ...
                        yDgoal2yUPrdout_wts,yDgoal2yDOrdout_wts,xDgoal2yUPrdout_wts,xDgoal2yDOrdout_wts, ...
                        yDcurr2yUPrdout_wts,yDcurr2yDOrdout_wts,xDcurr2yUPrdout_wts,xDcurr2yDOrdout_wts);
                end
                xIM_diff = x_diff;  % separate motor output from difference on grid
                yIM_diff = y_diff;
                if smallIMs
                    xIM_diff = round(x_diff*Cfac);  % scale motor gain
                    yIM_diff = round(y_diff*Cfac);
                end
                if Ninst == 1 && PLT
                    figure((i)*1000+no_reset);   % plotting 1
                    if count == 1
                        imagesc(IM);colormap(gray);axis square; hold on
                        plot(IMpts_bkp(1:9,1),IMpts_bkp(1:9,2),'co','LineWidth',2,'MarkerSize',15/(smallIMs+1));
                    end
                    hold on;
                    len = sqrt(x_diff^2+y_diff^2);
                    quiver(xIM,yIM,xIM_diff,yIM_diff,'Color','r','LineWidth',2,'MaxHeadSize',40.0/len,'AutoScale','off'); box off; axis off
                end
                vect2ID = [vect2ID; xIM yIM xIM_diff yIM_diff no_reset i q];   % everthing we need to redraw the saccades
                vectors = cat(3, vectors, [x x+x_diff; y y+y_diff]);
            end
            
            
            if ~RESET
                IDvT(:,pcount+1,i) = min(1,max(0,PRCs));   % plotting, start at 2 to have leading zeros for plot
            end
            if Ninst == 1 && PLT
                if (any(PRCs>dec_thresh))
                    ticklocs   = [];
                    ticklabels = [];
                    if any(resets(i,:))
                        for k = 1:length(find(resets(i,:)))
                            if k == 1
                                a = 1;
                                b = resets(i,k);
                                ticklocs   = [ticklocs a:b];
                                ticklabels = [ticklabels; (a:b)'-1];
                                IDtmp      = IDvT(:,a:b,i);
                                set(gca,'ColorOrderIndex',1)
                                figure(i*10000);plot(1:b,IDtmp','LineWidth',4);colormap(lines(N_IM)); hold on
                                quiver(b+0.2,0.1,0,0.2,'vk','LineWidth',3);
                            else
                                a  = b;
                                b  = resets(i,k)-(k-1);
                                ticklocs   = [ticklocs ((a+1):b)];
                                ticklabels = [ticklabels; (((a+1):b)-a)'];
                                IDtmp      = IDvT(:,a+(k-1):b+(k-1),i);
                                set(gca,'ColorOrderIndex',1)
                                figure(i*10000);plot(a:b,IDtmp','LineWidth',4);colormap(lines(N_IM)); hold on
                                quiver(b+0.2,0.1,0,0.2,'vk','LineWidth',3);
                            end 
                        end
                        b = max(ticklocs)+1;
                        ticklocs   = [ticklocs (b:pcount)];
                        ticklabels = [ticklabels; (1:length(b:pcount))'];
                        IDtmp      = IDvT(:,b+(k-1):pcount+1,i);
                        set(gca,'ColorOrderIndex',1)
                        figure(i*10000);hold on
                        p = plot((b-1):pcount-(k-1),IDtmp','LineWidth',4);colormap(lines(N_IM));
                        legend([p(i)],'Stimulus','Location','northwest');
                        xticks(ticklocs);
                        xticklabels(num2str(ticklabels));
                        axis square; box on
                        pcount = pcount - 1;
                    else
                        set(gca,'ColorOrderIndex',1)
                        figure(i*10000); hold on
                        p = plot(IDvT(:,1:pcount+1,i)','LineWidth',4);colormap(lines(N_IM)); hold on
                        legend([p(i)],'Stimulus','Location','northwest');
                        axis square; box on
                        %xticks(1:pcount+1);
                        %xticklabels(num2str((0:pcount)'));
                    end
                    set(gca,'ColorOrderIndex',1)
                    figure(i*10000);plot((1:pcount+1),ones(pcount+1,1)*0.9,'k:','LineWidth',4,'HandleVisibility','off');
                    xlabel('fixation #');
                    ylabel('r [a.u.]');
                    set(gca,'FontSize',22);
                    xlim([1 pcount+1]);
                end
            end
            

            if ~(any(PRCs>dec_thresh)) && ~ATTonly && ~RESET    
                y   = y+y_diff;                              % update postion on grid cell sheet
                x   = x+x_diff;
                yIM = yIM+yIM_diff;                          % update postion on image (scaled version of above
                xIM = xIM+xIM_diff;
                PTind = find(PLCs) - (ceil(find(PLCs)/N_PTSperIM)-1)*N_PTSperIM;
                if any(PLCs)                                 % only do this if we have an active PLC
                    microSx = abs(IMpts_bkp(PTind,1)-xIM);
                    microSy = abs(IMpts_bkp(PTind,2)-yIM);
                    if microSx<=tol && microSy<=tol    % allow 1% pixel tolerance, assumption is that microsaccades help align fovea
                        xIM = IMpts_bkp(PTind,1);
                        x   = IMpts_B(PTind,1);        % do the same for x,y because shift on image means chane to GC pop vector
                    end
                    if microSy<=tol && microSx<=tol
                        yIM = IMpts_bkp(PTind,2);
                        y   = IMpts_B(PTind,2);
                    end
                end
                if xIM>440-30 || xIM<31 || yIM>440-30 || yIM<31
                    disp('Saccade out of stimulus bounds -> reset');   % this can only happen if we are following a wrong hypothesis
                    RESET = 1;                                         % necessitates reset anyway
                end
            end
            
            
            if ~RESET
                count  = count+1;         % feature counter
                pcount = pcount+1;        % plot counter
            end
            
            
            if max(PRCs)>dec_thresh   % decision reached
                if find(PRCs==max(PRCs))==i
                    REC(i) = 1;
                    disp('Stimulus successfully recognized.'); 
                    OTtmp = squeeze(IDvT(:,:,i));
                    coff  = find(sum(OTtmp,1),1,'last');
                    OTtmp = OTtmp(:,2:coff);  
                    beg   = find(OTtmp(i,:)==0,1,'last')+1;
                    if isempty(beg)
                        beg=1;
                    end
                    OTtmp = OTtmp(:,beg:end);
                    Minds = find(OTtmp==repmat(max(OTtmp),N_IM,1));
                    Minds = Minds./([i:N_IM:length(Minds)*N_IM]');
                    if any(diff(Minds))
                        Otake(i) = 1;
                    end
                else
                    REC(i) = 0;
                    disp('RECOGNITION ERROR!');
                    %pause(0.25);
                end
                trackocc  = [trackocc occ];
                trackoccX = [trackoccX xocc];
                trackoccY = [trackoccY yocc];
            end
            
            if Ninst == 1 && PLT
                figure(56);bar(REC,'k');title('Reconition 1, Failure 0, Error -1, ');
            end
              
        end

        disp([q i]);
        disp('   ');
        
    end

    
    %%%%% END OF LOOP ACROSS STIMULI
    
    
    if q == 1   % for data analysis
        Hdata  = [];
        HdataA = [];
        HdataB = [];
    end
    
    
    for i = 1:N_IM
        if REC(i)==1
            lastnonzero = find(IDvT(i,:,i)>0,1,'last');
            allrecorded = IDvT(i,1:lastnonzero,i);
            N_nonzero   = length(find(allrecorded));
            Hdata = [Hdata N_nonzero];                       % total number of fixation, including resets   
            if any(find(allrecorded(2:end)==0))              % there was at least one reset
                lastreset = find(allrecorded==0,1,'last');            %
                nfix = lastnonzero-lastreset;
            else
                nfix = N_nonzero;
            end
            HdataA = [HdataA nfix];
        end
    end

    
    if Ninst > 1
        RECvP(:,q)      = REC;
        resetsP(:,:,q)  = resets;
        OtakevP(:,q)    = Otake;
    end
    
    Nresets = sum(resets>0,2);
    
    disp('   ');
    str1 = (['Successfully recognized ' int2str(sum(REC==1)) ' out of N_IM stimuli,']);   disp(str1);
    str0 = (['of which ' int2str(sum(Otake==1)) ' had an initially wrong hypothesis,']);   disp(str0);
    str2 = (['' int2str(sum(REC==0)) ' recognition errors.']);   disp(str2);
    str3 = (['' int2str(sum(REC==-1)) ' stimuli not recognized within ' int2str(maxN_reset) ' resets.']);   disp(str3);
    str4 = (['average TOTAL number of saccades to recognition threshold is ' int2str(mean(Hdata)) '.']);   disp(str4);
    str8 = (['AFTER LAST RESET average number of saccades to recognition threshold is ' int2str(mean(HdataA)) '.']);   disp(str8);
    str5 = (['median total number of saccades to recognition threshold is ' int2str(median(Hdata)) '.']);   disp(str5);
    str9 = (['AFTER LAST RESET median number of saccades to recognition threshold is ' int2str(median(HdataA)) '.']);   disp(str9);
    str6 = (['max TOTAL number of saccades to recognition threshold is ' int2str(max(Hdata)) '.']);   disp(str6);
    str7 = (['average number of resets for recognized stimuli is ' num2str(mean(Nresets(REC==1))) '.']);   disp(str7);
            
end



if Ninst > 1
    disp('   ');
    disp('DONE ALL PERMS');
    disp('   ');
    str1 = (['Successfully recognized ' num2str(sum(sum(RECvP==1))/Ninst) ' out of N_IM stimuli on average,']);   disp(str1);
    str0 = (['of which on average ' num2str(sum(sum(OtakevP==1))/Ninst) ' had an initially wrong hypothesis,']);   disp(str0);
    str2 = (['' num2str(sum(RECvP(:)==0)/Ninst) ' recognition errors on average.']);   disp(str2);
    str3 = (['On average ' num2str(sum(RECvP(:)==-1)/Ninst) ' stimuli not recognized within ' int2str(maxN_reset) ' resets.']);   disp(str3);
    str4 = (['average number of saccades to recognition threshold is ' num2str(mean(Hdata)) '.']);   disp(str4);
    str5 = (['median number of saccades to recognition threshold is ' num2str(median(Hdata)) '.']);   disp(str5);
    str6 = (['max number of saccades to recognition threshold is ' num2str(max(Hdata)) '.']);   disp(str6);
end



if savedata
    if Ninst == 1
        if cleanIMs && ~ATTonly && ~addedPTs
            save DOEdata_SensePred_1inst_DEFAULT IDvT vectors resets REC Hdata HdataA Otake vect2ID trackocc trackoccX trackoccY trackIMs
        end
        if smallIMs
            save DOEdata_SensePred_1inst_SMALL IDvT vectors resets REC Hdata HdataA Otake vect2ID trackocc trackoccX trackoccY trackIMs
        end
        if occluIMs
            save DOEdata_SensePred_1inst_OCCLUDED IDvT vectors resets REC Hdata HdataA Otake vect2ID trackocc trackoccX trackoccY trackIMs
        end
        if realoIMs && occfix1
            save DOEdata_SensePred_1inst_OCCFIX IDvT vectors resets REC Hdata HdataA Otake vect2ID trackocc trackoccX trackoccY trackIMs 
        end
        if realoIMs && ~occfix1
            save DOEdata_SensePred_1inst_realOCCL IDvT vectors resets REC Hdata HdataA Otake vect2ID trackocc trackoccX trackoccY trackIMs
        end
        if ATTonly && ~addedPTs
            save DOEdata_SensePred_1inst_ATTonly IDvT vectors resets REC Hdata HdataA Otake vect2ID trackocc trackoccX trackoccY trackIMs
        end
        if ATTonly && addedPTs
            save DOEdata_SensePred_1inst_addedPTs IDvT vectors resets REC Hdata HdataA Otake vect2ID trackocc trackoccX trackoccY trackIMs
        end
    else
        test = 1;
    end
end

toc;

% to increase accuracy, the GC rate maps can have a higher resolution than
% the image. then we scale down the vector to the IM and get a more accurate,
% potentially pixel perfect target?


