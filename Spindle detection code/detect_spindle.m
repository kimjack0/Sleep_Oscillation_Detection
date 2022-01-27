function spindles = detect_spindles(lfp,varargin)
%% [spindles] = detect_spindles(lfp)
%   looks for threshold crossings of the
%   envelope of a specific frequency band
%
%   INPUTS:
%       lfp - cell array, each cell is a session containing
%           - a vector of lfp data (samples x 1) in sleep1
%
%   OPTIONAL NAME/VALUE PAIRS:
%       'Fs' - sampling rate of LFP (default is 2.44140625 / 24)
%       'sleep_idx' - cell array, each cell contains a logical vector 
%                   - (1-sleep,0-awake) corresponding to that session's lfp
%       'artifact_idx' - cell array, each cell contains a logical vector 
%                   - (1-artifact,0-fine) corresponding to that session's lfp
%       'fpass' - for filtering lfp signal (default = [10,14])
%       'PLOT' - [0|1]
%       'single_thresh' - [0|1], if 1 uses one threshold for all sessions
%       'avg_flag' - [0|1], if 1 use average lfp across channels
%       'bad_channels' - vector of bad channels (ignored)
%       'sleep_classify' - [0|1], computing only for classified sleep (sleep_idx)
% 
%   OUTPUTS:
%       spindles - output structure for sleep1 (defined below)
% 
%       spindles.
%           pks - time peak of filtered lfp signal for each thresh crossing
%           start - upward thresh crossing
%           finish - downward thresh crossing
%           amp - peak amplitude of envelope
%           dur - duration (secs)
%

%% INPUTS
Fs = 2.44140625 / 24;
sleep_idx = cell(1,3);
fpass = [10,15];
PLOT = 0;
artifact_idx = cell(1,3);
single_thresh = 1;
avg_flag = 0;
bad_channels = [];
sleep_classify=1;
assignopts(who,varargin);
DEBUG = false;


%% SIZING INFO
NumSessions = length(lfp);
for s=1:NumSessions,
    N = length(lfp{s});
    dt = 1/Fs;
    time{s} = ((1:N)/Fs)' - dt;
    
    if sleep_classify==0, sleep_idx{s} = ones(size(lfp{s},1),1)==1; end
end

%% REMOVE ARTIFACT FROM LFP 
for s=1:NumSessions,
    lfp{s}(artifact_idx{s}) = mean(lfp{s}(~artifact_idx{s}));
end

%% AVERAGE LFP SIGNAL
if avg_flag,
    for s=1:NumSessions,
        good = setdiff(1:size(lfp{s},2),bad_channels);
        lfp{s} = mean(zscore(lfp{s}(:,good)),2);
    end
end

%% GET SPINDLE POWER

% FILTER SIGNAL
[bhigh,ahigh] = butter(6,fpass(1)/(Fs/2),'high');
[blow,alow] = butter(8,fpass(2)/(Fs/2),'low');
for s=1:NumSessions,
    flfp{s} = filtfilt(bhigh,ahigh,lfp{s});
    flfp{s} = filtfilt(blow,alow,flfp{s});
end

%% ZSCORE LFP
for s=1:NumSessions,
    zlfp{s} = flfp{s};
end

%% GET SIGNAL ENVELOPE
gwin = gausswin(round(200e-3*Fs));
gwin = gwin / sum(gwin);
for s=1:NumSessions,
    H{s} = abs(hilbert(zlfp{s}));
    H{s} = conv(H{s},gwin,'same'); % smooth w/ 200ms gaussian kernel 
    sleep_H{s} = nan(size(H{s}));
    sleep_H{s}(sleep_idx{s},:) = H{s}(sleep_idx{s},:);
    if sleep_classify==0, sleep_H{s} = H{s}; end
end

%% SET THRESHOLD
if single_thresh,
    total_sleep_H = [];
    for s=1:NumSessions,
%         total_sleep_H = vertcat(1,total_sleep_H,sleep_H{s}); %using envelop
        total_sleep_H = vertcat(total_sleep_H,zlfp{s}); %using LFP
    end
    sigma = nanstd(total_sleep_H);
    mu = nanmean(total_sleep_H);
    low = mu + 1.5*sigma; % 1.5
    high = mu + 2.5*sigma; % 2.5
    % can also try percentiles
%     low = prctile(total_sleep_H,85); 
%     high = prctile(total_sleep_H,95); 
    low_thresh = repmat(low,1,NumSessions);
    high_thresh = repmat(high,1,NumSessions);
else
    for s=1:NumSessions,
%         total_sleep = sleep_H{s}; %using envelop
        total_sleep = zlfp{s}; %using LFP
        sigma = nanstd(total_sleep);
        mu = nanmean(total_sleep);
        low = mu + 1.5*sigma; % 1.5
        high = mu + 2.5*sigma; % 2.5
        % can also try percentiles
%         low = prctile(H{s},85); 
%         high = prctile(H{s},95); 
        low_thresh(s) = low;
        high_thresh(s) = high;
    end
end

%% FIND SPINDLES
for s=1:NumSessions,

    % EACH SESSION
    mag = sleep_H{s};
    low = low_thresh(s);
    high = high_thresh(s);
    
    % LIMIT SEARCH
    % only search within 1s of start and finish of recording,
    % also make sure signal is below threshold at boundaries
    idx0 = round(1*Fs);
    idx0 = idx0 + find(mag(idx0:end)<low,1,'first')-1;
    idxf = length(mag) - round(1*Fs)+1;
    idxf = idx0 + find(mag(idx0:idxf)<low,1,'last')-1;
    idx = idx0:idxf;

    % THRESHOLD CROSSINGS
    upidx  = idx0 + find(mag(idx)<high & mag(idx+1)>=high)-1;
    startidx = [];
    finishidx = [];
    for i=1:length(upidx),
        startidx(i) = find(H{s}(1:upidx(i))<low,1,'last');
        finishidx(i) = upidx(i)-1 + find(H{s}(upidx(i):end)<low,1,'first');
    end
    start = time{s}(startidx);
    finish = time{s}(finishidx);

    % COMBINE SPINDLES < 1 CYCLE (1/10Hz ~ 100ms)
    [startidx,finishidx,start,finish] = ...
        combine_spindles(startidx,finishidx,start,finish,.3);
    
    % REMOVE SPINDLES WITH DURATION < 400MS (< 3 CYCLES)
    dur = finish - start;
    idx = dur>.5 & dur<2.5;
    startidx = startidx(idx);
    finishidx = finishidx(idx);
    dur = dur(idx);
    
    % FIND TIME OF PEAK OF FILTERED LFP SIGNAL
    pksidx = zeros(size(startidx));
    npksidx = zeros(size(startidx));
    for i=1:length(startidx),
        idx0 = startidx(i);
        idxf = finishidx(i);
        pksidx(i) = idx0 + find(flfp{s}(idx0:idxf)==max(flfp{s}(idx0:idxf)))-1;
        npksidx(i) = idx0 + find(flfp{s}(idx0:idxf)==min(flfp{s}(idx0:idxf)))-1;
    end

    % FIND PEAK OF ENVELOPE FILTERED LFP SIGNAL
    amp = zeros(size(startidx));
    for i=1:length(startidx),
        idx0 = startidx(i);
        idxf = finishidx(i);
        amp(i) = max(H{s}(idx0:idxf));
    end

    % CONVERT FROM INDICES TO TIME
    spindles{s}.pks = time{s}(pksidx);
    spindles{s}.npks = time{s}(npksidx);
    spindles{s}.start  = time{s}(startidx);
    spindles{s}.finish = time{s}(finishidx);
    spindles{s}.dur = dur;
    spindles{s}.amp = amp;
    spindles{s}.thr_low = low;
    spindles{s}.thr_high = high;
    spindlesidx{s} = pksidx;
end % session

%% plot
if PLOT,
    for s=1:NumSessions,
        low = low_thresh(s);
        high = high_thresh(s);
        
        figure;
        cc = get(gca,'ColorOrder');
        set(gcf,'Name',sprintf('Spindle_Oscillations'))
        set(gcf,'Position',[680,620,980,360])
        subplot(2,3,[1,4]), hold on
        [W,t] = triggered_lfp(zlfp{s},Fs,spindles{s}.pks,[1,1]);
        shadedErrorBar(t,mean(W{1},2),std(W{1},[],2),{'-','Color',cc(1,:)},1)
        xlim([t(1),t(end)])
        hline([low,high])
        
        ax(1)=subplot(2,3,2:3); hold on
        sig = mean(lfp{s},2);
        plot(time{s},sig,'Color',cc(1,:))
        plot(time{s}(sleep_idx{s}),repmat(1*max(sig),sum(sleep_idx{s}),1),'.','Color',cc(2,:))
        
        ax(2)=subplot(2,3,5:6); hold on
        plot(time{s},zlfp{s},'Color',cc(1,:))
        plot(spindles{s}.pks,zlfp{s}(spindlesidx{s}),'k*')
        hline([low,high])
        
        linkaxes(ax,'x')
    end
end % plot
end % detect_spindles

function zlfp = myzscore(lfp,sleep_idx)
    mu = mean(lfp(~sleep_idx,:));
    sigma = std(lfp(~sleep_idx,:));
    zlfp = bsxfun(@minus,lfp,mu);
    zlfp = bsxfun(@rdivide,zlfp,sigma);
    zlfp = mean(zlfp,2);
end

function varargout=shadedErrorBar(x,y,errBar,lineProps,transparent)
    % function H=shadedErrorBar(x,y,errBar,lineProps,transparent)
    %
    % Purpose 
    % Makes a 2-d line plot with a pretty shaded error bar made
    % using patch. Error bar color is chosen automatically.
    %
    % Inputs
    % x - vector of x values [optional, can be left empty]
    % y - vector of y values or a matrix of n observations by m cases
    %     where m has length(x);
    % errBar - if a vector we draw symmetric errorbars. If it has a size
    %          of [2,length(x)] then we draw asymmetric error bars with
    %          row 1 being the upper bar and row 2 being the lower bar
    %          (with respect to y). ** alternatively ** errBar can be a
    %          cellArray of two function handles. The first defines which
    %          statistic the line should be and the second defines the
    %          error bar.
    % lineProps - [optional,'-k' by default] defines the properties of
    %             the data line. e.g.:    
    %             'or-', or {'-or','markerfacecolor',[1,0.2,0.2]}
    % transparent - [optional, 0 by default] if ==1 the shaded error
    %               bar is made transparent, which forces the renderer
    %               to be openGl. However, if this is saved as .eps the
    %               resulting file will contain a raster not a vector
    %               image. 
    %
    % Outputs
    % H - a structure of handles to the generated plot objects.     
    %
    %
    % Examples
    % y=randn(30,80); x=1:size(y,2);
    % shadedErrorBar(x,mean(y,1),std(y),'g');
    % shadedErrorBar(x,y,{@median,@std},{'r-o','markerfacecolor','r'});    
    % shadedErrorBar([],y,{@median,@std},{'r-o','markerfacecolor','r'});    
    %
    % Overlay two transparent lines
    % y=randn(30,80)*10; x=(1:size(y,2))-40;
    % shadedErrorBar(x,y,{@mean,@std},'-r',1); 
    % hold on
    % y=ones(30,1)*x; y=y+0.06*y.^2+randn(size(y))*10;
    % shadedErrorBar(x,y,{@mean,@std},'-b',1); 
    % hold off
    %
    %
    % Rob Campbell - November 2009



    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Error checking    
    error(nargchk(3,5,nargin))


    %Process y using function handles if needed to make the error bar
    %dynamically
    if iscell(errBar) 
        fun1=errBar{1};
        fun2=errBar{2};
        errBar=fun2(y);
        y=fun1(y);
    else
        y=y(:)';
    end

    if isempty(x)
        x=1:length(y);
    else
        x=x(:)';
    end


    %Make upper and lower error bars if only one was specified
    if length(errBar)==length(errBar(:))
        errBar=repmat(errBar(:)',2,1);
    else
        s=size(errBar);
        f=find(s==2);
        if isempty(f), error('errBar has the wrong size'), end
        if f==2, errBar=errBar'; end
    end

    if length(x) ~= length(errBar)
        error('length(x) must equal length(errBar)')
    end

    %Set default options
    defaultProps={'-k'};
    if nargin<4, lineProps=defaultProps; end
    if isempty(lineProps), lineProps=defaultProps; end
    if ~iscell(lineProps), lineProps={lineProps}; end

    if nargin<5, transparent=0; end





    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    % Plot to get the parameters of the line 
    H.mainLine=plot(x,y,lineProps{:});


    % Work out the color of the shaded region and associated lines
    % Using alpha requires the render to be openGL and so you can't
    % save a vector image. On the other hand, you need alpha if you're
    % overlaying lines. There we have the option of choosing alpha or a
    % de-saturated solid colour for the patch surface .

    col=get(H.mainLine,'color');
    edgeColor=col+(1-col)*0.55;
    patchSaturation=0.15; %How de-saturated or transparent to make patch
    if transparent
        faceAlpha=patchSaturation;
        patchColor=col;
        set(gcf,'renderer','openGL')
    else
        faceAlpha=1;
        patchColor=col+(1-col)*(1-patchSaturation);
        set(gcf,'renderer','painters')
    end


    %Calculate the error bars
    uE=y+errBar(1,:);
    lE=y-errBar(2,:);


    %Add the patch error bar
    holdStatus=ishold;
    if ~holdStatus, hold on,  end


    %Make the patch
    yP=[lE,fliplr(uE)];
    xP=[x,fliplr(x)];

    %remove nans otherwise patch won't work
    xP(isnan(yP))=[];
    yP(isnan(yP))=[];


    H.patch=patch(xP,yP,1,'facecolor',patchColor,...
                  'edgecolor','none',...
                  'facealpha',faceAlpha);


    %Make pretty edges around the patch. 
    H.edge(1)=plot(x,lE,'-','color',edgeColor);
    H.edge(2)=plot(x,uE,'-','color',edgeColor);

    %Now replace the line (this avoids having to bugger about with z coordinates)
    delete(H.mainLine)
    H.mainLine=plot(x,y,lineProps{:});


    if ~holdStatus, hold off, end


    if nargout==1
        varargout{1}=H;
    end
end

function [waves,time] = triggered_lfp(lfp,Fs,events,win,varargin)
    %% triggered_lfp(lfp,Fs,events,win)
    %   collects waves triggered on times in events vector
    %       lfp should be a vector or matrix (samples x channels)
    %       Fs is the sample rate of lfp (default is 24414.0625 / 24)
    %       events is a cell array of event times for each channel
    %       win is a vector of length 2 (time window around events)
    % triggered_lfp(lfp,Fs,events,win,fpass)
    %       fpass is vector of freqs to filter at, default is [.1,200]
    % 
    % triggered_lfp(lfp,Fs,events,win,fpass,hilbert_flag)
    %       hilbert_flag - [0|1], if 1, do hilbert transform to filtered lfp
    %       signal before triggering
    % 
    % [waves,time] = triggered_lfp(lfp,...)
    %   returns cell array of waves (events x time) for each channel
    %   returns time vector

    %% deal with inputs
    narginchk(1,6)
    hilbert_flag = 0;
    if nargin==5,
        fpass = varargin{1};
        assert(isvector(fpass) && length(fpass)==2,...
            'fpass should be a vector of filter freqs')
    elseif nargin==6,
        fpass = varargin{1};
        assert(isvector(fpass) && length(fpass)==2,...
            'fpass should be a vector of filter freqs')
        hilbert_flag = varargin{2};
    else
        fpass = [];
    end

    %% filter 
    if ~isempty(fpass),
        [b,a] = butter(2,fpass/(Fs/2));
        lfp = filtfilt(b,a,lfp);
    end

    %% hilbert
    if hilbert_flag,
        lfp = abs(hilbert(lfp));
    end


    %% go through each channel and collect waves triggered on events
    waves = cell(size(lfp,2),1);
    for ch=1:size(lfp,2),
        waves{ch} = createdatamatc(lfp(:,ch),events,Fs,win);
        if isempty(waves{ch}),
            waves{ch} = nan(round((win(2)+win(1))*Fs));
        end
    end

    %% time vector
    time = linspace(-win(1),win(2),size(waves{1},1))';

end

function hhh=hline(y,in1,in2)
    % function h=hline(y, linetype, label)
    % 
    % Draws a horizontal line on the current axes at the location specified by 'y'.  Optional arguments are
    % 'linetype' (default is 'r:') and 'label', which applies a text label to the graph near the line.  The
    % label appears in the same color as the line.
    %
    % The line is held on the current axes, and after plotting the line, the function returns the axes to
    % its prior hold state.
    %
    % The HandleVisibility property of the line object is set to "off", so not only does it not appear on
    % legends, but it is not findable by using findobj.  Specifying an output argument causes the function to
    % return a handle to the line, so it can be manipulated or deleted.  Also, the HandleVisibility can be 
    % overridden by setting the root's ShowHiddenHandles property to on.
    %
    % h = hline(42,'g','The Answer')
    %
    % returns a handle to a green horizontal line on the current axes at y=42, and creates a text object on
    % the current axes, close to the line, which reads "The Answer".
    %
    % hline also supports vector inputs to draw multiple lines at once.  For example,
    %
    % hline([4 8 12],{'g','r','b'},{'l1','lab2','LABELC'})
    %
    % draws three lines with the appropriate labels and colors.
    % 
    % By Brandon Kuczenski for Kensington Labs.
    % brandon_kuczenski@kensingtonlabs.com
    % 8 November 2001

    if length(y)>1  % vector input
        for I=1:length(y)
            switch nargin
            case 1
                linetype='r:';
                label='';
            case 2
                if ~iscell(in1)
                    in1={in1};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                label='';
            case 3
                if ~iscell(in1)
                    in1={in1};
                end
                if ~iscell(in2)
                    in2={in2};
                end
                if I>length(in1)
                    linetype=in1{end};
                else
                    linetype=in1{I};
                end
                if I>length(in2)
                    label=in2{end};
                else
                    label=in2{I};
                end
            end
            h(I)=hline(y(I),linetype,label);
        end
    else
        switch nargin
        case 1
            linetype='r:';
            label='';
        case 2
            linetype=in1;
            label='';
        case 3
            linetype=in1;
            label=in2;
        end




        g=ishold(gca);
        hold on

        x=get(gca,'xlim');
        h=plot(x,[y y],linetype);
        if ~isempty(label)
            yy=get(gca,'ylim');
            yrange=yy(2)-yy(1);
            yunit=(y-yy(1))/yrange;
            if yunit<0.2
                text(x(1)+0.02*(x(2)-x(1)),y+0.02*yrange,label,'color',get(h,'color'))
            else
                text(x(1)+0.02*(x(2)-x(1)),y-0.02*yrange,label,'color',get(h,'color'))
            end
        end

        if g==0
        hold off
        end
        set(h,'tag','hline','handlevisibility','off') % this last part is so that it doesn't show up on legends
    end % else

    if nargout
        hhh=h;
    end
end
