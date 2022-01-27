function out = detect_so_delta(lfp,Fs,varargin)
%% [out] = detect_so_delta(lfp,Fs)
%   identifies times of down states during sleep block
%       1. averages activity across channels
%       2. filters average lfp signal for delta activity (.1-4Hz)
%       3. finds positive-to-negative zero crossings (zc) during sleep
%       4. for each pos-to-neg zc
%           4.1. get the prev peak (down state)
%           4.2. get the prev negative-to-positive zc
%           4.3. get the next negative-to-positive zc
%       5. returns slow oscillations that fit the following criteria
%           5.1. time btw neg-to-pos zcs is btw .8 and 2s
%           5.2. prev peak exceeds a threshold of 90% of all peaks
%       
%   Inputs:
%       lfp - lfp is a matrix (samples x channels) of lfp data from sleep
%       Fs - sample rate of lfp 
% 
%   Optional Name/Value Pairs:
%       'sleep_idx' - a logical vector 
%                   - (1-sleep,0-awake) corresponding to that session's lfp
%       'artifact_idx' - a logical vector 
%                   - (1-artifact,0-fine) corresponding to that session's lfp
%       'bad_channels' - vector of channels that are ignored (default = [])
%       'PLOT' - [0|1]
%       'sleep_classify' - [0|1], only for classified sleep (sleep_idx)
%       'mnl_parm' - [peak-thr trough-thr dur-min dur-max] or [peak-thr trough-thr dur-min dur-max high-pass low-pass] manual parameter setting
% 
%   Outputs:
%       out
%           .down_states - time of down states
%           .peak - peak of z-scored delta band lfp signal
%           .trough - trough of z-scored delta band lfp signal

%% INPUTS
bad_channels = [];
MinPeakHeight = 90;
sleep_idx = ones(size(lfp,1),1)==1;
PLOT = 1;
mnl_parm=[85 40 .15 .5];
artifact_idx = zeros(size(lfp,1),1)==1;
sleep_classify=1;
%%
assignopts(who,varargin);

%% ORGANIZE LFP
if sleep_classify==0, sleep_idx = ones(size(lfp,1),1)==1; end
% sizing
[N,num_ch] = size(lfp);
good = setdiff(1:num_ch,bad_channels);
dt = 1/Fs;
time = ((1:N)/Fs)' - dt;

%% REMOVE ARTIFACT FROM LFP 
lfp = lfp(:,good);
for ch=1:size(lfp,2),
    lfp(artifact_idx,:) = mean(lfp(~artifact_idx,ch),1);
end

%% AVERAGE ACROSS CHANNELS
lfp = nanmean(lfp,2);

%% FILTER LFP FOR DELTA
if length(mnl_parm)==4, fpass = [.1,4];%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif length(mnl_parm)==6, fpass=mnl_parm(5:6); end
lfp_delta = filter_delta(lfp,Fs,fpass);

%% POS-TO-NEG ZERO CROSSINGS
idx = round(Fs):(N-round(Fs));
ptnzc = round(Fs)-1 + find(lfp_delta(idx)>=0 & lfp_delta(idx+1)<0 & sleep_idx(idx)==1);

%% GET STATS ON ALL POS-TO-NEG ZERO CROSSINGS
prev_peak_idx = zeros(size(ptnzc));
prev_peak     = zeros(size(ptnzc));
next_trough   = zeros(size(ptnzc));
dur           = zeros(size(ptnzc));
isasleep      = zeros(size(ptnzc));
for i=2:length(ptnzc)-1,
    prev_ntpzc          = ptnzc(i-1)-1 + find(lfp_delta(ptnzc(i-1):ptnzc(i  ))<0,1,'last');
    next_ntpzc          = ptnzc(i)  -1 + find(lfp_delta(ptnzc(i  ):ptnzc(i+1))<0,1,'last');
    [prev_peak(i),idx]  = max(lfp_delta(ptnzc(i-1):ptnzc(i)));
    prev_peak_idx(i)    = ptnzc(i-1)-1 + idx;
    [next_trough(i),idx]= min(lfp_delta(ptnzc(i):next_ntpzc));
    next_trough_idx(i)  = ptnzc(i)-1 + idx;
    dur(i)              = time(next_trough_idx(i) - prev_peak_idx(i));
    isasleep(i)         = sleep_idx(prev_peak_idx(i));
end

%% APPLY CRITERIA FOR SLOW OSCILLATIONS
thresh_up = prctile(prev_peak,mnl_parm(1));%85%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresh_dwn = prctile(next_trough,mnl_parm(2));%40%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
thresh = [thresh_up,thresh_dwn];
idx = isasleep & prev_peak>=thresh_up & next_trough<thresh_dwn & dur>mnl_parm(3) & dur<mnl_parm(4);%.15%.5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
idx2 = isasleep & prev_peak<thresh_up & next_trough<thresh_dwn & dur>mnl_parm(3)-.05;

%% organize outputs
%SO
out.so_down_states = time(prev_peak_idx(idx));
out.so_up_states = time(next_trough_idx(idx));
out.so_peaks = prev_peak(idx);
out.so_troughs = next_trough(idx);
%delta-waves
out.delta_down_states = time(prev_peak_idx(idx2));
out.delta_up_states = time(next_trough_idx(idx2));
out.delta_peaks = prev_peak(idx2);
out.delta_troughs = next_trough(idx2);
 
%% plot
if PLOT,
    data_txt={'so','delta'};
    for data_i=1:2
        figure;
        cc = get(gca,'ColorOrder');

        clf
        set(gcf,'Name',sprintf(data_txt{data_i}))
        set(gcf,'Position',[680,620-460*(data_i-1),980,360])

        subplot(2,3,1), hold on
        ref_time=eval(['out.' data_txt{data_i} '_up_states']);
        ref_time(ref_time<2 | ref_time>time(end)-2) = [];
        [W,t] = triggered_lfp(lfp,Fs,ref_time,[2,2]);
        shadedErrorBar(t,mean(W{1},2),std(W{1},[],2),{'-','Color',cc(1,:)},1)
        xlim([t(1),t(end)])

        t = 0:.1:10;
        subplot(2,3,4), hold on
        isi = diff(eval(['out.' data_txt{data_i} '_up_states']));
        histogram(isi,0:.2:10,'Normalization','pdf')
        lam = lognfit(isi);
        plot(t,lognpdf(t,lam(1),lam(2)),'LineWidth',2)
        [mx,idx] = max(lognpdf(t,lam(1),lam(2)));
        plot(t(idx),mx,'k*','MarkerSize',10)
        text(t(idx)+1,mx,sprintf('Peak @ %.2fHz',1/t(idx)))
        xlabel('sec')
        ylabel('Autocorrelation')

        ax(1)=subplot(2,3,2:3); hold on
        plot(time,lfp,'Color',cc(1,:))
        plot(time(sleep_idx),repmat(1*max(lfp),sum(sleep_idx),1),'.','Color',cc(2,:))

        ax(2)=subplot(2,3,5:6); hold on
        plot(time,lfp_delta,'Color',cc(1,:))
        plot(eval(['out.' data_txt{data_i} '_up_states']),eval(['out.' data_txt{data_i} '_peaks']),'k*')
        hline(thresh,'k-')
        hline(0,'k--')    
        linkaxes(ax,'x')
    end
end
end 

function lfp_delta = filter_delta(lfp,Fs,fpass)
    % filter for delta
    [b,a] = butter(2,fpass(1)/(Fs/2),'high');
    lfp_delta = filtfilt(b,a,lfp);
    % lowpass
    [b,a] = butter(4,fpass(2)/(Fs/2),'low');
    lfp_delta = filtfilt(b,a,lfp_delta);
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
