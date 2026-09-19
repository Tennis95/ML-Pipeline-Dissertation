function MLstructureFinder()
% Contour viewer with cores/saddles overlay for Tecplot POINT zones.
% Header columns: 1:x  2:y  3:u  4:(u-uc)  5:v  6:w  7:sc1  8:drho
% Parses I=,J= per ZONE.
% Uses cell-centred fields (non-uniform aware) by default.

%% ----- shared state -----
S = struct();
S.inFile   = '';
S.nzones   = 0;
S.zoneMeta = struct('pos',{},'ni',{},'nj',{});
S.curZone  = 1;

S.opts = struct('minDelta',0.001, ...
                'tolSign',1e-8, ...
                'xTol',[], ...
                'jacWin',5, ...
                'x0',0, ...
                'xMin',0, ...
                'xMax',inf, ...
                'ratioMin',0.125, ...
                'ratioMax',8, ...
                'lxMin',0, ...
                'lxMax',inf);

S.plot     = struct('field','u-uc', ...
                    'filled',true, ...
                    'nlevels',25, ...
                    'useCLim',false, ...
                    'clim',[NaN NaN], ...
                    'showColorbar',true, ...
                    'cmap','parula', ...
                    'cellCentred',true);

S.cache    = struct('have',false);
S.smooth   = struct('enable',false,'kernel','binomial','passes',1);

S.featureStats = struct('data',[], 'meta',struct());

% 0 = uniform density (X_HS = sc1), 1 = He = HS stream, 2 = He = LS stream
S.hsMode = 0;

%% ----- UI -----
f = figure('Name','Structure Contour','NumberTitle','off', ...
           'Position',[80 80 1180 640]);

uicontrol(f,'Style','pushbutton','String','Open .dat file...', ...
    'Position',[20 600 110 28], ...
    'Callback',@onOpen);

uicontrol(f,'Style','pushbutton', ...
    'String','Load .mat State', ...
    'FontSize',10, ...
    'Position',[20 560 110 30], ...
    'Callback',@onLoadFlowVisState);

uicontrol(f,'Style','pushbutton', ...
    'String','Display feature information', ...
    'FontSize',9, ...
    'Position',[140 560 180 28], ...
    'Callback',@onShowFeatures);

uicontrol(f,'Style','pushbutton', ...
    'String','output features to file', ...
    'FontSize',9, ...
    'Position',[160 520 160 28], ...
    'Callback',@onOutputFeaturesToFile);

uicontrol(f,'Style','pushbutton', ...
    'String','Feature statistics', ...
    'FontSize',9, ...
    'Position',[20 520 120 28], ...
    'Callback',@onFeatureStats);

uicontrol(f,'Style','pushbutton', ...
    'String','Label triplets', ...
    'FontSize',9, ...
    'Position',[20 490 120 28], ...
    'Callback',@onLabelTriplets);

uicontrol(f,'Style','pushbutton', ...
    'String','Export triplets', ...
    'FontSize',9, ...
    'Position',[160 490 160 28], ...
    'Callback',@onExportLabelledTriplets);

uicontrol(f,'Style','text','String','Zone:', ...
    'Position',[150 604 40 18], ...
    'HorizontalAlignment','left');

S.sld = uicontrol(f,'Style','slider','Min',1,'Max',1,'Value',1, ...
    'SliderStep',[0.01 0.1], ...
    'Position',[190 602 280 20], ...
    'Callback',@(~,~)onZoneSlider());

S.lbl = uicontrol(f,'Style','text','String','/ 0', ...
    'Position',[480 604 60 18], ...
    'HorizontalAlignment','left');

uicontrol(f,'Style','pushbutton','String','Prev', ...
    'Position',[550 600 60 28], ...
    'Callback',@(~,~)bumpZone(-1));

uicontrol(f,'Style','pushbutton','String','Next', ...
    'Position',[615 600 60 28], ...
    'Callback',@(~,~)bumpZone(+1));

S.status = uicontrol(f,'Style','text','String','Select a file.', ...
    'Position',[690 604 470 18], ...
    'HorizontalAlignment','left');

% ----- smoothing -----
pSm = uipanel(f,'Title','Smoothing','Position',[0.015 0.60 0.25 0.18]);

cSmOn = uicontrol(pSm,'Style','checkbox','String','enable','Value',0, ...
    'Position',[10 60 80 22], ...
    'Callback',@(~,~)updateSmooth());

uicontrol(pSm,'Style','text','String','kernel', ...
    'Position',[100 62 50 18], ...
    'HorizontalAlignment','left');

popK = uicontrol(pSm,'Style','popupmenu','String',{'binomial','uniform'}, ...
    'Position',[150 60 90 22], ...
    'Callback',@(~,~)updateSmooth());

uicontrol(pSm,'Style','text','String','passes', ...
    'Position',[10 30 50 18], ...
    'HorizontalAlignment','left');

ePass = uicontrol(pSm,'Style','edit','String','1', ...
    'Position',[65 28 60 24], ...
    'Callback',@(~,~)updateSmooth());

% ----- detection -----
pDet = uipanel(f,'Title','Detection / Spacing','Position',[0.015 0.30 0.25 0.28]);

uicontrol(pDet,'Style','text','String','minDelta', ...
    'Position',[10 135 60 18], ...
    'HorizontalAlignment','left');

eMin = uicontrol(pDet,'Style','edit','String',num2str(S.opts.minDelta), ...
    'Position',[10 115 60 22], ...
    'Callback',@(~,~)updateDet());

uicontrol(pDet,'Style','text','String','tolSign', ...
    'Position',[90 135 60 18], ...
    'HorizontalAlignment','left');

eTol = uicontrol(pDet,'Style','edit','String',num2str(S.opts.tolSign), ...
    'Position',[90 115 60 22], ...
    'Callback',@(~,~)updateDet());

uicontrol(pDet,'Style','text','String','x_0', ...
    'Position',[10 90 25 18], ...
    'HorizontalAlignment','left');

eX0 = uicontrol(pDet,'Style','edit','String',num2str(S.opts.x0), ...
    'Position',[40 88 55 22], ...
    'Callback',@(~,~)updateDet());

uicontrol(pDet,'Style','text','String','x_{min}', ...
    'Position',[105 90 45 18], ...
    'HorizontalAlignment','left');

eXMin = uicontrol(pDet,'Style','edit','String',num2str(S.opts.xMin), ...
    'Position',[155 88 55 22], ...
    'Callback',@(~,~)updateDet());

uicontrol(pDet,'Style','text','String','x_{max}', ...
    'Position',[215 90 45 18], ...
    'HorizontalAlignment','left');

eXMax = uicontrol(pDet,'Style','edit','String',num2str(S.opts.xMax), ...
    'Position',[260 88 55 22], ...
    'Callback',@(~,~)updateDet());
uicontrol(pDet,'Style','text','String','(l/(x-x_0))min', ...
    'Position',[10 65 80 18], ...
    'HorizontalAlignment','left');

eLxMin = uicontrol(pDet,'Style','edit','String',num2str(S.opts.lxMin), ...
    'Position',[10 45 60 22], ...
    'Callback',@(~,~)updateDet());

uicontrol(pDet,'Style','text','String','(l/(x-x_0))max', ...
    'Position',[90 65 80 18], ...
    'HorizontalAlignment','left');

eLxMax = uicontrol(pDet,'Style','edit','String',num2str(S.opts.lxMax), ...
    'Position',[90 45 60 22], ...
    'Callback',@(~,~)updateDet());

uicontrol(pDet,'Style','text','String','b/a (min)', ...
    'Position',[170 65 120 18], ...
    'HorizontalAlignment','left');

eRatioMin = uicontrol(pDet,'Style','edit','String',num2str(S.opts.ratioMin), ...
    'Position',[170 45 60 22], ...
    'Callback',@(~,~)updateDet());

uicontrol(pDet,'Style','text','String','b/a (max)', ...
    'Position',[240 65 120 18], ...
    'HorizontalAlignment','left');

eRatioMax = uicontrol(pDet,'Style','edit','String',num2str(S.opts.ratioMax), ...
    'Position',[240 45 60 22], ...
    'Callback',@(~,~)updateDet());

uicontrol(pDet,'Style','text','String','HS species', ...
    'Position',[10 15 70 18], ...
    'HorizontalAlignment','left');

popHS = uicontrol(pDet,'Style','popupmenu', ...
    'String',{'uniform density','He = HS stream','He = LS stream'}, ...
    'Position',[85 15 180 22], ...
    'Callback',@(~,~)updateHS());

% ----- contour controls -----
pCtr = uipanel(f,'Title','Contour','Position',[0.015 0.02 0.25 0.26]);

uicontrol(pCtr,'Style','text','String','Field', ...
    'Position',[10 120 60 18], ...
    'HorizontalAlignment','left');

popField = uicontrol(pCtr,'Style','popupmenu', ...
    'String',{'u-uc','v','sc1','X_HS'}, ...
    'Position',[70 118 120 24], ...
    'Callback',@(~,~)updatePlot());

cFilled  = uicontrol(pCtr,'Style','checkbox','String','filled', ...
    'Value',double(S.plot.filled), ...
    'Position',[10 90 80 22], ...
    'Callback',@(~,~)updatePlot());

uicontrol(pCtr,'Style','text','String','#levels', ...
    'Position',[100 90 60 18], ...
    'HorizontalAlignment','left');

eNlev = uicontrol(pCtr,'Style','edit','String',num2str(S.plot.nlevels), ...
    'Position',[160 88 60 24], ...
    'Callback',@(~,~)updatePlot());

cCLim = uicontrol(pCtr,'Style','checkbox','String','manual CLim', ...
    'Value',double(S.plot.useCLim), ...
    'Position',[10 60 100 22], ...
    'Callback',@(~,~)updatePlot());

eCLo  = uicontrol(pCtr,'Style','edit','String','', ...
    'Position',[120 58 50 24], ...
    'Callback',@(~,~)updatePlot());

eCHi  = uicontrol(pCtr,'Style','edit','String','', ...
    'Position',[175 58 50 24], ...
    'Callback',@(~,~)updatePlot());

cBar  = uicontrol(pCtr,'Style','checkbox','String','colorbar', ...
    'Value',double(S.plot.showColorbar), ...
    'Position',[10 30 80 22], ...
    'Callback',@(~,~)updatePlot());

uicontrol(pCtr,'Style','text','String','colormap', ...
    'Position',[100 30 60 18], ...
    'HorizontalAlignment','left');

popMap = uicontrol(pCtr,'Style','popupmenu', ...
    'String',{'parula','turbo','jet','hot','gray'}, ...
    'Position',[165 28 60 24], ...
    'Callback',@(~,~)updatePlot());

cCell = uicontrol(pCtr,'Style','checkbox', ...
    'String','cell-centred (non-uniform aware)', ...
    'Value',double(S.plot.cellCentred), ...
    'Position',[10 5 220 22], ...
    'Callback',@(~,~)toggleCell());

S.ax = axes('Parent',f,'Position',[0.30 0.08 0.68 0.82]);
grid(S.ax,'on');
box(S.ax,'on');

guidata(f,S);
updateControls();

%% ----- callbacks -----
    function onLoadFlowVisState(~,~)
        [fn,fp] = uigetfile('*.mat','Select flowVisMainGUI state file');
        if isequal(fn,0)
            return;
        end

        fullp = fullfile(fp,fn);
        L = load(fullp);

        if ~isfield(L,'savedState')
            errordlg('Selected file does not contain "savedState" from flowVisMainGUI.', ...
                     'Invalid state file');
            return;
        end

        st = L.savedState;

        if isfield(st,'zoneData') && ~isempty(st.zoneData)
            zones = st.zoneData;
            srcName = 'flowviz.dat from savedState.zoneData';
        elseif isfield(st,'singflowvizData') && ~isempty(st.singflowvizData)
            zones = st.singflowvizData;
            srcName = 'singflowviz.dat from savedState.singflowvizData';
        else
            errordlg(['savedState does not contain "zoneData" or "singflowvizData".', ...
                      sprintf('\nNothing usable for structureFinder.')], ...
                     'No zones in state');
            return;
        end

        if ~iscell(zones)
            zones = {zones};
        end

        S = guidata(f);

        S.inFile   = fullp;
        S.nzones   = numel(zones);
        S.curZone  = 1;
        S.zoneMeta = struct('pos',cell(1,S.nzones), ...
                            'ni', cell(1,S.nzones), ...
                            'nj', cell(1,S.nzones));

        for k = 1:S.nzones
            z = zones{k};
            if isfield(z,'x') && isfield(z,'y')
                [ni,nj] = size(z.x);
            else
                error('Zone %d in state file has no x/y fields – cannot use in structureFinder.',k);
            end

            S.zoneMeta(k).ni  = ni;
            S.zoneMeta(k).nj  = nj;
            S.zoneMeta(k).pos = [];
        end

        S.zones = zones;
        S.cache.have = false;

        guidata(f,S);

        set(S.sld,'Min',1,'Max',max(1,S.nzones),'Value',1,'SliderStep',stepFor(S.nzones));
        set(S.lbl,'String',sprintf('/ %d',S.nzones));

        setStatus(sprintf('Loaded state (%d zones) from %s',S.nzones,fn));

        if S.nzones > 0
            drawZone();
        end

        msgbox(sprintf('Loaded %d zones from %s', S.nzones, srcName), ...
               'flowVis state loaded');
    end

    function onOpen(~,~)
        S = guidata(f);

        [fn, fp] = uigetfile({'*.dat;*.txt','Tecplot ASCII (*.dat, *.txt)';'*.*','All files'}, ...
                             'Select Tecplot POINT file');
        if isequal(fn,0)
            return;
        end

        S.inFile = fullfile(fp,fn);
        [nz, zmeta] = indexZones(S.inFile);

        S.nzones   = nz;
        S.zoneMeta = zmeta;
        S.curZone  = 1;
        S.cache.have = false;

        if isfield(S,'zones')
            S.zones = [];
        end

        guidata(f,S);

        set(S.sld,'Min',1,'Max',max(1,nz),'Value',1,'SliderStep',stepFor(nz));
        set(S.lbl,'String',sprintf('/ %d',nz));

        setStatus(sprintf('Loaded %s (%d zones)',fn,nz));

        if nz > 0
            drawZone();
        end
    end

    function onZoneSlider()
        S = guidata(f);
        if S.nzones < 1
            return;
        end

        S.curZone = max(1, min(S.nzones, round(get(S.sld,'Value'))));
        guidata(f,S);
        drawZone();
    end

    function bumpZone(dk)
        S = guidata(f);
        if S.nzones < 1
            return;
        end

        S.curZone = max(1, min(S.nzones, S.curZone + dk));
        set(S.sld,'Value',S.curZone);

        guidata(f,S);
        drawZone();
    end

function updateDet()
    S = guidata(f);

    S.opts.minDelta = str2double(get(eMin,'String'));
    S.opts.tolSign  = str2double(get(eTol,'String'));
    S.opts.x0       = str2double(get(eX0,'String'));
    S.opts.xMin     = str2double(get(eXMin,'String'));
    S.opts.xMax     = str2double(get(eXMax,'String'));

    S.opts.ratioMin = str2double(get(eRatioMin,'String'));
    S.opts.ratioMax = str2double(get(eRatioMax,'String'));
    S.opts.lxMin    = str2double(get(eLxMin,'String'));
    S.opts.lxMax    = str2double(get(eLxMax,'String'));

    if ~isfinite(S.opts.minDelta) || S.opts.minDelta < 0
        S.opts.minDelta = 0.001;
    end

    if ~isfinite(S.opts.tolSign) || S.opts.tolSign <= 0
        S.opts.tolSign = 1e-8;
    end

    if ~isfinite(S.opts.x0)
        S.opts.x0 = 0;
    end

    if ~isfinite(S.opts.xMin)
        S.opts.xMin = 0;
    end

    if ~isfinite(S.opts.xMax) || S.opts.xMax <= S.opts.xMin
        S.opts.xMax = inf;
    end

    if ~isfinite(S.opts.ratioMin) || S.opts.ratioMin <= 0
        S.opts.ratioMin = 0.125;
    end

    if ~isfinite(S.opts.ratioMax) || S.opts.ratioMax <= S.opts.ratioMin
        S.opts.ratioMax = 8;
    end

    if ~isfinite(S.opts.lxMin) || S.opts.lxMin < 0
        S.opts.lxMin = 0;
    end

    if ~isfinite(S.opts.lxMax) || S.opts.lxMax <= S.opts.lxMin
        S.opts.lxMax = inf;
    end

    set(eMin,'String',num2str(S.opts.minDelta));
    set(eTol,'String',num2str(S.opts.tolSign));
    set(eX0,'String',num2str(S.opts.x0));
    set(eXMin,'String',num2str(S.opts.xMin));

    if isfinite(S.opts.xMax)
        set(eXMax,'String',num2str(S.opts.xMax));
    else
        set(eXMax,'String','inf');
    end

    set(eRatioMin,'String',num2str(S.opts.ratioMin));
    set(eRatioMax,'String',num2str(S.opts.ratioMax));
    set(eLxMin,'String',num2str(S.opts.lxMin));

    if isfinite(S.opts.lxMax)
        set(eLxMax,'String',num2str(S.opts.lxMax));
    else
        set(eLxMax,'String','inf');
    end

    S.cache.have = false;

    guidata(f,S);
    drawZone();
end

    function updatePlot()
        S = guidata(f);

        strs = get(popField,'String');
        S.plot.field = strs{get(popField,'Value')};

        S.plot.filled   = logical(get(cFilled,'Value'));
        S.plot.nlevels  = max(2, round(str2double(get(eNlev,'String'))));
        S.plot.useCLim  = logical(get(cCLim,'Value'));

        lo = str2double(get(eCLo,'String'));
        hi = str2double(get(eCHi,'String'));

        if isfinite(lo) && isfinite(hi) && lo < hi
            S.plot.clim = [lo hi];
        else
            S.plot.clim = [NaN NaN];
        end

        maps = get(popMap,'String');
        S.plot.cmap = maps{get(popMap,'Value')};

        S.plot.showColorbar = logical(get(cBar,'Value'));

        guidata(f,S);
        drawZone();
    end

    function updateSmooth()
        S = guidata(f);

        S.smooth.enable = logical(get(cSmOn,'Value'));
        S.smooth.passes = max(1, round(str2double(get(ePass,'String'))));

        if ~isfinite(S.smooth.passes)
            S.smooth.passes = 1;
        end

        ks = get(popK,'String');
        S.smooth.kernel = ks{get(popK,'Value')};

        guidata(f,S);
        drawZone();
    end

    function updateHS()
        S = guidata(f);

        val = get(popHS,'Value');
        switch val
            case 1
                S.hsMode = 0;
            case 2
                S.hsMode = 1;
            case 3
                S.hsMode = 2;
            otherwise
                S.hsMode = 0;
        end

        guidata(f,S);
        drawZone();
    end

    function toggleCell()
        S = guidata(f);

        S.plot.cellCentred = logical(get(cCell,'Value'));
        S.cache.have = false;

        guidata(f,S);
        drawZone();
    end

    function drawZone()
        S = guidata(f);

        if S.nzones < 1
            cla(S.ax);
            return;
        end

        setStatus(sprintf('Zone %d / %d',S.curZone,S.nzones));

        [xN,yN,ucN,vN,sc1N,dxN] = getZoneXYUCVSc1(S);

        if S.plot.cellCentred
            [X,Y,UC,VV,SC] = toCellCentred(xN,yN,ucN,vN,sc1N);
            dxLoc = median(abs(diff(X,1,1)),'all','omitnan');
            if ~isfinite(dxLoc) || dxLoc <= 0
                dxLoc = max(eps, dxN);
            end
            xTol = dxLoc;
        else
            X = xN;
            Y = yN;
            UC = ucN;
            VV = vN;
            SC = sc1N;
            xTol = max(eps, dxN);
        end

        opts = S.opts;
        opts.xTol = xTol;

        [cores,saddles,combined] = detectStructuresRobust(X,Y,UC,VV,SC,S.curZone,opts);

        switch S.plot.field
            case 'u-uc'
                Z = UC;
            case 'v'
                Z = VV;
            case 'sc1'
                Z = SC;
            case 'X_HS'
                Z = highSpeedMoleFracField(SC, S.hsMode);
            otherwise
                Z = UC;
        end

        Zplot = Z;
        if S.smooth.enable
            Zplot = apply2DSmoothing(Zplot, S.smooth.kernel, S.smooth.passes);
        end

        cla(S.ax);

        if S.plot.filled
            contourf(S.ax, X, Y, Zplot, S.plot.nlevels);
        else
            contour(S.ax, X, Y, Zplot, S.plot.nlevels);
        end

        axis(S.ax,'tight');
        hold(S.ax,'on');

        colormap(S.ax, feval(S.plot.cmap));

        if S.plot.useCLim && all(isfinite(S.plot.clim))
            set(S.ax,'CLim',S.plot.clim);
        end

        if S.plot.showColorbar
            colorbar(S.ax);
        else
            cb = findobj(f,'Type','ColorBar');
            delete(cb);
        end

        try
            contour(S.ax, X, Y, UC, [0 0], 'k-', 'LineWidth', 1.2);
        catch
        end

        if ~isempty(cores)
            plot(S.ax, cores(:,1), cores(:,2), 'r.', 'MarkerSize',18);
        end

        if ~isempty(saddles)
            plot(S.ax, saddles(:,1), saddles(:,2), 'bx', 'MarkerSize',12,'LineWidth',2);
        end

        if ~isempty(combined)
            plot(S.ax, combined(:,1), combined(:,2), '.', 'Color',[0.5 0.5 0.5]);
        end

        hold(S.ax,'off');

        title(S.ax, sprintf('Zone %d — %s (%s)', S.curZone, S.plot.field, ...
            tern(S.plot.cellCentred,'cell-centred','node-based')));

        xlabel(S.ax,'x');
        ylabel(S.ax,'y');
    end

    function onShowFeatures(~,~)
        S = guidata(f);

        if S.nzones < 1
            errordlg('No zones loaded.','Feature information');
            return;
        end

        [xN,yN,ucN,vN,sc1N,dxN] = getZoneXYUCVSc1(S);

        if S.plot.cellCentred
            [X,Y,UC,VV,SC] = toCellCentred(xN,yN,ucN,vN,sc1N);
            dxLoc = median(abs(diff(X,1,1)),'all','omitnan');
            if ~isfinite(dxLoc) || dxLoc <= 0
                dxLoc = max(eps, dxN);
            end
            xTol = dxLoc;
        else
            X = xN;
            Y = yN;
            UC = ucN;
            VV = vN;
            SC = sc1N;
            xTol = max(eps, dxN);
        end

        opts = S.opts;
        opts.xTol = xTol;

        [cores,saddles,~] = detectStructuresRobust(X,Y,UC,VV,SC,S.curZone,opts);

        rows = buildFeatureRows(cores,saddles,S.opts);

        if isempty(rows)
            errordlg('No clean saddle–core–saddle structures found meeting the criteria.', ...
                     'Feature information');
            return;
        end

        colNames = { ...
            'feature', ...
            'x (upstream saddle)', ...
            'x (core)', ...
            'x (downstream saddle)', ...
            'l', ...
            'l / (x - x_0)', ...
            'b / a'};

        figInfo = figure('Name','Feature information','NumberTitle','off', ...
                         'Position',[300 200 900 400]);

        uitable('Parent',figInfo, ...
                'Data',rows, ...
                'ColumnName',colNames, ...
                'ColumnEditable',false(1,numel(colNames)), ...
                'Units','normalized', ...
                'Position',[0.02 0.02 0.96 0.96]);
    end

    function onOutputFeaturesToFile(~,~)
        S = guidata(f);

        if S.nzones < 1
            errordlg('No zones loaded.','Output features to file');
            return;
        end

        [fn, fp] = uiputfile({'*.dat','Tecplot data (*.dat)'}, ...
                             'Save feature data as');

        if isequal(fn,0)
            return;
        end

        answ = inputdlg({'Time interval between visualisation outputs:'}, ...
                        'Output time interval', 1, {'1'});

        if isempty(answ)
            return;
        end

        dtOut = str2double(strtrim(answ{1}));

        if ~isfinite(dtOut) || dtOut < 0
            errordlg('Time interval must be a finite non-negative number.', ...
                     'Output features to file');
            return;
        end

        outFile = fullfile(fp,fn);
        fid = fopen(outFile,'w');

        if fid < 0
            errordlg('Could not open output file for writing.', ...
                     'Output features to file');
            return;
        end

        fprintf(fid,'VARIABLES = "time" "x_upstream" "x_core" "x_downstream"\n');

        curZone0 = S.curZone;
        nZonesWritten = 0;
        nRowsWritten  = 0;

        fprintf('Writing feature data for %d zones to %s\n', S.nzones, outFile);
        fprintf('Using output time interval dt = %.15g\n', dtOut);

        for k = 1:S.nzones
            S.curZone = k;
            guidata(f,S);

            try
                [xN,yN,ucN,vN,sc1N,dxN] = getZoneXYUCVSc1(S);

                if S.plot.cellCentred
                    [X,Y,UC,VV,SC] = toCellCentred(xN,yN,ucN,vN,sc1N);
                    dxLoc = median(abs(diff(X,1,1)),'all','omitnan');
                    if ~isfinite(dxLoc) || dxLoc <= 0
                        dxLoc = max(eps, dxN);
                    end
                    xTol = dxLoc;
                else
                    X = xN;
                    Y = yN;
                    UC = ucN;
                    VV = vN;
                    SC = sc1N;
                    xTol = max(eps, dxN);
                end

                opts = S.opts;
                opts.xTol = xTol;

                [cores,saddles,~] = detectStructuresRobust(X,Y,UC,VV,SC,S.curZone,opts);
                rowsCell = buildFeatureRows(cores,saddles,S.opts);

                if isempty(rowsCell)
                    continue;
                end

                tZone = (k - 1) * dtOut;
                rows = zeros(size(rowsCell,1),4);

                for r = 1:size(rowsCell,1)
                    rows(r,:) = [tZone, rowsCell{r,2}, rowsCell{r,3}, rowsCell{r,4}];
                end

                fprintf(fid,'ZONE T="Zone %d", I=%d, J=1, F=POINT\n', k, size(rows,1));
                fprintf(fid,'%.15g %.15g %.15g %.15g\n', rows.');

                nZonesWritten = nZonesWritten + 1;
                nRowsWritten  = nRowsWritten + size(rows,1);

                if mod(k,100) == 0
                    fprintf('  Wrote zone %d / %d; cumulative structures: %d\n', ...
                            k, S.nzones, nRowsWritten);
                end

            catch ME
                fprintf('  Skipping zone %d due to error: %s\n', k, ME.message);
            end
        end

        fclose(fid);

        S = guidata(f);
        S.curZone = curZone0;
        guidata(f,S);
        drawZone();

        msg = sprintf('Wrote %d structures from %d zones to:\n%s', ...
                      nRowsWritten, nZonesWritten, outFile);

        fprintf('%s\n', msg);
        msgbox(msg, 'Feature export complete');
    end

    function onFeatureStats(~,~)
        S = guidata(f);

        if S.nzones < 1
            errordlg('No zones loaded.','Feature statistics');
            return;
        end

        fprintf('Scanning %d zones for saddle–core–saddle structures...\n', S.nzones);

        raw = zeros(0,8);
        curZone0 = S.curZone;

        for k = 1:S.nzones
            S.curZone = k;
            guidata(f,S);

            try
                [xN,yN,ucN,vN,sc1N,dxN] = getZoneXYUCVSc1(S);
            catch ME
                fprintf('  Skipping zone %d due to error: %s\n', k, ME.message);
                continue;
            end

            if S.plot.cellCentred
                [X,Y,UC,VV,SC] = toCellCentred(xN,yN,ucN,vN,sc1N);
                dxLoc = median(abs(diff(X,1,1)),'all','omitnan');
                if ~isfinite(dxLoc) || dxLoc <= 0
                    dxLoc = max(eps, dxN);
                end
                xTol = dxLoc;
            else
                X = xN;
                Y = yN;
                UC = ucN;
                VV = vN;
                SC = sc1N;
                xTol = max(eps, dxN);
            end

            opts = S.opts;
            opts.xTol = xTol;

            [cores,saddles,~] = detectStructuresRobust(X,Y,UC,VV,SC,S.curZone,opts);
            rawK = buildRawFeatureStats(cores,saddles,k);

            if ~isempty(rawK)
                raw = [raw; rawK]; %#ok<AGROW>
            end

            if mod(k,100) == 0
                fprintf('  Processed zone %d / %d; accumulated structures: %d\n', ...
                        k, S.nzones, size(raw,1));
            end
        end

        S.curZone = curZone0;
        guidata(f,S);
        drawZone();

        if isempty(raw)
            fprintf('No valid saddle–core–saddle structures found in any zone; no stats cached.\n');
            errordlg('No valid saddle–core–saddle structures found in any zone.', ...
                     'Feature statistics');
            return;
        end

        stats.meta = struct();
        stats.meta.nZones  = S.nzones;
        stats.meta.nStruct = size(raw,1);
        stats.meta.created = datestr(now);
        stats.data         = raw;

        S.featureStats = stats;
        S.tripletLabels.entries = struct( ...
    'zoneIdx',{}, ...
    'xL',{}, ...
    'xc',{}, ...
    'xR',{}, ...
    'l',{}, ...
    'lnorm',{}, ...
    'ratio',{}, ...
    'accepted',{});
        guidata(f,S);

        fprintf('Done. Cached %d saddle–core–saddle structures from %d zones.\n', ...
                stats.meta.nStruct, stats.meta.nZones);

        featureStatsGUI();
    end

    function featureStatsGUI()
        S = guidata(f);

        if ~isfield(S,'featureStats') || ...
           ~isfield(S.featureStats,'data') || ...
           isempty(S.featureStats.data)

            errordlg('No cached feature statistics. Run "Feature statistics" first.', ...
                     'Feature statistics');
            return;
        end

        data = S.featureStats.data;

        zoneIdxAll = data(:,1);
        l          = data(:,5);
        xref       = data(:,6);
        l_minus    = data(:,7);
        l_plus     = data(:,8);

        x0_default       = S.opts.x0;
        xMin_default     = S.opts.xMin;
        xMax_default     = S.opts.xMax;
        lxMin_default    = S.opts.lxMin;
        lxMax_default    = S.opts.lxMax;
        ratioMin_default = S.opts.ratioMin;
        ratioMax_default = S.opts.ratioMax;

        fsTick_default  = 10;
        fsLabel_default = 12;
        fsTitle_default = 12;

        nbinsLx_default      = 30;
        nbinsRatio_default   = 30;
        nbinsJointLx_default = 25;
        nbinsJointR_default  = 25;
        nbinsLxX_X_default   = 30;
        nbinsLxX_L_default   = 30;

        xPlotMin_default = xMin_default;

        cfg = struct();
        cfg.filters = struct('x0',x0_default, ...
                             'xMin',xMin_default, ...
                             'xMax',xMax_default, ...
                             'lxMin',lxMin_default, ...
                             'lxMax',lxMax_default, ...
                             'ratioMin',ratioMin_default, ...
                             'ratioMax',ratioMax_default);

        cfg.style = struct('fsTick',fsTick_default, ...
                           'fsLabel',fsLabel_default, ...
                           'fsTitle',fsTitle_default);

        cfg.lx    = struct('nbins',nbinsLx_default);
        cfg.ratio = struct('nbins',nbinsRatio_default);

        cfg.joint = struct('nbinsL',nbinsJointLx_default, ...
                           'nbinsR',nbinsJointR_default);

        cfg.lx_x  = struct('nbinsX',nbinsLxX_X_default, ...
                           'nbinsL',nbinsLxX_L_default, ...
                           'xPlotMin',xPlotMin_default, ...
                           'xPlotMax',NaN);

        figLx    = [];
        figRatio = [];
        figJoint = [];
        figLxX   = [];

        fStats = figure('Name','Feature statistics (controls)', ...
                        'NumberTitle','off', ...
                        'Position',[340 220 900 480]);

        pCtrl = uipanel('Parent',fStats, 'Title','Filters', ...
                        'Units','normalized', ...
                        'Position',[0.02 0.62 0.96 0.34]);

        uicontrol(pCtrl,'Style','text','String','x_0:', ...
                  'Units','normalized','Position',[0.02 0.65 0.06 0.22], ...
                  'HorizontalAlignment','left');

        hX0 = uicontrol(pCtrl,'Style','edit','String',num2str(x0_default), ...
                  'Units','normalized','Position',[0.08 0.67 0.08 0.24]);

        uicontrol(pCtrl,'Style','text','String','x_{min}:', ...
                  'Units','normalized','Position',[0.18 0.65 0.08 0.22], ...
                  'HorizontalAlignment','left');

        hXMin = uicontrol(pCtrl,'Style','edit','String',num2str(xMin_default), ...
                  'Units','normalized','Position',[0.26 0.67 0.08 0.24]);

        uicontrol(pCtrl,'Style','text','String','x_{max}:', ...
                  'Units','normalized','Position',[0.36 0.65 0.08 0.22], ...
                  'HorizontalAlignment','left');

        hXMax = uicontrol(pCtrl,'Style','edit','String',num2str(xMax_default), ...
                  'Units','normalized','Position',[0.44 0.67 0.08 0.24]);

        uicontrol(pCtrl,'Style','text','String','(l/(x-x_0))_{min}:', ...
                  'Units','normalized','Position',[0.02 0.18 0.18 0.22], ...
                  'HorizontalAlignment','left');

        hLxMin = uicontrol(pCtrl,'Style','edit','String',num2str(lxMin_default), ...
                  'Units','normalized','Position',[0.20 0.20 0.10 0.24]);

        uicontrol(pCtrl,'Style','text','String','(l/(x-x_0))_{max}:', ...
                  'Units','normalized','Position',[0.32 0.18 0.18 0.22], ...
                  'HorizontalAlignment','left');

        hLxMax = uicontrol(pCtrl,'Style','edit','String',num2str(lxMax_default), ...
                  'Units','normalized','Position',[0.50 0.20 0.10 0.24]);

        uicontrol(pCtrl,'Style','pushbutton','String','Update plots', ...
                  'Units','normalized','Position',[0.70 0.18 0.20 0.28], ...
                  'Callback',@(~,~)refreshPlots());

        uicontrol(pCtrl,'Style','text','String','b/a (min):', ...
                  'Units','normalized','Position',[0.02 0.02 0.18 0.22], ...
                  'HorizontalAlignment','left');

        hRMin = uicontrol(pCtrl,'Style','edit','String',num2str(ratioMin_default), ...
                  'Units','normalized','Position',[0.20 0.04 0.10 0.24]);

        uicontrol(pCtrl,'Style','text','String','b/a (max):', ...
                  'Units','normalized','Position',[0.32 0.02 0.18 0.22], ...
                  'HorizontalAlignment','left');

        hRMax = uicontrol(pCtrl,'Style','edit','String',num2str(ratioMax_default), ...
                  'Units','normalized','Position',[0.50 0.04 0.10 0.24]);

        pStyle = uipanel('Parent',fStats, 'Title','Style (fonts)', ...
                         'Units','normalized', ...
                         'Position',[0.02 0.32 0.46 0.26]);

        uicontrol(pStyle,'Style','text','String','Tick font size:', ...
                  'Units','normalized','Position',[0.05 0.60 0.30 0.25], ...
                  'HorizontalAlignment','left');

        hFsTick = uicontrol(pStyle,'Style','edit','String',num2str(fsTick_default), ...
                  'Units','normalized','Position',[0.36 0.62 0.20 0.25]);

        uicontrol(pStyle,'Style','text','String','Label font size:', ...
                  'Units','normalized','Position',[0.05 0.30 0.30 0.25], ...
                  'HorizontalAlignment','left');

        hFsLabel = uicontrol(pStyle,'Style','edit','String',num2str(fsLabel_default), ...
                  'Units','normalized','Position',[0.36 0.32 0.20 0.25]);

        uicontrol(pStyle,'Style','text','String','Title font size:', ...
                  'Units','normalized','Position',[0.05 0.02 0.30 0.25], ...
                  'HorizontalAlignment','left');

        hFsTitle = uicontrol(pStyle,'Style','edit','String',num2str(fsTitle_default), ...
                  'Units','normalized','Position',[0.36 0.04 0.20 0.25]);

        uicontrol(pStyle,'Style','pushbutton','String','Apply style to plots', ...
                  'Units','normalized','Position',[0.65 0.20 0.30 0.45], ...
                  'Callback',@(~,~)refreshPlots());

        pPlot = uipanel('Parent',fStats, 'Title','Plot settings (bins / x-range)', ...
                        'Units','normalized', ...
                        'Position',[0.50 0.32 0.48 0.26]);

        uicontrol(pPlot,'Style','text','String','Bins for l/(x-x_0):', ...
                  'Units','normalized','Position',[0.05 0.65 0.30 0.25], ...
                  'HorizontalAlignment','left');

        hNbinsLx = uicontrol(pPlot,'Style','edit','String',num2str(nbinsLx_default), ...
                  'Units','normalized','Position',[0.36 0.67 0.12 0.25]);

        uicontrol(pPlot,'Style','text','String','Bins for b/a:', ...
                  'Units','normalized','Position',[0.55 0.65 0.35 0.25], ...
                  'HorizontalAlignment','left');

        hNbinsRatio = uicontrol(pPlot,'Style','edit','String',num2str(nbinsRatio_default), ...
                  'Units','normalized','Position',[0.86 0.67 0.10 0.25]);

        uicontrol(pPlot,'Style','text','String','Joint l/(x-x_0) vs b/a bins [L,R]:', ...
                  'Units','normalized','Position',[0.05 0.35 0.55 0.25], ...
                  'HorizontalAlignment','left');

        hNbinsJointL = uicontrol(pPlot,'Style','edit','String',num2str(nbinsJointLx_default), ...
                  'Units','normalized','Position',[0.62 0.37 0.10 0.25]);

        hNbinsJointR = uicontrol(pPlot,'Style','edit','String',num2str(nbinsJointR_default), ...
                  'Units','normalized','Position',[0.75 0.37 0.10 0.25]);

        uicontrol(pPlot,'Style','text','String','l/(x-x_0) vs x bins [x,l/x]:', ...
                  'Units','normalized','Position',[0.05 0.05 0.40 0.25], ...
                  'HorizontalAlignment','left');

        hNbinsLxX_X = uicontrol(pPlot,'Style','edit','String',num2str(nbinsLxX_X_default), ...
                  'Units','normalized','Position',[0.46 0.07 0.10 0.25]);

        hNbinsLxX_L = uicontrol(pPlot,'Style','edit','String',num2str(nbinsLxX_L_default), ...
                  'Units','normalized','Position',[0.59 0.07 0.10 0.25]);

        uicontrol(pPlot,'Style','text','String','x-range for l/(x-x_0) vs x [xmin,xmax]:', ...
                  'Units','normalized','Position',[0.72 0.05 0.26 0.25], ...
                  'HorizontalAlignment','left');

        hXPlotMin = uicontrol(pPlot,'Style','edit','String',num2str(xPlotMin_default), ...
                  'Units','normalized','Position',[0.72 0.07 0.10 0.25]);

        hXPlotMax = uicontrol(pPlot,'Style','edit','String','', ...
                  'Units','normalized','Position',[0.85 0.07 0.10 0.25]);

        refreshPlots();

        function v = getVal(h, defaultVal)
            str = strtrim(get(h,'String'));
            if isempty(str)
                v = defaultVal;
                set(h,'String',num2str(defaultVal));
                return;
            end

            v = str2double(str);

            if ~isfinite(v)
                v = defaultVal;
                set(h,'String',num2str(defaultVal));
            end
        end

        function refreshPlots()
            cfg.filters.x0       = getVal(hX0,   x0_default);
            cfg.filters.xMin     = getVal(hXMin, xMin_default);
            cfg.filters.xMax     = getVal(hXMax, xMax_default);
            cfg.filters.lxMin    = getVal(hLxMin,lxMin_default);
            cfg.filters.lxMax    = getVal(hLxMax,lxMax_default);
            cfg.filters.ratioMin = getVal(hRMin, ratioMin_default);
            cfg.filters.ratioMax = getVal(hRMax, ratioMax_default);

            set(hX0,  'String',num2str(cfg.filters.x0));
            set(hXMin,'String',num2str(cfg.filters.xMin));

            if isfinite(cfg.filters.xMax)
                set(hXMax,'String',num2str(cfg.filters.xMax));
            else
                set(hXMax,'String','inf');
            end

            cfg.style.fsTick  = max(1, round(getVal(hFsTick,  fsTick_default)));
            cfg.style.fsLabel = max(1, round(getVal(hFsLabel, fsLabel_default)));
            cfg.style.fsTitle = max(1, round(getVal(hFsTitle, fsTitle_default)));

            set(hFsTick, 'String', num2str(cfg.style.fsTick));
            set(hFsLabel,'String', num2str(cfg.style.fsLabel));
            set(hFsTitle,'String', num2str(cfg.style.fsTitle));

            cfg.lx.nbins     = max(5, round(getVal(hNbinsLx,       nbinsLx_default)));
            cfg.ratio.nbins  = max(5, round(getVal(hNbinsRatio,    nbinsRatio_default)));
            cfg.joint.nbinsL = max(5, round(getVal(hNbinsJointL,   nbinsJointLx_default)));
            cfg.joint.nbinsR = max(5, round(getVal(hNbinsJointR,   nbinsJointR_default)));
            cfg.lx_x.nbinsX  = max(5, round(getVal(hNbinsLxX_X,    nbinsLxX_X_default)));
            cfg.lx_x.nbinsL  = max(5, round(getVal(hNbinsLxX_L,    nbinsLxX_L_default)));

            cfg.lx_x.xPlotMin = getVal(hXPlotMin, xPlotMin_default);

            strXMax = strtrim(get(hXPlotMax,'String'));
            if isempty(strXMax)
                cfg.lx_x.xPlotMax = NaN;
            else
                tmp = str2double(strXMax);
                if ~isfinite(tmp)
                    cfg.lx_x.xPlotMax = NaN;
                    set(hXPlotMax,'String','');
                else
                    cfg.lx_x.xPlotMax = tmp;
                end
            end

            set(hNbinsLx,     'String',num2str(cfg.lx.nbins));
            set(hNbinsRatio,  'String',num2str(cfg.ratio.nbins));
            set(hNbinsJointL, 'String',num2str(cfg.joint.nbinsL));
            set(hNbinsJointR, 'String',num2str(cfg.joint.nbinsR));
            set(hNbinsLxX_X,  'String',num2str(cfg.lx_x.nbinsX));
            set(hNbinsLxX_L,  'String',num2str(cfg.lx_x.nbinsL));

            denom = xref - cfg.filters.x0;
            lnorm = l ./ denom;
            ratio = l_plus ./ l_minus;

            mask = isfinite(lnorm) & isfinite(ratio) & isfinite(xref);

            if isfinite(cfg.filters.xMin)
                mask = mask & (xref >= cfg.filters.xMin);
            end
            if isfinite(cfg.filters.xMax)
                mask = mask & (xref <= cfg.filters.xMax);
            end
            if isfinite(cfg.filters.lxMin)
                mask = mask & (lnorm >= cfg.filters.lxMin);
            end
            if isfinite(cfg.filters.lxMax)
                mask = mask & (lnorm <= cfg.filters.lxMax);
            end
            if isfinite(cfg.filters.ratioMin)
                mask = mask & (ratio >= cfg.filters.ratioMin);
            end
            if isfinite(cfg.filters.ratioMax)
                mask = mask & (ratio <= cfg.filters.ratioMax);
            end

            lnormF = lnorm(mask);
            ratioF = ratio(mask);
            xrefF  = xref(mask);
            zoneF  = zoneIdxAll(mask);

            fprintf('Feature stats: %d / %d structures pass current thresholds.\n', ...
                    numel(lnormF), numel(lnorm));

            if ~isempty(lnormF)
                muL  = mean(lnormF);
                sigL = std(lnormF);

                [N_L, edges_L] = histcounts(lnormF, 30);
                if any(N_L)
                    [~, iMaxL] = max(N_L);
                    mpvL = 0.5 * (edges_L(iMaxL) + edges_L(iMaxL+1));
                else
                    mpvL = NaN;
                end

                fprintf('  l/(x - x0):   mean = %.4g, std = %.4g, most probable ≈ %.2f\n', ...
                        muL, sigL, mpvL);

                [maxL, idxMaxL] = max(lnormF);
                zoneMaxL = zoneF(idxMaxL);

                fprintf('  Max l/(x - x0) = %.4g in zone %d\n', maxL, zoneMaxL);
            else
                fprintf('  l/(x - x0):   no samples after thresholds.\n');
            end

            if ~isempty(ratioF)
                muR  = mean(ratioF);
                sigR = std(ratioF);

                [N_R, edges_R] = histcounts(ratioF, 30);
                if any(N_R)
                    [~, iMaxR] = max(N_R);
                    mpvR = 0.5 * (edges_R(iMaxR) + edges_R(iMaxR+1));
                else
                    mpvR = NaN;
                end

                fprintf('  b/a: mean = %.4g, std = %.4g, most probable ≈ %.2f\n', ...
                        muR, sigR, mpvR);

                [maxR, idxMaxR] = max(ratioF);
                zoneMaxR = zoneF(idxMaxR);

                fprintf('  Max b/a = %.4g in zone %d\n', maxR, zoneMaxR);
            else
                fprintf('  b/a: no samples after thresholds.\n');
            end

            if isempty(figLx) || ~isvalid(figLx)
                figLx = figure('Name','PDF of l/(x - x_0)','NumberTitle','off');
            end

            figure(figLx);
            cla;
            ax = gca;

            if ~isempty(lnormF)
                histogram(ax, lnormF, cfg.lx.nbins, 'Normalization','pdf');
                xlabel(ax,'l / (x - x_0)');
                ylabel(ax,'PDF');
                title(ax,'PDF of l / (x - x_0)');

                if isfinite(cfg.filters.lxMin) && isfinite(cfg.filters.lxMax)
                    xlim(ax,[cfg.filters.lxMin cfg.filters.lxMax]);
                end

                applyStyle(ax, cfg.style);
            else
                text(0.5,0.5,'No samples', ...
                     'HorizontalAlignment','center', ...
                     'VerticalAlignment','middle');
                axis(ax,'off');
            end

            if isempty(figRatio) || ~isvalid(figRatio)
                figRatio = figure('Name','PDF of b/a','NumberTitle','off');
            end

            figure(figRatio);
            cla;
            ax = gca;

            if ~isempty(ratioF)
                histogram(ax, ratioF, cfg.ratio.nbins, 'Normalization','pdf');
                xlabel(ax,'b/a');
                ylabel(ax,'PDF');
                title(ax,'PDF of b/a');

                if isfinite(cfg.filters.ratioMin) && isfinite(cfg.filters.ratioMax)
                    xlim(ax,[cfg.filters.ratioMin cfg.filters.ratioMax]);
                end

                applyStyle(ax, cfg.style);
            else
                text(0.5,0.5,'No samples', ...
                     'HorizontalAlignment','center', ...
                     'VerticalAlignment','middle');
                axis(ax,'off');
            end

            if isempty(figJoint) || ~isvalid(figJoint)
                figJoint = figure('Name','Joint PDF: l/(x-x_0) vs b/a','NumberTitle','off');
            end

            figure(figJoint);
            cla;
            ax = gca;

            if ~isempty(lnormF) && ~isempty(ratioF)
                [N2, edgesL2, edgesR2] = histcounts2(lnormF, ratioF, ...
                                                    [cfg.joint.nbinsL cfg.joint.nbinsR]);

                if any(N2(:))
                    centersL = 0.5*(edgesL2(1:end-1) + edgesL2(2:end));
                    centersR = 0.5*(edgesR2(1:end-1) + edgesR2(2:end));

                    [XL, YR] = meshgrid(centersL, centersR);

                    P2 = N2';
                    P2 = P2 / sum(P2(:));

                    surf(ax, XL, YR, P2);
                    shading(ax,'interp');

                    xlabel(ax,'l / (x - x_0)');
                    ylabel(ax,'b/a');
                    zlabel(ax,'Probability');
                    title(ax,'Joint PDF: l/(x-x_0) vs b/a');

                    view(ax,35,30);
                    applyStyle(ax, cfg.style);
                else
                    text(0.5,0.5,'No samples', ...
                         'HorizontalAlignment','center', ...
                         'VerticalAlignment','middle');
                    axis(ax,'off');
                end
            else
                text(0.5,0.5,'No samples', ...
                     'HorizontalAlignment','center', ...
                     'VerticalAlignment','middle');
                axis(ax,'off');
            end

            if isempty(figLxX) || ~isvalid(figLxX)
                figLxX = figure('Name','Joint PDF: l/(x-x_0) vs x','NumberTitle','off');
            end

            figure(figLxX);
            cla;
            ax = gca;

            if ~isempty(lnormF) && ~isempty(xrefF)
                xUse = xrefF;
                lUse = lnormF;

                if isfinite(cfg.lx_x.xPlotMin)
                    idx = xUse >= cfg.lx_x.xPlotMin;
                    xUse = xUse(idx);
                    lUse = lUse(idx);
                end

                if isfinite(cfg.lx_x.xPlotMax)
                    idx = xUse <= cfg.lx_x.xPlotMax;
                    xUse = xUse(idx);
                    lUse = lUse(idx);
                end

                if numel(xUse) > 0
                    [N2, edgesX, edgesL] = histcounts2(xUse, lUse, ...
                                                       [cfg.lx_x.nbinsX cfg.lx_x.nbinsL]);

                    if any(N2(:))
                        centersX = 0.5*(edgesX(1:end-1) + edgesX(2:end));
                        centersL = 0.5*(edgesL(1:end-1) + edgesL(2:end));

                        [XX, LL] = meshgrid(centersX, centersL);

                        P2 = N2';
                        P2 = P2 / sum(P2(:));

                        surf(ax, XX, LL, P2);
                        shading(ax,'interp');

                        xlabel(ax,'x');
                        ylabel(ax,'l / (x - x_0)');
                        zlabel(ax,'Probability');
                        title(ax,'Joint PDF: l/(x-x_0) vs x');

                        view(ax,35,30);
                        applyStyle(ax, cfg.style);
                    else
                        text(0.5,0.5,'No samples', ...
                             'HorizontalAlignment','center', ...
                             'VerticalAlignment','middle');
                        axis(ax,'off');
                    end
                else
                    text(0.5,0.5,'No samples in x-range', ...
                         'HorizontalAlignment','center', ...
                         'VerticalAlignment','middle');
                    axis(ax,'off');
                end
            else
                text(0.5,0.5,'No samples', ...
                     'HorizontalAlignment','center', ...
                     'VerticalAlignment','middle');
                axis(ax,'off');
            end
        end

        function applyStyle(ax, style)
            if ~ishandle(ax)
                return;
            end

            set(ax,'FontSize',style.fsTick);

            xl = get(ax,'XLabel');
            yl = get(ax,'YLabel');
            zl = get(ax,'ZLabel');
            tl = get(ax,'Title');

            if ~isempty(xl), set(xl,'FontSize',style.fsLabel); end
            if ~isempty(yl), set(yl,'FontSize',style.fsLabel); end
            if ~isempty(zl), set(zl,'FontSize',style.fsLabel); end
            if ~isempty(tl), set(tl,'FontSize',style.fsTitle); end
        end
    end

    function onLabelTriplets(~,~)
    S = guidata(f);

    if S.nzones < 1
        errordlg('No zones loaded.','Label triplets');
        return;
    end

    labelZone = S.curZone;

    [triplets,X,Y,Zplot] = computeDetailedTripletsForZone(labelZone);

    nTrip = numel(triplets);
    curTrip = 1;
    nAcceptedClicked = 0;
    nRejectedClicked = 0;

    figLab = figure('Name',sprintf('Label triplets — zone %d',labelZone), ...
                    'NumberTitle','off', ...
                    'Position',[180 100 1050 620]);

    axLab = axes('Parent',figLab,'Position',[0.08 0.18 0.72 0.76]);
    box(axLab,'on'); grid(axLab,'on');

    uicontrol(figLab,'Style','text','String','Zone:', ...
        'Units','normalized','Position',[0.82 0.93 0.05 0.035], ...
        'HorizontalAlignment','left');

    hZoneLabel = uicontrol(figLab,'Style','text', ...
        'String',sprintf('%d / %d',labelZone,S.nzones), ...
        'Units','normalized','Position',[0.88 0.93 0.10 0.035], ...
        'HorizontalAlignment','left');

hCount = uicontrol(figLab,'Style','text', ...
    'String','', ...
    'Units','normalized', ...
    'Position',[0.82 0.72 0.16 0.06], ...
    'HorizontalAlignment','left');

    hZoneSlider = uicontrol(figLab,'Style','slider', ...
        'Min',1,'Max',max(1,S.nzones),'Value',labelZone, ...
        'SliderStep',stepFor(S.nzones), ...
        'Units','normalized','Position',[0.82 0.89 0.16 0.035], ...
        'Callback',@(~,~)changeLabelZone());

    hInfo = uicontrol(figLab,'Style','text','String','', ...
        'Units','normalized','Position',[0.82 0.78 0.16 0.10], ...
        'HorizontalAlignment','left');

    uicontrol(figLab,'Style','pushbutton','String','Previous', ...
        'Units','normalized','Position',[0.82 0.64 0.16 0.06], ...
        'Callback',@(~,~)moveTriplet(-1));

    uicontrol(figLab,'Style','pushbutton','String','Next', ...
        'Units','normalized','Position',[0.82 0.58 0.16 0.06], ...
        'Callback',@(~,~)moveTriplet(+1));

    uicontrol(figLab,'Style','pushbutton','String','Accept triplet', ...
        'Units','normalized','Position',[0.82 0.46 0.16 0.07], ...
        'Callback',@(~,~)setTripletAccepted(true));

    uicontrol(figLab,'Style','pushbutton','String','Reject triplet', ...
        'Units','normalized','Position',[0.82 0.37 0.16 0.07], ...
        'Callback',@(~,~)setTripletAccepted(false));

    uicontrol(figLab,'Style','pushbutton','String','Save labels', ...
        'Units','normalized','Position',[0.82 0.24 0.16 0.07], ...
        'Callback',@(~,~)saveTripletLabels());

uicontrol(figLab,'Style','pushbutton','String','Close', ...
    'Units','normalized','Position',[0.82 0.10 0.16 0.07], ...
    'Callback',@(~,~)closeWithSave());

    drawTriplet();

    function updateTripletCounts()
        set(hCount,'String',sprintf('Accepted selected: %d\nRejected selected: %d', ...
            nAcceptedClicked,nRejectedClicked));
    end

function closeWithSave()
    saveTripletLabelsSilent();
    close(figLab);
end

    function changeLabelZone()
        saveTripletLabelsSilent();

        labelZone = max(1,min(S.nzones,round(get(hZoneSlider,'Value'))));
        set(hZoneSlider,'Value',labelZone);
        set(hZoneLabel,'String',sprintf('%d / %d',labelZone,S.nzones));

        [triplets,X,Y,Zplot] = computeDetailedTripletsForZone(labelZone);

        nTrip = numel(triplets);
        curTrip = 1;

        drawTriplet();
    end

    function moveTriplet(dk)
        if isempty(triplets), return; end
        curTrip = max(1,min(nTrip,curTrip + dk));
        drawTriplet();
    end

        function setTripletAccepted(tf)
    if isempty(triplets), return; end

    triplets(curTrip).labelled = true;
    triplets(curTrip).accepted = logical(tf);

    if tf
        nAcceptedClicked = nAcceptedClicked + 1;
    else
        nRejectedClicked = nRejectedClicked + 1;
    end

    saveTripletLabelsSilent();

    drawTriplet();
end

        function saveTripletLabels()
    saveTripletLabelsSilent();

    S = guidata(f);

    if isfield(S,'tripletLabels') && isfield(S.tripletLabels,'entries')
        nLab = numel(S.tripletLabels.entries);
    else
        nLab = 0;
    end

    fprintf('Saved %d labelled triplets globally across all zones.\n', nLab);
        end


function saveTripletLabelsSilent()
    if isempty(triplets)
        return;
    end

    S = guidata(f);

    template = makeTripletEntryTemplate_();
    emptyEntries = repmat(template,0,1);

    if ~isfield(S,'tripletLabels') || ~isfield(S.tripletLabels,'entries') || ...
       isempty(S.tripletLabels.entries)

        S.tripletLabels = struct();
        S.tripletLabels.entries = emptyEntries;

    else
        oldE = S.tripletLabels.entries;
        newE = emptyEntries;
        fns = fieldnames(template);

        for kk = 1:numel(oldE)
            e0 = template;

            for ff = 1:numel(fns)
                if isfield(oldE, fns{ff})
                    e0(1).(fns{ff}) = oldE(kk).(fns{ff});
                end
            end

            e0.accepted = logical(e0.accepted);
            newE(end+1) = e0; %#ok<AGROW>
        end

        S.tripletLabels.entries = newE;
    end

    keep = true(1,numel(S.tripletLabels.entries));
    for ii = 1:numel(S.tripletLabels.entries)
        if S.tripletLabels.entries(ii).zoneIdx == labelZone
            keep(ii) = false;
        end
    end
    S.tripletLabels.entries = S.tripletLabels.entries(keep);

    for ii = 1:numel(triplets)

        if ~isfield(triplets(ii),'labelled') || ~triplets(ii).labelled
            continue;
        end

        e = template;

        e.zoneIdx = labelZone;

        e.xL = triplets(ii).xL;
        e.yL = triplets(ii).yL;
        e.xc = triplets(ii).xc;
        e.yc = triplets(ii).yc;
        e.xR = triplets(ii).xR;
        e.yR = triplets(ii).yR;

        e.l     = triplets(ii).l;
        e.lnorm = triplets(ii).lnorm;
        e.ratio = triplets(ii).ratio;

        e.detJ_L     = getTripletField_(triplets(ii),'detJ_L');
        e.traceJ_L   = getTripletField_(triplets(ii),'traceJ_L');
        e.discJ_L    = getTripletField_(triplets(ii),'discJ_L');
        e.lambda1r_L = getTripletField_(triplets(ii),'lambda1r_L');
        e.lambda1i_L = getTripletField_(triplets(ii),'lambda1i_L');
        e.lambda2r_L = getTripletField_(triplets(ii),'lambda2r_L');
        e.lambda2i_L = getTripletField_(triplets(ii),'lambda2i_L');

        e.detJ_c     = getTripletField_(triplets(ii),'detJ_c');
        e.traceJ_c   = getTripletField_(triplets(ii),'traceJ_c');
        e.discJ_c    = getTripletField_(triplets(ii),'discJ_c');
        e.lambda1r_c = getTripletField_(triplets(ii),'lambda1r_c');
        e.lambda1i_c = getTripletField_(triplets(ii),'lambda1i_c');
        e.lambda2r_c = getTripletField_(triplets(ii),'lambda2r_c');
        e.lambda2i_c = getTripletField_(triplets(ii),'lambda2i_c');

        e.detJ_R     = getTripletField_(triplets(ii),'detJ_R');
        e.traceJ_R   = getTripletField_(triplets(ii),'traceJ_R');
        e.discJ_R    = getTripletField_(triplets(ii),'discJ_R');
        e.lambda1r_R = getTripletField_(triplets(ii),'lambda1r_R');
        e.lambda1i_R = getTripletField_(triplets(ii),'lambda1i_R');
        e.lambda2r_R = getTripletField_(triplets(ii),'lambda2r_R');
        e.lambda2i_R = getTripletField_(triplets(ii),'lambda2i_R');

        e.accepted = logical(triplets(ii).accepted);

        S.tripletLabels.entries(end+1) = e; %#ok<AGROW>
    end

    guidata(f,S);
end

    function drawTriplet()
        cla(axLab);

        if isempty(triplets)
            text(axLab,0.5,0.5,sprintf('No triplets in zone %d',labelZone), ...
                'HorizontalAlignment','center');
            axis(axLab,'off');
            set(hInfo,'String','No triplets');
            updateTripletCounts();
            return;
        end

        if S.plot.filled
            contourf(axLab,X,Y,Zplot,S.plot.nlevels);
        else
            contour(axLab,X,Y,Zplot,S.plot.nlevels);
        end

        hold(axLab,'on');

        tr = triplets(curTrip);

        plot(axLab,[tr.xL tr.xc tr.xR],[tr.yL tr.yc tr.yR], ...
            'ko-','MarkerFaceColor','y','LineWidth',2);

        plot(axLab,tr.xc,tr.yc,'ro','MarkerSize',10,'LineWidth',2);

        if isfield(tr,'labelled') && tr.labelled
            statusStr = tern(tr.accepted,'ACCEPTED','REJECTED');
        else
            statusStr = 'UNLABELLED';
        end

        title(axLab,sprintf('Zone %d — Triplet %d/%d (%s)', ...
            labelZone,curTrip,nTrip,statusStr));

        set(hInfo,'String',sprintf('l=%.4g\nl/(x-x0)=%.4g\nb/a=%.4g', ...
            tr.l,tr.lnorm,tr.ratio));

        hold(axLab,'off');

        updateTripletCounts();
    end

function v = getTripletField_(tr, name)
    if isfield(tr,name)
        v = tr.(name);
    else
        v = NaN;
    end
end

function e = makeTripletEntryTemplate_()
    e = struct( ...
        'zoneIdx',NaN, ...
        'xL',NaN,'yL',NaN, ...
        'xc',NaN,'yc',NaN, ...
        'xR',NaN,'yR',NaN, ...
        'l',NaN, ...
        'lnorm',NaN, ...
        'ratio',NaN, ...
        'detJ_L',NaN,'traceJ_L',NaN,'discJ_L',NaN, ...
        'lambda1r_L',NaN,'lambda1i_L',NaN,'lambda2r_L',NaN,'lambda2i_L',NaN, ...
        'detJ_c',NaN,'traceJ_c',NaN,'discJ_c',NaN, ...
        'lambda1r_c',NaN,'lambda1i_c',NaN,'lambda2r_c',NaN,'lambda2i_c',NaN, ...
        'detJ_R',NaN,'traceJ_R',NaN,'discJ_R',NaN, ...
        'lambda1r_R',NaN,'lambda1i_R',NaN,'lambda2r_R',NaN,'lambda2i_R',NaN, ...
        'accepted',false);
end

    end

function [triplets,X,Y,Zplot] = computeDetailedTripletsForZone(zoneIdx)
    S = guidata(f);

    curZone0 = S.curZone;
    S.curZone = zoneIdx;
    guidata(f,S);

    [xN,yN,ucN,vN,sc1N,dxN] = getZoneXYUCVSc1(S);

    if S.plot.cellCentred
        [X,Y,UC,VV,SC] = toCellCentred(xN,yN,ucN,vN,sc1N);

        dxLoc = median(abs(diff(X,1,1)),'all','omitnan');
        if ~isfinite(dxLoc) || dxLoc <= 0
            dxLoc = max(eps, dxN);
        end
        xTol = dxLoc;
    else
        X = xN;
        Y = yN;
        UC = ucN;
        VV = vN;
        SC = sc1N;
        xTol = max(eps, dxN);
    end

    opts = S.opts;
    opts.xTol = xTol;

    [cores,saddles,~] = detectStructuresRobust(X,Y,UC,VV,SC,zoneIdx,opts);

    switch S.plot.field
        case 'u-uc'
            Zplot = UC;
        case 'v'
            Zplot = VV;
        case 'sc1'
            Zplot = SC;
        case 'X_HS'
            Zplot = highSpeedMoleFracField(SC,S.hsMode);
        otherwise
            Zplot = UC;
    end

    if S.smooth.enable
        Zplot = apply2DSmoothing(Zplot,S.smooth.kernel,S.smooth.passes);
    end

    triplets = buildDetailedTriplets(cores,saddles,S.opts,S,zoneIdx);

    S = guidata(f);
    S.curZone = curZone0;
    guidata(f,S);
end

    function onExportLabelledTriplets(~,~)
    S = guidata(f);

    if S.nzones < 1
        errordlg('No zones loaded.','Export triplets');
        return;
    end

    [fn, fp] = uiputfile({'*.csv','CSV file (*.csv)'}, ...
                         'Choose base filename for triplet export');

    if isequal(fn,0)
        return;
    end

    [~,base,~] = fileparts(fn);

    allFile   = fullfile(fp,[base '_all_triplets.csv']);
    trainFile = fullfile(fp,[base '_labelled_training.csv']);

    allTriplets = collectAllTripletsForExport();

    if isempty(allTriplets)
        errordlg('No triplets found to export.','Export triplets');
        return;
    end

    writeTripletCSV(allFile, allTriplets, false);
    writeTripletCSV(trainFile, allTriplets, true);

    nLabelled = sum(arrayfun(@(q)isfield(q,'labelled') && q.labelled, allTriplets));

    fprintf('Exported %d total triplets to %s\n', numel(allTriplets), allFile);
    fprintf('Exported %d labelled training triplets to %s\n', nLabelled, trainFile);

    msgbox(sprintf(['Export complete.\n\nAll triplets:\n%s\n\n', ...
                    'Labelled training data:\n%s'], ...
                    allFile, trainFile), ...
           'Triplet export complete');
    end

function rows = computeTripletRowsForZone(zoneIdx)
    S = guidata(f);

    curZone0 = S.curZone;
    S.curZone = zoneIdx;
    guidata(f,S);

    [xN,yN,ucN,vN,sc1N,dxN] = getZoneXYUCVSc1(S);

    if S.plot.cellCentred
        [X,Y,UC,VV,SC] = toCellCentred(xN,yN,ucN,vN,sc1N);
        dxLoc = median(abs(diff(X,1,1)),'all','omitnan');
        if ~isfinite(dxLoc) || dxLoc <= 0
            dxLoc = max(eps, dxN);
        end
        xTol = dxLoc;
    else
        X = xN;
        Y = yN;
        UC = ucN;
        VV = vN;
        SC = sc1N;
        xTol = max(eps, dxN);
    end

    opts = S.opts;
    opts.xTol = xTol;

    [cores,saddles,~] = detectStructuresRobust(X,Y,UC,VV,SC,zoneIdx,opts);
    rows = buildFeatureRows(cores,saddles,S.opts);

    S = guidata(f);
    S.curZone = curZone0;
    guidata(f,S);
end

function [hasLabel, accepted] = lookupTripletLabelStatic(S,zoneIdx,xL,xc,xR)
hasLabel = false;
accepted = false;

if ~isfield(S,'tripletLabels') || ...
   ~isfield(S.tripletLabels,'entries') || ...
   isempty(S.tripletLabels.entries)
    return;
end

E = S.tripletLabels.entries;

tol = max(1e-10,1e-8*max(1,abs(xR - xL)));

for k = 1:numel(E)
    if E(k).zoneIdx ~= zoneIdx
        continue;
    end

    if abs(E(k).xL - xL) <= tol && ...
       abs(E(k).xc - xc) <= tol && ...
       abs(E(k).xR - xR) <= tol

        hasLabel = true;
        accepted = logical(E(k).accepted);
        return;
    end
end
end


function updateControls()
    S = guidata(f);

    strs = {'u-uc','v','sc1','X_HS'};
    idxField = find(strcmp(strs,S.plot.field),1,'first');
    if isempty(idxField)
        idxField = 1;
    end
    set(popField,'Value',idxField);

    set(cFilled,'Value',double(S.plot.filled));
    set(eNlev,'String',num2str(S.plot.nlevels));
    set(cCLim,'Value',double(S.plot.useCLim));

    if all(isfinite(S.plot.clim))
        set(eCLo,'String',num2str(S.plot.clim(1)));
        set(eCHi,'String',num2str(S.plot.clim(2)));
    else
        set(eCLo,'String','');
        set(eCHi,'String','');
    end

    set(cBar,'Value',double(S.plot.showColorbar));

    maps = {'parula','turbo','jet','hot','gray'};
    mapIdx = find(strcmp(maps,S.plot.cmap),1,'first');
    if isempty(mapIdx)
        mapIdx = 1;
    end
    set(popMap,'Value',mapIdx);

    set(cCell,'Value',double(S.plot.cellCentred));

    set(cSmOn,'Value',double(S.smooth.enable));
    set(ePass,'String',num2str(S.smooth.passes));

    kIdx = find(strcmp({'binomial','uniform'},S.smooth.kernel),1,'first');
    if isempty(kIdx)
        kIdx = 1;
    end
    set(popK,'Value',kIdx);

    set(eMin,'String',num2str(S.opts.minDelta));
    set(eTol,'String',num2str(S.opts.tolSign));
    set(eX0,'String',num2str(S.opts.x0));
    set(eXMin,'String',num2str(S.opts.xMin));

    if isfinite(S.opts.xMax)
        set(eXMax,'String',num2str(S.opts.xMax));
    else
        set(eXMax,'String','inf');
    end

    set(eRatioMin,'String',num2str(S.opts.ratioMin));
    set(eRatioMax,'String',num2str(S.opts.ratioMax));
    set(eLxMin,'String',num2str(S.opts.lxMin));

    if isfinite(S.opts.lxMax)
        set(eLxMax,'String',num2str(S.opts.lxMax));
    else
        set(eLxMax,'String','inf');
    end

    switch S.hsMode
        case 0
            set(popHS,'Value',1);
        case 1
            set(popHS,'Value',2);
        case 2
            set(popHS,'Value',3);
        otherwise
            set(popHS,'Value',1);
    end
end


    function setStatus(msg)
        S = guidata(f);
        if isfield(S,'status') && ishandle(S.status)
            set(S.status,'String',msg);
        end
    end

function allTriplets = collectAllTripletsForExport()
    S = guidata(f);

    allTriplets = struct([]);

    curZone0 = S.curZone;

    for k = 1:S.nzones
        try
            [tripK,~,~,~] = computeDetailedTripletsForZone(k);

            if ~isempty(tripK)
                if isempty(allTriplets)
                    allTriplets = tripK(:);
                else
                    allTriplets = [allTriplets; tripK(:)]; %#ok<AGROW>
                end
            end

            if mod(k,100) == 0
                fprintf('  Export scan processed zone %d / %d; accumulated triplets: %d\n', ...
                        k, S.nzones, numel(allTriplets));
            end

        catch ME
            fprintf('  Skipping zone %d during triplet export: %s\n', k, ME.message);
        end
    end

    S = guidata(f);
    S.curZone = curZone0;
    guidata(f,S);
end

function writeTripletCSV(outFile, triplets, labelledOnly)
    fid = fopen(outFile,'w');

    if fid < 0
        error('Could not open %s for writing.', outFile);
    end

    fprintf(fid,['zone,', ...
                 'x_upstream,y_upstream,x_core,y_core,x_downstream,y_downstream,', ...
                 'l,l_over_xminusx0,b_over_a,', ...
                 'detJ_upstream,traceJ_upstream,discJ_upstream,', ...
                 'lambda1r_upstream,lambda1i_upstream,lambda2r_upstream,lambda2i_upstream,', ...
                 'detJ_core,traceJ_core,discJ_core,', ...
                 'lambda1r_core,lambda1i_core,lambda2r_core,lambda2i_core,', ...
                 'detJ_downstream,traceJ_downstream,discJ_downstream,', ...
                 'lambda1r_downstream,lambda1i_downstream,lambda2r_downstream,lambda2i_downstream,', ...
                 'labelled,accepted (1=accept,0=reject,2=unlabelled)\n']);

    for k = 1:numel(triplets)
        tr = triplets(k);

        isLabelled = isfield(tr,'labelled') && tr.labelled;

        if labelledOnly && ~isLabelled
            continue;
        end

if isLabelled
    accVal = double(tr.accepted);   % 1 or 0
else
    accVal = 2;                     % not labelled
end

        fprintf(fid,['%d,', ...
                     '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,', ...
                     '%.15g,%.15g,%.15g,', ...
                     '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,', ...
                     '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,', ...
                     '%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,', ...
                     '%d,%.15g\n'], ...
                tr.zoneIdx, ...
                tr.xL,tr.yL,tr.xc,tr.yc,tr.xR,tr.yR, ...
                tr.l,tr.lnorm,tr.ratio, ...
                tr.detJ_L,tr.traceJ_L,tr.discJ_L, ...
                tr.lambda1r_L,tr.lambda1i_L,tr.lambda2r_L,tr.lambda2i_L, ...
                tr.detJ_c,tr.traceJ_c,tr.discJ_c, ...
                tr.lambda1r_c,tr.lambda1i_c,tr.lambda2r_c,tr.lambda2i_c, ...
                tr.detJ_R,tr.traceJ_R,tr.discJ_R, ...
                tr.lambda1r_R,tr.lambda1i_R,tr.lambda2r_R,tr.lambda2i_R, ...
                double(isLabelled),accVal);
    end

    fclose(fid);
end


end

%% ===== helpers (outside the GUI function) =====

function [nz, zmeta] = indexZones(inFile)
fid = fopen(inFile,'r');
assert(fid > 0,'Cannot open %s',inFile);

fgetl(fid);

zmeta = struct('pos',{},'ni',{},'nj',{});

while true
    L = fgetl(fid);
    if ~ischar(L)
        break;
    end

    if contains(upper(L),'ZONE')
        ni = nan;
        nj = nan;

        tI = regexp(L,'I\s*=\s*(\d+)','tokens','once');
        tJ = regexp(L,'J\s*=\s*(\d+)','tokens','once');

        if ~isempty(tI)
            ni = str2double(tI{1});
        end
        if ~isempty(tJ)
            nj = str2double(tJ{1});
        end

        zmeta(end+1) = struct('pos', ftell(fid), 'ni', ni, 'nj', nj); %#ok<AGROW>
    end
end

fclose(fid);
nz = numel(zmeta);
end


function [x,y,uc,v,sc1,dx] = getZoneXYUCVSc1(S)
if isfield(S,'zones') && ~isempty(S.zones)
    z = S.zones{S.curZone};

    if ~isfield(z,'x') || ~isfield(z,'y')
        error('Zone %d in state has no x/y fields.', S.curZone);
    end

    x = z.x;
    y = z.y;

    if isfield(z,'u_uc')
        uc = z.u_uc;
    elseif isfield(z,'u')
        uc = z.u;
    else
        error('Zone %d in state has neither u_uc nor u.', S.curZone);
    end

    if isfield(z,'v')
        v = z.v;
    else
        error('Zone %d in state has no v field.', S.curZone);
    end

    if isfield(z,'sc1')
        sc1 = z.sc1;
    else
        sc1 = 0.5*ones(size(x));
    end

    dx = median(abs(diff(x,1,1)),'all','omitnan');
    if ~isfinite(dx) || dx <= 0
        dx = 1e-9;
    end
else
    meta = S.zoneMeta(S.curZone);
    [x,y,uc,v,sc1,dx] = readZoneXYUCVSc1(S.inFile, meta);
end
end


function [x,y,uc,v,sc1,dx] = readZoneXYUCVSc1(inFile, meta)
fid = fopen(inFile,'r');
assert(fid > 0);

fseek(fid, meta.pos, 'bof');

cands = { {'%f',8}, {'%f',7}, {'%f',6} };

ok = false;
points = [];
nvars = 0;
ni = meta.ni;
nj = meta.nj;

if ~isnan(ni) && ~isnan(nj)
    for k = 1:numel(cands)
        fmt = cands{k}{1};
        nv  = cands{k}{2};
        buf = ni*nj*nv;
        p0  = ftell(fid);

        data = fscanf(fid, fmt, buf);

        if numel(data) == buf
            A = reshape(data', nv, ni, nj);
            points = shiftdim(A,1);
            nvars = nv;
            ok = true;
            break;
        else
            fseek(fid,p0,'bof');
        end
    end
else
    guesses = [1024 256; 768 256; 512 256; 2048 256; 1024 128; 1536 256];

    for g = 1:size(guesses,1)
        ni = guesses(g,1);
        nj = guesses(g,2);

        for k = 1:numel(cands)
            fmt = cands{k}{1};
            nv  = cands{k}{2};
            buf = ni*nj*nv;
            p0  = ftell(fid);

            data = fscanf(fid, fmt, buf);

            if numel(data) == buf
                A = reshape(data', nv, ni, nj);
                points = shiftdim(A,1);
                nvars = nv;
                ok = true;
                break;
            else
                fseek(fid,p0,'bof');
            end
        end

        if ok
            break;
        end
    end
end

fclose(fid);

if ~ok
    error('Failed to read zone at pos %d', meta.pos);
end

x  = points(:,:,1);
y  = points(:,:,2);

if nvars >= 4
    uc = points(:,:,4);
else
    error('Zone has fewer than 4 variables; cannot read u-uc.');
end

if nvars >= 5
    v = points(:,:,5);
else
    error('Zone has fewer than 5 variables; cannot read v.');
end

if nvars >= 7
    sc1 = points(:,:,7);
else
    sc1 = 0.5*ones(size(x));
end

dx = median(abs(diff(x,1,1)),'all','omitnan');

if ~isfinite(dx) || dx <= 0
    dx = 1e-9;
end
end


function [xC,yC,ucC,vC,sc1C] = toCellCentred(x,y,uc,v,sc1)
xC = 0.25*(x(1:end-1,1:end-1) + ...
           x(2:end,  1:end-1) + ...
           x(1:end-1,2:end)   + ...
           x(2:end,  2:end));

yC = 0.25*(y(1:end-1,1:end-1) + ...
           y(2:end,  1:end-1) + ...
           y(1:end-1,2:end)   + ...
           y(2:end,  2:end));

ucC = 0.25*(uc(1:end-1,1:end-1) + ...
            uc(2:end,  1:end-1) + ...
            uc(1:end-1,2:end)   + ...
            uc(2:end,  2:end));

vC = 0.25*(v(1:end-1,1:end-1) + ...
           v(2:end,  1:end-1) + ...
           v(1:end-1,2:end)   + ...
           v(2:end,  2:end));

sc1C = 0.25*(sc1(1:end-1,1:end-1) + ...
             sc1(2:end,  1:end-1) + ...
             sc1(1:end-1,2:end)   + ...
             sc1(2:end,  2:end));
end


function rows = buildFeatureRows(cores,saddles,opts)
rows = {};

if isempty(cores) || isempty(saddles)
    return;
end

coreX = cores(:,1);
sadX  = saddles(:,1);

sadX_sorted = sort(sadX);

nC = numel(coreX);
coreLIdx = zeros(nC,1);
coreRIdx = zeros(nC,1);

for ic = 1:nC
    xc = coreX(ic);

    iL = find(sadX_sorted < xc, 1, 'last');
    iR = find(sadX_sorted > xc, 1, 'first');

    if ~isempty(iL) && ~isempty(iR)
        coreLIdx(ic) = iL;
        coreRIdx(ic) = iR;
    end
end

valid = (coreLIdx > 0) & (coreRIdx > 0);

if ~any(valid)
    return;
end

pairs = [coreLIdx(valid), coreRIdx(valid)];
coreIdxValid = find(valid);

[uniqPairs,~,icGroup] = unique(pairs,'rows');

x0       = opts.x0;
xMin     = opts.xMin;
xMax     = opts.xMax;
ratioMin = opts.ratioMin;
ratioMax = opts.ratioMax;
lxMin    = opts.lxMin;
lxMax    = opts.lxMax;

for k = 1:size(uniqPairs,1)
    coresForPair = coreIdxValid(icGroup == k);

    if numel(coresForPair) ~= 1
        continue;
    end

    icore = coresForPair(1);

    iL = uniqPairs(k,1);
    iR = uniqPairs(k,2);

    xL = sadX_sorted(iL);
    xR = sadX_sorted(iR);
    xc = coreX(icore);

    l = xR - xL;

    if l <= 0
        continue;
    end

    xref = 0.5*(xL + xR);

    if xref < xMin
        continue;
    end

    if isfinite(xMax) && xref > xMax
        continue;
    end

    denom = xref - x0;

    if denom <= 0
        lnorm = NaN;
    else
        lnorm = l / denom;
    end

    l_minus = xc - xL;
    l_plus  = xR - xc;

    if l_minus <= 0 || l_plus <= 0
        continue;
    end

    ratio = l_plus / l_minus;

    if ~isfinite(lnorm) || lnorm < lxMin || lnorm > lxMax
        continue;
    end

    if ~isfinite(ratio) || ratio < ratioMin || ratio > ratioMax
        continue;
    end

    rows(end+1,:) = { ...
        'core', ...
        xL, ...
        xc, ...
        xR, ...
        l, ...
        lnorm, ...
        ratio}; %#ok<AGROW>
end
end


function raw = buildRawFeatureStats(cores,saddles,zoneIdx)
raw = zeros(0,8);

if isempty(cores) || isempty(saddles)
    return;
end

coreX = cores(:,1);
sadX  = saddles(:,1);

sadX_sorted = sort(sadX);

nC = numel(coreX);
coreLIdx = zeros(nC,1);
coreRIdx = zeros(nC,1);

for ic = 1:nC
    xc = coreX(ic);

    iL = find(sadX_sorted < xc, 1, 'last');
    iR = find(sadX_sorted > xc, 1, 'first');

    if ~isempty(iL) && ~isempty(iR)
        coreLIdx(ic) = iL;
        coreRIdx(ic) = iR;
    end
end

valid = (coreLIdx > 0) & (coreRIdx > 0);

if ~any(valid)
    return;
end

pairs = [coreLIdx(valid), coreRIdx(valid)];
coreIdxValid = find(valid);

[uniqPairs,~,icGroup] = unique(pairs,'rows');

for p = 1:size(uniqPairs,1)
    coresForPair = coreIdxValid(icGroup == p);

    if numel(coresForPair) ~= 1
        continue;
    end

    icore = coresForPair(1);

    iL = uniqPairs(p,1);
    iR = uniqPairs(p,2);

    xL = sadX_sorted(iL);
    xR = sadX_sorted(iR);
    xc = coreX(icore);

    l = xR - xL;

    if l <= 0
        continue;
    end

    xref    = 0.5*(xL + xR);
    l_minus = xc - xL;
    l_plus  = xR - xc;

    if l_minus <= 0 || l_plus <= 0
        continue;
    end

    raw(end+1,:) = [zoneIdx, xL, xc, xR, l, xref, l_minus, l_plus]; %#ok<AGROW>
end
end

function [cores, saddles, combined] = detectStructuresRobust(X,Y,UC,VV,SC,frame,opts)
candXY = zeroLineIntersections(X,Y,UC,VV);

cores   = zeros(0,11);
saddles = zeros(0,11);

if ~isempty(candXY)
    cls = classifyByJacobian(X,Y,UC,VV,candXY, max(3,opts.jacWin));

    if ~isempty(cls)
        isCore   = cls(:,10) == 1;
        isSaddle = cls(:,10) == 2;

        if any(isCore)
            cores = [ ...
                cls(isCore,1:2), ...
                repmat(frame,sum(isCore),1), ...
                repmat(1,sum(isCore),1), ...
                cls(isCore,3:9)];
        end

        if any(isSaddle)
            saddles = [ ...
                cls(isSaddle,1:2), ...
                repmat(frame,sum(isSaddle),1), ...
                repmat(2,sum(isSaddle),1), ...
                cls(isSaddle,3:9)];
        end
    end
end

if isempty(cores) && isempty(saddles)
    [cores,saddles] = detectBySignPattern(X,Y,UC,VV,SC,frame,opts);
end

combined = [cores; saddles];

if ~isempty(combined)
    combined = sortrows(combined,1);
    combined = proximityFilter(combined, opts.minDelta);

    combined = fillJacobianDiagnosticsForPoints(X,Y,UC,VV,combined,max(3,opts.jacWin));

    cores   = combined(combined(:,4) == 1,:);
    saddles = combined(combined(:,4) == 2,:);
end
end


function [cores, saddles] = detectBySignPattern(X,Y,UC,VV,SC,frame,opts) %#ok<INUSD>
minDelta = opts.minDelta;
tolSign  = opts.tolSign;

sUC = signT(UC, tolSign);
sV  = signT(VV, tolSign);

ucDiff = diff(sUC,1,2)/2;
vDiff  = diff(sV, 1,1);

ucDiff(:, end+1) = 0;
vDiff(end+1, :)  = 0;

poi = ucDiff - vDiff;

isCore   = (poi ==  3);
isSaddle = (poi == -1);

cores   = assemblePoints(X,Y,isCore,  frame,1);
saddles = assemblePoints(X,Y,isSaddle,frame,2);

cores   = collapseNearX(cores,  minDelta);
saddles = collapseNearX(saddles,minDelta);
end


function pts = zeroLineIntersections(X,Y,UC,VV)
pts = zeros(0,2);

try
    Cu = contourc(X, Y, UC, [0 0]);
    Cv = contourc(X, Y, VV, [0 0]);

    Su = segmentsFromC(Cu);
    Sv = segmentsFromC(Cv);

    if isempty(Su) || isempty(Sv)
        return;
    end

    pts = intersectSegments(Su,Sv);
catch
    pts = zeros(0,2);
end
end


function segs = segmentsFromC(C)
segs = zeros(0,4);

idx = 1;

while idx <= size(C,2)
    n = C(2,idx);
    idx = idx + 1;

    if idx+n-1 > size(C,2)
        break;
    end

    poly = C(:, idx:idx+n-1);
    idx = idx + n;

    if size(poly,2) < 2
        continue;
    end

    p1 = poly(:,1:end-1);
    p2 = poly(:,2:end);

    segs = [segs; [p1(1,:).', p1(2,:).', p2(1,:).', p2(2,:).']]; %#ok<AGROW>
end
end


function P = intersectSegments(Su,Sv)
P = zeros(0,2);

for i = 1:size(Su,1)
    x1 = Su(i,1);
    y1 = Su(i,2);
    x2 = Su(i,3);
    y2 = Su(i,4);

    for j = 1:size(Sv,1)
        x3 = Sv(j,1);
        y3 = Sv(j,2);
        x4 = Sv(j,3);
        y4 = Sv(j,4);

        [hit,xI,yI] = segInt(x1,y1,x2,y2, x3,y3,x4,y4);

        if hit
            P(end+1,:) = [xI,yI]; %#ok<AGROW>
        end
    end
end

if ~isempty(P)
    P = unique(round(P, 10), 'rows');
end
end


function [hit,xI,yI] = segInt(x1,y1,x2,y2, x3,y3,x4,y4)
xI = NaN;
yI = NaN;

den = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4);

if abs(den) < eps
    hit = false;
    return;
end

t = ((x1-x3)*(y3-y4) - (y1-y3)*(x3-x4)) / den;
u = ((x1-x3)*(y1-y2) - (y1-y3)*(x1-x2)) / den;

hit = (t >= 0) && (t <= 1) && (u >= 0) && (u <= 1);

if hit
    xI = x1 + t*(x2-x1);
    yI = y1 + t*(y2-y1);
end
end

function cls = classifyByJacobian(X,Y,UC,VV,XY,kWin)
cls = zeros(0,10);

if isempty(XY)
    return;
end

kWin = max(3, 2*floor(kWin/2)+1);

[ni,nj] = size(X);

for n = 1:size(XY,1)
    x0 = XY(n,1);
    y0 = XY(n,2);

    [i0,j0,ok] = nearestIJ(X,Y,x0,y0);

    if ~ok
        continue;
    end

    r = (kWin-1)/2;

    if i0-r < 1 || i0+r > ni || j0-r < 1 || j0+r > nj
        continue;
    end

    iR = (i0-r):(i0+r);
    jR = (j0-r):(j0+r);

    XX = X(iR,jR);
    YY = Y(iR,jR);
    UU = UC(iR,jR);
    VVloc = VV(iR,jR);

    [dux, duy, okU] = fitPlaneGrad(XX,YY,UU);
    [dvx, dvy, okV] = fitPlaneGrad(XX,YY,VVloc);

    if ~(okU && okV)
        continue;
    end

    J = [dux, duy; dvx, dvy];

    detJ   = det(J);
    traceJ = trace(J);
    discJ  = traceJ*traceJ - 4*detJ;

    if detJ < 0
        featureType = 2;   % saddle
    elseif detJ > 0 && discJ < 0
        featureType = 1;   % core / centre-like point
    else
        continue;
    end

    ev = eig(J);

    cls(end+1,:) = [ ...
        x0, y0, ...
        detJ, traceJ, discJ, ...
        real(ev(1)), imag(ev(1)), ...
        real(ev(2)), imag(ev(2)), ...
        featureType]; %#ok<AGROW>
end
end



function [i0,j0,ok] = nearestIJ(X,Y,x0,y0)
D = (X - x0).^2 + (Y - y0).^2;

[~,lin] = min(D(:));
[i0,j0] = ind2sub(size(X), lin);

ok = isfinite(X(i0,j0) + Y(i0,j0));
end


function [gx, gy, ok] = fitPlaneGrad(X,Y,Z)
xv = X(:);
yv = Y(:);
zv = Z(:);

A = [xv, yv, ones(size(xv))];

if any(~isfinite(A(:))) || any(~isfinite(zv))
    gx = NaN;
    gy = NaN;
    ok = false;
    return;
end

coef = A \ zv;

gx = coef(1);
gy = coef(2);

ok = all(isfinite(coef));
end


function S = signT(A,t)
S = zeros(size(A),'like',A);
S(A >  t) =  1;
S(A < -t) = -1;
end

function P = assemblePoints(x,y,mask,frame,type)
if ~any(mask(:))
    P = zeros(0,11);
    return;
end

xv = x(mask);
yv = y(mask);

n = numel(xv);

P = [ ...
    xv, ...
    yv, ...
    repmat(frame,n,1), ...
    repmat(type,n,1), ...
    nan(n,7)];

P = sortrows(P,1);
end


function P = collapseNearX(P,minDelta)
if isempty(P)
    return;
end

P = sortrows(P,1);

keep = true(size(P,1),1);
lastKept = 1;

for k = 2:size(P,1)
    if abs(P(k,1) - P(lastKept,1)) < minDelta
        keep(k) = false;
    else
        lastKept = k;
    end
end

P = P(keep,:);
end


function C = proximityFilter(C,minDelta)
if size(C,1) < 2
    return;
end

C = sortrows(C,1);

keep = true(size(C,1),1);
lastKept = 1;

for r = 2:size(C,1)
    if abs(C(r,1) - C(lastKept,1)) < minDelta
        keep(r) = false;
    else
        lastKept = r;
    end
end

C = C(keep,:);
end


function step = stepFor(n)
if n <= 1
    step = [1 1];
else
    step = [1/max(n-1,1), 10/max(n-1,1)];
end
end


function out = tern(cond,a,b)
if cond
    out = a;
else
    out = b;
end
end


function A = apply2DSmoothing(A, kernel, passes)
if ~isfinite(passes) || passes < 1
    return;
end

k = getSmoothingKernel(kernel);

for p = 1:passes
    A = smoothAlongX(A, k);
    A = smoothAlongY(A, k);
end
end


function k = getSmoothingKernel(name)
switch lower(name)
    case 'binomial'
        k = [1 2 1];
        k = k/sum(k);
    case 'uniform'
        k = [1 1 1]/3;
    otherwise
        k = [1 2 1];
        k = k/sum(k);
end
end


function A = smoothAlongX(A,k)
ni = size(A,1);
r  = floor((numel(k)-1)/2);

if r > 0
    A = [repmat(A(1,:),r,1); A; repmat(A(end,:),r,1)];
end

A = conv2(k(:), 1, A, 'same');

if r > 0
    A = A(r+1:r+ni, :);
end
end


function A = smoothAlongY(A,k)
nj = size(A,2);
r  = floor((numel(k)-1)/2);

if r > 0
    A = [repmat(A(:,1),1,r), A, repmat(A(:,end),1,r)];
end

A = conv2(1, k(:).', A, 'same');

if r > 0
    A = A(:, r+1:r+nj);
end
end


function XHS = highSpeedMoleFracField(sc1, hsMode)
XHS = sc1;

if isempty(sc1)
    return;
end

if hsMode == 0
    XHS = sc1;
    return;
end

M_He  = 4.0;
M_air = 29.0;

Y_HS = min(max(sc1,0),1);
Y_LS = 1 - Y_HS;

switch hsMode
    case 1
        Y_He  = Y_HS;
        Y_air = Y_LS;

        num = Y_He ./ M_He;
        den = num + Y_air ./ M_air;

        X = num ./ den;

    case 2
        Y_air = Y_HS;
        Y_He  = Y_LS;

        num = Y_air ./ M_air;
        den = num + Y_He ./ M_He;

        X = num ./ den;

    otherwise
        X = Y_HS;
end

X(~isfinite(X)) = NaN;
XHS = X;
end

function [hasLabel, accepted] = lookupTripletLabelStatic(S,zoneIdx,xL,xc,xR)
hasLabel = false;
accepted = false;

if ~isfield(S,'tripletLabels') || ...
   ~isfield(S.tripletLabels,'entries') || ...
   isempty(S.tripletLabels.entries)
    return;
end

E = S.tripletLabels.entries;

tol = max(1e-10,1e-8*max(1,abs(xR - xL)));

for k = 1:numel(E)
    if E(k).zoneIdx ~= zoneIdx
        continue;
    end

    if abs(E(k).xL - xL) <= tol && ...
       abs(E(k).xc - xc) <= tol && ...
       abs(E(k).xR - xR) <= tol

        hasLabel = true;
        accepted = logical(E(k).accepted);
        return;
    end
end
end

function triplets = buildDetailedTriplets(cores,saddles,opts,S,zoneIdx)
triplets = struct('zoneIdx',{}, ...
                  'xL',{},'yL',{}, ...
                  'xc',{},'yc',{}, ...
                  'xR',{},'yR',{}, ...
                  'l',{},'lnorm',{},'ratio',{}, ...
                  'detJ_L',{},'traceJ_L',{},'discJ_L',{}, ...
                  'lambda1r_L',{},'lambda1i_L',{},'lambda2r_L',{},'lambda2i_L',{}, ...
                  'detJ_c',{},'traceJ_c',{},'discJ_c',{}, ...
                  'lambda1r_c',{},'lambda1i_c',{},'lambda2r_c',{},'lambda2i_c',{}, ...
                  'detJ_R',{},'traceJ_R',{},'discJ_R',{}, ...
                  'lambda1r_R',{},'lambda1i_R',{},'lambda2r_R',{},'lambda2i_R',{}, ...
                  'labelled',{}, ...
                  'accepted',{});

if isempty(cores) || isempty(saddles)
    return;
end

coreX = cores(:,1);

sadX = saddles(:,1);

[sadX_sorted,idxS] = sort(sadX);
saddles_sorted = saddles(idxS,:);

nC = numel(coreX);
coreLIdx = zeros(nC,1);
coreRIdx = zeros(nC,1);

for ic = 1:nC
    xc = coreX(ic);

    iL = find(sadX_sorted < xc,1,'last');
    iR = find(sadX_sorted > xc,1,'first');

    if ~isempty(iL) && ~isempty(iR)
        coreLIdx(ic) = iL;
        coreRIdx(ic) = iR;
    end
end

valid = (coreLIdx > 0) & (coreRIdx > 0);

if ~any(valid)
    return;
end

pairs = [coreLIdx(valid), coreRIdx(valid)];
coreIdxValid = find(valid);

[uniqPairs,~,icGroup] = unique(pairs,'rows');

for k = 1:size(uniqPairs,1)
    coresForPair = coreIdxValid(icGroup == k);

    if numel(coresForPair) ~= 1
        continue;
    end

    icore = coresForPair(1);

    iL = uniqPairs(k,1);
    iR = uniqPairs(k,2);

    Lrow = saddles_sorted(iL,:);
    Crow = cores(icore,:);
    Rrow = saddles_sorted(iR,:);

    xL = Lrow(1); yL = Lrow(2);
    xc = Crow(1); yc = Crow(2);
    xR = Rrow(1); yR = Rrow(2);

    l = xR - xL;

    if l <= 0
        continue;
    end

    xref = 0.5*(xL + xR);

    if xref < opts.xMin
        continue;
    end

    if isfinite(opts.xMax) && xref > opts.xMax
        continue;
    end

    denom = xref - opts.x0;

    if denom <= 0
        lnorm = NaN;
    else
        lnorm = l / denom;
    end

    l_minus = xc - xL;
    l_plus  = xR - xc;

    if l_minus <= 0 || l_plus <= 0
        continue;
    end

    ratio = l_plus / l_minus;

    if ~isfinite(lnorm) || lnorm < opts.lxMin || lnorm > opts.lxMax
        continue;
    end

    if ~isfinite(ratio) || ratio < opts.ratioMin || ratio > opts.ratioMax
        continue;
    end

    [hasLabel, accepted] = lookupTripletLabelStatic(S, zoneIdx, xL, xc, xR);

    triplets(end+1) = struct( ...
        'zoneIdx',zoneIdx, ...
        'xL',xL,'yL',yL, ...
        'xc',xc,'yc',yc, ...
        'xR',xR,'yR',yR, ...
        'l',l, ...
        'lnorm',lnorm, ...
        'ratio',ratio, ...
        'detJ_L',Lrow(5),'traceJ_L',Lrow(6),'discJ_L',Lrow(7), ...
        'lambda1r_L',Lrow(8),'lambda1i_L',Lrow(9),'lambda2r_L',Lrow(10),'lambda2i_L',Lrow(11), ...
        'detJ_c',Crow(5),'traceJ_c',Crow(6),'discJ_c',Crow(7), ...
        'lambda1r_c',Crow(8),'lambda1i_c',Crow(9),'lambda2r_c',Crow(10),'lambda2i_c',Crow(11), ...
        'detJ_R',Rrow(5),'traceJ_R',Rrow(6),'discJ_R',Rrow(7), ...
        'lambda1r_R',Rrow(8),'lambda1i_R',Rrow(9),'lambda2r_R',Rrow(10),'lambda2i_R',Rrow(11), ...
        'labelled',hasLabel, ...
        'accepted',accepted); %#ok<AGROW>
end


end

function P = fillJacobianDiagnosticsForPoints(X,Y,UC,VV,P,kWin)
if isempty(P)
    return;
end

if size(P,2) < 11
    P(:,11) = NaN;
end

kWin = max(3, 2*floor(kWin/2)+1);
[ni,nj] = size(X);
r = (kWin-1)/2;

for n = 1:size(P,1)

    % Only fill missing diagnostics.
    if all(isfinite(P(n,5:11)))
        continue;
    end

    x0 = P(n,1);
    y0 = P(n,2);

    [i0,j0,ok] = nearestIJ(X,Y,x0,y0);

    if ~ok
        continue;
    end

    if i0-r < 1 || i0+r > ni || j0-r < 1 || j0+r > nj
        continue;
    end

    iR = (i0-r):(i0+r);
    jR = (j0-r):(j0+r);

    XX = X(iR,jR);
    YY = Y(iR,jR);
    UU = UC(iR,jR);
    VVloc = VV(iR,jR);

    [dux, duy, okU] = fitPlaneGrad(XX,YY,UU);
    [dvx, dvy, okV] = fitPlaneGrad(XX,YY,VVloc);

    if ~(okU && okV)
        continue;
    end

    J = [dux, duy; dvx, dvy];

    detJ   = det(J);
    traceJ = trace(J);
    discJ  = traceJ*traceJ - 4*detJ;

    ev = eig(J);

    P(n,5:11) = [ ...
        detJ, traceJ, discJ, ...
        real(ev(1)), imag(ev(1)), ...
        real(ev(2)), imag(ev(2))];
end
end
