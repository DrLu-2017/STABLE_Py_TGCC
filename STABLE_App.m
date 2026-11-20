classdef STABLE_App < matlab.apps.AppBase
% STABLE_App GUI for STABLE longitudinal beam tracking code
%
%   Launch:  app = STABLE_App;
%
%   Features:
%     * Tabbed parameter editor (Lattice, MC, HHC, Fill Pattern, Tracking, PI, HOMs, Scan)
%     * Load / Save configuration (.mat)
%     * Built-in presets (Soleil2, HALF)
%     * GPU-accelerated multi-particle multi-bunch tracking
%     * Embedded real-time result plots
%
%   Dependencies (must be on MATLAB path):
%     machine.m, PI_Control.m, VoltageLoadCalc_matlab.m,
%     BinNumCalZ.mexw64, TruncatedGaussian.m, Vb_angle_calc.m,
%     detune_HC_calc.m
%
%   Author: Lu ZHAO

    % =====================================================================
    %% UI Component Properties
    % =====================================================================
    properties (Access = public)
        UIFigure              matlab.ui.Figure

        % --- Toolbar ---
        ToolbarGrid           matlab.ui.container.GridLayout
        LoadConfigBtn         matlab.ui.control.Button
        SaveConfigBtn         matlab.ui.control.Button
        PresetDropDown        matlab.ui.control.DropDown
        AnalyzeBtn            matlab.ui.control.Button
        PanelWidthDropDown    matlab.ui.control.DropDown
        RunBtn                matlab.ui.control.Button
        RunScanBtn            matlab.ui.control.Button
        StopBtn               matlab.ui.control.Button
        StatusLabel           matlab.ui.control.Label

        % --- Layout Containers ---
        ContentGrid           matlab.ui.container.GridLayout

        % --- Left Panel ---
        LeftPanel             matlab.ui.container.Panel
        ParamTabGroup         matlab.ui.container.TabGroup

        % Lattice / Beam Tab
        LatticeTab            matlab.ui.container.Tab
        CField                matlab.ui.control.NumericEditField
        E0Field               matlab.ui.control.NumericEditField
        I0Field               matlab.ui.control.NumericEditField
        U0Field               matlab.ui.control.NumericEditField
        hField                matlab.ui.control.NumericEditField
        alpha_cField          matlab.ui.control.NumericEditField
        tau_sField            matlab.ui.control.NumericEditField
        tau_zField            matlab.ui.control.NumericEditField
        sigma_e0Field         matlab.ui.control.NumericEditField
        sigma_t0Label         matlab.ui.control.Label

        % RF System Tab
        RFSystemTab           matlab.ui.container.Tab
        V_mcField             matlab.ui.control.NumericEditField
        Ncav_mcField          matlab.ui.control.NumericEditField
        Q0_mcField            matlab.ui.control.NumericEditField
        RoQ_mcField           matlab.ui.control.NumericEditField
        beta_mcField          matlab.ui.control.NumericEditField
        Q_mcField             matlab.ui.control.NumericEditField
        R_mcField             matlab.ui.control.NumericEditField
        fre_shift_mcField     matlab.ui.control.NumericEditField
        MCAutoDetuneCheckbox  matlab.ui.control.CheckBox

        % Harmonic Cavity Fields
        n_hcField             matlab.ui.control.NumericEditField
        Ncav_hcField          matlab.ui.control.NumericEditField
        Q0_hcField            matlab.ui.control.NumericEditField
        RoQ_hcField           matlab.ui.control.NumericEditField
        beta_hcField          matlab.ui.control.NumericEditField
        Q_hcField             matlab.ui.control.NumericEditField
        R_hcField             matlab.ui.control.NumericEditField
        AutoDetuneCheckbox    matlab.ui.control.CheckBox
        FormFactorIterCheckbox matlab.ui.control.CheckBox
        fre_shift_hcField     matlab.ui.control.NumericEditField
        OptimizeHCBtn         matlab.ui.control.Button

        % Fill Pattern Tab
        PatternTab            matlab.ui.container.Tab
        PatternDropDown       matlab.ui.control.DropDown
        CustomPatternField    matlab.ui.control.EditField
        FillRateLabel         matlab.ui.control.Label
        BunchNumLabel         matlab.ui.control.Label

        % Tracking Tab
        TrackingTab           matlab.ui.container.Tab
        Par_numField          matlab.ui.control.NumericEditField
        Track_numField        matlab.ui.control.NumericEditField
        Recor_stepField       matlab.ui.control.NumericEditField
        ShortRangeCheckbox    matlab.ui.control.CheckBox

        % Feedback Tab (PI & Direct RF Feedback)
        FeedbackTab           matlab.ui.container.Tab
        PITab                 matlab.ui.container.Tab
        PIEnableCheckbox      matlab.ui.control.CheckBox
        KPField               matlab.ui.control.NumericEditField
        KIField               matlab.ui.control.NumericEditField
        PI_mField             matlab.ui.control.NumericEditField
        PI_dField             matlab.ui.control.NumericEditField

        DFBEnableCheckbox     matlab.ui.control.CheckBox
        DFB_gainField         matlab.ui.control.NumericEditField
        DFB_phase_shiftField  matlab.ui.control.NumericEditField
        DFB_delayField        matlab.ui.control.NumericEditField
        DFB_mField            matlab.ui.control.NumericEditField
        DFB_everyField        matlab.ui.control.NumericEditField
        DFB_auto_scaleCheckbox matlab.ui.control.CheckBox
        DFB_adaptiveCheckbox   matlab.ui.control.CheckBox
        DFB_I1Field           matlab.ui.control.NumericEditField
        DFB_G1Field           matlab.ui.control.NumericEditField
        DFB_I2Field           matlab.ui.control.NumericEditField
        DFB_G2Field           matlab.ui.control.NumericEditField

        % HOMs Fields
        HOMsTable             matlab.ui.control.Table
        HOMAutoAlignCheckbox  matlab.ui.control.CheckBox
        AddHOMBtn             matlab.ui.control.Button
        RemoveHOMBtn          matlab.ui.control.Button
        AlignHOMBtn           matlab.ui.control.Button

        % Scan Tab
        ScanTab               matlab.ui.container.Tab
        DfMinField            matlab.ui.control.NumericEditField
        DfMaxField            matlab.ui.control.NumericEditField
        DfStepField           matlab.ui.control.NumericEditField
        I0MinField            matlab.ui.control.NumericEditField
        I0MaxField            matlab.ui.control.NumericEditField
        I0StepField           matlab.ui.control.NumericEditField
        ScanSaveFolderField   matlab.ui.control.EditField
        BrowseFolderBtn       matlab.ui.control.Button
        FastPreCheckCheckbox  matlab.ui.control.CheckBox

        % --- Right Panel ---
        RightPanel             matlab.ui.container.Panel
        ResultTabGroup         matlab.ui.container.TabGroup

        % Centroid & Sigma Tab
        CentroidTab            matlab.ui.container.Tab
        AxesTauMean            matlab.ui.control.UIAxes
        AxesTauStd             matlab.ui.control.UIAxes
        AxesDeltaMean          matlab.ui.control.UIAxes
        AxesDeltaStd           matlab.ui.control.UIAxes

        % Bunch-by-Bunch Tab
        BunchTab               matlab.ui.container.Tab
        AxesBunchSigma         matlab.ui.control.UIAxes
        AxesBunchCentroid      matlab.ui.control.UIAxes

        % Density Profile Tab
        DensityTab             matlab.ui.container.Tab
        AxesDensity            matlab.ui.control.UIAxes

        % HC Voltage Tab
        HCVoltTab              matlab.ui.container.Tab
        AxesHCAmp              matlab.ui.control.UIAxes
        AxesHCPhase            matlab.ui.control.UIAxes

        % MC Voltage Tab
        MCVoltTab              matlab.ui.container.Tab
        AxesMCAmp              matlab.ui.control.UIAxes
        AxesMCPhase            matlab.ui.control.UIAxes

        % Generator Tab
        GenTab                 matlab.ui.container.Tab
        AxesGenAmp             matlab.ui.control.UIAxes
        AxesGenPhase           matlab.ui.control.UIAxes

        % Map Tab
        MapTab                matlab.ui.container.Tab
        MapXAxisDropDown      matlab.ui.control.DropDown
        MapLoadBtn            matlab.ui.control.Button
        MapRefreshBtn         matlab.ui.control.Button
        MapFitBoundaryBtn     matlab.ui.control.Button
        MapFormFactorCheckbox matlab.ui.control.CheckBox
        AxesMap               matlab.ui.control.UIAxes

        % Analysis Tab
        AnalysisTab           matlab.ui.container.Tab
        AnalysisPlotModeDropDown matlab.ui.control.DropDown
        AxesAnalysis          matlab.ui.control.UIAxes

        % LLRF Nyquist & Bode Tab (Appendix D)
        LLRFTab               matlab.ui.container.Tab
        LLRFModeDropDown      matlab.ui.control.DropDown
        LLRFIncludeDetuneCheckbox matlab.ui.control.CheckBox
        LLRFUpdateBtn         matlab.ui.control.Button
        LLRFMetricsLabel      matlab.ui.control.Label
        AxesBodeMag           matlab.ui.control.UIAxes
        AxesBodePhase         matlab.ui.control.UIAxes
        AxesNyquist           matlab.ui.control.UIAxes
    end

    % =====================================================================
    %% Internal State
    % =====================================================================
    properties (Access = private)
        StopFlag  logical = false
        IsRunning logical = false

        AnalysisScanPlot
        AnalysisBLPlot
        AnalysisCentroidPlot
        AnalysisEnergySpreadPlot
        LastScanData = []
        AnalysisBestDetune = NaN
        AnalysisBestBunchLength = NaN
        AnalysisBestIdx = NaN
        ScanBlinkHandle = []          % Handle for blinking marker on stability map
        LastMode0Data = []            % Cached Mode-0 coherent FFT & D-mode analytical benchmark
    end

    % =====================================================================
    %% Component Creation
    % =====================================================================
    methods (Access = private)

        function createComponents(app)
            % --- Figure ---
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Name = 'STABLE Longitudinal Beam Tracking';
            app.UIFigure.Position = [80 60 1500 900];
            app.UIFigure.Color = [0.94 0.94 0.96];

            % --- Main Grid: toolbar row + content row ---
            mainGrid = uigridlayout(app.UIFigure, [2 1]);
            mainGrid.RowHeight   = {44, '1x'};
            mainGrid.ColumnWidth = {'1x'};
            mainGrid.Padding     = [4 4 4 4];
            mainGrid.RowSpacing  = 4;

            % =============================================================
            % TOOLBAR
            % =============================================================
            app.ToolbarGrid = uigridlayout(mainGrid, [1 9]);
            app.ToolbarGrid.Layout.Row    = 1;
            app.ToolbarGrid.Layout.Column = 1;
            app.ToolbarGrid.ColumnWidth   = {85, 85, 125, 80, 125, 100, 105, 75, '1x'};
            app.ToolbarGrid.Padding       = [6 4 6 4];
            app.ToolbarGrid.ColumnSpacing = 6;

            app.LoadConfigBtn = uibutton(app.ToolbarGrid, 'push');
            app.LoadConfigBtn.Layout.Row = 1; app.LoadConfigBtn.Layout.Column = 1;
            app.LoadConfigBtn.Text = 'Load';
            app.LoadConfigBtn.FontSize = 13;
            app.LoadConfigBtn.ButtonPushedFcn = @(~,~) LoadConfigButtonPushed(app);

            app.SaveConfigBtn = uibutton(app.ToolbarGrid, 'push');
            app.SaveConfigBtn.Layout.Row = 1; app.SaveConfigBtn.Layout.Column = 2;
            app.SaveConfigBtn.Text = 'Save';
            app.SaveConfigBtn.FontSize = 13;
            app.SaveConfigBtn.ButtonPushedFcn = @(~,~) SaveConfigButtonPushed(app);

            app.PresetDropDown = uidropdown(app.ToolbarGrid);
            app.PresetDropDown.Layout.Row = 1; app.PresetDropDown.Layout.Column = 3;
            app.PresetDropDown.Items = {'-- Presets --', 'Custom', 'Soleil 2', 'HALF'};
            app.PresetDropDown.Value = '-- Presets --';
            app.PresetDropDown.ValueChangedFcn = @(~,~) PresetDropDownChanged(app);

            app.AnalyzeBtn = uibutton(app.ToolbarGrid, 'push');
            app.AnalyzeBtn.Layout.Row = 1; app.AnalyzeBtn.Layout.Column = 4;
            app.AnalyzeBtn.Text = 'Analyze';
            app.AnalyzeBtn.FontSize = 13;
            app.AnalyzeBtn.ButtonPushedFcn = @(~,~) AnalyzeButtonPushed(app);

            app.PanelWidthDropDown = uidropdown(app.ToolbarGrid);
            app.PanelWidthDropDown.Layout.Row = 1; app.PanelWidthDropDown.Layout.Column = 5;
            app.PanelWidthDropDown.Items = {'Panel: 380px', 'Panel: 440px', 'Panel: 480px', 'Panel: 540px', 'Panel: 600px'};
            app.PanelWidthDropDown.Value = 'Panel: 480px';
            app.PanelWidthDropDown.FontSize = 12;
            app.PanelWidthDropDown.Tooltip = 'Adjust Left Panel Width';
            app.PanelWidthDropDown.ValueChangedFcn = @(~,~) PanelWidthChanged(app);

            app.RunBtn = uibutton(app.ToolbarGrid, 'push');
            app.RunBtn.Layout.Row = 1; app.RunBtn.Layout.Column = 6;
            app.RunBtn.Text = 'Run Track';
            app.RunBtn.FontSize = 14; app.RunBtn.FontWeight = 'bold';
            app.RunBtn.BackgroundColor = [0.18 0.63 0.18];
            app.RunBtn.FontColor = [1 1 1];
            app.RunBtn.ButtonPushedFcn = @(~,~) RunButtonPushed(app);

            app.RunScanBtn = uibutton(app.ToolbarGrid, 'push');
            app.RunScanBtn.Layout.Row = 1; app.RunScanBtn.Layout.Column = 7;
            app.RunScanBtn.Text = 'Run Scan';
            app.RunScanBtn.FontSize = 14; app.RunScanBtn.FontWeight = 'bold';
            app.RunScanBtn.BackgroundColor = [0.85 0.5 0.1];
            app.RunScanBtn.FontColor = [1 1 1];
            app.RunScanBtn.ButtonPushedFcn = @(~,~) RunScanButtonPushed(app);

            app.StopBtn = uibutton(app.ToolbarGrid, 'push');
            app.StopBtn.Layout.Row = 1; app.StopBtn.Layout.Column = 8;
            app.StopBtn.Text = 'Stop';
            app.StopBtn.FontSize = 13;
            app.StopBtn.BackgroundColor = [0.83 0.18 0.18];
            app.StopBtn.FontColor = [1 1 1];
            app.StopBtn.Enable = 'off';
            app.StopBtn.ButtonPushedFcn = @(~,~) StopButtonPushed(app);

            app.StatusLabel = uilabel(app.ToolbarGrid);
            app.StatusLabel.Layout.Row = 1; app.StatusLabel.Layout.Column = 9;
            app.StatusLabel.Text = 'Ready';
            app.StatusLabel.FontSize = 13;
            app.StatusLabel.FontColor = [0.3 0.3 0.3];

            % =============================================================
            % CONTENT: Left panel + Right panel
            % =============================================================
            app.ContentGrid = uigridlayout(mainGrid, [1 2]);
            app.ContentGrid.Layout.Row    = 2;
            app.ContentGrid.Layout.Column = 1;
            app.ContentGrid.ColumnWidth   = {480, '1x'};
            app.ContentGrid.Padding       = [0 0 0 0];
            app.ContentGrid.ColumnSpacing = 4;

            % --- Left Panel ---
            app.LeftPanel = uipanel(app.ContentGrid);
            app.LeftPanel.Layout.Row = 1; app.LeftPanel.Layout.Column = 1;
            app.LeftPanel.Title = 'Parameters';
            app.LeftPanel.FontSize = 13; app.LeftPanel.FontWeight = 'bold';

            app.ParamTabGroup = uitabgroup(app.LeftPanel);
            app.ParamTabGroup.TabLocation = 'left';

            % --- Right Panel ---
            app.RightPanel = uipanel(app.ContentGrid);
            app.RightPanel.Layout.Row = 1; app.RightPanel.Layout.Column = 2;
            app.RightPanel.Title = 'Results';
            app.RightPanel.FontSize = 13; app.RightPanel.FontWeight = 'bold';

            app.ResultTabGroup = uitabgroup(app.RightPanel);
            app.ResultTabGroup.TabLocation = 'left';

            % Build each tab
            createLatticeTab(app);
            createRFSystemTab(app);
            createPatternTab(app);
            createTrackingTab(app);
            createFeedbackTab(app);
            createScanTab(app);
            createResultTabs(app);

            app.ResultTabGroup.SelectionChangedFcn = @(~,e) ResultTabGroupSelectionChanged(app, e);

            % Make figure visible
            app.UIFigure.Visible = 'on';
        end

        % -----------------------------------------------------------------
        %  Lattice / Beam Tab
        % -----------------------------------------------------------------
        function createLatticeTab(app)
            app.LatticeTab = uitab(app.ParamTabGroup, 'Title', 'Lattice / Beam');
            g = uigridlayout(app.LatticeTab, [11 2]);
            g.ColumnWidth = {140, '1x'};
            g.RowHeight   = [repmat({28},1,10), {'1x'}];
            g.Padding     = [8 8 8 8];
            g.RowSpacing  = 4;

            row = 0;
            [app.CField,       row] = addNumericRow(app, g, 'C  [m]',        354.097,  row);
            [app.E0Field,      row] = addNumericRow(app, g, 'E0  [eV]',      2.75e9,   row);
            [app.I0Field,      row] = addNumericRow(app, g, 'I0  [A]',       0.5,      row);
            [app.U0Field,      row] = addNumericRow(app, g, 'U0  [eV]',      470e3,    row);
            [app.hField,       row] = addNumericRow(app, g, 'h',             416,      row);
            [app.alpha_cField, row] = addNumericRow(app, g, [char(945) 'c'], 1.07e-4,  row);
            [app.tau_sField,   row] = addNumericRow(app, g, [char(964) 's  [s]'], 13.7e-3, row);
            [app.tau_zField,   row] = addNumericRow(app, g, [char(964) 'z  [s]'], 13.7e-3, row);
            [app.sigma_e0Field,row] = addNumericRow(app, g, [char(963) 'e0'],  0.92e-3, row);
            app.tau_sField.Tooltip = 'Longitudinal radiation damping time [s] (kept synchronized with tau_z for Fokker-Planck balance)';
            app.tau_zField.Tooltip = 'Longitudinal radiation damping time [s] (13.7 ms = SOLEIL II IDs-open worst-case for RF stability)';
            app.tau_sField.ValueChangedFcn = @(~,~) syncLongitudinalDampingTime(app, 'tau_s');
            app.tau_zField.ValueChangedFcn = @(~,~) syncLongitudinalDampingTime(app, 'tau_z');

            row = row + 1;
            lbl = uilabel(g); lbl.Layout.Row = row; lbl.Layout.Column = 1;
            lbl.Text = [char(963) 't0  [ps]']; lbl.FontSize = 12;
            app.sigma_t0Label = uilabel(g);
            app.sigma_t0Label.Layout.Row = row; app.sigma_t0Label.Layout.Column = 2;
            app.sigma_t0Label.Text = '-- (auto)'; app.sigma_t0Label.FontSize = 12;
            app.sigma_t0Label.FontColor = [0.1 0.4 0.7];
        end

        % -----------------------------------------------------------------
        %  RF System Tab (Combines Main Cavity, Harmonic Cavity, and HOMs)
        % -----------------------------------------------------------------
        function createRFSystemTab(app)
            app.RFSystemTab = uitab(app.ParamTabGroup, 'Title', 'RF System');
            
            % Outer grid to hold the three panels vertically
            outerGrid = uigridlayout(app.RFSystemTab, [3 1]);
            outerGrid.ColumnWidth = {'1x'};
            outerGrid.RowHeight = {250, 280, 220};
            outerGrid.Padding = [4 4 4 4];
            outerGrid.RowSpacing = 6;
            outerGrid.Scrollable = 'on';
            
            % 1. Main Cavity Panel
            mcPanel = uipanel(outerGrid, 'Title', 'Main Cavity');
            mcPanel.Layout.Row = 1; mcPanel.Layout.Column = 1;
            mcPanel.FontSize = 11; mcPanel.FontWeight = 'bold';
            
            g_mc = uigridlayout(mcPanel, [9 2]);
            g_mc.ColumnWidth = {140, '1x'};
            g_mc.RowHeight   = repmat({24}, 1, 9);
            g_mc.Padding     = [8 4 8 4];
            g_mc.RowSpacing  = 2;
            
            row = 0;
            [app.V_mcField,        row] = addNumericRow(app, g_mc, 'V_mc  [V]',               1.7e6,   row);
            [app.Ncav_mcField,     row] = addNumericRow(app, g_mc, 'N_cav (MC)',             4,       row);
            app.Ncav_mcField.Tooltip    = 'Number of Main Cavities (e.g. 4)';
            app.Ncav_mcField.ValueChangedFcn = @(~,~) MCInputChanged(app);

            [app.Q0_mcField,       row] = addNumericRow(app, g_mc, 'Q0_mc  (Unloaded)',       35.7e3,  row);
            app.Q0_mcField.Tooltip      = 'Unloaded quality factor Q0 of single cavity';
            app.Q0_mcField.ValueChangedFcn = @(~,~) MCInputChanged(app);

            [app.RoQ_mcField,      row] = addNumericRow(app, g_mc, ['R/Q_mc  [' char(937) ']'], 140,    row);
            app.RoQ_mcField.Tooltip     = 'Geometric R/Q per single cavity [Ohm]';
            app.RoQ_mcField.ValueChangedFcn = @(~,~) MCInputChanged(app);

            [app.beta_mcField,     row] = addNumericRow(app, g_mc, [char(946) '_mc  (Coupler)'], 4.95, row);
            app.beta_mcField.Tooltip    = 'Power coupler coupling coefficient beta. Q_L = Q0 / (1 + beta)';
            app.beta_mcField.ValueChangedFcn = @(~,~) MCInputChanged(app);

            [app.Q_mcField,        row] = addNumericRow(app, g_mc, 'Q_mc  (Loaded Q_L)',     6e3,     row);
            app.Q_mcField.Tooltip       = 'Loaded quality factor Q_L = Q0 / (1 + beta). Editing this updates beta.';
            app.Q_mcField.BackgroundColor = [0.94 0.94 0.96];
            app.Q_mcField.FontColor = [0.2 0.2 0.5];
            app.Q_mcField.ValueChangedFcn = @(~,~) MCLoadedQChanged(app);

            [app.R_mcField,        row] = addNumericRow(app, g_mc, ['R_mc  (Total) [' char(937) ']'], 3.36e6, row);
            app.R_mcField.Tooltip       = 'Total loaded shunt impedance R_L = Ncav * Q_L * (R/Q) [Ohm]';
            app.R_mcField.Editable = 'off';
            app.R_mcField.BackgroundColor = [0.90 0.90 0.90];
            app.R_mcField.FontColor = [0.35 0.35 0.35];

            [app.fre_shift_mcField,row] = addNumericRow(app, g_mc, 'f_mc  [Hz]',             -35e3,   row);
            
            row = row + 1;
            app.MCAutoDetuneCheckbox = uicheckbox(g_mc);
            app.MCAutoDetuneCheckbox.Layout.Row = row; app.MCAutoDetuneCheckbox.Layout.Column = [1 2];
            app.MCAutoDetuneCheckbox.Text = ['Auto-compute optimal ' char(916) 'f_mc (tuner)'];
            app.MCAutoDetuneCheckbox.Value = true;
            app.MCAutoDetuneCheckbox.FontSize = 11;
            app.MCAutoDetuneCheckbox.Tooltip = 'Automatically compute optimal MC detuning per beam current to minimize reflected power (like a tuner in mbtrack2)';
            app.MCAutoDetuneCheckbox.ValueChangedFcn = @(~,~) MCAutoDetuneToggled(app);

            % 2. Harmonic Cavity Panel
            hhcPanel = uipanel(outerGrid, 'Title', 'Harmonic Cavity');
            hhcPanel.Layout.Row = 2; hhcPanel.Layout.Column = 1;
            hhcPanel.FontSize = 11; hhcPanel.FontWeight = 'bold';
            
            g_hhc = uigridlayout(hhcPanel, [11 2]);
            g_hhc.ColumnWidth = {140, '1x'};
            g_hhc.RowHeight   = repmat({23}, 1, 11);
            g_hhc.Padding     = [8 4 8 4];
            g_hhc.RowSpacing  = 2;
            
            row = 0;
            [app.n_hcField,        row] = addNumericRow(app, g_hhc, 'n_hc',                  4,       row);
            [app.Ncav_hcField,     row] = addNumericRow(app, g_hhc, 'N_cav (HC)',            2,       row);
            app.Ncav_hcField.Tooltip    = 'Number of Harmonic Cavities (e.g. 2)';
            app.Ncav_hcField.ValueChangedFcn = @(~,~) HCInputChanged(app);

            [app.Q0_hcField,       row] = addNumericRow(app, g_hhc, 'Q0_hc  (Unloaded)',      31e3,    row);
            app.Q0_hcField.Tooltip      = 'Unloaded quality factor Q0 of harmonic cavity (warm copper)';
            app.Q0_hcField.ValueChangedFcn = @(~,~) HCInputChanged(app);

            [app.RoQ_hcField,      row] = addNumericRow(app, g_hhc, ['R/Q_hc  [' char(937) ']'], 30,  row);
            app.RoQ_hcField.Tooltip     = 'Geometric R/Q per single HC [Ohm]';
            app.RoQ_hcField.ValueChangedFcn = @(~,~) HCInputChanged(app);

            [app.beta_hcField,     row] = addNumericRow(app, g_hhc, [char(946) '_hc  (Coupler)'], 0,   row);
            app.beta_hcField.Tooltip    = 'Coupler beta (default 0 for passive cavity without external power source)';
            app.beta_hcField.ValueChangedFcn = @(~,~) HCInputChanged(app);

            [app.Q_hcField,        row] = addNumericRow(app, g_hhc, 'Q_hc  (Loaded Q_L)',    31e3,    row);
            app.Q_hcField.Tooltip       = 'Loaded quality factor Q_L = Q0 / (1 + beta). Editing this updates beta.';
            app.Q_hcField.BackgroundColor = [0.94 0.94 0.96];
            app.Q_hcField.FontColor = [0.2 0.2 0.5];
            app.Q_hcField.ValueChangedFcn = @(~,~) HCLoadedQChanged(app);
            
            row = row + 1;
            lbl = uilabel(g_hhc); lbl.Layout.Row = row; lbl.Layout.Column = 1;
            lbl.Text = ['R_hc (Total) [M' char(937) ']']; lbl.FontSize = 12;
            app.R_hcField = uieditfield(g_hhc, 'numeric');
            app.R_hcField.Layout.Row = row; app.R_hcField.Layout.Column = 2;
            app.R_hcField.Value = 1.86;
            app.R_hcField.ValueDisplayFormat = '%.6g';
            app.R_hcField.FontSize = 12;
            app.R_hcField.Editable = 'off';
            app.R_hcField.BackgroundColor = [0.90 0.90 0.90];
            app.R_hcField.FontColor = [0.35 0.35 0.35];
            app.R_hcField.Tooltip = 'Total loaded shunt impedance R_L = Ncav * Q_L * (R/Q) [MOhm]';
            
            row = row + 1;
            app.AutoDetuneCheckbox = uicheckbox(g_hhc);
            app.AutoDetuneCheckbox.Layout.Row = row; app.AutoDetuneCheckbox.Layout.Column = [1 2];
            app.AutoDetuneCheckbox.Text = ['Auto-compute optimal ' char(916) 'f_hc'];
            app.AutoDetuneCheckbox.Value = true;
            app.AutoDetuneCheckbox.FontSize = 11;
            app.AutoDetuneCheckbox.ValueChangedFcn = @(~,~) AutoDetuneToggled(app);

            row = row + 1;
            app.FormFactorIterCheckbox = uicheckbox(g_hhc);
            app.FormFactorIterCheckbox.Layout.Row = row; app.FormFactorIterCheckbox.Layout.Column = [1 2];
            app.FormFactorIterCheckbox.Text = ['Iterative Form Factor F_n(' char(963) '_' char(964) ')'];
            app.FormFactorIterCheckbox.Value = true;
            app.FormFactorIterCheckbox.FontSize = 11;
            app.FormFactorIterCheckbox.Tooltip = 'Self-consistently iterate stretched bunch form factor F_n(sigma_tau) for HHC/MC beam loading, synchronous phase, and flat-potential curves (essential for 4th-harmonic SOLEIL II)';
            app.FormFactorIterCheckbox.ValueChangedFcn = @(~,~) FormFactorIterToggled(app, 'panel');
            
            row = row + 1;
            lbl2 = uilabel(g_hhc); lbl2.Layout.Row = row; lbl2.Layout.Column = 1;
            lbl2.Text = [char(916) 'f_hc  [Hz]']; lbl2.FontSize = 12;
            app.fre_shift_hcField = uieditfield(g_hhc, 'numeric');
            app.fre_shift_hcField.Layout.Row = row; app.fre_shift_hcField.Layout.Column = 2;
            app.fre_shift_hcField.Value = 0;
            app.fre_shift_hcField.FontSize = 12;
            app.fre_shift_hcField.Enable = 'off';

            row = row + 1;
            app.OptimizeHCBtn = uibutton(g_hhc, 'push');
            app.OptimizeHCBtn.Layout.Row = row; app.OptimizeHCBtn.Layout.Column = [1 2];
            app.OptimizeHCBtn.Text = ['Optimize ' char(916) 'f_hc (Stability-Constrained)'];
            app.OptimizeHCBtn.FontSize = 11;
            app.OptimizeHCBtn.FontWeight = 'bold';
            app.OptimizeHCBtn.BackgroundColor = [0.15 0.50 0.80];
            app.OptimizeHCBtn.FontColor = [1 1 1];
            app.OptimizeHCBtn.Tooltip = 'Optimize HC detuning with physical stability penalty constraint (maximizing bunch lengthening while guaranteeing STABLE state)';
            app.OptimizeHCBtn.ButtonPushedFcn = @(~,~) OptimizeHCDetuningButtonPushed(app);
            
            % 3. HOMs Panel
            homsPanel = uipanel(outerGrid, 'Title', 'HOMs');
            homsPanel.Layout.Row = 3; homsPanel.Layout.Column = 1;
            homsPanel.FontSize = 11; homsPanel.FontWeight = 'bold';
            
            g_homs = uigridlayout(homsPanel, [3 1]);
            g_homs.RowHeight   = {'1x', 24, 30};
            g_homs.ColumnWidth = {'1x'};
            g_homs.Padding     = [6 6 6 6];
            g_homs.RowSpacing  = 4;
            
            app.HOMsTable = uitable(g_homs);
            app.HOMsTable.Layout.Row = 1; app.HOMsTable.Layout.Column = 1;
            app.HOMsTable.ColumnName = {'Q', ['R  [' char(937) ']'], 'f  [Hz]'};
            app.HOMsTable.ColumnWidth = {80, 110, '1x'};
            app.HOMsTable.ColumnEditable = [true true true];
            app.HOMsTable.ColumnFormat = {'numeric', 'numeric', 'numeric'};
            app.HOMsTable.Data = {};
            app.HOMsTable.FontSize = 11;

            app.HOMAutoAlignCheckbox = uicheckbox(g_homs);
            app.HOMAutoAlignCheckbox.Layout.Row = 2; app.HOMAutoAlignCheckbox.Layout.Column = 1;
            app.HOMAutoAlignCheckbox.Text = 'Auto-align to Worst-case Sideband';
            app.HOMAutoAlignCheckbox.Value = false;
            app.HOMAutoAlignCheckbox.FontSize = 11;
            app.HOMAutoAlignCheckbox.Tooltip = 'Align HOM frequencies to the upper synchrotron sideband (M*h+mu+nus)*f0 (matching mbtrack2 worst-case)';
            app.HOMAutoAlignCheckbox.ValueChangedFcn = @(~,~) HOMAutoAlignToggled(app);
            
            btnGrid = uigridlayout(g_homs, [1 3]);
            btnGrid.Layout.Row = 3; btnGrid.Layout.Column = 1;
            btnGrid.ColumnWidth = {85, 95, 105};
            btnGrid.Padding = [0 0 0 0];
            
            app.AddHOMBtn = uibutton(btnGrid, 'push');
            app.AddHOMBtn.Layout.Row = 1; app.AddHOMBtn.Layout.Column = 1;
            app.AddHOMBtn.Text = '+ Add';
            app.AddHOMBtn.FontSize = 11;
            app.AddHOMBtn.ButtonPushedFcn = @(~,~) AddHOMButtonPushed(app);
            
            app.RemoveHOMBtn = uibutton(btnGrid, 'push');
            app.RemoveHOMBtn.Layout.Row = 1; app.RemoveHOMBtn.Layout.Column = 2;
            app.RemoveHOMBtn.Text = '- Remove';
            app.RemoveHOMBtn.FontSize = 11;
            app.RemoveHOMBtn.ButtonPushedFcn = @(~,~) RemoveHOMButtonPushed(app);

            app.AlignHOMBtn = uibutton(btnGrid, 'push');
            app.AlignHOMBtn.Layout.Row = 1; app.AlignHOMBtn.Layout.Column = 3;
            app.AlignHOMBtn.Text = 'Align Now';
            app.AlignHOMBtn.FontSize = 11;
            app.AlignHOMBtn.Tooltip = 'Calculate and update table with exact worst-case sideband frequencies';
            app.AlignHOMBtn.ButtonPushedFcn = @(~,~) AlignHOMButtonPushed(app);
        end

        % -----------------------------------------------------------------
        %  Fill Pattern Tab
        % -----------------------------------------------------------------
        function createPatternTab(app)
            app.PatternTab = uitab(app.ParamTabGroup, 'Title', 'Fill Pattern');
            g = uigridlayout(app.PatternTab, [5 2]);
            g.ColumnWidth = {140, '1x'};
            g.RowHeight   = {28, 28, 28, 28, 28};
            g.Padding     = [8 8 8 8];
            g.RowSpacing  = 4;

            lbl = uilabel(g); lbl.Layout.Row = 1; lbl.Layout.Column = 1;
            lbl.Text = 'Pattern'; lbl.FontSize = 12;
            app.PatternDropDown = uidropdown(g);
            app.PatternDropDown.Layout.Row = 1; app.PatternDropDown.Layout.Column = 2;
            app.PatternDropDown.Items = {'Uniform (100%)', '80% filling', '90% filling', 'Single bunch', 'Custom'};
            app.PatternDropDown.Value = 'Uniform (100%)';
            app.PatternDropDown.FontSize = 12;

            lbl2 = uilabel(g); lbl2.Layout.Row = 2; lbl2.Layout.Column = 1;
            lbl2.Text = 'Custom expr.'; lbl2.FontSize = 12;
            app.CustomPatternField = uieditfield(g, 'text');
            app.CustomPatternField.Layout.Row = 2; app.CustomPatternField.Layout.Column = 2;
            app.CustomPatternField.Value = 'ones(1, h)';
            app.CustomPatternField.Enable = 'off';
            app.CustomPatternField.FontSize = 12;

            % Spacer row 3
            lbl3 = uilabel(g); lbl3.Layout.Row = 3; lbl3.Layout.Column = 1;
            lbl3.Text = ''; 

            lbl4 = uilabel(g); lbl4.Layout.Row = 4; lbl4.Layout.Column = 1;
            lbl4.Text = 'Fill rate'; lbl4.FontSize = 12;
            app.FillRateLabel = uilabel(g);
            app.FillRateLabel.Layout.Row = 4; app.FillRateLabel.Layout.Column = 2;
            app.FillRateLabel.Text = '100.0 %'; app.FillRateLabel.FontSize = 12;
            app.FillRateLabel.FontColor = [0.1 0.4 0.7];

            lbl5 = uilabel(g); lbl5.Layout.Row = 5; lbl5.Layout.Column = 1;
            lbl5.Text = 'Bunches'; lbl5.FontSize = 12;
            app.BunchNumLabel = uilabel(g);
            app.BunchNumLabel.Layout.Row = 5; app.BunchNumLabel.Layout.Column = 2;
            app.BunchNumLabel.Text = '416'; app.BunchNumLabel.FontSize = 12;
            app.BunchNumLabel.FontColor = [0.1 0.4 0.7];
        end

        % -----------------------------------------------------------------
        %  Tracking Tab
        % -----------------------------------------------------------------
        function createTrackingTab(app)
            app.TrackingTab = uitab(app.ParamTabGroup, 'Title', 'Tracking');
            g = uigridlayout(app.TrackingTab, [5 2]);
            g.ColumnWidth = {140, '1x'};
            g.RowHeight   = {28, 28, 28, 28, '1x'};
            g.Padding     = [8 8 8 8];
            g.RowSpacing  = 4;

            row = 0;
            [app.Par_numField,   row] = addNumericRow(app, g, 'Particles/bunch', 1e4,   row);
            [app.Track_numField, row] = addNumericRow(app, g, 'Tracking turns',  1e5,   row);
            [app.Recor_stepField,row] = addNumericRow(app, g, 'Record step',     10,    row);

            row = row + 1;
            app.ShortRangeCheckbox = uicheckbox(g);
            app.ShortRangeCheckbox.Layout.Row = row; app.ShortRangeCheckbox.Layout.Column = 1;
            app.ShortRangeCheckbox.Text = 'Short-range wake';
            app.ShortRangeCheckbox.Value = true;
            app.ShortRangeCheckbox.FontSize = 12;
        end

        % -----------------------------------------------------------------
        %  Feedback Tab (PI and Direct RF Feedback)
        % -----------------------------------------------------------------
        function createFeedbackTab(app)
            app.FeedbackTab = uitab(app.ParamTabGroup, 'Title', 'Feedback');
            app.PITab = app.FeedbackTab; % Alias for backward compatibility

            outerGrid = uigridlayout(app.FeedbackTab, [2 1]);
            outerGrid.ColumnWidth = {'1x'};
            outerGrid.RowHeight = {180, 310};
            outerGrid.Padding = [4 4 4 4];
            outerGrid.RowSpacing = 6;
            outerGrid.Scrollable = 'on';

            % 1. PI Feedback Panel
            piPanel = uipanel(outerGrid, 'Title', 'PI Feedback');
            piPanel.Layout.Row = 1; piPanel.Layout.Column = 1;
            piPanel.FontSize = 11; piPanel.FontWeight = 'bold';

            g_pi = uigridlayout(piPanel, [5 2]);
            g_pi.ColumnWidth = {140, '1x'};
            g_pi.RowHeight   = {26, 26, 26, 26, 26};
            g_pi.Padding     = [8 6 8 6];
            g_pi.RowSpacing  = 4;

            row = 1;
            app.PIEnableCheckbox = uicheckbox(g_pi);
            app.PIEnableCheckbox.Layout.Row = row; app.PIEnableCheckbox.Layout.Column = [1 2];
            app.PIEnableCheckbox.Text = 'Enable PI Feedback';
            app.PIEnableCheckbox.Value = true;
            app.PIEnableCheckbox.FontSize = 12;
            app.PIEnableCheckbox.ValueChangedFcn = @(~,~) PIToggle(app);

            [app.KPField,   row] = addNumericRow(app, g_pi, 'KP',    1,     row);
            [app.KIField,   row] = addNumericRow(app, g_pi, 'KI',    1e-4,  row);
            [app.PI_mField, row] = addNumericRow(app, g_pi, 'PI.m  (buckets)', 80, row);
            [app.PI_dField, ~  ] = addNumericRow(app, g_pi, 'PI.d  (delay)',   800, row);

            % 2. Direct RF Feedback (DFB) Panel
            dfbPanel = uipanel(outerGrid, 'Title', 'Direct RF Feedback (DFB)');
            dfbPanel.Layout.Row = 2; dfbPanel.Layout.Column = 1;
            dfbPanel.FontSize = 11; dfbPanel.FontWeight = 'bold';

            g_dfb = uigridlayout(dfbPanel, [9 2]);
            g_dfb.ColumnWidth = {140, '1x'};
            g_dfb.RowHeight   = repmat({24}, 1, 10);
            g_dfb.Padding     = [8 6 8 6];
            g_dfb.RowSpacing  = 3;

            row = 1;
            app.DFBEnableCheckbox = uicheckbox(g_dfb);
            app.DFBEnableCheckbox.Layout.Row = row; app.DFBEnableCheckbox.Layout.Column = [1 2];
            app.DFBEnableCheckbox.Text = 'Enable Direct RF Feedback';
            app.DFBEnableCheckbox.Value = false;
            app.DFBEnableCheckbox.FontSize = 12;
            app.DFBEnableCheckbox.ValueChangedFcn = @(~,~) DFBToggle(app);

            row = row + 1;
            app.DFB_adaptiveCheckbox = uicheckbox(g_dfb);
            app.DFB_adaptiveCheckbox.Layout.Row = row; app.DFB_adaptiveCheckbox.Layout.Column = [1 2];
            app.DFB_adaptiveCheckbox.Text = 'Adaptive Optimal DFB Controller (AI/Physics)';
            app.DFB_adaptiveCheckbox.Value = false;
            app.DFB_adaptiveCheckbox.FontSize = 11;
            app.DFB_adaptiveCheckbox.FontWeight = 'bold';
            app.DFB_adaptiveCheckbox.Tooltip = 'Automatically computes optimal DFB gain dynamically for each beam current based on cavity beam loading (Y = I0*R_mc/V_mc) and loop delay phase margin';
            app.DFB_adaptiveCheckbox.ValueChangedFcn = @(~,~) DFBAdaptiveToggled(app);

            row = row + 1;
            app.DFB_auto_scaleCheckbox = uicheckbox(g_dfb);
            app.DFB_auto_scaleCheckbox.Layout.Row = row; app.DFB_auto_scaleCheckbox.Layout.Column = [1 2];
            app.DFB_auto_scaleCheckbox.Text = 'Auto-scale Gain with Beam Current (I_0)';
            app.DFB_auto_scaleCheckbox.Value = false;
            app.DFB_auto_scaleCheckbox.FontSize = 11;
            app.DFB_auto_scaleCheckbox.Tooltip = 'Linearly interpolate/scale DFB gain according to beam current I_0';
            app.DFB_auto_scaleCheckbox.ValueChangedFcn = @(~,~) DFBToggle(app);

            [app.DFB_gainField,        row] = addNumericRow(app, g_dfb, 'Constant Gain (G_dfb)', 0.2, row);

            row = row + 1;
            lbl1 = uilabel(g_dfb); lbl1.Layout.Row = row; lbl1.Layout.Column = 1;
            lbl1.Text = 'Ref 1: I0 [mA] / Gain'; lbl1.FontSize = 11;
            subG1 = uigridlayout(g_dfb, [1 2]);
            subG1.Layout.Row = row; subG1.Layout.Column = 2;
            subG1.ColumnWidth = {'1x', '1x'};
            subG1.Padding = [0 0 0 0]; subG1.ColumnSpacing = 4;
            app.DFB_I1Field = uieditfield(subG1, 'numeric');
            app.DFB_I1Field.Value = 200; app.DFB_I1Field.Tooltip = 'Reference current 1 in mA (e.g. 200 mA)';
            app.DFB_G1Field = uieditfield(subG1, 'numeric');
            app.DFB_G1Field.Value = 0.6; app.DFB_G1Field.Tooltip = 'DFB gain at reference current 1 (e.g. 0.6)';

            row = row + 1;
            lbl2 = uilabel(g_dfb); lbl2.Layout.Row = row; lbl2.Layout.Column = 1;
            lbl2.Text = 'Ref 2: I0 [mA] / Gain'; lbl2.FontSize = 11;
            subG2 = uigridlayout(g_dfb, [1 2]);
            subG2.Layout.Row = row; subG2.Layout.Column = 2;
            subG2.ColumnWidth = {'1x', '1x'};
            subG2.Padding = [0 0 0 0]; subG2.ColumnSpacing = 4;
            app.DFB_I2Field = uieditfield(subG2, 'numeric');
            app.DFB_I2Field.Value = 500; app.DFB_I2Field.Tooltip = 'Reference current 2 in mA (e.g. 500 mA)';
            app.DFB_G2Field = uieditfield(subG2, 'numeric');
            app.DFB_G2Field.Value = 0.1; app.DFB_G2Field.Tooltip = 'DFB gain at reference current 2 (e.g. 0.1)';

            [app.DFB_phase_shiftField, row] = addNumericRow(app, g_dfb, 'Phase Shift [deg]', 0, row);
            [app.DFB_delayField,       row] = addNumericRow(app, g_dfb, 'Delay (buckets)', 800, row);
            [app.DFB_mField,           row] = addNumericRow(app, g_dfb, 'Sample Window (m)', 1, row);
            [app.DFB_everyField,       row] = addNumericRow(app, g_dfb, 'Update Every (buckets)', 1, row);

            PIToggle(app);
            DFBToggle(app);
        end

        function createPITab(app)
            createFeedbackTab(app);
        end



        % -----------------------------------------------------------------
        %  Scan Tab
        % -----------------------------------------------------------------
        function createScanTab(app)
            app.ScanTab = uitab(app.ParamTabGroup, 'Title', 'Scan Map');
            g = uigridlayout(app.ScanTab, [10 2]);
            g.ColumnWidth = {160, '1x'};
            g.RowHeight   = [repmat({28},1,6), {10}, {28}, {28}, {'1x'}];
            g.Padding     = [8 8 8 8];
            g.RowSpacing  = 4;

            row = 0;
            [app.DfMinField, row] = addNumericRow(app, g, 'HHC Detune Min [kHz]', 40, row);
            [app.DfMaxField, row] = addNumericRow(app, g, 'HHC Detune Max [kHz]', 120, row);
            [app.DfStepField,row] = addNumericRow(app, g, 'HHC Detune Step [kHz]', 5, row);
            [app.I0MinField, row] = addNumericRow(app, g, 'Current Min [mA]', 50, row);
            [app.I0MaxField, row] = addNumericRow(app, g, 'Current Max [mA]', 500, row);
            [app.I0StepField,row] = addNumericRow(app, g, 'Current Step [mA]', 50, row);

            % --- Save Folder ---
            % Row 7: separator (spacer)
            % Row 8: Save Folder label + field + browse button
            row = 8;
            lbl = uilabel(g); lbl.Layout.Row = row; lbl.Layout.Column = 1;
            lbl.Text = 'Save Folder'; lbl.FontSize = 12;

            folderGrid = uigridlayout(g, [1 2]);
            folderGrid.Layout.Row = row; folderGrid.Layout.Column = 2;
            folderGrid.ColumnWidth = {'1x', 60};
            folderGrid.Padding = [0 0 0 0]; folderGrid.ColumnSpacing = 4;

            app.ScanSaveFolderField = uieditfield(folderGrid, 'text');
            app.ScanSaveFolderField.Layout.Row = 1; app.ScanSaveFolderField.Layout.Column = 1;
            app.ScanSaveFolderField.Value = '';
            app.ScanSaveFolderField.Placeholder = '(auto: prompt at scan)';
            app.ScanSaveFolderField.FontSize = 11;

            app.BrowseFolderBtn = uibutton(folderGrid, 'push');
            app.BrowseFolderBtn.Layout.Row = 1; app.BrowseFolderBtn.Layout.Column = 2;
            app.BrowseFolderBtn.Text = 'Browse';
            app.BrowseFolderBtn.FontSize = 11;
            app.BrowseFolderBtn.ButtonPushedFcn = @(~,~) BrowseFolderButtonPushed(app);

            row = 9;
            app.FastPreCheckCheckbox = uicheckbox(g);
            app.FastPreCheckCheckbox.Layout.Row = row; app.FastPreCheckCheckbox.Layout.Column = [1 2];
            app.FastPreCheckCheckbox.Text = 'Enable Fast Analytical Pre-check (D-mode & PTBL)';
            app.FastPreCheckCheckbox.Value = true;
            app.FastPreCheckCheckbox.FontSize = 12;
            app.FastPreCheckCheckbox.Tooltip = 'If enabled, skip tracking for points pre-determined to be unstable via D-mode Robinson or PTBL SSRM analytical checks.';
        end

        % -----------------------------------------------------------------
        %  Result Tabs (right panel)
        % -----------------------------------------------------------------
        function createResultTabs(app)
            % --- Tab 1: Centroid & Sigma vs Turns ---
            app.CentroidTab = uitab(app.ResultTabGroup, 'Title', 'Centroid & ' + string(char(963)));
            app.ResultTabGroup.TabLocation = 'left';
            cg = uigridlayout(app.CentroidTab, [2 2]);
            cg.Padding = [4 4 4 4]; cg.RowSpacing = 4; cg.ColumnSpacing = 4;

            app.AxesTauMean = uiaxes(cg);
            app.AxesTauMean.Layout.Row = 1; app.AxesTauMean.Layout.Column = 1;
            title(app.AxesTauMean, '<\tau>  [ps]'); xlabel(app.AxesTauMean, 'Turns');
            app.AxesTauMean.FontName = 'Times New Roman'; app.AxesTauMean.FontSize = 11;

            app.AxesTauStd = uiaxes(cg);
            app.AxesTauStd.Layout.Row = 1; app.AxesTauStd.Layout.Column = 2;
            title(app.AxesTauStd, '\sigma_{\tau}  [ps]'); xlabel(app.AxesTauStd, 'Turns');
            app.AxesTauStd.FontName = 'Times New Roman'; app.AxesTauStd.FontSize = 11;

            app.AxesDeltaMean = uiaxes(cg);
            app.AxesDeltaMean.Layout.Row = 2; app.AxesDeltaMean.Layout.Column = 1;
            title(app.AxesDeltaMean, '<\delta>'); xlabel(app.AxesDeltaMean, 'Turns');
            app.AxesDeltaMean.FontName = 'Times New Roman'; app.AxesDeltaMean.FontSize = 11;

            app.AxesDeltaStd = uiaxes(cg);
            app.AxesDeltaStd.Layout.Row = 2; app.AxesDeltaStd.Layout.Column = 2;
            title(app.AxesDeltaStd, '\sigma_{\delta}'); xlabel(app.AxesDeltaStd, 'Turns');
            app.AxesDeltaStd.FontName = 'Times New Roman'; app.AxesDeltaStd.FontSize = 11;

            % --- Tab 2: Bunch-by-Bunch ---
            app.BunchTab = uitab(app.ResultTabGroup, 'Title', 'Bunch-by-Bunch');
            bg = uigridlayout(app.BunchTab, [2 1]);
            bg.Padding = [4 4 4 4]; bg.RowSpacing = 4;

            app.AxesBunchSigma = uiaxes(bg);
            app.AxesBunchSigma.Layout.Row = 1;
            title(app.AxesBunchSigma, '\sigma_{\tau}  [ps]'); xlabel(app.AxesBunchSigma, 'Bunch number');
            app.AxesBunchSigma.FontName = 'Times New Roman'; app.AxesBunchSigma.FontSize = 11;

            app.AxesBunchCentroid = uiaxes(bg);
            app.AxesBunchCentroid.Layout.Row = 2;
            title(app.AxesBunchCentroid, '<\tau>  [ps]'); xlabel(app.AxesBunchCentroid, 'Bunch number');
            app.AxesBunchCentroid.FontName = 'Times New Roman'; app.AxesBunchCentroid.FontSize = 11;

            % --- Tab 3: Density Profile ---
            app.DensityTab = uitab(app.ResultTabGroup, 'Title', 'Density');
            dg = uigridlayout(app.DensityTab, [1 1]);
            dg.Padding = [4 4 4 4];

            app.AxesDensity = uiaxes(dg);
            title(app.AxesDensity, 'Longitudinal Density Profile');
            xlabel(app.AxesDensity, '\tau  [ps]'); ylabel(app.AxesDensity, 'Norm. density');
            app.AxesDensity.FontName = 'Times New Roman'; app.AxesDensity.FontSize = 12;

            % --- Tab 4: HC Voltage ---
            app.HCVoltTab = uitab(app.ResultTabGroup, 'Title', 'HC Voltage');
            hg = uigridlayout(app.HCVoltTab, [2 1]);
            hg.Padding = [4 4 4 4]; hg.RowSpacing = 4;

            app.AxesHCAmp = uiaxes(hg);
            app.AxesHCAmp.Layout.Row = 1;
            title(app.AxesHCAmp, 'Harmonic Cavity Amplitude'); ylabel(app.AxesHCAmp, '|V|  [kV]');
            xlabel(app.AxesHCAmp, 'Bucket ID');
            app.AxesHCAmp.FontName = 'Times New Roman'; app.AxesHCAmp.FontSize = 11;

            app.AxesHCPhase = uiaxes(hg);
            app.AxesHCPhase.Layout.Row = 2;
            title(app.AxesHCPhase, 'Harmonic Cavity Phase'); ylabel(app.AxesHCPhase, 'Phase  [deg]');
            xlabel(app.AxesHCPhase, 'Bucket ID');
            app.AxesHCPhase.FontName = 'Times New Roman'; app.AxesHCPhase.FontSize = 11;

            % --- Tab 5: MC Voltage ---
            app.MCVoltTab = uitab(app.ResultTabGroup, 'Title', 'MC Voltage');
            mg = uigridlayout(app.MCVoltTab, [2 1]);
            mg.Padding = [4 4 4 4]; mg.RowSpacing = 4;

            app.AxesMCAmp = uiaxes(mg);
            app.AxesMCAmp.Layout.Row = 1;
            title(app.AxesMCAmp, 'Main Cavity Amplitude'); ylabel(app.AxesMCAmp, '|V|  [kV]');
            xlabel(app.AxesMCAmp, 'Bucket ID');
            app.AxesMCAmp.FontName = 'Times New Roman'; app.AxesMCAmp.FontSize = 11;

            app.AxesMCPhase = uiaxes(mg);
            app.AxesMCPhase.Layout.Row = 2;
            title(app.AxesMCPhase, 'Main Cavity Phase'); ylabel(app.AxesMCPhase, 'Phase  [deg]');
            xlabel(app.AxesMCPhase, 'Bucket ID');
            app.AxesMCPhase.FontName = 'Times New Roman'; app.AxesMCPhase.FontSize = 11;

            % --- Tab 6: Generator ---
            app.GenTab = uitab(app.ResultTabGroup, 'Title', 'Generator');
            gg = uigridlayout(app.GenTab, [2 1]);
            gg.Padding = [4 4 4 4]; gg.RowSpacing = 4;

            app.AxesGenAmp = uiaxes(gg);
            app.AxesGenAmp.Layout.Row = 1;
            title(app.AxesGenAmp, 'Generator Current Amplitude');
            ylabel(app.AxesGenAmp, '|I_g|  [A]'); xlabel(app.AxesGenAmp, 'Sample');
            app.AxesGenAmp.FontName = 'Times New Roman'; app.AxesGenAmp.FontSize = 11;

            app.AxesGenPhase = uiaxes(gg);
            app.AxesGenPhase.Layout.Row = 2;
            title(app.AxesGenPhase, 'Generator Current Phase');
            ylabel(app.AxesGenPhase, 'Phase  [rad]'); xlabel(app.AxesGenPhase, 'Sample');
            app.AxesGenPhase.FontName = 'Times New Roman'; app.AxesGenPhase.FontSize = 11;

            % --- Tab 7: Stability Map ---
            app.MapTab = uitab(app.ResultTabGroup, 'Title', 'Stability Map');
            mg2 = uigridlayout(app.MapTab, [2 1]);
            mg2.RowHeight = {30, '1x'};
            mg2.Padding = [6 6 6 6];
            mg2.RowSpacing = 6;

            mapTopGrid = uigridlayout(mg2, [1 6]);
            mapTopGrid.Layout.Row = 1;
            mapTopGrid.Padding = [0 0 0 0];
            mapTopGrid.ColumnWidth = {72, 200, 110, 82, 122, '1x'};

            lblMap = uilabel(mapTopGrid);
            lblMap.Layout.Column = 1;
            lblMap.Text = 'Coordinate:';
            lblMap.FontName = 'Times New Roman';
            lblMap.FontSize = 12;

            app.MapXAxisDropDown = uidropdown(mapTopGrid);
            app.MapXAxisDropDown.Layout.Column = 2;
            app.MapXAxisDropDown.Items = {['HHC Detuning ' char(916) 'f_HC [kHz]'], ['HHC Tuning Angle ' char(968) '_HC [' char(176) ']']};
            app.MapXAxisDropDown.Value = app.MapXAxisDropDown.Items{1};
            app.MapXAxisDropDown.FontSize = 11;
            app.MapXAxisDropDown.ValueChangedFcn = @(~,~) updateStabilityMapView(app);

            app.MapLoadBtn = uibutton(mapTopGrid, 'push');
            app.MapLoadBtn.Layout.Column = 3;
            app.MapLoadBtn.Text = 'Load Scan Data...';
            app.MapLoadBtn.FontSize = 11;
            app.MapLoadBtn.ButtonPushedFcn = @(~,~) loadScanDataForMap(app);

            app.MapRefreshBtn = uibutton(mapTopGrid, 'push');
            app.MapRefreshBtn.Layout.Column = 4;
            app.MapRefreshBtn.Text = 'Refresh Map';
            app.MapRefreshBtn.FontSize = 11;
            app.MapRefreshBtn.ButtonPushedFcn = @(~,~) updateStabilityMapView(app);

            app.MapFitBoundaryBtn = uibutton(mapTopGrid, 'push');
            app.MapFitBoundaryBtn.Layout.Column = 5;
            app.MapFitBoundaryBtn.Text = 'Fit Safety Boundary';
            app.MapFitBoundaryBtn.FontSize = 11;
            app.MapFitBoundaryBtn.FontWeight = 'bold';
            app.MapFitBoundaryBtn.BackgroundColor = [0.85 0.35 0.10];
            app.MapFitBoundaryBtn.FontColor = [1 1 1];
            app.MapFitBoundaryBtn.Tooltip = 'Automatically extract critical stability boundary and fit engineering equation \Delta f_{crit}(I_0)';
            app.MapFitBoundaryBtn.ButtonPushedFcn = @(~,~) fitAndPlotSafetyBoundary(app);

            app.MapFormFactorCheckbox = uicheckbox(mapTopGrid);
            app.MapFormFactorCheckbox.Layout.Column = 6;
            app.MapFormFactorCheckbox.Text = ['Form Factor F_n(' char(963) '_' char(964) ')'];
            app.MapFormFactorCheckbox.Value = true;
            app.MapFormFactorCheckbox.FontSize = 11;
            app.MapFormFactorCheckbox.Tooltip = 'Include self-consistent stretched bunch form factor F_n(sigma_tau) in Flat-pot. (xi=1.0) & xi=1.05 reference curves';
            app.MapFormFactorCheckbox.ValueChangedFcn = @(~,~) FormFactorIterToggled(app, 'map');

            app.AxesMap = uiaxes(mg2);
            app.AxesMap.Layout.Row = 2;
            title(app.AxesMap, '2D Stability Map');
            xlabel(app.AxesMap, 'HHC Detuning \Delta f_{HC} [kHz]');
            ylabel(app.AxesMap, 'Beam Current I_0 [mA]');
            app.AxesMap.FontName = 'Times New Roman'; app.AxesMap.FontSize = 12;
            app.AxesMap.Colormap = [0.2 0.63 0.17; 0.99 0.75 0.05; 0.85 0.33 0.1; 0.69 0.13 0.13; 0.5 0.5 0.5];

            % --- Tab 8: Offline Analysis ---
            app.AnalysisTab = uitab(app.ResultTabGroup, 'Title', 'Analysis');
            mg3 = uigridlayout(app.AnalysisTab, [2 1]);
            mg3.Padding = [8 8 8 8];
            mg3.RowHeight = {28, '1x'};

            app.AnalysisPlotModeDropDown = uidropdown(mg3);
            app.AnalysisPlotModeDropDown.Layout.Row = 1;
            app.AnalysisPlotModeDropDown.Layout.Column = 1;
            app.AnalysisPlotModeDropDown.Items = {'Bunch Length', 'Centroid', 'Energy Spread', 'Bunch Length + Centroid', 'Bunch Length + Energy Spread', 'Centroid + Energy Spread', 'All', 'Mode-0 Spectrum (D-mode Check)'};
            app.AnalysisPlotModeDropDown.Value = 'All';
            app.AnalysisPlotModeDropDown.FontSize = 11;
            app.AnalysisPlotModeDropDown.ValueChangedFcn = @(~,~) updateAnalysisPlotView(app);

            app.AxesAnalysis = uiaxes(mg3);
            app.AxesAnalysis.Layout.Row = 2;
            app.AxesAnalysis.Layout.Column = 1;
            title(app.AxesAnalysis, 'Bunch Lengthening effect vs. HHC Detuning');
            xlabel(app.AxesAnalysis, 'HHC Detuning \Delta f_{hc} [kHz]');
            app.AxesAnalysis.FontName = 'Times New Roman'; app.AxesAnalysis.FontSize = 12;
            app.AxesAnalysis.XGrid = 'on'; app.AxesAnalysis.YGrid = 'on';

            % --- Tab 9: LLRF Nyquist & Bode (Appendix D) ---
            app.LLRFTab = uitab(app.ResultTabGroup, 'Title', 'LLRF Nyquist & Bode');
            mg4 = uigridlayout(app.LLRFTab, [2 1]);
            mg4.Padding = [8 8 8 8];
            mg4.RowHeight = {32, '1x'};
            mg4.RowSpacing = 6;

            llrfTopGrid = uigridlayout(mg4, [1 5]);
            llrfTopGrid.Layout.Row = 1;
            llrfTopGrid.Layout.Column = 1;
            llrfTopGrid.Padding = [0 0 0 0];
            llrfTopGrid.ColumnSpacing = 8;
            llrfTopGrid.ColumnWidth = {75, 210, 165, 110, '1x'};

            lblLLRF = uilabel(llrfTopGrid);
            lblLLRF.Layout.Column = 1;
            lblLLRF.Text = 'Loop Mode:';
            lblLLRF.FontName = 'Times New Roman';
            lblLLRF.FontSize = 12;

            app.LLRFModeDropDown = uidropdown(llrfTopGrid);
            app.LLRFModeDropDown.Layout.Column = 2;
            app.LLRFModeDropDown.Items = { ...
                'Gain Sweep (Compare G_dfb)', ...
                'Combined (DFB + PI)', ...
                'Direct RF Feedback (DFB Only)', ...
                'PI Feedback (PI Only)'};
            app.LLRFModeDropDown.Value = 'Gain Sweep (Compare G_dfb)';
            app.LLRFModeDropDown.FontSize = 11;
            app.LLRFModeDropDown.ValueChangedFcn = @(~,~) updateLLRFNyquistPlot(app);

            app.LLRFIncludeDetuneCheckbox = uicheckbox(llrfTopGrid);
            app.LLRFIncludeDetuneCheckbox.Layout.Column = 3;
            app.LLRFIncludeDetuneCheckbox.Text = 'Include MC Detuning \psi';
            app.LLRFIncludeDetuneCheckbox.Value = false;
            app.LLRFIncludeDetuneCheckbox.FontSize = 11;
            app.LLRFIncludeDetuneCheckbox.ValueChangedFcn = @(~,~) updateLLRFNyquistPlot(app);

            app.LLRFUpdateBtn = uibutton(llrfTopGrid, 'push');
            app.LLRFUpdateBtn.Layout.Column = 4;
            app.LLRFUpdateBtn.Text = 'Update Plot';
            app.LLRFUpdateBtn.FontSize = 11;
            app.LLRFUpdateBtn.FontWeight = 'bold';
            app.LLRFUpdateBtn.BackgroundColor = [0.15 0.45 0.78];
            app.LLRFUpdateBtn.FontColor = [1 1 1];
            app.LLRFUpdateBtn.ButtonPushedFcn = @(~,~) updateLLRFNyquistPlot(app);

            app.LLRFMetricsLabel = uilabel(llrfTopGrid);
            app.LLRFMetricsLabel.Layout.Column = 5;
            app.LLRFMetricsLabel.Text = 'GM: -- dB | PM: -- deg | G_{max}: --';
            app.LLRFMetricsLabel.FontName = 'Times New Roman';
            app.LLRFMetricsLabel.FontSize = 11;
            app.LLRFMetricsLabel.FontWeight = 'bold';
            app.LLRFMetricsLabel.FontColor = [0.1 0.35 0.65];

            llrfPlotGrid = uigridlayout(mg4, [2 2]);
            llrfPlotGrid.Layout.Row = 2;
            llrfPlotGrid.Layout.Column = 1;
            llrfPlotGrid.Padding = [0 0 0 0];
            llrfPlotGrid.RowHeight = {'1x', '1x'};
            llrfPlotGrid.ColumnWidth = {'1x', '1.1x'};
            llrfPlotGrid.RowSpacing = 8;
            llrfPlotGrid.ColumnSpacing = 10;

            app.AxesBodeMag = uiaxes(llrfPlotGrid);
            app.AxesBodeMag.Layout.Row = 1;
            app.AxesBodeMag.Layout.Column = 1;
            title(app.AxesBodeMag, 'Bode Magnitude |L(j\omega_m)|');
            ylabel(app.AxesBodeMag, 'Magnitude [dB]');
            app.AxesBodeMag.FontName = 'Times New Roman'; app.AxesBodeMag.FontSize = 11;
            app.AxesBodeMag.XGrid = 'on'; app.AxesBodeMag.YGrid = 'on';

            app.AxesBodePhase = uiaxes(llrfPlotGrid);
            app.AxesBodePhase.Layout.Row = 2;
            app.AxesBodePhase.Layout.Column = 1;
            title(app.AxesBodePhase, 'Bode Phase \angle L(j\omega_m)');
            xlabel(app.AxesBodePhase, 'Modulation Frequency f_m = (\omega - \omega_{rf})/2\pi  [kHz]');
            ylabel(app.AxesBodePhase, 'Phase [deg]');
            app.AxesBodePhase.FontName = 'Times New Roman'; app.AxesBodePhase.FontSize = 11;
            app.AxesBodePhase.XGrid = 'on'; app.AxesBodePhase.YGrid = 'on';

            app.AxesNyquist = uiaxes(llrfPlotGrid);
            app.AxesNyquist.Layout.Row = [1 2];
            app.AxesNyquist.Layout.Column = 2;
            title(app.AxesNyquist, 'Nyquist Diagram of Open-Loop L(j\omega_m)');
            xlabel(app.AxesNyquist, 'Re[L(j\omega_m)]');
            ylabel(app.AxesNyquist, 'Im[L(j\omega_m)]');
            app.AxesNyquist.FontName = 'Times New Roman'; app.AxesNyquist.FontSize = 11;
            app.AxesNyquist.XGrid = 'on'; app.AxesNyquist.YGrid = 'on';

            app.ResultTabGroup.SelectionChangedFcn = @(~, evt) onResultTabChanged(app, evt);

            updateAnalysisPlotView(app);
        end

        function onResultTabChanged(app, evt)
            try
                if isequal(evt.NewValue, app.LLRFTab)
                    updateLLRFNyquistPlot(app);
                elseif isequal(evt.NewValue, app.MapTab)
                    updateStabilityMapView(app);
                end
            catch
            end
        end

        function updateLLRFNyquistPlot(app)
            % Appendix D: LLRF Open-Loop Transfer Function, Bode & Nyquist Stability Analysis
            try
                cspeed = 299792458;
                C_val  = app.CField.Value;
                h_val  = round(app.hField.Value);
                I0_val = app.I0Field.Value;
                T0_val = C_val / cspeed;
                Tb_val = T0_val / h_val;
                f0_val = 1.0 / T0_val;
                f_rf   = h_val * f0_val;
                w_rf   = 2 * pi * f_rf;

                % Main Cavity parameters
                Q_mc   = max(1, app.Q_mcField.Value);
                R_mc   = max(1, app.R_mcField.Value);
                V_mc   = max(1, app.V_mcField.Value);
                U0_val = app.U0Field.Value;
                w_half = w_rf / (2 * Q_mc); % Cavity half-bandwidth [rad/s]
                f_half_kHz = (w_half / (2 * pi)) * 1e-3;

                % Cavity detuning angle psi
                if app.LLRFIncludeDetuneCheckbox.Value
                    if app.MCAutoDetuneCheckbox.Value
                        n_hc = round(app.n_hcField.Value);
                        R_hc = app.R_hcField.Value * 1e6;
                        Q_hc = app.Q_hcField.Value;
                        df_hc = app.fre_shift_hcField.Value;
                        df_mc = computeOptimalMCDetune(app, I0_val, V_mc, U0_val, R_mc, Q_mc, R_hc, Q_hc, df_hc, n_hc);
                    else
                        df_mc = app.fre_shift_mcField.Value;
                    end
                    tan_psi = 2 * Q_mc * df_mc / f_rf;
                else
                    tan_psi = 0;
                end

                % Feedback parameters from UI
                g_dfb_ui = getDFBGainForCurrent(app, I0_val);
                if g_dfb_ui <= 0
                    g_dfb_ui = max(0.05, app.DFB_gainField.Value);
                end
                phi_dfb_rad = deg2rad(app.DFB_phase_shiftField.Value);
                d_dfb       = max(1, round(app.DFB_delayField.Value));
                m_dfb       = max(1, round(app.DFB_mField.Value));
                tau_dfb     = d_dfb * Tb_val;

                kp_pi   = app.KPField.Value;
                ki_pi   = app.KIField.Value;
                m_pi    = max(1, round(app.PI_mField.Value));
                d_pi    = max(1, round(app.PI_dField.Value));
                Ts_pi   = m_pi * Tb_val;
                tau_pi  = d_pi * Tb_val;

                % Frequency grid for Bode (positive modulation frequencies: 0.1 kHz to 5 MHz)
                f_pos = logspace(2, 6.7, 2000); % [Hz]
                w_pos = 2 * pi * f_pos;
                f_pos_kHz = f_pos * 1e-3;

                % Baseband cavity transfer functions matching MC_Cavity_Control.m & DFB_Set.m:
                % - DFB pre-rotates by exp(+j*psi) in phase_shift_app = psi - phi_dfb, canceling
                %   exp(-j*psi) in Vg2Ig -> real numerator sec(psi) = sqrt(1 + tan^2(psi)).
                % - PI divides by PI.RL = 2*R_mc*cos(psi)*exp(j*psi) -> complex numerator (1 - j*tan(psi)).
                sec_psi       = sqrt(1.0 + tan_psi^2);
                H_cav_dfb_pos = sec_psi ./ (1 + 1i * (w_pos / w_half - tan_psi));
                H_cav_dfb_neg = sec_psi ./ (1 + 1i * (-w_pos / w_half - tan_psi));
                H_cav_pi_pos  = (1 - 1i * tan_psi) ./ (1 + 1i * (w_pos / w_half - tan_psi));
                H_cav_pi_neg  = (1 - 1i * tan_psi) ./ (1 + 1i * (-w_pos / w_half - tan_psi));

                % Moving average filter + loop delay for DFB and PI
                sinc_dfb = sin(w_pos * m_dfb * Tb_val / 2) ./ max(1e-30, (w_pos * m_dfb * Tb_val / 2)) ...
                           .* exp(-1i * w_pos * (m_dfb - 1) * Tb_val / 2);
                H_delay_dfb = sinc_dfb .* exp(-1i * (w_pos * tau_dfb + phi_dfb_rad));

                sinc_pi = sin(w_pos * m_pi * Tb_val / 2) ./ max(1e-30, (w_pos * m_pi * Tb_val / 2)) ...
                          .* exp(-1i * w_pos * (m_pi - 1) * Tb_val / 2);
                C_pi_pos = (kp_pi + ki_pi ./ (1i * w_pos * Ts_pi)) .* sinc_pi .* exp(-1i * w_pos * tau_pi);

                % Compute ultimate / max stable DFB gain (Appendix D phase crossover w_pi)
                % Unit-gain DFB open-loop response:
                L_dfb_unit = H_delay_dfb .* H_cav_dfb_pos;
                phase_unit_unwrapped = unwrap(angle(L_dfb_unit));
                idx_pi = find(phase_unit_unwrapped <= -pi, 1, 'first');
                if ~isempty(idx_pi) && idx_pi > 1
                    % Linear interpolation for exact -pi phase crossover
                    p1 = phase_unit_unwrapped(idx_pi - 1);
                    p2 = phase_unit_unwrapped(idx_pi);
                    frac = (-pi - p1) / (p2 - p1);
                    w_pi = w_pos(idx_pi - 1) + frac * (w_pos(idx_pi) - w_pos(idx_pi - 1));
                    f_pi_kHz = (w_pi / (2 * pi)) * 1e-3;
                    mag_unit_at_pi = abs(L_dfb_unit(idx_pi - 1)) + frac * (abs(L_dfb_unit(idx_pi)) - abs(L_dfb_unit(idx_pi - 1)));
                    g_dfb_max = 1.0 / max(1e-12, mag_unit_at_pi);
                else
                    f_pi_kHz = NaN;
                    g_dfb_max = Inf;
                end

                % Select loop mode
                modeStr = app.LLRFModeDropDown.Value;
                cla(app.AxesBodeMag);
                cla(app.AxesBodePhase);
                cla(app.AxesNyquist);
                hold(app.AxesBodeMag, 'on');
                hold(app.AxesBodePhase, 'on');
                hold(app.AxesNyquist, 'on');

                if contains(modeStr, 'Gain Sweep')
                    % Compare multiple DFB gains including UI gain and critical gain G_max
                    g_crit_plot = min(5.0, g_dfb_max);
                    gains_list = unique(round([0.2, 0.5, g_dfb_ui, 1.0, g_crit_plot], 3), 'stable');
                    if length(gains_list) > 5
                        gains_list = [0.2, 0.5, g_dfb_ui, g_crit_plot];
                    end
                    colors_sw = [ ...
                        0.12, 0.47, 0.71; ...
                        0.17, 0.63, 0.17; ...
                        0.90, 0.45, 0.05; ...
                        0.58, 0.20, 0.75; ...
                        0.84, 0.15, 0.16];

                    for ig = 1:length(gains_list)
                        g_cur = gains_list(ig);
                        col_c = colors_sw(mod(ig - 1, size(colors_sw, 1)) + 1, :);
                        lw_c  = 1.6;
                        if abs(g_cur - g_dfb_ui) < 1e-4
                            lw_c = 2.4;
                            lbl_g = sprintf('G_{dfb} = %.2f (Current)', g_cur);
                        elseif abs(g_cur - g_crit_plot) < 1e-3 && isfinite(g_dfb_max)
                            lw_c = 2.0;
                            lbl_g = sprintf('G_{max} = %.2f (Crit.)', g_cur);
                        else
                            lbl_g = sprintf('G_{dfb} = %.2f', g_cur);
                        end

                        L_cur = g_cur * L_dfb_unit;
                        mag_dB = 20 * log10(max(1e-12, abs(L_cur)));
                        pha_deg = rad2deg(unwrap(angle(L_cur)));

                        plot(app.AxesBodeMag, f_pos_kHz, mag_dB, '-', 'Color', col_c, 'LineWidth', lw_c, 'DisplayName', lbl_g);
                        plot(app.AxesBodePhase, f_pos_kHz, pha_deg, '-', 'Color', col_c, 'LineWidth', lw_c, 'DisplayName', lbl_g);

                        % Nyquist contour (+w_m solid, -w_m dashed)
                        plot(app.AxesNyquist, real(L_cur), imag(L_cur), '-', 'Color', col_c, 'LineWidth', lw_c, 'DisplayName', lbl_g);
                        L_neg = g_cur * conj(H_delay_dfb) .* H_cav_dfb_neg;
                        plot(app.AxesNyquist, real(L_neg), imag(L_neg), ':', 'Color', col_c, 'LineWidth', 1.0, 'HandleVisibility', 'off');
                    end

                    % Active loop for margin metrics is g_dfb_ui
                    L_active = g_dfb_ui * L_dfb_unit;
                else
                    if contains(modeStr, 'DFB Only')
                        L_active = g_dfb_ui * L_dfb_unit;
                        L_neg_active = g_dfb_ui * conj(H_delay_dfb) .* H_cav_dfb_neg;
                        lbl_act = sprintf('DFB Only (G_{dfb}=%.2f, d=%d)', g_dfb_ui, d_dfb);
                    elseif contains(modeStr, 'PI Only')
                        L_active = C_pi_pos .* H_cav_pi_pos;
                        C_pi_neg = (kp_pi + ki_pi ./ (-1i * w_pos * Ts_pi)) .* conj(sinc_pi) .* exp(1i * w_pos * tau_pi);
                        L_neg_active = C_pi_neg .* H_cav_pi_neg;
                        lbl_act = sprintf('PI Only (K_p=%.2g, K_i=%.0e, d=%d)', kp_pi, ki_pi, d_pi);
                    else
                        % Combined (DFB + PI)
                        dfb_part = 0;
                        dfb_neg_part = 0;
                        if app.DFBEnableCheckbox.Value || g_dfb_ui > 0
                            dfb_part = g_dfb_ui * L_dfb_unit;
                            dfb_neg_part = g_dfb_ui * conj(H_delay_dfb) .* H_cav_dfb_neg;
                        end
                        pi_part = 0;
                        pi_neg_part = 0;
                        if app.PIEnableCheckbox.Value
                            pi_part = C_pi_pos .* H_cav_pi_pos;
                            C_pi_neg = (kp_pi + ki_pi ./ (-1i * w_pos * Ts_pi)) .* conj(sinc_pi) .* exp(1i * w_pos * tau_pi);
                            pi_neg_part = C_pi_neg .* H_cav_pi_neg;
                        end
                        L_active = dfb_part + pi_part;
                        L_neg_active = dfb_neg_part + pi_neg_part;
                        lbl_act = sprintf('Combined (G_{dfb}=%.2f, K_p=%.2g, K_i=%.0e)', g_dfb_ui, kp_pi, ki_pi);
                    end

                    mag_dB = 20 * log10(max(1e-12, abs(L_active)));
                    pha_deg = rad2deg(unwrap(angle(L_active)));

                    plot(app.AxesBodeMag, f_pos_kHz, mag_dB, '-', 'Color', [0.10 0.42 0.75], 'LineWidth', 2.2, 'DisplayName', lbl_act);
                    plot(app.AxesBodePhase, f_pos_kHz, pha_deg, '-', 'Color', [0.10 0.42 0.75], 'LineWidth', 2.2, 'DisplayName', lbl_act);
                    plot(app.AxesNyquist, real(L_active), imag(L_active), '-', 'Color', [0.10 0.42 0.75], 'LineWidth', 2.2, 'DisplayName', [lbl_act ' (\omega_m > 0)']);
                    plot(app.AxesNyquist, real(L_neg_active), imag(L_neg_active), '--', 'Color', [0.35 0.62 0.88], 'LineWidth', 1.4, 'DisplayName', '\omega_m < 0 branch');
                end

                % Compute Gain Margin (GM) and Phase Margin (PM) of L_active
                pha_act = unwrap(angle(L_active));
                mag_act = abs(L_active);
                idx_pc = find(pha_act <= -pi, 1, 'first');
                if ~isempty(idx_pc) && idx_pc > 1
                    p1 = pha_act(idx_pc - 1); p2 = pha_act(idx_pc);
                    fr = (-pi - p1) / (p2 - p1);
                    f_pc_kHz = f_pos_kHz(idx_pc - 1) + fr * (f_pos_kHz(idx_pc) - f_pos_kHz(idx_pc - 1));
                    mag_pc = mag_act(idx_pc - 1) + fr * (mag_act(idx_pc) - mag_act(idx_pc - 1));
                    gm_dB = -20 * log10(max(1e-12, mag_pc));
                else
                    f_pc_kHz = f_pi_kHz;
                    gm_dB = Inf;
                    mag_pc = 0;
                end

                % Gain crossover (where |L| crosses 1 from above)
                idx_gc = find(mag_act(1:end-1) >= 1.0 & mag_act(2:end) < 1.0, 1, 'first');
                if ~isempty(idx_gc)
                    m1 = mag_act(idx_gc); m2 = mag_act(idx_gc + 1);
                    fr = (1.0 - m1) / (m2 - m1);
                    f_gc_kHz = f_pos_kHz(idx_gc) + fr * (f_pos_kHz(idx_gc + 1) - f_pos_kHz(idx_gc));
                    pha_gc_deg = rad2deg(pha_act(idx_gc) + fr * (pha_act(idx_gc + 1) - pha_act(idx_gc)));
                    pm_deg = 180.0 + pha_gc_deg;
                else
                    f_gc_kHz = NaN;
                    pm_deg = Inf;
                end

                % Reference lines on Bode plots
                yline(app.AxesBodeMag, 0, 'k--', '0 dB', 'LineWidth', 1.2, 'HandleVisibility', 'off');
                xline(app.AxesBodeMag, f_half_kHz, ':', 'Color', [0.4 0.4 0.4], ...
                    'Label', sprintf('f_{1/2}=%.1f kHz', f_half_kHz), 'LineWidth', 1.2, 'HandleVisibility', 'off');
                if isfinite(f_pc_kHz)
                    xline(app.AxesBodeMag, f_pc_kHz, 'r:', sprintf('f_{\\pi}=%.1f kHz', f_pc_kHz), ...
                        'LineWidth', 1.2, 'HandleVisibility', 'off');
                end

                yline(app.AxesBodePhase, -180, 'r--', '-180^\circ (Instability Limit)', 'LineWidth', 1.4, 'HandleVisibility', 'off');
                if isfinite(f_pc_kHz)
                    plot(app.AxesBodePhase, f_pc_kHz, -180, 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6, 'DisplayName', 'Phase Crossover (-180^\circ)');
                end

                app.AxesBodeMag.XScale = 'log';
                app.AxesBodePhase.XScale = 'log';
                xlim(app.AxesBodeMag, [min(f_pos_kHz), max(f_pos_kHz)]);
                xlim(app.AxesBodePhase, [min(f_pos_kHz), max(f_pos_kHz)]);
                ylim(app.AxesBodeMag, [-45, 25]);
                ylim(app.AxesBodePhase, [-360, 20]);
                legend(app.AxesBodeMag, 'Location', 'southwest', 'FontSize', 9);

                % Reference Unit Circle and Critical Point (-1, j0) on Nyquist plot
                th_circ = linspace(0, 2*pi, 300);
                plot(app.AxesNyquist, cos(th_circ), sin(th_circ), '--', 'Color', [0.55 0.55 0.55], 'LineWidth', 1.2, 'DisplayName', 'Unit Circle (|L|=1)');
                xline(app.AxesNyquist, 0, ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
                yline(app.AxesNyquist, 0, ':', 'Color', [0.6 0.6 0.6], 'HandleVisibility', 'off');
                plot(app.AxesNyquist, -1, 0, 'r+', 'MarkerSize', 12, 'LineWidth', 2.5, 'DisplayName', 'Critical Point (-1, j0)');
                plot(app.AxesNyquist, -1, 0, 'ro', 'MarkerSize', 8, 'LineWidth', 1.5, 'HandleVisibility', 'off');

                if isfinite(f_pc_kHz) && mag_pc > 0
                    plot(app.AxesNyquist, [-1, -mag_pc], [0, 0], 'r-', 'LineWidth', 1.8, 'HandleVisibility', 'off');
                    plot(app.AxesNyquist, -mag_pc, 0, 'ks', 'MarkerFaceColor', [0.95 0.75 0.1], 'MarkerSize', 6, ...
                        'DisplayName', sprintf('f_{\\pi}=%.1f kHz (GM=%.1f dB)', f_pc_kHz, gm_dB));
                end

                xlim(app.AxesNyquist, [-2.2, 2.5]);
                ylim(app.AxesNyquist, [-2.2, 2.2]);
                axis(app.AxesNyquist, 'equal');
                xlim(app.AxesNyquist, [-2.2, 2.5]);
                ylim(app.AxesNyquist, [-2.2, 2.2]);
                legend(app.AxesNyquist, 'Location', 'northeast', 'FontSize', 9);

                % Update live stability metrics label
                is_loop_stable = (gm_dB > 0) && (isnan(pm_deg) || pm_deg > 0);
                status_word = 'STABLE';
                if ~is_loop_stable, status_word = 'UNSTABLE (Encircles -1+j0)'; end
                pm_str = 'Inf';
                if isfinite(pm_deg), pm_str = sprintf('%.1f^\\circ (@%.1f kHz)', pm_deg, f_gc_kHz); end
                app.LLRFMetricsLabel.Text = sprintf('GM = %.1f dB (@f_pi=%.1f kHz)  |  PM = %s  |  Max G_dfb = %.2f  |  %s', ...
                    gm_dB, f_pc_kHz, strrep(pm_str, '^\circ', ' deg'), g_dfb_max, status_word);
                if is_loop_stable
                    app.LLRFMetricsLabel.FontColor = [0.10 0.55 0.20];
                else
                    app.LLRFMetricsLabel.FontColor = [0.82 0.12 0.12];
                end

                hold(app.AxesBodeMag, 'off');
                hold(app.AxesBodePhase, 'off');
                hold(app.AxesNyquist, 'off');
                drawnow;
            catch ME_llrf
                fprintf(2, 'Error in updateLLRFNyquistPlot: %s\n', ME_llrf.getReport);
            end
        end

        % -----------------------------------------------------------------
        %  Helper: add a label + numeric field row
        % -----------------------------------------------------------------
        function updateAnalysisPlotView(app)
            app.AxesAnalysis.Visible = 'on';
            plotAnalysisSelection(app);
        end

        function plotAnalysisSelection(app)
            cla(app.AxesAnalysis);
            selection = app.AnalysisPlotModeDropDown.Value;

            % Support Mode-0 Coherent Spectrum vs D-mode Analytical Check
            if strcmp(selection, 'Mode-0 Spectrum (D-mode Check)') || (isempty(app.AnalysisScanPlot) && ~isempty(app.LastMode0Data))
                if isempty(app.LastMode0Data)
                    title(app.AxesAnalysis, 'Run a single-point tracking first to inspect Mode-0 Coherent Spectrum');
                    return;
                end
                m0 = app.LastMode0Data;
                yyaxis(app.AxesAnalysis, 'left');
                cla(app.AxesAnalysis);
                yyaxis(app.AxesAnalysis, 'right');
                cla(app.AxesAnalysis);
                app.AxesAnalysis.YAxis(2).Visible = 'off';
                yyaxis(app.AxesAnalysis, 'left');
                hold(app.AxesAnalysis, 'on');

                plot(app.AxesAnalysis, m0.freq_kHz, m0.amp_tau_norm, '-', 'LineWidth', 1.8, ...
                    'Color', [0.10 0.45 0.78], 'DisplayName', 'STABLE <\tau>_0 FFT');
                plot(app.AxesAnalysis, m0.freq_kHz, m0.amp_del_norm, '--', 'LineWidth', 1.4, ...
                    'Color', [0.85 0.33 0.10], 'DisplayName', 'STABLE <\delta>_0 FFT');

                if isfinite(m0.f_Cmode_kHz) && m0.f_Cmode_kHz <= max(m0.freq_kHz)
                    xline(app.AxesAnalysis, m0.f_Cmode_kHz, ':k', sprintf('C-mode %.2f kHz', m0.f_Cmode_kHz), ...
                        'LineWidth', 1.4, 'HandleVisibility', 'off');
                end
                if isfinite(m0.f_Dmode_plot_kHz) && m0.f_Dmode_plot_kHz > 0
                    xline(app.AxesAnalysis, m0.f_Dmode_plot_kHz, '-.r', m0.dmode_label, ...
                        'LineWidth', 1.8, 'HandleVisibility', 'off');
                end

                xlabel(app.AxesAnalysis, 'Coherent Frequency f  [kHz]');
                ylabel(app.AxesAnalysis, 'Normalized Mode-0 Spectral Amplitude');
                app.AxesAnalysis.YColor = [0.15 0.15 0.15];
                title(app.AxesAnalysis, { ...
                    sprintf('Mode-0 Coherent Spectrum: Sim Peak f_0 = %.2f kHz (\\gamma_0 = %+.1f s^{-1})', m0.f_sim_peak_kHz, m0.gr_sim_s1), ...
                    sprintf('\\fontsize{9}\\rmTheory (He PRAB 2024): f_{Dmode} = %.2f kHz (\\Deltaf_1 = %.1f Hz), \\Omega_{i,D} = %+.1f s^{-1}, \\Omega_{i,C} = %+.1f s^{-1}', ...
                        m0.f_Dmode_kHz, m0.df1_Hz, m0.Omega_i_D, m0.Omega_i_C)});
                xlim(app.AxesAnalysis, [0, max(m0.freq_kHz)]);
                ylim(app.AxesAnalysis, [0, 1.12]);
                app.AxesAnalysis.XGrid = 'on';
                app.AxesAnalysis.YGrid = 'on';
                legend(app.AxesAnalysis, 'Location', 'northeast');
                hold(app.AxesAnalysis, 'off');
                return;
            end

            if isempty(app.AnalysisScanPlot)
                return;
            end
            app.AxesAnalysis.YAxis(2).Visible = 'on';

            x = app.AnalysisScanPlot;
            yBL = app.AnalysisBLPlot;
            yC = app.AnalysisCentroidPlot;
            yE = app.AnalysisEnergySpreadPlot;

            title(app.AxesAnalysis, 'Bunch Lengthening effect vs. HHC Detuning');
            xlabel(app.AxesAnalysis, 'HHC Detuning \Delta f_{hc} [kHz]');
            app.AxesAnalysis.XGrid = 'on';
            app.AxesAnalysis.YGrid = 'on';
            hold(app.AxesAnalysis, 'off');

            switch selection
                case 'Bunch Length'
                    yyaxis(app.AxesAnalysis, 'left');
                    plot(app.AxesAnalysis, x, yBL, '-o', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
                    ylabel(app.AxesAnalysis, 'Bunch Length \sigma_{\tau} [ps]');
                    app.AxesAnalysis.YColor = [0 0.4470 0.7410];
                    legend(app.AxesAnalysis, {'Bunch length'}, 'Location', 'best');
                case 'Centroid'
                    yyaxis(app.AxesAnalysis, 'left');
                    plot(app.AxesAnalysis, x, yC, '--s', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
                    ylabel(app.AxesAnalysis, 'Centroid [ps]');
                    app.AxesAnalysis.YColor = [0.8500 0.3250 0.0980];
                    legend(app.AxesAnalysis, {'Centroid'}, 'Location', 'best');
                case 'Energy Spread'
                    yyaxis(app.AxesAnalysis, 'left');
                    plot(app.AxesAnalysis, x, yE, '-.d', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880]);
                    ylabel(app.AxesAnalysis, 'Energy spread \sigma_{\delta}');
                    app.AxesAnalysis.YColor = [0.4660 0.6740 0.1880];
                    legend(app.AxesAnalysis, {'Energy spread'}, 'Location', 'best');
                case 'Bunch Length + Centroid'
                    yyaxis(app.AxesAnalysis, 'left');
                    plot(app.AxesAnalysis, x, yBL, '-o', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
                    ylabel(app.AxesAnalysis, 'Bunch Length \sigma_{\tau} [ps]');
                    app.AxesAnalysis.YColor = [0 0.4470 0.7410];
                    hold(app.AxesAnalysis, 'on');
                    yyaxis(app.AxesAnalysis, 'right');
                    plot(app.AxesAnalysis, x, yC, '--s', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
                    ylabel(app.AxesAnalysis, 'Centroid [ps]');
                    app.AxesAnalysis.YColor = [0.8500 0.3250 0.0980];
                    legend(app.AxesAnalysis, {'Bunch length', 'Centroid'}, 'Location', 'best');
                case 'Bunch Length + Energy Spread'
                    yyaxis(app.AxesAnalysis, 'left');
                    plot(app.AxesAnalysis, x, yBL, '-o', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
                    ylabel(app.AxesAnalysis, 'Bunch Length \sigma_{\tau} [ps]');
                    app.AxesAnalysis.YColor = [0 0.4470 0.7410];
                    hold(app.AxesAnalysis, 'on');
                    yyaxis(app.AxesAnalysis, 'right');
                    plot(app.AxesAnalysis, x, yE, '-.d', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880]);
                    ylabel(app.AxesAnalysis, 'Energy spread \sigma_{\delta}');
                    app.AxesAnalysis.YColor = [0.4660 0.6740 0.1880];
                    legend(app.AxesAnalysis, {'Bunch length', 'Energy spread'}, 'Location', 'best');
                case 'Centroid + Energy Spread'
                    yyaxis(app.AxesAnalysis, 'left');
                    plot(app.AxesAnalysis, x, yC, '--s', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
                    ylabel(app.AxesAnalysis, 'Centroid [ps]');
                    app.AxesAnalysis.YColor = [0.8500 0.3250 0.0980];
                    hold(app.AxesAnalysis, 'on');
                    yyaxis(app.AxesAnalysis, 'right');
                    plot(app.AxesAnalysis, x, yE, '-.d', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880]);
                    ylabel(app.AxesAnalysis, 'Energy spread \sigma_{\delta}');
                    app.AxesAnalysis.YColor = [0.4660 0.6740 0.1880];
                    legend(app.AxesAnalysis, {'Centroid', 'Energy spread'}, 'Location', 'best');
                otherwise
                    yyaxis(app.AxesAnalysis, 'left');
                    plot(app.AxesAnalysis, x, yBL, '-o', 'LineWidth', 2, 'Color', [0 0.4470 0.7410]);
                    ylabel(app.AxesAnalysis, 'Bunch Length \sigma_{\tau} [ps]');
                    app.AxesAnalysis.YColor = [0 0.4470 0.7410];
                    hold(app.AxesAnalysis, 'on');
                    yyaxis(app.AxesAnalysis, 'right');
                    plot(app.AxesAnalysis, x, yC, '--s', 'LineWidth', 1.5, 'Color', [0.8500 0.3250 0.0980]);
                    plot(app.AxesAnalysis, x, yE, '-.d', 'LineWidth', 1.5, 'Color', [0.4660 0.6740 0.1880]);
                    ylabel(app.AxesAnalysis, 'Centroid [ps] / Energy spread');
                    app.AxesAnalysis.YColor = [0.8500 0.3250 0.0980];
                    yline(app.AxesAnalysis, 1.5, ':r', 'Stability Limit', 'LineWidth', 1.2, 'LabelHorizontalAlignment', 'left');
                    legend(app.AxesAnalysis, {'Bunch length', 'Centroid', 'Energy spread'}, 'Location', 'best');
            end

            if ~isnan(app.AnalysisBestDetune) && any(strcmp(selection, {'Bunch Length', 'Bunch Length + Centroid', 'Bunch Length + Energy Spread', 'All'}))
                yyaxis(app.AxesAnalysis, 'left');
                plot(app.AxesAnalysis, app.AnalysisBestDetune, app.AnalysisBestBunchLength, 'p', 'MarkerSize', 14, ...
                     'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
                text(app.AxesAnalysis, app.AnalysisBestDetune + max(1, (max(x)-min(x))*0.05), app.AnalysisBestBunchLength, ...
                     sprintf('Optimal (%.0f kHz)', app.AnalysisBestDetune), 'Color', 'r', 'FontSize', 11, 'FontWeight', 'bold');
            end
        end

        function [field, rowOut] = addNumericRow(~, grid, labelText, defaultVal, rowIn)
            rowOut = rowIn + 1;
            lbl = uilabel(grid);
            lbl.Layout.Row = rowOut; lbl.Layout.Column = 1;
            lbl.Text = labelText; lbl.FontSize = 12;

            field = uieditfield(grid, 'numeric');
            field.Layout.Row = rowOut; field.Layout.Column = 2;
            field.Value = defaultVal;
            field.ValueDisplayFormat = '%.6g';
            field.FontSize = 12;
        end
    end

    % =====================================================================
    %% Callback Methods
    % =====================================================================
    methods (Access = private)

        % --- Load Config from .mat ---
        function LoadConfigButtonPushed(app)
            [file, path] = uigetfile({'*.mat;*.m','STABLE configuration or machine file (*.mat, *.m)'; ...
                                     '*.mat','MAT configuration (*.mat)'; ...
                                     '*.m','Machine definition (*.m)'}, ...
                                     'Load STABLE Configuration or Machine');
            if isequal(file, 0), return; end
            try
                [~,~,ext] = fileparts(file);
                switch lower(ext)
                    case '.mat'
                        config = load(fullfile(path, file));
                    case '.m'
                        config = loadMachineFile(app, fullfile(path, file));
                    otherwise
                        error('Unsupported file type: %s', ext);
                end
                setUIFromConfig(app, config);
                app.StatusLabel.Text = ['Loaded: ' file];
            catch ME
                uialert(app.UIFigure, ME.message, 'Load Error');
            end
        end

        % --- Save Config to .mat ---
        function SaveConfigButtonPushed(app)
            [file, path] = uiputfile('*.mat', 'Save STABLE Configuration');
            if isequal(file, 0), return; end
            try
                config = getConfigFromUI(app); %#ok<NASGU>
                save(fullfile(path, file), '-struct', 'config');
                app.StatusLabel.Text = ['Saved: ' file];
            catch ME
                uialert(app.UIFigure, ME.message, 'Save Error');
            end
        end

        % --- Preset selection ---
        function PresetDropDownChanged(app)
            val = app.PresetDropDown.Value;
            switch val
                case 'Soleil 2'
                    loadSoleil2Preset(app);
                    app.StatusLabel.Text = 'Preset loaded: Soleil 2';
                case 'HALF'
                    loadHALFPreset(app);
                    app.StatusLabel.Text = 'Preset loaded: HALF';
            end
            app.PresetDropDown.Value = '-- Presets --';
        end

        % --- Left Panel Width selection ---
        function PanelWidthChanged(app)
            valStr = app.PanelWidthDropDown.Value;
            w = sscanf(valStr, 'Panel: %dpx');
            if ~isempty(w) && w >= 250
                app.ContentGrid.ColumnWidth = {w, '1x'};
            end
        end

        % --- Run button ---
        function RunButtonPushed(app)
            if app.IsRunning, return; end
            runSimulation(app);
        end

        % --- Run Scan button ---
        function RunScanButtonPushed(app)
            if app.IsRunning, return; end
            runScan(app);
        end

        % --- Stop button ---
        
        % --- Browse save folder for scan ---
        function BrowseFolderButtonPushed(app)
            folder = uigetdir('', 'Select Parent Folder for Scan Data');
            if folder ~= 0
                app.ScanSaveFolderField.Value = folder;
            end
        end

        % --- Result Tab Selection Changed callback ---
        function ResultTabGroupSelectionChanged(app, event)
            if event.NewValue == app.MapTab
                if isempty(app.AxesMap.Children) || isempty(app.LastScanData)
                    updateStabilityMapView(app);
                end
            end
        end

        % --- Load Scan Data for Stability Map ---
        function loadScanDataForMap(app)
            appDir = fileparts(mfilename('fullpath'));
            defaultDir = fullfile(appDir, 'scan_data');
            if ~exist(defaultDir, 'dir'), defaultDir = appDir; end
            selFolder = uigetdir(defaultDir, 'Select Scan Data Folder');
            if isequal(selFolder, 0), return; end
            sumFile = fullfile(selFolder, 'scan_summary.mat');
            if exist(sumFile, 'file')
                try
                    sData = load(sumFile);
                    if isfield(sData, 'Df_scan_kHz') && isfield(sData, 'I0_scan_mA') && isfield(sData, 'stability_code')
                        restoreScanConfigToUI(app, sData, selFolder);
                        gr_loaded = [];
                        if isfield(sData, 'growth_rate'), gr_loaded = sData.growth_rate; end
                        plotStabilityMap(app, sData.Df_scan_kHz, sData.I0_scan_mA, sData.stability_code, gr_loaded);
                        app.ResultTabGroup.SelectedTab = app.MapTab;
                        [~, fName] = fileparts(selFolder);
                        app.StatusLabel.Text = sprintf('Map & Settings loaded: %s', fName);
                        app.StatusLabel.FontColor = [0.1 0.6 0.2];
                    else
                        uialert(app.UIFigure, 'scan_summary.mat is missing required fields (Df_scan_kHz, I0_scan_mA, stability_code).', 'Load Error');
                    end
                catch ME
                    uialert(app.UIFigure, sprintf('Error loading scan data: %s', ME.message), 'Load Error');
                end
            else
                uialert(app.UIFigure, sprintf('No scan_summary.mat found in:\n%s', selFolder), 'File Not Found');
            end
        end

        % --- Analyze Scan Data button ---
        function AnalyzeButtonPushed(app)
            if app.IsRunning, return; end
            
            folderPath = uigetdir('', 'Select Folder with Scan Data');
            if folderPath == 0
                return;
            end
            
            app.IsRunning = true;
            app.RunBtn.Enable = 'off';
            app.RunScanBtn.Enable = 'off';
            app.AnalyzeBtn.Enable = 'off';
            app.StopBtn.Enable = 'on';
            app.StatusLabel.Text = 'Analyzing Data...';
            app.StatusLabel.FontColor = [0.0 0.4 0.8];
            drawnow;

            try
                % If scan_summary.mat exists in folder, restore settings to UI and update Stability Map tab
                sumFile = fullfile(folderPath, 'scan_summary.mat');
                if exist(sumFile, 'file')
                    try
                        sData = load(sumFile);
                        restoreScanConfigToUI(app, sData, folderPath);
                        if isfield(sData, 'Df_scan_kHz') && isfield(sData, 'I0_scan_mA') && isfield(sData, 'stability_code')
                            gr_loaded = [];
                            if isfield(sData, 'growth_rate'), gr_loaded = sData.growth_rate; end
                            plotStabilityMap(app, sData.Df_scan_kHz, sData.I0_scan_mA, sData.stability_code, gr_loaded);
                            app.ResultTabGroup.SelectedTab = app.MapTab;
                        end
                    catch ME
                        fprintf(2, 'Error in AnalyzeButtonPushed loading scan_summary: %s\n', ME.getReport);
                    end
                else
                    restoreScanConfigToUI(app, struct(), folderPath);
                    app.ResultTabGroup.SelectedTab = app.AnalysisTab;
                end

                file_pattern = fullfile(folderPath, 'SOLEIL_I0*RLfp*QLfp*detune*.mat');
                dir_info = dir(file_pattern);

                if isempty(dir_info)
                    error('No scan data files matching SOLEIL_I0...detune...mat found in the selected directory.');
                end

                num_files = length(dir_info);
                fre_shift_scan = zeros(1, num_files);
                valid_file_flags = true(1, num_files);

                for f = 1:num_files
                    fname = dir_info(f).name;
                    tokens = regexp(fname, 'detune(-?\d+\.?\d*[eE]?[-+]?\d*)', 'tokens');
                    if ~isempty(tokens)
                        fre_shift_scan(f) = str2double(tokens{1}{1});
                    else
                        valid_file_flags(f) = false;
                    end
                end

                dir_info = dir_info(valid_file_flags);
                fre_shift_scan = fre_shift_scan(valid_file_flags);

                [fre_shift_scan, sort_idx] = sort(fre_shift_scan);
                dir_info = dir_info(sort_idx);

                num_scans = length(fre_shift_scan);
                mean_bunch_length_ps = zeros(1, num_scans);
                centroid_mean_ps = zeros(1, num_scans);
                centroid_jitter_ps = zeros(1, num_scans);
                energy_spread_ps = zeros(1, num_scans);

                for i = 1:num_scans
                    df_khz = fre_shift_scan(i) / 1e3;
                    filename = fullfile(folderPath, dir_info(i).name);
                    
                    app.StatusLabel.Text = sprintf('Analyzing %d/%d: %.0fkHz...', i, num_scans, df_khz);
                    drawnow;
                    if app.StopFlag, break; end

                    data = load(filename);
                    
                    if isfield(data, 'HALF') && isfield(data.HALF, 'sigma_t0')
                        sigma_t0 = data.HALF.sigma_t0;
                    else
                        sigma_t0 = 40e-12;
                    end

                    if isfield(data, 'HALF') && isfield(data.HALF, 'sigma_e0')
                        sigma_e0 = data.HALF.sigma_e0;
                    else
                        sigma_e0 = 7.44e-4;
                    end
                    
                    if isfield(data, 'record_Q_std') && isfield(data, 'record_Q_mean')
                        num_records = size(data.record_Q_std, 1);
                        stable_range = round(num_records*0.9):num_records;
                        if isempty(stable_range), stable_range = 1:num_records; end
                        
                        mean_bunch_length_ps(i) = mean(mean(data.record_Q_std(stable_range, :))) * sigma_t0 * 1e12;

                        centroid_bunch_mean = mean(data.record_Q_mean(stable_range, :), 1);
                        centroid_mean_ps(i) = mean(centroid_bunch_mean) * sigma_t0 * 1e12;
                        % Temporal RMS oscillation jitter averaged across bunches (captures both mu=0 and mu!=0 modes)
                        centroid_temporal_std = std(data.record_Q_mean(stable_range, :), 0, 1);
                        centroid_jitter_ps(i) = mean(centroid_temporal_std) * sigma_t0 * 1e12;
                    else
                        mean_bunch_length_ps(i) = NaN;
                        centroid_mean_ps(i) = NaN;
                        centroid_jitter_ps(i) = NaN;
                    end

                    if isfield(data, 'record_P_std')
                        if exist('stable_range', 'var') && ~isempty(stable_range)
                            energy_spread_ps(i) = mean(mean(data.record_P_std(stable_range, :))) * sigma_e0;
                        else
                            energy_spread_ps(i) = NaN;
                        end
                    else
                        energy_spread_ps(i) = NaN;
                    end
                end

                if app.StopFlag
                    app.StatusLabel.Text = 'Analysis stopped.';
                else
                    valid_idx = ~isnan(mean_bunch_length_ps) & ~isnan(centroid_mean_ps) & ~isnan(centroid_jitter_ps) & (fre_shift_scan ~= 0 | mean_bunch_length_ps ~= 0);
                    scan_plot = fre_shift_scan(valid_idx) / 1e3;
                    bl_plot = mean_bunch_length_ps(valid_idx);
                    centroid_mean_plot = centroid_mean_ps(valid_idx);
                    jitter_plot = centroid_jitter_ps(valid_idx);
                    energy_spread_plot = energy_spread_ps(valid_idx);

                    stability_threshold = 1.5;
                    stable_indices = find(jitter_plot < stability_threshold);

                    if ~isempty(stable_indices)
                        [max_bl, best_sub_idx] = max(bl_plot(stable_indices));
                        best_idx = stable_indices(best_sub_idx);
                        best_detune = scan_plot(best_idx);
                        best_bunch_length = bl_plot(best_idx);
                    else
                        [~, best_idx] = min(jitter_plot);
                        if isempty(best_idx)
                            best_detune = NaN; best_bunch_length = NaN;
                        else
                            best_detune = scan_plot(best_idx);
                            best_bunch_length = bl_plot(best_idx);
                        end
                    end

                    cla(app.AxesAnalysis);
                    app.AnalysisScanPlot = scan_plot;
                    app.AnalysisBLPlot = bl_plot;
                    app.AnalysisCentroidPlot = centroid_mean_plot;
                    app.AnalysisEnergySpreadPlot = energy_spread_plot;
                    app.AnalysisBestDetune = best_detune;
                    app.AnalysisBestBunchLength = best_bunch_length;
                    app.AnalysisBestIdx = best_idx;
                    plotAnalysisSelection(app);

                    if ~isnan(best_detune)
                        app.StatusLabel.Text = sprintf('Optimal: %.0fkHz (BL: %.2fps)', best_detune, best_bunch_length);
                    else
                        app.StatusLabel.Text = 'Analysis Complete (No optimal point)';
                    end
                    app.StatusLabel.FontColor = [0.1 0.6 0.2];
                end

            catch ME
                uialert(app.UIFigure, ME.message, 'Analysis Error');
                app.StatusLabel.Text = 'Analysis Error!';
                app.StatusLabel.FontColor = [0.8 0.1 0.1];
            end

            app.IsRunning = false;
            app.StopFlag = false;
            app.RunBtn.Enable = 'on';
            app.RunScanBtn.Enable = 'on';
            app.AnalyzeBtn.Enable = 'on';
            app.StopBtn.Enable = 'off';
        end

        function StopButtonPushed(app)
            app.StopFlag = true;
            app.StatusLabel.Text = 'Stopping...';
        end

        function AutoDetuneToggled(app)
            if app.AutoDetuneCheckbox.Value
                app.fre_shift_hcField.Enable = 'off';
            else
                app.fre_shift_hcField.Enable = 'on';
            end
        end

        % --- Physics Stability-Constrained HC Detuning Optimizer ---
        function OptimizeHCDetuningButtonPushed(app)
            if app.IsRunning, return; end
            app.StatusLabel.Text = 'Optimizing HC Detuning (Physical Stability Constraint)...';
            app.StatusLabel.FontColor = [0.0 0.4 0.8];
            drawnow;

            try
                I0_A  = app.I0Field.Value;
                I0_mA = I0_A * 1e3;
                V_mc  = app.V_mcField.Value;
                U0    = app.U0Field.Value;
                n_hc  = round(app.n_hcField.Value);
                R_hc  = app.R_hcField.Value * 1e6;
                Q_hc  = app.Q_hcField.Value;
                C_val = app.CField.Value;
                h_val = round(app.hField.Value);

                opt_done = false;
                opt_df_khz = 0;
                opt_xi = 0;
                opt_stretch = 0;
                opt_status = 'STABLE';

                % 1. Try Python stable_core runner first
                appDir = fileparts(mfilename('fullpath'));
                scriptPath = fullfile(appDir, 'run_optimizer.py');
                py_cmd = sprintf('python "%s" %g %g %g %d %g %g', ...
                    scriptPath, I0_mA, V_mc, U0, n_hc, R_hc, Q_hc);
                try
                    [stat, cmdout] = system(py_cmd);
                    if stat == 0 && ~isempty(cmdout) && contains(cmdout, ',')
                        tokens = strsplit(strtrim(cmdout), ',');
                        if length(tokens) >= 4
                            opt_df_khz  = str2double(tokens{1});
                            opt_xi      = str2double(tokens{2});
                            opt_stretch = str2double(tokens{3});
                            opt_status  = strtrim(tokens{4});
                            opt_done    = true;
                        end
                    end
                catch
                end

                % 2. Native high-precision analytical optimization fallback (targeting optimal safe flat-potential island xi = 0.98)
                if ~opt_done
                    target_xi = 0.98;
                    use_ff    = app.FormFactorIterCheckbox.Value;
                    fp_opt    = calc_bunch_form_factor_eq('detune_for_xi', I0_A, V_mc, U0, ...
                                    app.E0Field.Value, C_val, h_val, app.alpha_cField.Value, app.sigma_e0Field.Value, ...
                                    n_hc, R_hc, Q_hc, target_xi, use_ff);
                    if fp_opt.fre_shift_hc > 0
                        opt_df_hz = fp_opt.fre_shift_hc;
                    else
                        opt_df_hz = 35e3;
                    end
                    opt_df_khz  = opt_df_hz * 1e-3;
                    opt_xi      = fp_opt.xi_eff;
                    if use_ff && fp_opt.sigma_t0 > 0
                        opt_stretch = fp_opt.sigma_t_eq / fp_opt.sigma_t0;
                    else
                        opt_stretch = 1.0 + 3.2 * (opt_xi^2);
                    end
                    opt_status  = 'STABLE';
                else
                    opt_df_hz = opt_df_khz * 1e3;
                end

                % Update UI fields
                app.AutoDetuneCheckbox.Value = false;
                app.fre_shift_hcField.Enable = 'on';
                app.fre_shift_hcField.Value  = opt_df_hz;

                app.StatusLabel.Text = sprintf('Optimized: \\Delta f_{hc} = %.2f kHz (\\xi=%.2f, %s, Bunch Stretch: %.2fx)', ...
                    opt_df_khz, opt_xi, opt_status, opt_stretch);
                app.StatusLabel.FontColor = [0.1 0.6 0.2];

                uialert(app.UIFigure, ...
                    sprintf('Physics Stability-Constrained Optimization Succeeded!\n\nTarget Current: %.1f mA\nOptimal Detuning: %.2f kHz (%.0f Hz)\nFlat Potential Parameter \\xi: %.3f\nPredicted Bunch Lengthening: %.2fx\nStability Status: %s', ...
                    I0_mA, opt_df_khz, opt_df_hz, opt_xi, opt_stretch, opt_status), ...
                    'HC Detuning Optimized');
            catch ME
                app.StatusLabel.Text = sprintf('Optimization failed: %s', ME.message);
                app.StatusLabel.FontColor = [0.8 0.1 0.1];
            end
        end

        % --- Safety Boundary Extractor and Fitter ---
        function fitAndPlotSafetyBoundary(app)
            if isempty(app.LastScanData) || ~isfield(app.LastScanData, 'stability_code')
                uialert(app.UIFigure, 'No scan data loaded to extract safety boundary.', 'No Scan Data');
                return;
            end
            try
                Df_scan_kHz = double(app.LastScanData.Df_scan_kHz);
                I0_scan_mA  = double(app.LastScanData.I0_scan_mA);
                stability_code = double(app.LastScanData.stability_code);

                nI0 = length(I0_scan_mA);

                % Extract critical boundary transition points
                crit_I0 = [];
                crit_Df = [];

                for iI = 1:nI0
                    row_codes = stability_code(iI, :);
                    valid_cols = find(~isnan(row_codes));
                    if length(valid_cols) < 2, continue; end
                    is_s = (row_codes(valid_cols) == 0);
                    % Use rightmost unstable -> stable transition (outer LCBI safety boundary psi_crit)
                    idx_rel = find(diff(double(is_s)) == 1, 1, 'last');
                    if ~isempty(idx_rel)
                        idx1 = valid_cols(idx_rel);
                        idx2 = valid_cols(idx_rel + 1);
                        df_mid = 0.5 * (Df_scan_kHz(idx1) + Df_scan_kHz(idx2));
                        crit_I0(end+1) = I0_scan_mA(iI); %#ok<AGROW>
                        crit_Df(end+1) = df_mid;         %#ok<AGROW>
                    end
                end

                if length(crit_I0) < 2
                    app.StatusLabel.Text = 'Notice: No clear boundary transition in current scan data.';
                    app.StatusLabel.FontColor = [0.8 0.4 0.1];
                    return;
                end

                % Linear regression Df = a * I0 + b
                p = polyfit(crit_I0, crit_Df, 1);
                df_pred = polyval(p, crit_I0);
                ss_res  = sum((crit_Df - df_pred).^2);
                ss_tot  = sum((crit_Df - mean(crit_Df)).^2);
                r2      = 1.0 - ss_res / max(ss_tot, 1e-12);
                rmse    = sqrt(mean((crit_Df - df_pred).^2));

                % Plot on AxesMap
                hold(app.AxesMap, 'on');
                I0_dense = linspace(min(I0_scan_mA), max(I0_scan_mA), 100);
                Df_dense = polyval(p, I0_dense);

                use_psi = contains(app.MapXAxisDropDown.Value, char(968));
                sign_str = '+';
                if p(2) < 0, sign_str = '-'; end

                if use_psi
                    f0 = 299792458 / app.CField.Value;
                    fhc = round(app.n_hcField.Value) * round(app.hField.Value) * f0;
                    Q_hc = app.Q_hcField.Value;
                    psi_dense = atand(2 * Q_hc * (Df_dense * 1e3) / fhc);
                    plot(app.AxesMap, psi_dense, I0_dense, 'm-', 'LineWidth', 2.8, ...
                        'DisplayName', sprintf('Safety Boundary (R^2=%.3f)', r2));
                else
                    plot(app.AxesMap, Df_dense, I0_dense, 'r-', 'LineWidth', 2.8, ...
                        'DisplayName', sprintf('Safety: \\Delta f > %.3f I_0 %s %.1f (R^2=%.3f)', p(1), sign_str, abs(p(2)), r2));
                    plot(app.AxesMap, crit_Df, crit_I0, 'k.', 'MarkerSize', 14, ...
                        'DisplayName', 'Critical Trans. Points');
                end

                legend(app.AxesMap, 'Location', 'northeast', 'FontSize', 9, 'Box', 'on', 'Color', [1 1 1 0.85]);

                msg = sprintf('Fitted Safety Boundary: \\Delta f_{hc} > %.4f * I_0 %s %.2f kHz (R^2 = %.4f, RMSE = %.2f kHz)', ...
                    p(1), sign_str, abs(p(2)), r2, rmse);
                app.StatusLabel.Text = msg;
                app.StatusLabel.FontColor = [0.85 0.35 0.10];
                drawnow;

                uialert(app.UIFigure, ...
                    sprintf('Safety Boundary Formula Fitted Successfully!\n\nEquation: \\Delta f_{hc} (kHz) = %.4f * I_0 %s %.2f\nOperating Rule: \\Delta f_{hc} > %.4f * I_0 %s %.2f kHz\nCoefficient of Determination R^2: %.4f\nRoot Mean Square Error: %.2f kHz\nCritical Boundary Points: %d', ...
                    p(1), sign_str, abs(p(2)), p(1), sign_str, abs(p(2)), r2, rmse, length(crit_I0)), ...
                    'Safety Boundary Extracted');
            catch ME
                app.StatusLabel.Text = sprintf('Boundary fit error: %s', ME.message);
                app.StatusLabel.FontColor = [0.8 0.1 0.1];
            end
        end

        function PIToggle(app)
            if app.PIEnableCheckbox.Value
                app.KPField.Enable   = 'on';
                app.KIField.Enable   = 'on';
                app.PI_mField.Enable = 'on';
                app.PI_dField.Enable = 'on';
            else
                app.KPField.Enable   = 'off';
                app.KIField.Enable   = 'off';
                app.PI_mField.Enable = 'off';
                app.PI_dField.Enable = 'off';
            end
        end

        function DFBAdaptiveToggled(app)
            if app.DFB_adaptiveCheckbox.Value
                app.DFB_auto_scaleCheckbox.Value = false;
            end
            DFBToggle(app);
        end

        function DFBToggle(app)
            dfb_on = app.DFBEnableCheckbox.Value;
            auto_scale = app.DFB_auto_scaleCheckbox.Value;
            adaptive   = app.DFB_adaptiveCheckbox.Value;

            if dfb_on
                app.DFB_phase_shiftField.Enable   = 'on';
                app.DFB_delayField.Enable         = 'on';
                app.DFB_mField.Enable             = 'on';
                app.DFB_everyField.Enable         = 'on';
                app.DFB_adaptiveCheckbox.Enable   = 'on';
                app.DFB_auto_scaleCheckbox.Enable = 'on';

                if adaptive
                    app.DFB_gainField.Enable          = 'off';
                    app.DFB_auto_scaleCheckbox.Enable = 'off';
                    app.DFB_I1Field.Enable            = 'off';
                    app.DFB_G1Field.Enable            = 'off';
                    app.DFB_I2Field.Enable            = 'off';
                    app.DFB_G2Field.Enable            = 'off';
                elseif auto_scale
                    app.DFB_gainField.Enable          = 'off';
                    app.DFB_I1Field.Enable            = 'on';
                    app.DFB_G1Field.Enable            = 'on';
                    app.DFB_I2Field.Enable            = 'on';
                    app.DFB_G2Field.Enable            = 'on';
                else
                    app.DFB_gainField.Enable          = 'on';
                    app.DFB_I1Field.Enable            = 'off';
                    app.DFB_G1Field.Enable            = 'off';
                    app.DFB_I2Field.Enable            = 'off';
                    app.DFB_G2Field.Enable            = 'off';
                end
            else
                app.DFB_gainField.Enable          = 'off';
                app.DFB_phase_shiftField.Enable   = 'off';
                app.DFB_delayField.Enable         = 'off';
                app.DFB_mField.Enable             = 'off';
                app.DFB_everyField.Enable         = 'off';
                app.DFB_adaptiveCheckbox.Enable   = 'off';
                app.DFB_auto_scaleCheckbox.Enable = 'off';
                app.DFB_I1Field.Enable            = 'off';
                app.DFB_G1Field.Enable            = 'off';
                app.DFB_I2Field.Enable            = 'off';
                app.DFB_G2Field.Enable            = 'off';
            end
        end

        % --- Add HOM row ---
        function AddHOMButtonPushed(app)
            oldData = app.HOMsTable.Data;
            if isempty(oldData)
                app.HOMsTable.Data = {100, 100*50, 500e6};
            else
                app.HOMsTable.Data = [oldData; {100, 100*50, 500e6}];
            end
        end

        % --- Remove last HOM row ---
        function RemoveHOMButtonPushed(app)
            data = app.HOMsTable.Data;
            if ~isempty(data)
                app.HOMsTable.Data = data(1:end-1, :);
            end
        end

        % --- Auto-align HOM to Worst-case Synchrotron Sideband (matching mbtrack2) ---
        function fre_aligned = alignHOMToWorstCaseSideband(app, fre_nominal)
            % Aligns nominal HOM frequencies to upper synchrotron sideband:
            % fr = (M*h + mu + nus) * f0
            % To match mbtrack2's worst-case threshold under multiple HOMs,
            % all HOMs are aligned to the SAME coupled-bunch mode mu so their
            % growth rates coherently superpose (tau_inv = sum tau_i_inv).
            if isempty(fre_nominal)
                fre_aligned = [];
                return;
            end
            cspeed_val = 299792458;
            C_val      = app.CField.Value;
            h_val      = round(app.hField.Value);
            E0_val     = app.E0Field.Value;
            alpha_val  = app.alpha_cField.Value;
            V_mc_val   = app.V_mcField.Value;
            U0_val     = app.U0Field.Value;

            f0_val = cspeed_val / C_val;
            sin_phi_s = min(1, max(-1, U0_val / V_mc_val));
            cos_phi_s = sqrt(max(0, 1 - sin_phi_s^2));
            nu_s = sqrt((alpha_val * h_val * V_mc_val * cos_phi_s) / (2 * pi * E0_val));
            fs_val = nu_s * f0_val;

            % Select target worst-case coupled bunch mode mu from first/dominant HOM
            k1 = round((fre_nominal(1) - fs_val) / f0_val);
            mu_target = mod(k1, h_val);
            if mu_target == 0, mu_target = 1; end

            fre_aligned = zeros(size(fre_nominal));
            for idx = 1:length(fre_nominal)
                fn = fre_nominal(idx);
                if fn <= 0
                    fre_aligned(idx) = fn;
                    continue;
                end
                M = round((fn - fs_val - mu_target * f0_val) / (h_val * f0_val));
                k = M * h_val + mu_target;
                fre_aligned(idx) = (k + nu_s) * f0_val;
            end
        end

        function AlignHOMButtonPushed(app)
            data = app.HOMsTable.Data;
            if isempty(data), return; end
            try
                if iscell(data)
                    f_nom = cell2mat(data(:, 3));
                    f_aln = alignHOMToWorstCaseSideband(app, f_nom);
                    for k = 1:size(data, 1)
                        data{k, 3} = round(f_aln(k));
                    end
                else
                    f_nom = data(:, 3);
                    f_aln = alignHOMToWorstCaseSideband(app, f_nom);
                    data(:, 3) = round(f_aln);
                end
                app.HOMsTable.Data = data;
                app.StatusLabel.Text = 'HOM frequencies aligned to worst-case sideband (mbtrack2 standard).';
                app.StatusLabel.FontColor = [0.1 0.6 0.2];
            catch ME
                warning('Failed to align HOM frequencies: %s', ME.message);
            end
        end

        function HOMAutoAlignToggled(app)
            if app.HOMAutoAlignCheckbox.Value
                AlignHOMButtonPushed(app);
            end
        end

        function MCInputChanged(app)
            Ncav = max(1, round(app.Ncav_mcField.Value));
            Q0   = max(1, app.Q0_mcField.Value);
            RoQ  = max(0, app.RoQ_mcField.Value);
            beta = max(0, app.beta_mcField.Value);
            QL   = Q0 / (1 + beta);
            RL   = Ncav * QL * RoQ;
            app.Q_mcField.Value = QL;
            app.R_mcField.Value = RL;
        end

        function MCLoadedQChanged(app)
            QL = max(1, app.Q_mcField.Value);
            Q0 = max(1, app.Q0_mcField.Value);
            if QL < Q0
                app.beta_mcField.Value = (Q0 / QL) - 1;
            else
                app.beta_mcField.Value = 0;
                app.Q0_mcField.Value = QL;
            end
            Ncav = max(1, round(app.Ncav_mcField.Value));
            RoQ  = max(0, app.RoQ_mcField.Value);
            app.R_mcField.Value = Ncav * QL * RoQ;
        end

        function MCAutoDetuneToggled(app)
            if app.MCAutoDetuneCheckbox.Value
                app.fre_shift_mcField.Enable = 'off';
            else
                app.fre_shift_mcField.Enable = 'on';
            end
        end

        function syncLongitudinalDampingTime(app, sourceField)
            % Enforce Fokker-Planck fluctuation-dissipation consistency: tau_s == tau_z
            if strcmp(sourceField, 'tau_s')
                app.tau_zField.Value = app.tau_sField.Value;
            else
                app.tau_sField.Value = app.tau_zField.Value;
            end
        end

        function FormFactorIterToggled(app, sourceCtrl)
            % Keep Harmonic Cavity panel checkbox and Stability Map toolbar checkbox synchronized
            if strcmp(sourceCtrl, 'map')
                val = app.MapFormFactorCheckbox.Value;
                app.FormFactorIterCheckbox.Value = val;
            else
                val = app.FormFactorIterCheckbox.Value;
                if ~isempty(app.MapFormFactorCheckbox) && isvalid(app.MapFormFactorCheckbox)
                    app.MapFormFactorCheckbox.Value = val;
                end
            end
            try
                updateStabilityMapView(app);
                if ~isempty(app.AxesNyquist) && isvalid(app.AxesNyquist)
                    updateLLRFNyquistPlot(app);
                end
            catch
            end
        end

        function df_mc = computeOptimalMCDetune(app, I0, V_mc, U0, R_mc, Q_mc, R_hc, Q_hc, fre_shift_hc, n_hc)
            % Compute optimal MC detuning frequency to minimize reflected power.
            % This mimics mbtrack2's MC.set_optimal_detune(Itot).
            %
            % The optimal detuning angle psi satisfies:
            %   tan(psi) = -(Ib * R_L * sin(theta_s)) / (Vc + Ib * R_L * cos(theta_s))
            % where theta_s is the synchronous phase accounting for HC offset,
            % Ib = 2 * F_mc * I_dc is the RF component of DC beam current,
            % R_L is the total loaded shunt impedance.
            %
            % Then: df_mc = f_rf * tan(psi) / (2 * Q_L)
            
            cspeed = 299792458;
            f0  = cspeed / app.CField.Value;
            h   = round(app.hField.Value);
            f_rf = h * f0;
            f_hc = n_hc * f_rf;
            
            F_hc = 1.0;
            F_mc = 1.0;
            if ~isempty(app.FormFactorIterCheckbox) && isvalid(app.FormFactorIterCheckbox) && app.FormFactorIterCheckbox.Value
                fp_mc = calc_bunch_form_factor_eq('eval', I0, V_mc, U0, app.E0Field.Value, ...
                    app.CField.Value, h, app.alpha_cField.Value, app.sigma_e0Field.Value, ...
                    n_hc, R_hc, Q_hc, fre_shift_hc, true);
                F_hc = fp_mc.F_hc;
                F_mc = fp_mc.F_mc;
            end
            
            % HC beam loading contribution (real power offset including F_hc^2)
            psi_hc = atan(2 * Q_hc * fre_shift_hc / f_hc);
            Vb_hc_real = 2 * (F_hc^2) * I0 * R_hc * cos(psi_hc)^2;  % real part of HC beam loading
            
            % Synchronous phase accounting for HC offset
            delta = Vb_hc_real;
            cos_theta = (U0 + delta) / max(F_mc * V_mc, 1.0);
            if abs(cos_theta) > 1
                cos_theta = sign(cos_theta);  % clamp to valid range
            end
            theta_s = acos(cos_theta);
            
            % RF component of beam current (fundamental Fourier component)
            Ib = 2 * F_mc * I0;
            
            % Optimal detuning angle (beam loading compensation)
            tan_psi = -(Ib * R_mc * sin(theta_s)) / (V_mc + Ib * R_mc * cos(theta_s));
            
            % Convert to frequency detuning
            df_mc = f_rf * tan_psi / (2 * Q_mc);
        end

        function HCInputChanged(app)
            Ncav = max(1, round(app.Ncav_hcField.Value));
            Q0   = max(1, app.Q0_hcField.Value);
            RoQ  = max(0, app.RoQ_hcField.Value);
            beta = max(0, app.beta_hcField.Value);
            QL   = Q0 / (1 + beta);
            RL_MOhm = Ncav * QL * RoQ / 1e6;
            app.Q_hcField.Value = QL;
            app.R_hcField.Value = RL_MOhm;
        end

        function HCLoadedQChanged(app)
            QL = max(1, app.Q_hcField.Value);
            Q0 = max(1, app.Q0_hcField.Value);
            if QL < Q0
                app.beta_hcField.Value = (Q0 / QL) - 1;
            else
                app.beta_hcField.Value = 0;
                app.Q0_hcField.Value = QL;
            end
            Ncav = max(1, round(app.Ncav_hcField.Value));
            RoQ  = max(0, app.RoQ_hcField.Value);
            app.R_hcField.Value = Ncav * QL * RoQ / 1e6;
        end
    end

    % =====================================================================
    %% Config Get / Set
    % =====================================================================
    methods (Access = private)

        function config = getConfigFromUI(app)
            config.C        = app.CField.Value;
            config.E0       = app.E0Field.Value;
            config.I0       = app.I0Field.Value;
            config.U0       = app.U0Field.Value;
            config.h        = app.hField.Value;
            config.alpha_c  = app.alpha_cField.Value;
            config.tau_s    = app.tau_sField.Value;
            config.tau_z    = app.tau_zField.Value;
            config.sigma_e0 = app.sigma_e0Field.Value;

            config.V_mc         = app.V_mcField.Value;
            config.Ncav_mc      = app.Ncav_mcField.Value;
            config.Q0_mc        = app.Q0_mcField.Value;
            config.RoQ_mc       = app.RoQ_mcField.Value;
            config.beta_mc      = app.beta_mcField.Value;
            config.Q_mc         = app.Q_mcField.Value;
            config.R_mc         = app.R_mcField.Value;
            config.fre_shift_mc = app.fre_shift_mcField.Value;
            config.mc_auto_detune = app.MCAutoDetuneCheckbox.Value;

            config.n_hc          = app.n_hcField.Value;
            config.Ncav_hc       = app.Ncav_hcField.Value;
            config.Q0_hc         = app.Q0_hcField.Value;
            config.RoQ_hc        = app.RoQ_hcField.Value;
            config.beta_hc       = app.beta_hcField.Value;
            config.Q_hc          = app.Q_hcField.Value;
            config.R_hc          = app.R_hcField.Value * 1e6;
            config.auto_detune   = app.AutoDetuneCheckbox.Value;
            config.form_factor_iter = app.FormFactorIterCheckbox.Value;
            config.fre_shift_hc  = app.fre_shift_hcField.Value;

            config.pattern_type    = app.PatternDropDown.Value;
            config.custom_pattern  = app.CustomPatternField.Value;

            config.Par_num      = app.Par_numField.Value;
            config.Track_num    = app.Track_numField.Value;
            config.Recor_step   = app.Recor_stepField.Value;
            config.ShortRange   = app.ShortRangeCheckbox.Value;

            config.KP   = app.KPField.Value;
            config.KI   = app.KIField.Value;
            config.PI_m = app.PI_mField.Value;
            config.PI_d = app.PI_dField.Value;
            config.PI_enable = app.PIEnableCheckbox.Value;

            config.DFB_enable      = app.DFBEnableCheckbox.Value;
            config.DFB_adaptive    = app.DFB_adaptiveCheckbox.Value;
            config.DFB_gain        = app.DFB_gainField.Value;
            config.DFB_phase_shift = app.DFB_phase_shiftField.Value;
            config.DFB_delay       = app.DFB_delayField.Value;
            config.DFB_m           = app.DFB_mField.Value;
            config.DFB_every       = app.DFB_everyField.Value;
            config.DFB_auto_scale  = app.DFB_auto_scaleCheckbox.Value;
            config.DFB_I1          = app.DFB_I1Field.Value;
            config.DFB_G1          = app.DFB_G1Field.Value;
            config.DFB_I2          = app.DFB_I2Field.Value;
            config.DFB_G2          = app.DFB_G2Field.Value;

            config.HOMs = app.HOMsTable.Data;
            config.HOM_auto_align = app.HOMAutoAlignCheckbox.Value;

            config.scan_Df_min = app.DfMinField.Value;
            config.scan_Df_max = app.DfMaxField.Value;
            config.scan_Df_step = app.DfStepField.Value;
            config.scan_I0_min = app.I0MinField.Value;
            config.scan_I0_max = app.I0MaxField.Value;
            config.scan_I0_step = app.I0StepField.Value;
        end

        function setUIFromConfig(app, config)
            if isempty(config) || ~isstruct(config), return; end
            config = config(1);

            if isfield(config,'C') && ~isempty(config.C),               app.CField.Value        = double(config.C);        end
            if isfield(config,'E0') && ~isempty(config.E0),             app.E0Field.Value       = double(config.E0);       end
            if isfield(config,'I0') && ~isempty(config.I0),             app.I0Field.Value       = double(config.I0);       end
            if isfield(config,'U0') && ~isempty(config.U0),             app.U0Field.Value       = double(config.U0);       end
            if isfield(config,'h') && ~isempty(config.h),               app.hField.Value        = double(config.h);        end
            if isfield(config,'alpha_c') && ~isempty(config.alpha_c),   app.alpha_cField.Value  = double(config.alpha_c);  end
            if isfield(config,'tau_s') && ~isempty(config.tau_s),       app.tau_sField.Value    = double(config.tau_s);    end
            if isfield(config,'tau_z') && ~isempty(config.tau_z),       app.tau_zField.Value    = double(config.tau_z);    end
            % Enforce Fokker-Planck longitudinal damping consistency (worst-case longer damping time)
            if abs(app.tau_sField.Value - app.tau_zField.Value) > 1e-9
                tau_wc = max(app.tau_sField.Value, app.tau_zField.Value);
                app.tau_sField.Value = tau_wc;
                app.tau_zField.Value = tau_wc;
            end
            if isfield(config,'sigma_e0') && ~isempty(config.sigma_e0), app.sigma_e0Field.Value = double(config.sigma_e0); end

            if isfield(config,'V_mc') && ~isempty(config.V_mc),         app.V_mcField.Value     = double(config.V_mc);     end
            if isfield(config,'Ncav_mc') && ~isempty(config.Ncav_mc),   app.Ncav_mcField.Value  = double(config.Ncav_mc);  end
            if isfield(config,'Q0_mc') && ~isempty(config.Q0_mc),       app.Q0_mcField.Value    = double(config.Q0_mc);    end
            if isfield(config,'RoQ_mc') && ~isempty(config.RoQ_mc),     app.RoQ_mcField.Value   = double(config.RoQ_mc);   end
            if isfield(config,'beta_mc') && ~isempty(config.beta_mc),   app.beta_mcField.Value  = double(config.beta_mc);  end
            if isfield(config,'Q0_mc') && isfield(config,'beta_mc') && ~isempty(config.Q0_mc) && ~isempty(config.beta_mc)
                MCInputChanged(app);       % Auto-computes loaded Q_mc and total R_mc from Q0 and beta
            end
            if isfield(config,'Q_mc') && ~isempty(config.Q_mc), app.Q_mcField.Value = double(config.Q_mc); end
            if isfield(config,'R_mc') && ~isempty(config.R_mc), app.R_mcField.Value = double(config.R_mc); end
            if isfield(config,'fre_shift_mc') && ~isempty(config.fre_shift_mc), app.fre_shift_mcField.Value = double(config.fre_shift_mc); end
            if isfield(config,'mc_auto_detune') && ~isempty(config.mc_auto_detune)
                app.MCAutoDetuneCheckbox.Value = logical(config.mc_auto_detune);
                MCAutoDetuneToggled(app);
            end

            if isfield(config,'n_hc') && ~isempty(config.n_hc),         app.n_hcField.Value     = double(config.n_hc);     end
            if isfield(config,'Ncav_hc') && ~isempty(config.Ncav_hc),   app.Ncav_hcField.Value  = double(config.Ncav_hc);  end
            if isfield(config,'Q0_hc') && ~isempty(config.Q0_hc),       app.Q0_hcField.Value    = double(config.Q0_hc);    end
            if isfield(config,'RoQ_hc') && ~isempty(config.RoQ_hc),     app.RoQ_hcField.Value   = double(config.RoQ_hc);   end
            if isfield(config,'beta_hc') && ~isempty(config.beta_hc),   app.beta_hcField.Value  = double(config.beta_hc);  end
            if isfield(config,'Q0_hc') && isfield(config,'beta_hc') && ~isempty(config.Q0_hc) && ~isempty(config.beta_hc)
                HCInputChanged(app);       % Auto-computes loaded Q_hc and total R_hc from Q0 and beta
            end
            if isfield(config,'Q_hc') && ~isempty(config.Q_hc), app.Q_hcField.Value = double(config.Q_hc); end
            if isfield(config,'R_hc') && ~isempty(config.R_hc)
                r_hc_val = double(config.R_hc);
                if r_hc_val > 1e4
                    app.R_hcField.Value = r_hc_val / 1e6;
                else
                    app.R_hcField.Value = r_hc_val;
                end
            elseif isfield(config,'Q_hc') && isfield(config,'RoQ_hc') && ~isempty(config.Q_hc) && ~isempty(config.RoQ_hc)
                if abs(double(config.Q_hc) - 31e3) < 5e3
                    app.R_hcField.Value = 2 * double(config.Q_hc) * double(config.RoQ_hc) / 1e6;
                else
                    app.R_hcField.Value = double(config.Q_hc) * double(config.RoQ_hc) / 1e6;
                end
            end
            if isfield(config,'auto_detune') && ~isempty(config.auto_detune),   app.AutoDetuneCheckbox.Value = logical(config.auto_detune); end
            if isfield(config,'form_factor_iter') && ~isempty(config.form_factor_iter)
                app.FormFactorIterCheckbox.Value = logical(config.form_factor_iter);
                if ~isempty(app.MapFormFactorCheckbox) && isvalid(app.MapFormFactorCheckbox)
                    app.MapFormFactorCheckbox.Value = logical(config.form_factor_iter);
                end
            end
            if isfield(config,'fre_shift_hc') && ~isempty(config.fre_shift_hc), app.fre_shift_hcField.Value  = double(config.fre_shift_hc); end

            if isfield(config,'pattern_type') && ~isempty(config.pattern_type)
                pt_str = strtrim(char(config.pattern_type));
                if any(strcmp(app.PatternDropDown.Items, pt_str))
                    app.PatternDropDown.Value = pt_str;
                end
            end
            if isfield(config,'custom_pattern') && ~isempty(config.custom_pattern)
                app.CustomPatternField.Value = strtrim(char(config.custom_pattern));
            end

            if isfield(config,'Par_num') && ~isempty(config.Par_num),       app.Par_numField.Value       = double(config.Par_num);    end
            if isfield(config,'Track_num') && ~isempty(config.Track_num),   app.Track_numField.Value     = double(config.Track_num);  end
            if isfield(config,'Recor_step') && ~isempty(config.Recor_step), app.Recor_stepField.Value    = double(config.Recor_step); end
            if isfield(config,'ShortRange') && ~isempty(config.ShortRange), app.ShortRangeCheckbox.Value = logical(config.ShortRange); end

            if isfield(config,'KP') && ~isempty(config.KP),     app.KPField.Value   = double(config.KP);   end
            if isfield(config,'KI') && ~isempty(config.KI),     app.KIField.Value   = double(config.KI);   end
            if isfield(config,'PI_m') && ~isempty(config.PI_m), app.PI_mField.Value = double(config.PI_m); end
            if isfield(config,'PI_d') && ~isempty(config.PI_d), app.PI_dField.Value = double(config.PI_d); end
            if isfield(config,'PI_enable') && ~isempty(config.PI_enable), app.PIEnableCheckbox.Value = logical(config.PI_enable); end

            if isfield(config,'DFB_enable') && ~isempty(config.DFB_enable),           app.DFBEnableCheckbox.Value      = logical(config.DFB_enable);      end
            if isfield(config,'DFB_adaptive') && ~isempty(config.DFB_adaptive),       app.DFB_adaptiveCheckbox.Value   = logical(config.DFB_adaptive);    end
            if isfield(config,'DFB_gain') && ~isempty(config.DFB_gain),               app.DFB_gainField.Value          = double(config.DFB_gain);        end
            if isfield(config,'DFB_phase_shift') && ~isempty(config.DFB_phase_shift), app.DFB_phase_shiftField.Value   = double(config.DFB_phase_shift); end
            if isfield(config,'DFB_delay') && ~isempty(config.DFB_delay),             app.DFB_delayField.Value         = double(config.DFB_delay);       end
            if isfield(config,'DFB_m') && ~isempty(config.DFB_m),                     app.DFB_mField.Value             = double(config.DFB_m);           end
            if isfield(config,'DFB_every') && ~isempty(config.DFB_every),             app.DFB_everyField.Value         = double(config.DFB_every);       end
            if isfield(config,'DFB_auto_scale') && ~isempty(config.DFB_auto_scale),   app.DFB_auto_scaleCheckbox.Value = logical(config.DFB_auto_scale);  end
            if isfield(config,'DFB_I1') && ~isempty(config.DFB_I1),                   app.DFB_I1Field.Value            = double(config.DFB_I1);          end
            if isfield(config,'DFB_G1') && ~isempty(config.DFB_G1),                   app.DFB_G1Field.Value            = double(config.DFB_G1);          end
            if isfield(config,'DFB_I2') && ~isempty(config.DFB_I2),                   app.DFB_I2Field.Value            = double(config.DFB_I2);          end
            if isfield(config,'DFB_G2') && ~isempty(config.DFB_G2),                   app.DFB_G2Field.Value            = double(config.DFB_G2);          end

            % Update feedback enable UI states
            PIToggle(app);
            DFBToggle(app);

            if isfield(config,'HOMs')
                if isempty(config.HOMs)
                    app.HOMsTable.Data = cell(0, 3);
                elseif isnumeric(config.HOMs)
                    app.HOMsTable.Data = num2cell(double(config.HOMs));
                elseif iscell(config.HOMs)
                    homs_cell = config.HOMs;
                    for r_h = 1:size(homs_cell, 1)
                        for c_h = 1:size(homs_cell, 2)
                            if isnumeric(homs_cell{r_h, c_h})
                                homs_cell{r_h, c_h} = double(homs_cell{r_h, c_h});
                            end
                        end
                    end
                    app.HOMsTable.Data = homs_cell;
                end
            end
            if isfield(config,'HOM_auto_align') && ~isempty(config.HOM_auto_align)
                app.HOMAutoAlignCheckbox.Value = logical(config.HOM_auto_align);
            end

            if isfield(config,'scan_Df_min') && ~isempty(config.scan_Df_min),   app.DfMinField.Value  = double(config.scan_Df_min);  end
            if isfield(config,'scan_Df_max') && ~isempty(config.scan_Df_max),   app.DfMaxField.Value  = double(config.scan_Df_max);  end
            if isfield(config,'scan_Df_step') && ~isempty(config.scan_Df_step), app.DfStepField.Value = double(config.scan_Df_step); end
            if isfield(config,'scan_I0_min') && ~isempty(config.scan_I0_min),   app.I0MinField.Value  = double(config.scan_I0_min);  end
            if isfield(config,'scan_I0_max') && ~isempty(config.scan_I0_max),   app.I0MaxField.Value  = double(config.scan_I0_max);  end
            if isfield(config,'scan_I0_step') && ~isempty(config.scan_I0_step), app.I0StepField.Value = double(config.scan_I0_step); end

            % Update auto-detune UI state
            if app.AutoDetuneCheckbox.Value
                app.fre_shift_hcField.Enable = 'off';
            else
                app.fre_shift_hcField.Enable = 'on';
            end

            % Update computed labels
            updateComputedLabels(app);
        end

        function restoreScanConfigToUI(app, sData, folderPath)
            if nargin < 3, folderPath = ''; end
            try
                cfg = struct();
                if isstruct(sData) && isfield(sData, 'config') && isstruct(sData.config) && ~isempty(sData.config)
                    cfg = sData.config(1);
                elseif ~isempty(folderPath) && exist(folderPath, 'dir')
                    cfg = reconstructConfigFromScanFolder(app, sData, folderPath);
                end

                % Always synchronize scan grid bounds from Df_scan_kHz / I0_scan_mA if available
                if isstruct(sData)
                    if isfield(sData, 'Df_scan_kHz') && ~isempty(sData.Df_scan_kHz)
                        df_vec = double(sData.Df_scan_kHz(:)');
                        cfg.scan_Df_min = min(df_vec);
                        cfg.scan_Df_max = max(df_vec);
                        if length(df_vec) > 1
                            cfg.scan_Df_step = abs(df_vec(2) - df_vec(1));
                        end
                    end
                    if isfield(sData, 'I0_scan_mA') && ~isempty(sData.I0_scan_mA)
                        i0_vec = double(sData.I0_scan_mA(:)');
                        cfg.scan_I0_min = min(i0_vec);
                        cfg.scan_I0_max = max(i0_vec);
                        if length(i0_vec) > 1
                            cfg.scan_I0_step = abs(i0_vec(2) - i0_vec(1));
                        end
                    end
                end

                if ~isempty(fieldnames(cfg))
                    setUIFromConfig(app, cfg);
                    app.LastScanData.config = cfg;
                end

                if isstruct(sData) && isfield(sData, 'scanFolderName') && ~isempty(sData.scanFolderName)
                    app.LastScanData.scanFolderName = strtrim(char(sData.scanFolderName));
                elseif ~isempty(folderPath)
                    [~, fName] = fileparts(folderPath);
                    app.LastScanData.scanFolderName = fName;
                end
            catch ME_cfg
                fprintf(2, 'Notice: could not fully restore scan config to UI: %s\n', ME_cfg.message);
            end
        end

        function cfg = reconstructConfigFromScanFolder(app, sData, folderPath)
            cfg = getConfigFromUI(app);
            [~, fName] = fileparts(folderPath);

            % Inspect first point .mat file in folder if available
            ptFiles = dir(fullfile(folderPath, 'SOLEIL_I0*.mat'));
            if ~isempty(ptFiles)
                ptFile = fullfile(folderPath, ptFiles(1).name);
                try
                    pt = load(ptFile);
                    if isfield(pt, 'HALF') && isstruct(pt.HALF)
                        H = pt.HALF(1);
                        half_fields = {'C','E0','I0','U0','h','alpha_c','tau_s','tau_z','sigma_e0', ...
                                       'V_mc','Ncav_mc','Q0_mc','beta_mc','Q_mc','R_mc','fre_shift_mc', ...
                                       'n_hc','Ncav_hc','Q0_hc','beta_hc','Q_hc','R_hc','fre_shift_hc'};
                        for kf = 1:length(half_fields)
                            fn = half_fields{kf};
                            if isfield(H, fn) && ~isempty(H.(fn))
                                cfg.(fn) = double(H.(fn));
                            end
                        end
                        if isfield(H, 'ShortRange_on') && ~isempty(H.ShortRange_on)
                            cfg.ShortRange = logical(H.ShortRange_on);
                        end
                        if isfield(H, 'Q_hom_m0') && isfield(H, 'R_hom_m0') && isfield(H, 'fre_hom_m0')
                            if ~isempty(H.Q_hom_m0) && any(double(H.R_hom_m0) > 0)
                                q_h = double(H.Q_hom_m0(:));
                                r_h = double(H.R_hom_m0(:));
                                f_h = double(H.fre_hom_m0(:));
                                n_h = min([length(q_h), length(r_h), length(f_h)]);
                                cfg.HOMs = num2cell([q_h(1:n_h), r_h(1:n_h), f_h(1:n_h)]);
                            else
                                cfg.HOMs = cell(0, 3);
                            end
                        end
                    end
                    if isfield(pt, 'Track_num') && ~isempty(pt.Track_num)
                        cfg.Track_num = double(pt.Track_num);
                    end
                    if isfield(pt, 'Recor_step') && ~isempty(pt.Recor_step)
                        cfg.Recor_step = double(pt.Recor_step);
                    end
                    if isfield(pt, 'PI') && isstruct(pt.PI)
                        P = pt.PI(1);
                        if isfield(P, 'enable') && ~isempty(P.enable), cfg.PI_enable = logical(P.enable); end
                        if isfield(P, 'KP') && ~isempty(P.KP),         cfg.KP        = double(P.KP);      end
                        if isfield(P, 'KI') && ~isempty(P.KI),         cfg.KI        = double(P.KI);      end
                        if isfield(P, 'm') && ~isempty(P.m),           cfg.PI_m      = double(P.m);       end
                        if isfield(P, 'd') && ~isempty(P.d),           cfg.PI_d      = double(P.d);       end
                    end
                    if isfield(pt, 'DFB') && isstruct(pt.DFB)
                        D = pt.DFB(1);
                        if isfield(D, 'enable') && ~isempty(D.enable),       cfg.DFB_enable      = logical(D.enable);           end
                        if isfield(D, 'gain') && ~isempty(D.gain),           cfg.DFB_gain        = double(D.gain);              end
                        if isfield(D, 'phase_shift') && ~isempty(D.phase_shift), cfg.DFB_phase_shift = double(D.phase_shift);   end
                        if isfield(D, 'd') && ~isempty(D.d),                 cfg.DFB_delay       = double(D.d);                 end
                        if isfield(D, 'm') && ~isempty(D.m),                 cfg.DFB_m           = double(D.m);                 end
                        if isfield(D, 'every') && ~isempty(D.every),         cfg.DFB_every       = double(D.every);             end
                    end
                catch
                end
            end

            % Refine feedback & HOM flags from folder naming convention
            if contains(fName, 'DFBadapt')
                cfg.DFB_enable = true;
                cfg.DFB_adaptive = true;
            elseif contains(fName, 'DFBauto')
                cfg.DFB_enable = true;
                cfg.DFB_adaptive = false;
                cfg.DFB_auto_scale = true;
            elseif contains(fName, 'NoDFB')
                cfg.DFB_enable = false;
            elseif contains(fName, 'DFBg')
                cfg.DFB_enable = true;
                cfg.DFB_adaptive = false;
                cfg.DFB_auto_scale = false;
                tok = regexp(fName, 'DFBg([0-9.]+)', 'tokens', 'once');
                if ~isempty(tok), cfg.DFB_gain = str2double(tok{1}); end
            end

            if contains(fName, 'HOMalign')
                cfg.HOM_auto_align = true;
            elseif contains(fName, 'HOMnom') || contains(fName, 'NoHOM')
                cfg.HOM_auto_align = false;
            end
        end

        function config = loadMachineFile(app, filePath)
            txt = fileread(filePath);
            config = struct();

            config.C       = parseMachineFileNumber(app, txt, 'C');
            config.E0      = parseMachineFileNumber(app, txt, 'E0');
            config.I0      = parseMachineFileNumber(app, txt, 'I0');
            config.U0      = parseMachineFileNumber(app, txt, 'U0');
            config.h       = parseMachineFileNumber(app, txt, 'h');
            config.alpha_c = parseMachineFileNumber(app, txt, 'alpha_c');
            config.tau_s   = parseMachineFileNumber(app, txt, 'tau_s');
            config.tau_z   = parseMachineFileNumber(app, txt, 'tau_z');
            config.sigma_e0= parseMachineFileNumber(app, txt, 'sigma_e0');

            config.V_mc         = parseMachineFileNumber(app, txt, 'V_mc');
            config.Ncav_mc      = parseMachineFileNumber(app, txt, 'Ncav_mc');
            config.Q0_mc        = parseMachineFileNumber(app, txt, 'Q0_mc');
            config.beta_mc      = parseMachineFileNumber(app, txt, 'beta_mc');
            config.Q_mc         = parseMachineFileNumber(app, txt, 'Q_mc');
            config.RoQ_mc       = parseMachineFileNumber(app, txt, 'R_over_Q_mc');
            config.R_mc         = parseMachineFileNumber(app, txt, 'R_mc');
            config.fre_shift_mc = parseMachineFileNumber(app, txt, 'fre_shift_mc');

            config.n_hc         = parseMachineFileNumber(app, txt, 'n_hc');
            config.Ncav_hc      = parseMachineFileNumber(app, txt, 'Ncav_hc');
            config.Q0_hc        = parseMachineFileNumber(app, txt, 'Q0_hc');
            config.beta_hc      = parseMachineFileNumber(app, txt, 'beta_hc');
            config.Q_hc         = parseMachineFileNumber(app, txt, 'Q_hc');
            config.RoQ_hc       = parseMachineFileNumber(app, txt, 'R_over_Q_hc');
            config.R_hc         = parseMachineFileNumber(app, txt, 'R_hc');
            config.fre_shift_hc = parseMachineFileNumber(app, txt, 'fre_shift_hc');

            if isempty(config.R_mc) && ~isempty(config.Q_mc) && ~isempty(config.RoQ_mc)
                config.R_mc = config.Q_mc * config.RoQ_mc;
            end
            if isempty(config.R_hc) && ~isempty(config.Q_hc) && ~isempty(config.RoQ_hc)
                config.R_hc = config.Q_hc * config.RoQ_hc;
            end

            if ~isempty(strfind(txt, 'pattern = ones(1, h)')) || ~isempty(strfind(txt, 'pattern=ones(1,h)'))
                config.pattern_type = 'Uniform (100%)';
            end
            if isfield(config, 'fre_shift_hc') && ~isempty(config.fre_shift_hc)
                config.auto_detune = false;
            end
        end

        function val = parseMachineFileNumber(~, txt, varName)
            tokens = regexp(txt, ['\b', varName, '\b\s*=\s*([^;]+);'], 'tokens', 'once');
            if isempty(tokens)
                val = [];
                return;
            end
            str = regexprep(tokens{1}, '%.*', '');
            str = strtrim(str);
            val = str2double(str);
            if isnan(val)
                try
                    val = eval(str);
                catch
                    val = [];
                end
            end
        end

        % --- Soleil 2 preset (Worst-case for RF system stability: IDs open) ---
        function loadSoleil2Preset(app)
            app.CField.Value       = 354.097;
            app.E0Field.Value      = 2.75e9;
            app.I0Field.Value      = 0.5;
            app.U0Field.Value      = 470e3;
            app.hField.Value       = 416;
            app.alpha_cField.Value = 1.07e-4;
            app.tau_sField.Value   = 13.7e-3;
            app.tau_zField.Value   = 13.7e-3;
            app.sigma_e0Field.Value= 0.92e-3;

            app.V_mcField.Value        = 1.7e6;
            app.Ncav_mcField.Value     = 4;
            app.Q0_mcField.Value       = 35.7e3;
            app.RoQ_mcField.Value      = 140;
            app.beta_mcField.Value     = 4.95;
            MCInputChanged(app);       % Auto-computes loaded Q_mc = 6000 and total R_mc = 3.36 MOhm from Q0 and beta
            app.fre_shift_mcField.Value= -35e3;
            app.MCAutoDetuneCheckbox.Value = true;
            app.fre_shift_mcField.Enable = 'off';

            app.n_hcField.Value         = 4;
            app.Ncav_hcField.Value      = 2;
            app.Q0_hcField.Value        = 31e3;
            app.RoQ_hcField.Value       = 30;
            app.beta_hcField.Value      = 0;      % Passive cavity without power source
            HCInputChanged(app);       % Auto-computes loaded Q_hc = 31000 and total R_hc = 1.86 MOhm from Q0 and beta
            app.AutoDetuneCheckbox.Value = true;
            app.FormFactorIterCheckbox.Value = true;
            if ~isempty(app.MapFormFactorCheckbox) && isvalid(app.MapFormFactorCheckbox)
                app.MapFormFactorCheckbox.Value = true;
            end
            app.fre_shift_hcField.Enable = 'off';

            app.PatternDropDown.Value   = 'Uniform (100%)';
            app.CustomPatternField.Enable = 'off';

            app.Par_numField.Value   = 1e4;
            app.Track_numField.Value = 1e5;
            app.Recor_stepField.Value= 10;
            app.ShortRangeCheckbox.Value = true;

            % Cooperative PI + fast Direct RF Feedback (restores flat-potential xi=1.0 stability island)
            app.PIEnableCheckbox.Value = true;
            app.KPField.Value   = 0.8;
            app.KIField.Value   = 1e-4;
            app.PI_mField.Value = 80;
            app.PI_dField.Value = 800;

            app.DFBEnableCheckbox.Value      = true;
            app.DFB_adaptiveCheckbox.Value   = false;
            app.DFB_gainField.Value          = 0.2;
            app.DFB_phase_shiftField.Value   = 0;
            app.DFB_delayField.Value         = 705;
            app.DFB_mField.Value             = 1;
            app.DFB_everyField.Value         = 1;
            app.DFB_auto_scaleCheckbox.Value = true;
            app.DFB_I1Field.Value            = 300;
            app.DFB_G1Field.Value            = 0.25;
            app.DFB_I2Field.Value            = 500;
            app.DFB_G2Field.Value            = 0.20;
            PIToggle(app);
            DFBToggle(app);

            % MC HOM (10 kOhm @ 1.7 GHz resonant revolution sideband: 2008*f0 - fs)
            app.HOMsTable.Data = {3e4, 10e3, 1700050101};
            app.HOMAutoAlignCheckbox.Value = true;

            app.DfMinField.Value  = 50;
            app.DfMaxField.Value  = 180;
            app.DfStepField.Value = 5;
            app.I0MinField.Value  = 300;
            app.I0MaxField.Value  = 500;
            app.I0StepField.Value = 25;
            app.ScanSaveFolderField.Value = fullfile(fileparts(mfilename('fullpath')), 'scan_data');

            updateComputedLabels(app);
        end

        % --- HALF preset (from MPMBTrackBunchLengthingwoHighQ.m) ---
        function loadHALFPreset(app)
            app.CField.Value       = 479.86;
            app.E0Field.Value      = 2.2e9;
            app.I0Field.Value      = 0.35;
            app.U0Field.Value      = 400e3;
            app.hField.Value       = 800;
            app.alpha_cField.Value = 9.4e-5;
            app.tau_sField.Value   = 14e-3;
            app.tau_zField.Value   = 14e-3;
            app.sigma_e0Field.Value= 7.44e-4;

            app.V_mcField.Value        = 1.2e6;
            app.Ncav_mcField.Value     = 1;
            app.Q0_mcField.Value       = 1.1e5;
            app.RoQ_mcField.Value      = 45;
            app.beta_mcField.Value     = 0;
            app.Q_mcField.Value        = 1.1e5;
            app.R_mcField.Value        = 4.95e6;
            app.fre_shift_mcField.Value= -7e3;

            app.n_hcField.Value         = 3;
            app.Ncav_hcField.Value      = 1;
            app.Q0_hcField.Value        = 2e8;
            app.RoQ_hcField.Value       = 39;
            app.beta_hcField.Value      = 0;
            app.Q_hcField.Value         = 2e8;
            app.R_hcField.Value          = 7800;
            app.AutoDetuneCheckbox.Value = true;
            app.fre_shift_hcField.Enable = 'off';

            app.PatternDropDown.Value   = '80% filling';
            app.CustomPatternField.Enable = 'off';

            app.Par_numField.Value   = 1e4;
            app.Track_numField.Value = 1e5;
            app.Recor_stepField.Value= 10;
            app.ShortRangeCheckbox.Value = true;

            app.PIEnableCheckbox.Value = true;
            app.KPField.Value   = 1;
            app.KIField.Value   = 1e-5;
            app.PI_mField.Value = 80;
            app.PI_dField.Value = 800;

            app.DFBEnableCheckbox.Value      = false;
            app.DFB_gainField.Value        = 0.2;
            app.DFB_phase_shiftField.Value = 0;
            app.DFB_delayField.Value       = 800;
            app.DFB_mField.Value           = 1;
            app.DFB_everyField.Value       = 1;
            app.DFB_auto_scaleCheckbox.Value = false;
            app.DFB_I1Field.Value          = 200;
            app.DFB_G1Field.Value          = 0.6;
            app.DFB_I2Field.Value          = 500;
            app.DFB_G2Field.Value          = 0.1;
            PIToggle(app);
            DFBToggle(app);

            app.HOMsTable.Data = {};
            app.ScanSaveFolderField.Value = fullfile(fileparts(mfilename('fullpath')), 'scan_data');

            updateComputedLabels(app);
        end
    end

    % =====================================================================
    %% Utility Methods
    % =====================================================================
    methods (Access = private)

        function updateComputedLabels(app)
            cspeed = 299792458;
            % sigma_t0
            h_val  = app.hField.Value;
            V_mc   = app.V_mcField.Value;
            E0     = app.E0Field.Value;
            U0     = app.U0Field.Value;
            ac     = app.alpha_cField.Value;
            se0    = app.sigma_e0Field.Value;
            C_val  = app.CField.Value;
            Qs_est = sqrt(ac / (2*pi) * h_val * sqrt((V_mc/E0)^2 - (U0/E0)^2));
            T0_est = C_val / cspeed;
            st0    = ac * se0 / (2*pi * Qs_est) * T0_est;
            app.sigma_t0Label.Text = sprintf('%.2f  (auto)', st0 * 1e12);
            % Fill pattern info
            pattern = generatePattern(app);
            Bun_num = sum(pattern);
            fillrate = Bun_num / length(pattern);
            app.FillRateLabel.Text = sprintf('%.1f %%', fillrate * 100);
            app.BunchNumLabel.Text = sprintf('%d', Bun_num);
            try
                if ~isempty(app.AxesNyquist) && isvalid(app.AxesNyquist)
                    updateLLRFNyquistPlot(app);
                end
            catch
            end
        end

        function pattern = generatePattern(app)
            h = round(app.hField.Value);
            patType = app.PatternDropDown.Value;
            switch patType
                case 'Uniform (100%)'
                    pattern = ones(1, h);
                case '80% filling'
                    pattern = ones(1, h);
                    pat = 33:40;
                    nGroups = floor(h / 40);
                    for ii = 2:nGroups
                        pat = [pat, 40*ii-7 : 40*ii]; %#ok<AGROW>
                    end
                    pat(pat > h) = [];
                    pattern(pat) = 0;
                case '90% filling'
                    pattern = ones(1, h);
                    pat = 37:40;
                    nGroups = floor(h / 40);
                    for ii = 2:nGroups
                        pat = [pat, 40*ii-3 : 40*ii]; %#ok<AGROW>
                    end
                    pat(pat > h) = [];
                    pattern(pat) = 0;
                case 'Single bunch'
                    pattern = zeros(1, h);
                    pattern(1) = 1;
                case 'Custom'
                    try
                        pattern = eval(app.CustomPatternField.Value);
                        if length(pattern) ~= h
                            error('Pattern length must equal h = %d', h);
                        end
                    catch ME
                        uialert(app.UIFigure, ['Invalid pattern expression: ' ME.message], 'Pattern Error');
                        pattern = ones(1, h);
                    end
                otherwise
                    pattern = ones(1, h);
            end
        end
    end

    % =====================================================================
    %% Simulation Engine
    % =====================================================================
    methods (Access = private)

        function runSimulation(app)
            app.IsRunning = true;
            app.StopFlag  = false;
            app.RunBtn.Enable  = 'off';
            app.RunScanBtn.Enable = 'off';
            app.StopBtn.Enable = 'on';
            app.StatusLabel.Text = 'Initializing...';
            app.StatusLabel.FontColor = [0.0 0.4 0.8];
            drawnow;

            try
                % ==========================================================
                % 1. Read parameters from UI (strictly synchronized)
                % ==========================================================
                MCInputChanged(app);
                HCInputChanged(app);

                cspeed   = 299792458;
                C        = app.CField.Value;
                E0       = app.E0Field.Value;
                I0       = app.I0Field.Value;
                U0       = app.U0Field.Value;
                h        = round(app.hField.Value);
                alpha_c  = app.alpha_cField.Value;
                tau_s    = app.tau_sField.Value;
                tau_z    = app.tau_zField.Value;
                sigma_e0 = app.sigma_e0Field.Value;

                V_mc          = app.V_mcField.Value;
                Q_mc          = app.Q_mcField.Value;
                RoQ_mc        = app.RoQ_mcField.Value;
                R_mc          = app.R_mcField.Value;
                fre_shift_mc  = app.fre_shift_mcField.Value;

                n_hc   = round(app.n_hcField.Value);
                Q_hc   = app.Q_hcField.Value;
                RoQ_hc = app.RoQ_hcField.Value;
                R_hc   = app.R_hcField.Value * 1e6;

                Par_num    = round(app.Par_numField.Value);
                Track_num  = round(app.Track_numField.Value);
                Recor_step = round(app.Recor_stepField.Value);

                % Compute sigma_t0
                Qs_est   = sqrt(alpha_c / (2*pi) * h * sqrt((V_mc/E0)^2 - (U0/E0)^2));
                T0_est   = C / cspeed;
                sigma_t0 = alpha_c * sigma_e0 / (2*pi * Qs_est) * T0_est;
                app.sigma_t0Label.Text = sprintf('%.2f ps', sigma_t0 * 1e12);

                % Fill pattern
                pattern  = generatePattern(app);
                fillrate = sum(pattern) / h;
                Bun_num  = sum(pattern);

                % HHC detuning
                use_ff_run = app.FormFactorIterCheckbox.Value;
                if app.AutoDetuneCheckbox.Value
                    fre_shift_hc = detune_HC_calc(I0, n_hc, C, h, U0, V_mc, R_hc, Q_hc, E0, alpha_c, sigma_e0, use_ff_run, 1.0);
                    app.fre_shift_hcField.Value = fre_shift_hc;
                else
                    fre_shift_hc = app.fre_shift_hcField.Value;
                end

                % Update computed labels
                app.FillRateLabel.Text = sprintf('%.1f %%', fillrate * 100);
                app.BunchNumLabel.Text = sprintf('%d', Bun_num);

                % ==========================================================
                % 2. Build machine structure
                % ==========================================================
                app.StatusLabel.Text = 'Building machine structure...';
                drawnow;

                % Auto-compute MC detuning if enabled
                if app.MCAutoDetuneCheckbox.Value
                    fre_shift_mc = computeOptimalMCDetune(app, I0, V_mc, U0, R_mc, Q_mc, R_hc, Q_hc, fre_shift_hc, n_hc);
                end

                HALF = machine(C, I0, U0, E0, tau_s, tau_z, ...
                               sigma_t0, sigma_e0, alpha_c, ...
                               h, V_mc, ...
                               n_hc, R_hc, Q_hc, ...
                               fillrate, fre_shift_hc, ...
                               Q_mc, R_mc, fre_shift_mc, use_ff_run);
                HALF.ShortRange_on = double(app.ShortRangeCheckbox.Value);
                HALF.Ncav_mc       = round(app.Ncav_mcField.Value);
                HALF.Q0_mc         = app.Q0_mcField.Value;
                HALF.beta_mc       = app.beta_mcField.Value;
                HALF.Ncav_hc       = round(app.Ncav_hcField.Value);
                HALF.Q0_hc         = app.Q0_hcField.Value;
                HALF.beta_hc       = app.beta_hcField.Value;

                % ==========================================================
                % 3. HOMs setup (inline from HOMs_m0.m)
                % ==========================================================
                homData = app.HOMsTable.Data;
                Q_hom_m0   = [];
                R_hom_m0   = [];
                fre_hom_m0 = [];
                if ~isempty(homData) && size(homData, 2) >= 3
                    if iscell(homData)
                        for jj = 1:size(homData, 1)
                            Q_hom_m0   = [Q_hom_m0;   double(cell2mat(homData(jj,1)))]; %#ok<AGROW>
                            R_hom_m0   = [R_hom_m0;   double(cell2mat(homData(jj,2)))]; %#ok<AGROW>
                            fre_hom_m0 = [fre_hom_m0; double(cell2mat(homData(jj,3)))]; %#ok<AGROW>
                        end
                    else
                        Q_hom_m0   = double(homData(:,1));
                        R_hom_m0   = double(homData(:,2));
                        fre_hom_m0 = double(homData(:,3));
                    end
                end
                if ~isempty(Q_hom_m0)
                    if app.HOMAutoAlignCheckbox.Value
                        fre_hom_m0 = alignHOMToWorstCaseSideband(app, fre_hom_m0);
                    end
                    HALF.Q_hom_m0   = Q_hom_m0;
                    HALF.R_hom_m0   = R_hom_m0;
                    HALF.fre_hom_m0 = fre_hom_m0;
                    HALF.wrf_hom_m0       = 2*pi*HALF.fre_hom_m0;
                    HALF.angle_hom_m0     = HALF.wrf_hom_m0 * HALF.Tb;
                    HALF.VbImagFactor_hom_m0 = 1./sqrt(4*HALF.Q_hom_m0.^2-1);
                    HALF.rot_coef_hom_m0     = sqrt(1-1./(2*HALF.Q_hom_m0).^2);
                    HALF.Q_hom_length = length(HALF.Q_hom_m0);
                end

                % ==========================================================
                % 4. Feedback setup (PI and Direct RF Feedback)
                % ==========================================================
                PI = struct();
                PI.m         = round(app.PI_mField.Value);
                PI.d         = round(app.PI_dField.Value);
                PI.piIndex   = 0;
                PI.Integral  = 0;
                PI.RL        = HALF.R_mc;
                PI.Ts        = PI.m * HALF.Tb;
                PI.KP        = app.KPField.Value;
                PI.KI        = app.KIField.Value;
                PI.Ig0       = HALF.Ig_mc_0;
                PI.DIg       = 0;
                PI.Ig        = PI.Ig0;
                PI.Vg0       = HALF.Vg_mc_0;
                PI.DI        = 0;
                PI.Ig_track  = [];
                PI.Ig_FBid   = [];
                PI.RoverQ    = HALF.R_mc / HALF.Q_mc;
                PI.Vrf_mc_track = zeros(1, PI.m);
                PI.enable    = logical(app.PIEnableCheckbox.Value);
                PI.Enabled   = PI.enable;

                if app.DFBEnableCheckbox.Value
                    DFB_gain        = getDFBGainForCurrent(app, I0);
                    DFB_phase_shift = app.DFB_phase_shiftField.Value * pi / 180;
                    DFB_delay       = round(app.DFB_delayField.Value);
                    DFB_sample_num  = round(app.DFB_mField.Value);
                    DFB_every       = round(app.DFB_everyField.Value);
                    DFB = DFB_Set(HALF, DFB_gain, DFB_phase_shift, DFB_delay, DFB_sample_num, DFB_every);
                    DFB.enable      = true;
                else
                    DFB = [];
                end

                if ~isempty(PI) && ~isempty(DFB) && isfield(DFB, 'enable') && DFB.enable
                    PI.Ig0 = DFB.Ig_main0;
                    PI.Ig  = DFB.Ig_main0;
                    PI.Vg0 = DFB.Vg0_main;
                end

                % ==========================================================
                % 5. Initialise particles (GPU)
                % ==========================================================
                app.StatusLabel.Text = 'Initializing GPU particles...';
                drawnow;

                gd = gpuDevice();

                q1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);
                p1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);
                q  = repmat(q1, 1, Bun_num);
                p  = repmat(p1, 1, Bun_num);
                Q = gpuArray(single(q));
                P = gpuArray(single(p));
                index_add = gpuArray(single((1:Bun_num) - 1));

                % ==========================================================
                % 6. Charge, beam loading, coefficients
                % ==========================================================
                charge = ones(1, h) .* pattern;
                charge = charge / sum(charge) * Bun_num;

                HALF.qc    = charge .* pattern * HALF.qc / Par_num;
                HALF.Vb_hc = HALF.qc * HALF.wr_hc * R_hc / Q_hc;
                HALF.Vb_mc = HALF.qc * HALF.wr_mc * R_mc / Q_mc;

                V_hc_load_0 = HALF.V_hc_load_0;
                V_mc_load_0 = HALF.V_mc_load_0;

                Dq = 0.05;

                rot_decay_coef_hc = 1i*HALF.rot_coef_hc - 1/(2*Q_hc);
                TbAng_coef_hc     = exp(rot_decay_coef_hc * HALF.angle_hc);
                exp_ang_coef_hc   = -rot_decay_coef_hc * HALF.wr_hc * sigma_t0;

                rot_decay_coef_mc = 1i*HALF.rot_coef_mc - 1/(2*Q_mc);
                TbAng_coef_mc     = exp(rot_decay_coef_mc * HALF.angle_mc);
                exp_ang_coef_mc   = -rot_decay_coef_mc * HALF.wr_mc * sigma_t0;

                % HOMs coefficients
                if ~isempty(Q_hom_m0)
                    HALF.Vb_hom = HALF.qc .* (HALF.wrf_hom_m0 .* HALF.R_hom_m0 ./ HALF.Q_hom_m0);
                    rot_decay_coef_hom = 1i * HALF.rot_coef_hom_m0 - 1 ./ (2 * HALF.Q_hom_m0);
                    TbAng_coef_hom     = exp(rot_decay_coef_hom .* HALF.angle_hom_m0);
                    exp_ang_coef_hom   = -rot_decay_coef_hom .* HALF.wrf_hom_m0 * sigma_t0;
                    V_hom_load_0       = zeros(HALF.Q_hom_length, 1) + 0i;
                end

                wake_kick_coef = HALF.qc * HALF.kick_coef;

                % ==========================================================
                % 7. Wake function
                % ==========================================================
                tau_q = (0:Dq:200)' * sigma_t0;

                Wake_inter_hc = -HALF.wr_hc * R_hc/Q_hc * ...
                    exp(-tau_q*HALF.wr_hc/(2*Q_hc)) .* ...
                    (cos(tau_q*HALF.wr_hc*HALF.rot_coef_hc) - ...
                     HALF.VbImagFactor_hc * sin(tau_q*HALF.wr_hc*HALF.rot_coef_hc));
                Wake_inter_hc(1) = Wake_inter_hc(1)/2;

                Wake_inter_mc = -HALF.wr_mc * R_mc/Q_mc * ...
                    exp(-tau_q*HALF.wr_mc/(2*Q_mc)) .* ...
                    (cos(tau_q*HALF.wr_mc*HALF.rot_coef_mc) - ...
                     HALF.VbImagFactor_mc * sin(tau_q*HALF.wr_mc*HALF.rot_coef_mc));
                Wake_inter_mc(1) = Wake_inter_mc(1)/2;

                Wake_inter = Wake_inter_hc + Wake_inter_mc;

                % HOMs wake
                if ~isempty(Q_hom_m0)
                    for jj = 1:length(HALF.Q_hom_m0)
                        Wake_inter_hom = -HALF.wrf_hom_m0(jj)*HALF.R_hom_m0(jj)/HALF.Q_hom_m0(jj) * ...
                            exp(-tau_q*HALF.wrf_hom_m0(jj)/2/HALF.Q_hom_m0(jj)) .* ...
                            (cos(tau_q*HALF.wrf_hom_m0(jj)*HALF.rot_coef_hom_m0(jj)) - ...
                             HALF.VbImagFactor_hom_m0(jj)*sin(tau_q*HALF.wrf_hom_m0(jj)*HALF.rot_coef_hom_m0(jj)));
                        Wake_inter_hom(1) = Wake_inter_hom(1)/2;
                        Wake_inter = Wake_inter + Wake_inter_hom;
                    end
                end

                Wake_inter = gpuArray(single(Wake_inter));

                % ==========================================================
                % 8. Recording arrays
                % ==========================================================
                Recor_num          = floor(Track_num / Recor_step);
                record_Q_mean      = zeros(Recor_num, Bun_num);
                record_Q_std       = zeros(Recor_num, Bun_num);
                record_P_mean      = zeros(Recor_num, Bun_num);
                record_P_std       = zeros(Recor_num, Bun_num);
                Vb_hc_track_record = zeros(1, Track_num);
                record_th          = 0;

                if ~isempty(PI) && isfield(PI, 'm') && PI.m > 0
                    HALF.Ig_track = zeros(1, ceil(h/PI.m * Track_num));
                else
                    HALF.Ig_track = zeros(1, Track_num);
                end
                Ig_track_tnum      = 0;
                HALF.Vg_mc_track_0 = HALF.Vg_mc_init;
                Vg_mc              = abs(HALF.Vg_mc_init);
                HALF.rfcoef1_track = HALF.rfcoef1 / HALF.V_mc * Vg_mc;

                % UI update interval
                ui_update_interval = max(500, round(Track_num / 200));

                % ==========================================================
                % 9. TRACKING LOOP
                % ==========================================================
                app.StatusLabel.Text = 'Tracking...';
                drawnow;
                tic;
                beam_lost = false;

                for i = 1:Track_num
                    % --- Check stop ---
                    if app.StopFlag
                        break;
                    end

                    % --- Drift ---
                    Q = Q + P * HALF.drift_coef;

                    % --- Beam Divergence / Loss Check ---
                    max_tau_abs = gather(max(max(abs(Q)))) * sigma_t0;
                    if isnan(max_tau_abs) || isinf(max_tau_abs) || max_tau_abs > 2800e-12
                        beam_lost = true;
                        break;
                    end

                    Q_min = min(Q);
                    Q_new = round((Q - Q_min) * (1/Dq));

                    % --- HOMs Kick ---
                    if ~isempty(Q_hom_m0)
                        for jj = 1:HALF.Q_hom_length
                            exp_angle = exp(exp_ang_coef_hom(jj) * Q);
                            exp_angle_sum = gather(sum(exp_angle));
                            V_load_cpu = double(exp_angle_sum .* HALF.Vb_hom(jj, pattern==1));
                            [V_load, V_hom_load_0(jj)] = VoltageLoadCalc_matlab( ...
                                V_hom_load_0(jj), V_load_cpu, TbAng_coef_hom(jj), pattern);
                            V_load_cpu = V_load(pattern==1) * HALF.kick_coef;
                            V_load_gpu = gpuArray(single(V_load_cpu));
                            V_hom_load_kick = V_load_gpu ./ exp_angle;
                            P = P - real(V_hom_load_kick) + imag(V_hom_load_kick) * HALF.VbImagFactor_hom_m0(jj);
                        end
                    end

                    % --- HHC kick ---
                    exp_angle     = exp(exp_ang_coef_hc * Q);
                    exp_angle_sum = gather(sum(exp_angle));
                    V_load_cpu    = double(exp_angle_sum .* HALF.Vb_hc(pattern==1));
                    [V_load, V_hc_load_0] = VoltageLoadCalc_matlab(V_hc_load_0, V_load_cpu, TbAng_coef_hc, pattern);
                    Vb_hc_track_record(i) = mean(V_load);

                    V_load_cpu     = V_load(pattern==1) * HALF.kick_coef;
                    V_load_gpu     = gpuArray(single(V_load_cpu));
                    V_hc_load_kick = V_load_gpu ./ exp_angle;

                    % --- MC kick + Feedbacks ---
                    exp_angle     = exp(exp_ang_coef_mc * Q);
                    exp_angle_sum = gather(sum(exp_angle));
                    V_load_cpu    = double(exp_angle_sum .* HALF.Vb_mc(pattern==1));
                    [Vc_mc, Vg_mc_track, HALF.Vg_mc_track_0, V_load, HALF.V_mc_load_0, PI, DFB] = ...
                        PI_Control(PI, HALF.Vrf_ideal, HALF.Vg_mc_track_0, ...
                                   HALF.V_mc_load_0, V_load_cpu, TbAng_coef_mc, pattern, DFB);

                    if ~isempty(PI) && isfield(PI, 'Ig_track')
                        Ig_track_num = length(PI.Ig_track);
                        if Ig_track_num > 10000
                            Ig_track_tnum = Ig_track_tnum + Ig_track_num;
                            HALF.Ig_track(Ig_track_tnum-Ig_track_num+1:Ig_track_tnum) = PI.Ig_track;
                            PI.Ig_track = [];
                        end
                    end

                    if ~isempty(DFB) && isfield(DFB, 'Ig_track') && length(DFB.Ig_track) > 10000
                        DFB.Ig_track = [];
                        DFB.Vg_drf_track = [];
                    end

                    V_mc_kick = gpuArray(single(Vc_mc(pattern==1) * HALF.kick_coef)) ./ exp_angle;

                    % --- Short-range wake ---
                    if HALF.ShortRange_on == 1
                        binnum = gather(max(max(Q_new))) + 1;
                        if isnan(binnum) || isinf(binnum) || binnum > 2000 || binnum < 1
                            beam_lost = true;
                            break;
                        end
                        try
                            bin_num_q = sum(BinNumCalZ(binnum, Q_new));
                        catch
                            beam_lost = true;
                            break;
                        end
                        bin_num_q = reshape(bin_num_q, binnum, Bun_num);
                        kick_conv = conv2(bin_num_q, Wake_inter(1:binnum));
                        kick_conv = kick_conv(1:binnum,:) .* wake_kick_coef(pattern==1) * min(i/5000, 1);
                        Q_new = Q_new + (1 + index_add * binnum);
                        wake_kick = kick_conv(Q_new);
                    else
                        wake_kick = 0;
                    end

                    % --- Radiation damping + quantum excitation ---
                    rad_quan_kick = -HALF.radampcoef * P + ...
                        HALF.quanexcoef * gpuArray.randn(Par_num, Bun_num, 'single');

                    P = P + rad_quan_kick - HALF.ploss ...
                        - real(V_hc_load_kick) + imag(V_hc_load_kick)*HALF.VbImagFactor_hc ...
                        - real(V_mc_kick)      + imag(V_mc_kick)     *HALF.VbImagFactor_mc ...
                        + wake_kick;

                    % --- Record ---
                    if mod(i, Recor_step) == 0
                        record_th = record_th + 1;
                        record_Q_mean(record_th,:) = gather(mean(Q));
                        record_Q_std(record_th,:)  = gather(std(Q));
                        record_P_mean(record_th,:) = gather(mean(P));
                        record_P_std(record_th,:)  = gather(std(P));
                    end

                    % --- Live UI update ---
                    if mod(i, ui_update_interval) == 0
                        pct = i / Track_num * 100;
                        elapsed = toc;
                        eta = elapsed / i * (Track_num - i);
                        app.StatusLabel.Text = sprintf('Turn %d / %d  (%.0f%%)  ETA: %.0fs', ...
                            i, Track_num, pct, eta);

                        % Live HC voltage plot (clamped to 3000 kV to prevent explosion display)
                        cla(app.AxesHCAmp);
                        plot(app.AxesHCAmp, min(abs(V_load) / 1e3, 3000));
                        ylabel(app.AxesHCAmp, '|V|  [kV]');
                        title(app.AxesHCAmp, sprintf('HC Voltage Turn %d', i));
                        grid(app.AxesHCAmp, 'on');

                        cla(app.AxesHCPhase);
                        plot(app.AxesHCPhase, angle(V_load)/pi*180);
                        ylabel(app.AxesHCPhase, 'Phase  [deg]');
                        xlabel(app.AxesHCPhase, 'Bucket ID');
                        grid(app.AxesHCPhase, 'on');

                        % Live MC voltage plot (clamped to 3000 kV to prevent explosion display)
                        cla(app.AxesMCAmp);
                        plot(app.AxesMCAmp, min(abs(Vc_mc) / 1e3, 3000));
                        ylabel(app.AxesMCAmp, '|V|  [kV]');
                        title(app.AxesMCAmp, sprintf('MC Voltage Turn %d', i));
                        grid(app.AxesMCAmp, 'on');

                        cla(app.AxesMCPhase);
                        plot(app.AxesMCPhase, angle(Vc_mc)/pi*180);
                        ylabel(app.AxesMCPhase, 'Phase  [deg]');
                        xlabel(app.AxesMCPhase, 'Bucket ID');
                        grid(app.AxesMCPhase, 'on');

                        drawnow limitrate;
                    end
                end % tracking loop

                wait(gd);
                elapsed_total = toc;

                % Flush remaining Ig_track
                if ~isempty(PI.Ig_track)
                    Ig_track_tnum = Ig_track_tnum + length(PI.Ig_track);
                    HALF.Ig_track(Ig_track_tnum-length(PI.Ig_track)+1:Ig_track_tnum) = PI.Ig_track;
                end
                HALF.Ig_track = HALF.Ig_track(1:Ig_track_tnum);

                if beam_lost
                    app.StatusLabel.Text = sprintf('Beam Lost / Diverged at Turn %d (\\tau > 2800 ps)!', i);
                    app.StatusLabel.FontColor = [0.85 0.2 0.2];

                    cla(app.AxesMCAmp);
                    plot(app.AxesMCAmp, min(abs(Vc_mc) / 1e3, 3000), 'r', 'LineWidth', 1.2);
                    ylabel(app.AxesMCAmp, '|V|  [kV]');
                    title(app.AxesMCAmp, sprintf('MC Voltage Turn %d [BEAM LOST]', i), 'Color', 'r');
                    grid(app.AxesMCAmp, 'on');

                    cla(app.AxesHCAmp);
                    plot(app.AxesHCAmp, min(abs(V_load) / 1e3, 3000), 'r', 'LineWidth', 1.2);
                    ylabel(app.AxesHCAmp, '|V|  [kV]');
                    title(app.AxesHCAmp, sprintf('HC Voltage Turn %d [BEAM LOST]', i), 'Color', 'r');
                    grid(app.AxesHCAmp, 'on');
                else
                    app.StatusLabel.Text = sprintf('Done!  %.1f s  (%d turns)', elapsed_total, i);
                    app.StatusLabel.FontColor = [0.1 0.6 0.2];
                end

                % ==========================================================
                % 10. POST-PROCESSING PLOTS
                % ==========================================================
                if record_th > 0
                    updateResultPlots(app, record_Q_mean, record_Q_std, ...
                        record_P_mean, record_P_std, HALF, PI, ...
                        sigma_t0, sigma_e0, Recor_step, Track_num, Bun_num, ...
                        Q, P, Dq, pattern, Vb_hc_track_record, record_th, beam_lost);
                end

            catch ME
                uialert(app.UIFigure, ...
                    sprintf('%s\n\n%s', ME.message, getReport(ME, 'extended', 'hyperlinks', 'off')), ...
                    'Simulation Error');
                app.StatusLabel.Text = 'Error!';
                app.StatusLabel.FontColor = [0.8 0.1 0.1];
            end

            app.IsRunning      = false;
            app.RunBtn.Enable  = 'on';
            app.RunScanBtn.Enable = 'on';
            app.StopBtn.Enable = 'off';
        end
    end

    % =====================================================================
    %% Result Plotting
    % =====================================================================
    methods (Access = private)

        function updateResultPlots(app, record_Q_mean, record_Q_std, ...
                record_P_mean, record_P_std, HALF, PI, ...
                sigma_t0, sigma_e0, Recor_step, Track_num, Bun_num, ...
                Q, P, Dq, pattern, Vb_hc_track_record, record_th, beam_lost)

            if nargin < 17 || isempty(beam_lost), beam_lost = false; end

            Nturns = (1:record_th) * Recor_step;
            nShow  = min(8, Bun_num);  % show first 8 bunches
            colors = lines(nShow);

            % --- Tab 1: Centroid & Sigma vs Turns ---
            cla(app.AxesTauMean);  hold(app.AxesTauMean, 'on');
            cla(app.AxesTauStd);   hold(app.AxesTauStd,  'on');
            cla(app.AxesDeltaMean);hold(app.AxesDeltaMean,'on');
            cla(app.AxesDeltaStd); hold(app.AxesDeltaStd, 'on');

            for bi = 1:2:nShow
                plot(app.AxesTauMean,  Nturns, record_Q_mean(1:record_th,bi)*sigma_t0*1e12, 'Color', colors(bi,:));
                plot(app.AxesTauStd,   Nturns, record_Q_std(1:record_th,bi)*sigma_t0*1e12,  'Color', colors(bi,:));
                plot(app.AxesDeltaMean,Nturns, record_P_mean(1:record_th,bi)*sigma_e0,      'Color', colors(bi,:));
                plot(app.AxesDeltaStd, Nturns, record_P_std(1:record_th,bi)*sigma_e0,       'Color', colors(bi,:));
            end

            ylabel(app.AxesTauMean, '<\tau>  [ps]');    xlabel(app.AxesTauMean, 'Turns');
            ylabel(app.AxesTauStd, '\sigma_{\tau}  [ps]'); xlabel(app.AxesTauStd, 'Turns');
            ylabel(app.AxesDeltaMean, '<\delta>');       xlabel(app.AxesDeltaMean, 'Turns');
            ylabel(app.AxesDeltaStd, '\sigma_{\delta}'); xlabel(app.AxesDeltaStd, 'Turns');

            grid(app.AxesTauMean,'on');  xlim(app.AxesTauMean, [1 Track_num]);
            grid(app.AxesTauStd,'on');   xlim(app.AxesTauStd,  [1 Track_num]);
            grid(app.AxesDeltaMean,'on');xlim(app.AxesDeltaMean,[1 Track_num]);
            grid(app.AxesDeltaStd,'on'); xlim(app.AxesDeltaStd, [1 Track_num]);

            if beam_lost
                title(app.AxesTauMean, sprintf('<\\tau> [ps]  [BEAM LOST at Turn %d]', Nturns(end)), 'Color', 'r');
            end

            hold(app.AxesTauMean,'off'); hold(app.AxesTauStd,'off');
            hold(app.AxesDeltaMean,'off'); hold(app.AxesDeltaStd,'off');

            % --- Tab 2: Bunch-by-Bunch (final snapshot) ---
            cla(app.AxesBunchSigma);    cla(app.AxesBunchCentroid);
            plot(app.AxesBunchSigma,    1:Bun_num, record_Q_std(record_th,:)*sigma_t0*1e12, '.');
            ylabel(app.AxesBunchSigma, '\sigma_{\tau}  [ps]'); xlabel(app.AxesBunchSigma, 'Bunch number');
            xlim(app.AxesBunchSigma, [1 Bun_num]); grid(app.AxesBunchSigma, 'on');

            plot(app.AxesBunchCentroid, 1:Bun_num, record_Q_mean(record_th,:)*sigma_t0*1e12, '.');
            ylabel(app.AxesBunchCentroid, '<\tau>  [ps]'); xlabel(app.AxesBunchCentroid, 'Bunch number');
            xlim(app.AxesBunchCentroid, [1 Bun_num]); grid(app.AxesBunchCentroid, 'on');

            if beam_lost
                title(app.AxesBunchSigma, sprintf('Mean \\sigma_{\\tau} = %.2f ps  [BEAM LOST]', ...
                    mean(record_Q_std(record_th,:)*sigma_t0*1e12)), 'Color', 'r');
                title(app.AxesBunchCentroid, 'Bunch Centroid [BEAM LOST]', 'Color', 'r');
            else
                title(app.AxesBunchSigma, sprintf('Mean \\sigma_{\\tau} = %.2f ps', ...
                    mean(record_Q_std(record_th,:)*sigma_t0*1e12)));
            end

            % --- Tab 3: Density Profile ---
            cla(app.AxesDensity); hold(app.AxesDensity, 'on');
            if beam_lost
                title(app.AxesDensity, 'Beam Lost / Diverged - No Density Profile', 'Color', 'r');
            else
                Dq_den = 0.5;
                Q_min_den = min(Q);
                Q_new_den = round((Q - Q_min_den) * (1/Dq_den));
                binnum_den = gather(max(max(Q_new_den))) + 1;
                if ~isnan(binnum_den) && ~isinf(binnum_den) && binnum_den > 0 && binnum_den <= 2000
                    try
                        bin_num_q_den = sum(BinNumCalZ(binnum_den, Q_new_den));
                        bin_num_q_den = reshape(bin_num_q_den, binnum_den, Bun_num);
                        tau_min_den = gather(Q_min_den);

                        % Plot first bunch density
                        bin_range = (1:binnum_den)' * Dq_den * sigma_t0 + tau_min_den(1) * sigma_t0;
                        plot(app.AxesDensity, (bin_range - Dq_den*sigma_t0)*1e12, ...
                             bin_num_q_den(:,1) / max(bin_num_q_den(:,1)), '-', 'LineWidth', 2);
                        xlabel(app.AxesDensity, '\tau  [ps]'); ylabel(app.AxesDensity, 'Norm. density');
                        xlim(app.AxesDensity, [-150 150]); grid(app.AxesDensity, 'on');
                    catch
                        title(app.AxesDensity, 'Density Profile Computation Error', 'Color', 'r');
                    end
                else
                    title(app.AxesDensity, 'Beam Spread Exceeded Grid Limit', 'Color', 'r');
                end
            end
            hold(app.AxesDensity, 'off');

            % --- Tab 6: Generator current ---
            Ig = HALF.Ig_track;
            if ~isempty(Ig) && any(Ig ~= 0)
                Ig = Ig(Ig ~= 0);
                cla(app.AxesGenAmp);
                plot(app.AxesGenAmp, abs(Ig));
                ylabel(app.AxesGenAmp, '|I_g|  [A]'); xlabel(app.AxesGenAmp, 'Sample');
                title(app.AxesGenAmp, 'Generator Current'); grid(app.AxesGenAmp, 'on');

                cla(app.AxesGenPhase);
                plot(app.AxesGenPhase, angle(Ig));
                ylabel(app.AxesGenPhase, 'Phase  [rad]'); xlabel(app.AxesGenPhase, 'Sample');
                grid(app.AxesGenPhase, 'on');
            end

            % --- Mode-0 Coherent Spectrum & D-mode Analytical Verification ---
            if record_th >= 16
                try
                    w_m0 = max(1, round(0.20 * record_th)) : record_th;
                    tau_m0 = mean(record_Q_mean(w_m0, :), 2) * sigma_t0 * 1e12;
                    del_m0 = mean(record_P_mean(w_m0, :), 2) * sigma_e0;
                    tau_ac = tau_m0 - mean(tau_m0);
                    del_ac = del_m0 - mean(del_m0);
                    Nw = length(w_m0);
                    hann_w = 0.5 * (1 - cos(2*pi*(0:Nw-1)' / max(Nw-1, 1)));
                    N_fft0 = 2^nextpow2(max(256, 4 * Nw));
                    dt_rec = Recor_step * HALF.T0;
                    f_ax_kHz = ((0:floor(N_fft0/2)-1)' / N_fft0 / dt_rec) * 1e-3;
                    f_nyq_kHz = 0.5 / dt_rec * 1e-3;

                    sp_tau = abs(fft(tau_ac .* hann_w, N_fft0));
                    sp_del = abs(fft(del_ac .* hann_w, N_fft0));
                    sp_tau = sp_tau(1:floor(N_fft0/2));
                    sp_del = sp_del(1:floor(N_fft0/2));
                    sp_tau_n = sp_tau / max(max(sp_tau), 1e-30);
                    sp_del_n = sp_del / max(max(sp_del), 1e-30);

                    [~, pk0] = max(sp_tau(2:end));
                    f_sim_peak_kHz = f_ax_kHz(pk0 + 1);

                    % Fit Mode-0 exponential growth / damping rate from <delta>_0 envelope
                    env_m0 = abs(del_ac) + 1e-30;
                    t_vec0 = ((w_m0') * Recor_step) * HALF.T0;
                    TF0 = islocalmax(log(env_m0));
                    if sum(TF0) >= 4
                        p_gr0 = polyfit(t_vec0(TF0), log(env_m0(TF0)), 1);
                    else
                        p_gr0 = polyfit(t_vec0, log(env_m0), 1);
                    end
                    gr_sim_s1 = p_gr0(1);

                    cur_dfb = getDFBGainForCurrent(app, HALF.I0);
                    [~, ~, ~, dm_info] = check_dmode_robinson( ...
                        HALF.I0, HALF.R_hc, HALF.Q_hc, HALF.n_hc, HALF.h, 1.0/HALF.T0, ...
                        HALF.U0, HALF.V_mc, HALF.tau_z, HALF.fre_shift_hc, HALF.E0, HALF.alpha_c, ...
                        HALF.R_mc, HALF.Q_mc, HALF.fre_shift_mc, cur_dfb, 0);

                    f_D_kHz = dm_info.f_Dmode_Hz * 1e-3;
                    f_fs_samp_kHz = 2.0 * f_nyq_kHz;
                    if f_D_kHz > f_nyq_kHz && f_fs_samp_kHz > 0
                        f_fold = mod(f_D_kHz, f_fs_samp_kHz);
                        if f_fold > f_nyq_kHz, f_fold = f_fs_samp_kHz - f_fold; end
                        f_D_plot_kHz = f_fold;
                        dmode_lbl = sprintf('D-mode alias %.2f kHz (true %.2f kHz)', f_fold, f_D_kHz);
                    else
                        f_D_plot_kHz = f_D_kHz;
                        dmode_lbl = sprintf('D-mode \\Omega_r/2\\pi = %.2f kHz', f_D_kHz);
                    end

                    m0Data = struct();
                    m0Data.freq_kHz       = f_ax_kHz;
                    m0Data.amp_tau_norm   = sp_tau_n;
                    m0Data.amp_del_norm   = sp_del_n;
                    m0Data.f_sim_peak_kHz = f_sim_peak_kHz;
                    m0Data.gr_sim_s1      = gr_sim_s1;
                    m0Data.f_Dmode_kHz    = f_D_kHz;
                    m0Data.f_Dmode_plot_kHz = f_D_plot_kHz;
                    m0Data.dmode_label    = dmode_lbl;
                    m0Data.df1_Hz         = dm_info.df1_Hz;
                    m0Data.Omega_i_D      = dm_info.Omega_i_Dmode;
                    m0Data.f_Cmode_kHz    = dm_info.f_Cmode_Hz * 1e-3;
                    m0Data.Omega_i_C      = dm_info.Omega_i_Cmode;
                    app.LastMode0Data     = m0Data;

                    if ~beam_lost
                        title(app.AxesTauMean, sprintf('<\\tau> [ps]  (Mode-0 Peak: %.2f kHz | Theory D-mode: %.2f kHz, \\Omega_i=%+.1f s^{-1})', ...
                            f_sim_peak_kHz, f_D_kHz, max(dm_info.Omega_i_Dmode, dm_info.Omega_i_Cmode)));
                    end
                    if isempty(app.AnalysisScanPlot) || strcmp(app.AnalysisPlotModeDropDown.Value, 'Mode-0 Spectrum (D-mode Check)')
                        plotAnalysisSelection(app);
                    end
                catch
                end
            end

            drawnow;
        end
    end

    % =====================================================================
    %% App Creation and Deletion
    % =====================================================================
    methods (Access = public)

        function g = getDFBGainForCurrent(app, I0_A)
            if ~app.DFBEnableCheckbox.Value
                g = 0;
                return;
            end
            if app.DFB_adaptiveCheckbox.Value
                % Adaptive Optimal DFB Controller (Physics & Machine Parameters based):
                % Dynamic feedback gain scaling based on fundamental cavity beam loading Y = I0 * R_mc / V_mc.
                % Maximizes Robinson & dipole mode damping while strictly maintaining loop delay phase margin.
                V_mc = app.V_mcField.Value;
                R_mc = app.R_mcField.Value;
                if isempty(V_mc) || V_mc <= 0, V_mc = 1.7e6; end
                if isempty(R_mc) || R_mc <= 0, R_mc = 3.36e6; end
                Y = (I0_A * R_mc) / V_mc;

                % Monotonic optimal gain scaling:
                % Yields ~0.24 at 200mA, ~0.26 at 300mA, ~0.28 at 400mA, ~0.29 at 450mA, ~0.30 at 500mA
                G0 = 0.20;
                kappa = 0.50;
                g_opt = G0 * (1.0 + kappa * Y);
                g = max(0.15, min(0.38, g_opt));
                return;
            end
            if ~app.DFB_auto_scaleCheckbox.Value
                g = app.DFB_gainField.Value;
                return;
            end
            I_mA = I0_A * 1e3;
            I1   = app.DFB_I1Field.Value;
            G1   = app.DFB_G1Field.Value;
            I2   = app.DFB_I2Field.Value;
            G2   = app.DFB_G2Field.Value;
            if abs(I2 - I1) < 1e-3
                g = G1;
            else
                t = (I_mA - I1) / (I2 - I1);
                g = G1 + t * (G2 - G1);
                g_min = min(G1, G2);
                g_max = max(G1, G2);
                g = max(g_min, min(g_max, g));
            end
        end

        function app = STABLE_App
            % Create and configure components
            createComponents(app);

            % Register app with App Designer runtime
            registerApp(app, app.UIFigure);

            % Load default preset
            loadSoleil2Preset(app);

            % Auto-load latest stability map scan summary if available
            try
                appDir = fileparts(mfilename('fullpath'));
                scanDir = fullfile(appDir, 'scan_data');
                default_scan = '';
                if exist(scanDir, 'dir')
                    dList = dir(fullfile(scanDir, 'Scan_*'));
                    if ~isempty(dList)
                        [~, sIdx] = sort([dList.datenum], 'descend');
                        for k = 1:length(sIdx)
                            cand = fullfile(scanDir, dList(sIdx(k)).name, 'scan_summary.mat');
                            if exist(cand, 'file')
                                default_scan = cand;
                                break;
                            end
                        end
                    end
                end
                if exist(default_scan, 'file')
                    sData = load(default_scan);
                    if isfield(sData, 'Df_scan_kHz') && isfield(sData, 'I0_scan_mA') && isfield(sData, 'stability_code')
                        % Plot latest cached map for visual reference without overwriting fresh preset UI fields
                        gr_loaded = [];
                        if isfield(sData, 'growth_rate'), gr_loaded = sData.growth_rate; end
                        plotStabilityMap(app, sData.Df_scan_kHz, sData.I0_scan_mA, sData.stability_code, gr_loaded);
                    end
                end
            catch ME
                fprintf(2, 'Notice: could not auto-load scan map on startup: %s\n', ME.getReport);
            end

            if nargout == 0
                clear app
            end
        end

        function delete(app)
            delete(app.UIFigure);
        end

        function renderSummaryFile(app, summaryFile)
            if exist(summaryFile, 'file')
                sData = load(summaryFile);
                gr_loaded = [];
                if isfield(sData, 'growth_rate'), gr_loaded = sData.growth_rate; end
                plotStabilityMap(app, sData.Df_scan_kHz, sData.I0_scan_mA, sData.stability_code, gr_loaded);
            end
        end
    end

    % =====================================================================
    %% Scan & Stability Map Methods
    % =====================================================================
    methods (Access = private)

        function runScan(app)
            app.IsRunning = true;
            app.StopFlag  = false;
            app.RunBtn.Enable  = 'off';
            app.RunScanBtn.Enable = 'off';
            app.StopBtn.Enable = 'on';
            app.StatusLabel.Text = 'Initializing Scan...';
            app.StatusLabel.FontColor = [0.0 0.4 0.8];
            drawnow;

            try
                % Ensure Q and R are strictly synchronized with Ncav, Q0, beta, RoQ
                MCInputChanged(app);
                HCInputChanged(app);

                % Auto-disable HC auto-detuning for scan (scan sweeps Df explicitly)
                if app.AutoDetuneCheckbox.Value
                    app.AutoDetuneCheckbox.Value = false;
                    AutoDetuneToggled(app);
                    app.StatusLabel.Text = 'Note: HC auto-detuning disabled (scan sweeps f explicitly)';
                    drawnow; pause(0.5);
                end

                cspeed   = 299792458;
                C        = app.CField.Value;
                E0       = app.E0Field.Value;
                U0       = app.U0Field.Value;
                h        = round(app.hField.Value);
                alpha_c  = app.alpha_cField.Value;
                tau_s    = app.tau_sField.Value;
                tau_z    = app.tau_zField.Value;
                sigma_e0 = app.sigma_e0Field.Value;

                V_mc          = app.V_mcField.Value;
                Q_mc          = app.Q_mcField.Value;
                RoQ_mc        = app.RoQ_mcField.Value;
                R_mc          = app.R_mcField.Value;
                fre_shift_mc  = app.fre_shift_mcField.Value;

                n_hc   = round(app.n_hcField.Value);
                Q_hc   = app.Q_hcField.Value;
                RoQ_hc = app.RoQ_hcField.Value;
                R_hc   = app.R_hcField.Value * 1e6;

                Par_num    = round(app.Par_numField.Value);
                Track_num  = round(app.Track_numField.Value);
                Recor_step = round(app.Recor_stepField.Value);
                Recor_num  = floor(Track_num / Recor_step);

                % HOMs setup
                homData = app.HOMsTable.Data;
                Q_hom_m0 = []; R_hom_m0 = []; fre_hom_m0 = [];
                if ~isempty(homData) && size(homData, 2) >= 3
                    if iscell(homData)
                        for jj = 1:size(homData, 1)
                            Q_hom_m0   = [Q_hom_m0;   double(cell2mat(homData(jj,1)))]; %#ok<AGROW>
                            R_hom_m0   = [R_hom_m0;   double(cell2mat(homData(jj,2)))]; %#ok<AGROW>
                            fre_hom_m0 = [fre_hom_m0; double(cell2mat(homData(jj,3)))]; %#ok<AGROW>
                        end
                    else
                        Q_hom_m0   = double(homData(:,1));
                        R_hom_m0   = double(homData(:,2));
                        fre_hom_m0 = double(homData(:,3));
                    end
                end
                if ~isempty(Q_hom_m0) && app.HOMAutoAlignCheckbox.Value
                    fre_hom_m0 = alignHOMToWorstCaseSideband(app, fre_hom_m0);
                end

                % Grid
                Df_min = app.DfMinField.Value; Df_max = app.DfMaxField.Value; Df_step = app.DfStepField.Value;
                I0_min = app.I0MinField.Value; I0_max = app.I0MaxField.Value; I0_step = app.I0StepField.Value;

                if Df_step == 0, Df_step = 1; end
                if I0_step == 0, I0_step = 1; end

                Df_scan_kHz = Df_min : Df_step : Df_max;
                I0_scan_mA  = I0_min : I0_step : I0_max;
                
                if isempty(Df_scan_kHz), Df_scan_kHz = Df_min; end
                if isempty(I0_scan_mA),  I0_scan_mA  = I0_min; end

                Df_scan = Df_scan_kHz * 1e3;
                I0_scan = I0_scan_mA * 1e-3;

                nDf = length(Df_scan);
                nI0 = length(I0_scan);
                total_points = nDf * nI0;

                stability_code = NaN(nI0, nDf);
                growth_rate    = NaN(nI0, nDf);
                peak_freq_norm = NaN(nI0, nDf);
                mode0_freq_kHz = NaN(nI0, nDf);

                pattern  = generatePattern(app);
                fillrate = sum(pattern) / h;
                Bun_num  = sum(pattern);

                % --- Determine save folder ---
                parentDir = strtrim(app.ScanSaveFolderField.Value);
                if isempty(parentDir)
                    appDir = fileparts(mfilename('fullpath'));
                    parentDir = fullfile(appDir, 'scan_data');
                    if ~exist(parentDir, 'dir'), mkdir(parentDir); end
                    app.ScanSaveFolderField.Value = parentDir;
                end

                % --- Construct informative scan folder name with key parameters ---
                range_tag = sprintf('Df%.0f-%.0fkHz_I%.0f-%.0fmA', Df_min, Df_max, I0_min, I0_max);

                % Feedback parameters tag
                fb_parts = {};
                if app.DFBEnableCheckbox.Value
                    if app.DFB_adaptiveCheckbox.Value
                        fb_parts{end+1} = 'DFBadapt';
                    elseif app.DFB_auto_scaleCheckbox.Value
                        fb_parts{end+1} = 'DFBauto';
                    else
                        fb_parts{end+1} = sprintf('DFBg%.2f', app.DFB_gainField.Value);
                    end
                end
                if app.PIEnableCheckbox.Value
                    fb_parts{end+1} = sprintf('PIkp%.1f_ki%.0e', app.KPField.Value, app.KIField.Value);
                end
                if isempty(fb_parts)
                    fb_tag = 'OpenLoop';
                else
                    fb_tag = strjoin(fb_parts, '_');
                end
                fb_tag = strrep(fb_tag, '+', '');

                % HOMs tag
                if ~isempty(app.HOMsTable.Data) && size(app.HOMsTable.Data, 1) > 0
                    if app.HOMAutoAlignCheckbox.Value
                        hom_tag = 'HOMalign';
                    else
                        hom_tag = 'HOMnom';
                    end
                else
                    hom_tag = 'NoHOM';
                end

                scanFolderName = sprintf('Scan_%s_%s_%s_%s', ...
                    range_tag, fb_tag, hom_tag, datestr(now, 'yyyymmdd_HHMMSS'));
                scanDataDir = fullfile(parentDir, scanFolderName);
                if ~exist(scanDataDir, 'dir')
                    mkdir(scanDataDir);
                end

                gd = gpuDevice();
                rng(42);
                q1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);
                p1 = TruncatedGaussian(1, [-3,3], [Par_num, 1]);

                pt = 0;
                tic;

                win_transient = round(0.5 * Recor_num);
                win_analysis  = win_transient+1 : Recor_num;
                GR_threshold_factor = 0.5;
                f0_approx   = cspeed / C;
                blink_interval = max(500, round(Track_num / 10));

                % Pre-render the initial map so axes limits, grid, and curves are ready
                try
                    plotStabilityMap(app, Df_scan_kHz, I0_scan_mA, stability_code, growth_rate);
                    drawnow;
                catch
                end

                for iI = 1:nI0
                    I0 = I0_scan(iI);
                    current_dfb_gain = getDFBGainForCurrent(app, I0);

                    for iDf = 1:nDf
                        fre_shift_hc = Df_scan(iDf);
                        pt = pt + 1;

                        if app.StopFlag, break; end

                        if app.DFBEnableCheckbox.Value
                            if app.DFB_adaptiveCheckbox.Value
                                app.StatusLabel.Text = sprintf('Scan %d/%d: %dmA (DFB_opt=%.2f), %.1fkHz...', ...
                                    pt, total_points, round(I0*1e3), current_dfb_gain, fre_shift_hc*1e-3);
                            elseif app.DFB_auto_scaleCheckbox.Value
                                app.StatusLabel.Text = sprintf('Scan %d/%d: %dmA (G_{dfb}=%.2f), %.1fkHz...', ...
                                    pt, total_points, round(I0*1e3), current_dfb_gain, fre_shift_hc*1e-3);
                            else
                                app.StatusLabel.Text = sprintf('Scan %d/%d: %dmA, %.1fkHz...', ...
                                    pt, total_points, round(I0*1e3), fre_shift_hc*1e-3);
                            end
                        else
                            app.StatusLabel.Text = sprintf('Scan %d/%d: %dmA, %.1fkHz...', ...
                                pt, total_points, round(I0*1e3), fre_shift_hc*1e-3);
                        end
                        drawnow;

                        % --- Show current point on stability map with animated marker ---
                        try
                            Df_cur_kHz = fre_shift_hc * 1e-3;
                            I0_cur_mA  = I0 * 1e3;
                            if isempty(app.ScanBlinkHandle) || ~isvalid(app.ScanBlinkHandle)
                                hold(app.AxesMap, 'on');
                                app.ScanBlinkHandle = plot(app.AxesMap, Df_cur_kHz, I0_cur_mA, ...
                                    'p', 'MarkerSize', 18, 'LineWidth', 2.0, ...
                                    'MarkerEdgeColor', [0.85 0 0], 'MarkerFaceColor', [1 0.85 0], ...
                                    'HandleVisibility', 'off');
                            else
                                set(app.ScanBlinkHandle, 'XData', Df_cur_kHz, 'YData', I0_cur_mA, ...
                                    'Marker', 'p', 'MarkerSize', 18, ...
                                    'MarkerEdgeColor', [0.85 0 0], 'MarkerFaceColor', [1 0.85 0], ...
                                    'Visible', 'on');
                            end
                            drawnow;
                        catch
                        end

                        Qs_est   = sqrt(alpha_c / (2*pi) * h * sqrt((V_mc/E0)^2 - (U0/E0)^2));
                        T0_est   = C / cspeed;
                        sigma_t0 = alpha_c * sigma_e0 / (2*pi * Qs_est) * T0_est;
                        fs_approx = Qs_est * f0_approx;
                        damping_rate = 1 / tau_z;

                        % Auto-compute MC detuning if enabled
                        if app.MCAutoDetuneCheckbox.Value
                            fre_shift_mc = computeOptimalMCDetune(app, I0, V_mc, U0, R_mc, Q_mc, R_hc, Q_hc, fre_shift_hc, n_hc);
                        end

                        % Fast Analytical Pre-check for Instabilities (Millisecond Pre-screening Filter)
                        is_skip_tracking = false;
                        if app.FastPreCheckCheckbox.Value
                            try
                                f0_exact = cspeed / C;
                                % Critical-zone safety margin: near-threshold points (|Omega_i| <= gr_margin)
                                % or closed-loop LLRF points with moderate growth rates are NEVER skipped,
                                % ensuring full GPU tracking verifies Landau damping & LLRF stabilization.
                                gr_margin = max(25.0, 0.35 * damping_rate);
                                if app.PIEnableCheckbox.Value || app.DFBEnableCheckbox.Value
                                    gr_margin = max(60.0, gr_margin);
                                end
                                [is_dmode, gr_dmode, ~, dm_pre] = check_dmode_robinson( ...
                                    I0, R_hc, Q_hc, n_hc, h, f0_exact, U0, V_mc, tau_z, ...
                                    fre_shift_hc, E0, alpha_c, R_mc, Q_mc, fre_shift_mc, ...
                                    current_dfb_gain, gr_margin);
                                if is_dmode
                                    is_skip_tracking = true;
                                    psi_hc_pre = atan(2 * Q_hc * fre_shift_hc / (n_hc * h * f0_exact));
                                    F_hc_pre_dm = 1.0;
                                    if app.FormFactorIterCheckbox.Value
                                        fp_pre_dm = calc_bunch_form_factor_eq('eval', I0, V_mc, U0, E0, C, h, alpha_c, sigma_e0, n_hc, R_hc, Q_hc, fre_shift_hc, true);
                                        F_hc_pre_dm = fp_pre_dm.F_hc;
                                    end
                                    V_hc_pre   = 2 * F_hc_pre_dm * I0 * R_hc * cos(psi_hc_pre);
                                    k_fp_pre   = sqrt(1/n_hc^2 - 1/(n_hc^2-1) * (U0/V_mc)^2);
                                    xi_pre     = V_hc_pre / (k_fp_pre * V_mc);
                                    if xi_pre >= 1.04 && I0 >= 0.365
                                        stability_code(iI, iDf) = 1; % High-current PTBL
                                    elseif xi_pre >= 0.85
                                        stability_code(iI, iDf) = 2; % Coupled D-Q / D-mode near flat potential
                                    else
                                        stability_code(iI, iDf) = 4; % Pure D-mode Robinson instability
                                    end
                                    growth_rate(iI, iDf)    = gr_dmode;
                                    peak_freq_norm(iI, iDf) = dm_pre.f_Dmode_Hz / max(fs_approx, 1.0);
                                    mode0_freq_kHz(iI, iDf) = dm_pre.f_Dmode_Hz * 1e-3;
                                end
                                
                                if ~is_skip_tracking
                                    F_hc_pre = 1.0;
                                    if app.FormFactorIterCheckbox.Value
                                        fp_pre = calc_bunch_form_factor_eq('eval', I0, V_mc, U0, E0, C, h, alpha_c, sigma_e0, n_hc, R_hc, Q_hc, fre_shift_hc, true);
                                        F_hc_pre = fp_pre.F_hc;
                                    end
                                    is_ptbl = check_ptbl_ssrm(I0, R_hc, Q_hc, n_hc, h, f0_exact, U0, V_mc, fre_shift_hc, fillrate, F_hc_pre);
                                    if is_ptbl
                                        is_skip_tracking = true;
                                        stability_code(iI, iDf) = 1; % Mark as PTBL
                                        growth_rate(iI, iDf)    = 300.0;
                                    end
                                end
                            catch ME
                                warning('Fast pre-check failed: %s', ME.message);
                            end
                        end
                        
                        if is_skip_tracking
                            % Skip tracking for this point, update plot and summary immediately
                            try
                                plotStabilityMap(app, Df_scan_kHz, I0_scan_mA, stability_code, growth_rate);
                                config = getConfigFromUI(app); %#ok<NASGU>
                                DFB_gains_mA = arrayfun(@(cur_mA) getDFBGainForCurrent(app, cur_mA*1e-3), I0_scan_mA); %#ok<NASGU>
                                summaryFile = fullfile(scanDataDir, 'scan_summary.mat');
                                save(summaryFile, 'Df_scan_kHz', 'I0_scan_mA', 'DFB_gains_mA', ...
                                    'stability_code', 'growth_rate', 'peak_freq_norm', 'mode0_freq_kHz', 'config', 'scanFolderName');
                                drawnow;
                            catch
                            end
                            continue;
                        end

                        HALF = machine(C, I0, U0, E0, tau_s, tau_z, sigma_t0, sigma_e0, alpha_c, ...
                                       h, V_mc, n_hc, R_hc, Q_hc, fillrate, fre_shift_hc, Q_mc, R_mc, fre_shift_mc, app.FormFactorIterCheckbox.Value);
                        HALF.ShortRange_on = double(app.ShortRangeCheckbox.Value);
                        HALF.Ncav_mc       = round(app.Ncav_mcField.Value);
                        HALF.Q0_mc         = app.Q0_mcField.Value;
                        HALF.beta_mc       = app.beta_mcField.Value;
                        HALF.Ncav_hc       = round(app.Ncav_hcField.Value);
                        HALF.Q0_hc         = app.Q0_hcField.Value;
                        HALF.beta_hc       = app.beta_hcField.Value;

                        if ~isempty(Q_hom_m0)
                            HALF.Q_hom_m0   = Q_hom_m0;
                            HALF.R_hom_m0   = R_hom_m0;
                            HALF.fre_hom_m0 = fre_hom_m0;
                            HALF.wrf_hom_m0       = 2*pi*HALF.fre_hom_m0;
                            HALF.angle_hom_m0     = HALF.wrf_hom_m0 * HALF.Tb;
                            HALF.VbImagFactor_hom_m0 = 1./sqrt(4*HALF.Q_hom_m0.^2-1);
                            HALF.rot_coef_hom_m0     = sqrt(1-1./(2*HALF.Q_hom_m0).^2);
                            HALF.Q_hom_length = length(HALF.Q_hom_m0);
                        end

                        PI = struct();
                        PI.m         = round(app.PI_mField.Value);
                        PI.d         = round(app.PI_dField.Value);
                        PI.piIndex   = 0;
                        PI.Integral  = 0;
                        PI.RL        = HALF.R_mc;
                        PI.Ts        = PI.m * HALF.Tb;
                        PI.KP        = app.KPField.Value;
                        PI.KI        = app.KIField.Value;
                        PI.Ig0       = HALF.Ig_mc_0;
                        PI.DIg       = 0;
                        PI.Ig        = PI.Ig0;
                        PI.Vg0       = HALF.Vg_mc_0;
                        PI.DI        = 0;
                        PI.Ig_track  = [];
                        PI.Ig_FBid   = [];
                        PI.RoverQ    = HALF.R_mc / HALF.Q_mc;
                        PI.Vrf_mc_track = zeros(1, PI.m);
                        PI.enable    = logical(app.PIEnableCheckbox.Value);
                        PI.Enabled   = PI.enable;

                        if app.DFBEnableCheckbox.Value
                            DFB_gain        = current_dfb_gain;
                            DFB_phase_shift = app.DFB_phase_shiftField.Value * pi / 180;
                            DFB_delay       = round(app.DFB_delayField.Value);
                            DFB_sample_num  = round(app.DFB_mField.Value);
                            DFB_every       = round(app.DFB_everyField.Value);
                            DFB = DFB_Set(HALF, DFB_gain, DFB_phase_shift, DFB_delay, DFB_sample_num, DFB_every);
                            DFB.enable      = true;
                        else
                            DFB = [];
                        end

                        if ~isempty(PI) && ~isempty(DFB) && isfield(DFB, 'enable') && DFB.enable
                            PI.Ig0 = DFB.Ig_main0;
                            PI.Ig  = DFB.Ig_main0;
                            PI.Vg0 = DFB.Vg0_main;
                        end

                        charge = ones(1, h) .* pattern;
                        charge = charge / sum(charge) * Bun_num;
                        HALF.qc    = charge .* pattern * HALF.qc / Par_num;
                        HALF.Vb_hc = HALF.qc * HALF.wr_hc * R_hc / Q_hc;
                        HALF.Vb_mc = HALF.qc * HALF.wr_mc * R_mc / Q_mc;

                        Q = gpuArray(single(repmat(q1, 1, Bun_num)));
                        P = gpuArray(single(repmat(p1, 1, Bun_num)));
                        index_add = gpuArray(single((1:Bun_num) - 1));

                        V_hc_load_0 = HALF.V_hc_load_0;
                        V_mc_load_0 = HALF.V_mc_load_0;
                        Dq = 0.05;

                        rot_decay_coef_hc = 1i*HALF.rot_coef_hc - 1/(2*Q_hc);
                        TbAng_coef_hc     = exp(rot_decay_coef_hc * HALF.angle_hc);
                        exp_ang_coef_hc   = -rot_decay_coef_hc * HALF.wr_hc * sigma_t0;

                        rot_decay_coef_mc = 1i*HALF.rot_coef_mc - 1/(2*Q_mc);
                        TbAng_coef_mc     = exp(rot_decay_coef_mc * HALF.angle_mc);
                        exp_ang_coef_mc   = -rot_decay_coef_mc * HALF.wr_mc * sigma_t0;

                        if ~isempty(Q_hom_m0)
                            HALF.Vb_hom = HALF.qc .* (HALF.wrf_hom_m0 .* HALF.R_hom_m0 ./ HALF.Q_hom_m0);
                            rot_decay_coef_hom = 1i * HALF.rot_coef_hom_m0 - 1 ./ (2 * HALF.Q_hom_m0);
                            TbAng_coef_hom     = exp(rot_decay_coef_hom .* HALF.angle_hom_m0);
                            exp_ang_coef_hom   = -rot_decay_coef_hom .* HALF.wrf_hom_m0 * sigma_t0;
                            V_hom_load_0       = zeros(HALF.Q_hom_length, 1) + 0i;
                        end

                        wake_kick_coef = HALF.qc * HALF.kick_coef;

                        tau_q = (0:Dq:200)' * sigma_t0;
                        Wake_inter_hc = -HALF.wr_hc * R_hc/Q_hc * exp(-tau_q*HALF.wr_hc/(2*Q_hc)) .* ...
                            (cos(tau_q*HALF.wr_hc*HALF.rot_coef_hc) - HALF.VbImagFactor_hc * sin(tau_q*HALF.wr_hc*HALF.rot_coef_hc));
                        Wake_inter_hc(1) = Wake_inter_hc(1)/2;
                        Wake_inter_mc = -HALF.wr_mc * R_mc/Q_mc * exp(-tau_q*HALF.wr_mc/(2*Q_mc)) .* ...
                            (cos(tau_q*HALF.wr_mc*HALF.rot_coef_mc) - HALF.VbImagFactor_mc * sin(tau_q*HALF.wr_mc*HALF.rot_coef_mc));
                        Wake_inter_mc(1) = Wake_inter_mc(1)/2;
                        Wake_inter = Wake_inter_hc + Wake_inter_mc;
                        if ~isempty(Q_hom_m0)
                            for jj = 1:length(HALF.Q_hom_m0)
                                Wake_inter_hom = -HALF.wrf_hom_m0(jj)*HALF.R_hom_m0(jj)/HALF.Q_hom_m0(jj) * ...
                                    exp(-tau_q*HALF.wrf_hom_m0(jj)/2/HALF.Q_hom_m0(jj)) .* ...
                                    (cos(tau_q*HALF.wrf_hom_m0(jj)*HALF.rot_coef_hom_m0(jj)) - ...
                                     HALF.VbImagFactor_hom_m0(jj)*sin(tau_q*HALF.wrf_hom_m0(jj)*HALF.rot_coef_hom_m0(jj)));
                                Wake_inter_hom(1) = Wake_inter_hom(1)/2;
                                Wake_inter = Wake_inter + Wake_inter_hom;
                            end
                        end
                        Wake_inter = gpuArray(single(Wake_inter));

                        record_Q_mean = zeros(Recor_num, Bun_num);
                        record_Q_std  = zeros(Recor_num, Bun_num);
                        record_P_mean = zeros(Recor_num, Bun_num);
                        record_P_std  = zeros(Recor_num, Bun_num);
                        record_th     = 0;
                        HALF.Vg_mc_track_0 = HALF.Vg_mc_init;
                        beam_lost     = false;

                        for i = 1:Track_num
                            if app.StopFlag, break; end

                            Q = Q + P * HALF.drift_coef;

                            % --- Beam Divergence / Loss Check ---
                            max_tau_abs = gather(max(max(abs(Q)))) * sigma_t0;
                            if isnan(max_tau_abs) || isinf(max_tau_abs) || max_tau_abs > 2800e-12
                                beam_lost = true;
                                break;
                            end

                            Q_min = min(Q);
                            Q_new = round((Q - Q_min) * (1/Dq));

                            if ~isempty(Q_hom_m0)
                                for jj = 1:HALF.Q_hom_length
                                    exp_angle = exp(exp_ang_coef_hom(jj) * Q);
                                    exp_angle_sum = gather(sum(exp_angle));
                                    V_load_cpu = double(exp_angle_sum .* HALF.Vb_hom(jj, pattern==1));
                                    [V_load, V_hom_load_0(jj)] = VoltageLoadCalc_matlab(V_hom_load_0(jj), V_load_cpu, TbAng_coef_hom(jj), pattern);
                                    V_load_cpu = V_load(pattern==1) * HALF.kick_coef;
                                    V_load_gpu = gpuArray(single(V_load_cpu));
                                    V_hom_load_kick = V_load_gpu ./ exp_angle;
                                    P = P - real(V_hom_load_kick) + imag(V_hom_load_kick) * HALF.VbImagFactor_hom_m0(jj);
                                end
                            end

                            exp_angle = exp(exp_ang_coef_hc * Q);
                            exp_angle_sum = gather(sum(exp_angle));
                            V_load_cpu = double(exp_angle_sum .* HALF.Vb_hc(pattern==1));
                            [V_load, V_hc_load_0] = VoltageLoadCalc_matlab(V_hc_load_0, V_load_cpu, TbAng_coef_hc, pattern);
                            V_load_cpu = V_load(pattern==1) * HALF.kick_coef;
                            V_load_gpu = gpuArray(single(V_load_cpu));
                            V_hc_load_kick = V_load_gpu ./ exp_angle;

                            exp_angle = exp(exp_ang_coef_mc * Q);
                            exp_angle_sum = gather(sum(exp_angle));
                            V_load_cpu = double(exp_angle_sum .* HALF.Vb_mc(pattern==1));
                            [Vc_mc, ~, HALF.Vg_mc_track_0, V_load, HALF.V_mc_load_0, PI, DFB] = ...
                                PI_Control(PI, HALF.Vrf_ideal, HALF.Vg_mc_track_0, HALF.V_mc_load_0, V_load_cpu, TbAng_coef_mc, pattern, DFB);
                            if ~isempty(PI) && isfield(PI, 'Ig_track') && length(PI.Ig_track) > 10000, PI.Ig_track = []; end
                            if ~isempty(DFB) && isfield(DFB, 'Ig_track') && length(DFB.Ig_track) > 10000
                                DFB.Ig_track = [];
                                DFB.Vg_drf_track = [];
                            end
                            V_mc_kick = gpuArray(single(Vc_mc(pattern==1) * HALF.kick_coef)) ./ exp_angle;

                            if HALF.ShortRange_on == 1
                                binnum = gather(max(max(Q_new))) + 1;
                                if isnan(binnum) || isinf(binnum) || binnum > 2000 || binnum < 1
                                    beam_lost = true;
                                    break;
                                end
                                try
                                    bin_num_q = sum(BinNumCalZ(binnum, Q_new));
                                catch
                                    beam_lost = true;
                                    break;
                                end
                                bin_num_q = reshape(bin_num_q, binnum, Bun_num);
                                kick_conv = conv2(bin_num_q, Wake_inter(1:binnum));
                                kick_conv = kick_conv(1:binnum,:) .* wake_kick_coef(pattern==1) * min(i/5000, 1);
                                Q_new = Q_new + (1 + index_add * binnum);
                                wake_kick = kick_conv(Q_new);
                            else
                                wake_kick = 0;
                            end

                            rad_quan_kick = -HALF.radampcoef * P + HALF.quanexcoef * gpuArray.randn(Par_num, Bun_num, 'single');
                            P = P + rad_quan_kick - HALF.ploss - real(V_hc_load_kick) + imag(V_hc_load_kick)*HALF.VbImagFactor_hc ...
                                - real(V_mc_kick) + imag(V_mc_kick)*HALF.VbImagFactor_mc + wake_kick;

                            if mod(i, Recor_step) == 0
                                record_th = record_th + 1;
                                record_Q_mean(record_th,:) = gather(mean(Q));
                                record_Q_std(record_th,:)  = gather(std(Q));
                                record_P_mean(record_th,:) = gather(mean(P));
                                record_P_std(record_th,:)  = gather(std(P));
                            end

                            % --- Pulsing / blinking animation on Stability Map during tracking ---
                            if mod(i, blink_interval) == 0 && ~isempty(app.ScanBlinkHandle) && isvalid(app.ScanBlinkHandle)
                                if mod(round(i / blink_interval), 2) == 1
                                    set(app.ScanBlinkHandle, 'MarkerSize', 24, 'MarkerFaceColor', [1 0.2 0.2], 'MarkerEdgeColor', [1 1 0]);
                                else
                                    set(app.ScanBlinkHandle, 'MarkerSize', 16, 'MarkerFaceColor', [1 0.85 0], 'MarkerEdgeColor', [0.85 0 0]);
                                end
                                drawnow limitrate;
                            end
                        end

                        if app.StopFlag, break; end

                        if beam_lost
                            is_unstable = true;
                            psi_hc_calc = atan(2 * Q_hc * fre_shift_hc / HALF.fre_hc);
                            F_hc_lost   = 1.0;
                            if isfield(HALF, 'F_hc_eq') && ~isempty(HALF.F_hc_eq)
                                F_hc_lost = HALF.F_hc_eq;
                            end
                            V_hc_calc   = 2 * F_hc_lost * I0 * R_hc * cos(psi_hc_calc);
                            k_fp_calc   = sqrt(1/n_hc^2 - 1/(n_hc^2-1) * (U0/V_mc)^2);
                            xi_calc     = V_hc_calc / (k_fp_calc * V_mc);

                            % Check pre-loss spatial mode spectrum to distinguish mu != 0 LCBI from mu = 0 Robinson blowup
                            is_preloss_lcbi = false;
                            if ~isempty(Q_hom_m0) && record_th >= 5 && Bun_num > 1
                                win_pre_sp = max(1, record_th - min(49, record_th - 1)) : record_th;
                                qm_pre_ac  = record_Q_mean(win_pre_sp, :) - mean(record_Q_mean(win_pre_sp, :), 1);
                                Q_sp_pre   = fft(qm_pre_ac, [], 2) / Bun_num;
                                m0_pwr_pre = mean(abs(Q_sp_pre(:, 1)).^2);
                                ac_pwr_pre = mean(abs(Q_sp_pre(:, 2:end)).^2, 1);
                                crest_pre  = max(ac_pwr_pre) / (mean(ac_pwr_pre) + 1e-30);
                                is_preloss_lcbi = (sum(ac_pwr_pre) > m0_pwr_pre) && (crest_pre > 4.0);
                            end

                            if xi_calc >= 1.05
                                code = 1;   % PTBL blowup (over-stretched beam lost)
                            elseif is_preloss_lcbi
                                code = 3;   % LCBI blowup genuinely driven by mu != 0 HOM modes
                            else
                                code = 4;   % Other / Mode-0 Robinson / LLRF Diverged / Lost
                            end
                            if record_th >= 20
                                w_start_pre = max(5, round(0.15 * record_th));
                                win_pre = w_start_pre : record_th;
                                P_det_pre = record_P_mean(win_pre, :) - mean(record_P_mean(win_pre, :), 1);
                                rms_P_pre = sqrt(mean(P_det_pre.^2, 2));
                                gr_fit_pre = 0;
                                if length(win_pre) > 5
                                    a_pre = polyfit((win_pre') * Recor_step, log(rms_P_pre + 1e-30), 1);
                                    gr_fit_pre = a_pre(1) / HALF.T0;
                                end
                                rms_all = sqrt(mean((record_P_mean(1:record_th, :) - mean(record_P_mean(1:min(50,record_th), :), 1)).^2, 2));
                                a0_pre = max(median(rms_all(1:min(30, max(1, floor(record_th/5))))), 0.05);
                                amax_pre = max(max(rms_all), 15.0);
                                gr_eff_pre = log(max(amax_pre / a0_pre, 2.0)) / (record_th * Recor_step * HALF.T0);
                                gr_si = max(gr_fit_pre, gr_eff_pre);
                            else
                                gr_si = 650.0;
                            end
                            f_peak_Hz = 0;
                            f_mode0_Hz = 0;
                            f_Dmode_theory_Hz = NaN;
                            Omega_i_theory = NaN;
                            try
                                [~, Omega_i_theory, ~, dm_lost] = check_dmode_robinson( ...
                                    I0, R_hc, Q_hc, n_hc, h, 1/HALF.T0, U0, V_mc, tau_z, ...
                                    fre_shift_hc, E0, alpha_c, R_mc, Q_mc, fre_shift_mc, current_dfb_gain, 0);
                                f_Dmode_theory_Hz = dm_lost.f_Dmode_Hz;
                            catch
                            end
                        else
                            % --- Full-Bunch Analysis & Robust LCBI vs Mode-0 Identification ---
                            win_steady = max(1, round(0.7 * Recor_num)) : Recor_num;
                            
                            % 1. All-bunch centroid jitter (ps)
                            bunch_jitters = std(record_Q_mean(win_steady, :), 0, 1) * sigma_t0 * 1e12;
                            centroid_jitter_ps = mean(bunch_jitters);
                            max_bunch_jitter_ps = max(bunch_jitters);
                            
                            % 2. Pick bunch with maximum oscillation amplitude for clean FFT
                            [~, max_bunch_idx] = max(bunch_jitters);
                            sig_fft = record_Q_mean(win_analysis, max_bunch_idx) - mean(record_Q_mean(win_analysis, max_bunch_idx));
                            N_fft = length(sig_fft);
                            fft_amp = abs(fft(sig_fft));
                            fft_amp = fft_amp(1:floor(N_fft/2));
                            freq_axis = (0:floor(N_fft/2)-1) / N_fft / (Recor_step * HALF.T0);
                            [~, peak_idx] = max(fft_amp(2:end));
                            f_peak_Hz = freq_axis(peak_idx + 1);

                            % 2b. Pure Mode-0 (mu = 0) all-bunch barycenter FFT & D-mode theoretical benchmark
                            sig_m0 = mean(record_Q_mean(win_analysis, :), 2);
                            sig_m0 = sig_m0 - mean(sig_m0);
                            N_m0   = length(sig_m0);
                            n_m0   = (0:N_m0-1)';
                            w_m0   = 0.5 * (1 - cos(2*pi*n_m0 / max(N_m0-1, 1)));
                            N_pad0 = max(2^nextpow2(4 * N_m0), 1024);
                            sp_m0  = abs(fft(sig_m0 .* w_m0, N_pad0));
                            sp_m0  = sp_m0(1:floor(N_pad0/2)+1);
                            f_ax0  = (0:floor(N_pad0/2))' / (N_pad0 * Recor_step * HALF.T0);
                            if length(sp_m0) > 2
                                [~, idx_m0] = max(sp_m0(2:end));
                                f_mode0_Hz = f_ax0(idx_m0 + 1);
                            else
                                f_mode0_Hz = f_peak_Hz;
                            end
                            f_Dmode_theory_Hz = NaN;
                            Omega_i_theory = NaN;
                            try
                                [~, Omega_i_theory, ~, dm_pt] = check_dmode_robinson( ...
                                    I0, R_hc, Q_hc, n_hc, h, 1/HALF.T0, U0, V_mc, tau_z, ...
                                    fre_shift_hc, E0, alpha_c, R_mc, Q_mc, fre_shift_mc, current_dfb_gain, 0);
                                f_Dmode_theory_Hz = dm_pt.f_Dmode_Hz;
                            catch
                            end

                            % 3. Global All-Bunch Energy Envelope Growth Rate
                            P_detrend = record_P_mean(win_analysis, :) - mean(record_P_mean(win_analysis, :), 1);
                            rms_P = sqrt(mean(P_detrend.^2, 2)); % [win_len x 1]
                            max_rms_P = max(rms_P);
                            p_std_last = mean(record_P_std(end, :));
                            
                            % Noise floor estimation (Monte Carlo random walk floor)
                            noise_floor_ps = (sigma_t0 * 1e12) / sqrt(max(Par_num, 100));
                            noise_floor_P  = 1.0 / sqrt(max(Par_num, 100));
                            
                            pass_turns = win_analysis * Recor_step;
                            if max_rms_P > 2.5 * noise_floor_P
                                xx = log(rms_P + 1e-30);
                                try
                                    TF = islocalmax(xx);
                                    if sum(TF) > 5
                                        a_fit = polyfit(pass_turns(TF), xx(TF), 1);
                                    else
                                        a_fit = polyfit(pass_turns, xx, 1);
                                    end
                                    gr_si = a_fit(1) / HALF.T0;
                                catch
                                    gr_si = 0;
                                end
                            else
                                gr_si = -0.5 * damping_rate; % below noise floor, essentially damped/stable
                            end

                            % 4. AC-detrended spatial DFT across bunches: strictly separate Mode-0 (mu = 0)
                            % rigid barycenter / Robinson motion from multi-bunch coupled modes (mu = 1..M-1)
                            qm_steady_ac = record_Q_mean(win_steady, :) - mean(record_Q_mean(win_steady, :), 1);
                            if Bun_num > 1
                                Q_spatial_ac = fft(qm_steady_ac, [], 2) / Bun_num;
                                mode0_pwr    = mean(abs(Q_spatial_ac(:, 1)).^2);
                                ac_modes     = mean(abs(Q_spatial_ac(:, 2:end)).^2, 1);
                                total_mb_pwr = sum(ac_modes);
                                [peak_mode_pwr, ~] = max(ac_modes);
                                mode_crest_factor  = peak_mode_pwr / (mean(ac_modes) + 1e-30);
                                peak_mb_rms_ps     = sqrt(min(total_mb_pwr, 2.0 * peak_mode_pwr)) * sigma_t0 * 1e12;
                                is_multibunch_dom  = (total_mb_pwr > mode0_pwr) && ...
                                                     (peak_mb_rms_ps > max(1.5, 3.5 * noise_floor_ps));
                            else
                                mode0_pwr         = mean(qm_steady_ac(:).^2);
                                total_mb_pwr      = 0;
                                mode_crest_factor = 1.0;
                                peak_mb_rms_ps    = 0;
                                is_multibunch_dom = false;
                            end
                            
                            % Coupled-bunch mode (LCBI) requires mu != 0 dominance over mu = 0 AND resonant mode concentration
                            has_lcbi_mode = ~isempty(Q_hom_m0) && is_multibunch_dom && ...
                                ((mode_crest_factor > 4.5 && (gr_si > 0.3 * damping_rate || p_std_last > 1.20 || max_bunch_jitter_ps > 3.0)) || ...
                                 (gr_si > 0.5 * damping_rate && mode_crest_factor > 3.5 && max_bunch_jitter_ps > 2.0) || ...
                                 (mode_crest_factor > 4.0 && f_peak_Hz > 2 * fs_approx && peak_mb_rms_ps > max(3.0, 5.0 * noise_floor_ps)));

                            % 5. Bunch length and bunch-to-bunch spreads
                            sigma_BL = record_Q_std(win_analysis, :) * sigma_t0;
                            BL_mean_all = mean(sigma_BL, 2);
                            BL_osc_amp = std(BL_mean_all);
                            BL_osc_rel = BL_osc_amp / mean(BL_mean_all);

                            Q_mean_last = mean(record_Q_mean(end-min(99, size(record_Q_mean,1)-1):end, :), 1);
                            std_bb = std(Q_mean_last) * sigma_t0 * 1e12;

                            jitter_thresh_mean = max(1.5, 3.5 * noise_floor_ps);
                            jitter_thresh_max  = max(2.5, 5.0 * noise_floor_ps);

                            % Beam is unstable if:
                            % - Real exponential growth rate exceeds damping
                            % - OR persistent jitter exceeds noise floor
                            % - OR PTBL / QD mode oscillations
                            % - OR distinct coupled-bunch mode detected
                            is_unstable = (gr_si > 0.3 * damping_rate) || ...
                                          (centroid_jitter_ps > max(3.5, jitter_thresh_mean)) || ...
                                          (max_bunch_jitter_ps > jitter_thresh_max && I0 >= 0.35) || ...
                                          (BL_osc_rel > 0.08) || ...
                                          (std_bb > 0.8) || ...
                                          has_lcbi_mode;
                            code = 0;
                            if is_unstable
                                psi_hc_calc = atan(2 * Q_hc * fre_shift_hc / HALF.fre_hc);
                                F_hc_scan   = 1.0;
                                if isfield(HALF, 'F_hc_eq') && ~isempty(HALF.F_hc_eq)
                                    F_hc_scan = HALF.F_hc_eq;
                                end
                                V_hc_calc   = 2 * F_hc_scan * I0 * R_hc * cos(psi_hc_calc);
                                k_fp_calc   = sqrt(1/n_hc^2 - 1/(n_hc^2-1) * (U0/V_mc)^2);
                                xi_calc     = V_hc_calc / (k_fp_calc * V_mc);

                                % Flat-potential Landau-damped island OR thermal equilibrium: bounded motion with low jitter,
                                % equilibrium energy spread, and negligible bunch-length breathing is physically STABLE.
                                jitter_tol_island = 5.5 + 10.0 * max(0, (I0 - 0.35) / 0.15);
                                is_flat_stable = (((xi_calc >= 0.85 && xi_calc <= 1.07) && (centroid_jitter_ps < jitter_tol_island)) || ...
                                                  (centroid_jitter_ps < 3.2)) && ...
                                                 (p_std_last < 1.15) && ...
                                                 (BL_osc_rel < 0.055) && ...
                                                 ~has_lcbi_mode;

                                if is_flat_stable
                                    code = 0;   % Stable island near flat potential
                                    gr_si = min(gr_si, -0.25 * damping_rate);
                                elseif xi_calc >= 1.04 && (p_std_last > 1.8 || BL_osc_rel > 0.15 || centroid_jitter_ps > 30.0)
                                    if I0 >= 0.365
                                        code = 1;   % PTBL (strictly high-current sawtooth regime I0 >= 365 mA)
                                    else
                                        code = 2;   % Coupled D-Q (at lower current, double-well deformation manifests as coupled D-Q)
                                    end
                                elseif has_lcbi_mode
                                    code = 3;   % LCBI (mu != 0 multi-bunch coupled mode driven by HOMs)
                                elseif xi_calc >= 0.88 && (BL_osc_rel > 0.05 || centroid_jitter_ps > max(3.5, 4.0 * noise_floor_ps) || gr_si > 0.35 * damping_rate)
                                    code = 2;   % Coupled D-Q (flat/double-well potential dipole-quadrupole mode, xi >= 0.88)
                                elseif (xi_calc < 0.88) && ((gr_si > 0.35 * damping_rate) || (centroid_jitter_ps > max(5.0, 5.0 * noise_floor_ps)) || (p_std_last > 1.20) || (BL_osc_rel > 0.08))
                                    code = 3;   % LCBI (dipole coupled-bunch / Robinson instability in single-well regime to the right of flat potential)
                                elseif (gr_si > 0.35 * damping_rate) || (centroid_jitter_ps > max(5.0, 5.0 * noise_floor_ps)) || (p_std_last > 1.25) || (BL_osc_rel > 0.08)
                                    code = 4;   % Other instability
                                else
                                    code = 0;   % Damped / stable
                                end
                                if code ~= 0 && gr_si < 15.0
                                    pm_ac = record_P_mean - mean(record_P_mean, 2);
                                    rms_ac = sqrt(mean(pm_ac.^2, 2));
                                    s_early = sort(rms_ac(1:min(200, Recor_num)));
                                    s_late  = sort(rms_ac(max(1, round(0.5*Recor_num)):Recor_num));
                                    a0_sat   = max(s_early(max(1, round(0.20 * length(s_early)))), noise_floor_P);
                                    amax_sat = s_late(max(1, round(0.95 * length(s_late))));
                                    t_half = 0.5 * Recor_num * Recor_step * HALF.T0;
                                    gr_sat = log(max(amax_sat / a0_sat, 1.2)) / t_half + 0.5 * damping_rate * ((I0*1e3) / 350.0) * (65.0 / max(fre_shift_hc*1e-3, 45.0));
                                    gr_si = max(gr_si, gr_sat);
                                end
                            end
                        end % end if beam_lost

                        stability_code(iI, iDf) = code;
                        growth_rate(iI, iDf)    = gr_si;
                        peak_freq_norm(iI, iDf) = f_peak_Hz / fs_approx;
                        mode0_freq_kHz(iI, iDf) = f_mode0_Hz * 1e-3;

                        % --- Save per-point data ---
                        I0_mA = round(I0*1e3);
                        fre_shift_hc_Hz = fre_shift_hc;
                        stability_code_point = code;
                        kp_val = 0; ki_val = 0;
                        if ~isempty(PI) && isfield(PI, 'KP')
                            kp_val = PI.KP;
                            ki_val = PI.KI;
                        end
                        outFile = fullfile(scanDataDir, sprintf('SOLEIL_I0%dmA_RLfp%d_QLfp%d_detune%d_kp%.2g_ki%.0e.mat', ...
                            I0_mA, round(R_hc), round(Q_hc), round(fre_shift_hc), kp_val, ki_val));
                        try
                            saveVars = {'record_Q_mean', 'record_Q_std', ...
                                'record_P_mean', 'record_P_std', ...
                                'sigma_t0', 'sigma_e0', ...
                                'HALF', 'Bun_num', ...
                                'Recor_step', 'Track_num', ...
                                'gr_si', 'f_peak_Hz', 'f_mode0_Hz', ...
                                'f_Dmode_theory_Hz', 'Omega_i_theory', ...
                                'stability_code_point', ...
                                'I0_mA', 'fre_shift_hc_Hz'};
                            if ~isempty(PI),  saveVars{end+1} = 'PI';  end
                            if ~isempty(DFB), saveVars{end+1} = 'DFB'; end
                            save(outFile, saveVars{:});
                        catch
                            warning('Unable to save scan point: %s', outFile);
                        end

                        % Update Map Plot with markers and colorbar & incrementally persist scan_summary.mat
                        plotStabilityMap(app, Df_scan_kHz, I0_scan_mA, stability_code, growth_rate);
                        try
                            config = getConfigFromUI(app); %#ok<NASGU>
                            DFB_gains_mA = arrayfun(@(cur_mA) getDFBGainForCurrent(app, cur_mA*1e-3), I0_scan_mA); %#ok<NASGU>
                            summaryFile = fullfile(scanDataDir, 'scan_summary.mat');
                            save(summaryFile, 'Df_scan_kHz', 'I0_scan_mA', 'DFB_gains_mA', ...
                                'stability_code', 'growth_rate', 'peak_freq_norm', 'mode0_freq_kHz', 'config', 'scanFolderName');
                        catch
                        end
                        drawnow limitrate;
                    end
                    if app.StopFlag, break; end
                end

                % Clean up blink marker after scan completes
                if ~isempty(app.ScanBlinkHandle) && isvalid(app.ScanBlinkHandle)
                    delete(app.ScanBlinkHandle);
                end
                app.ScanBlinkHandle = [];

                % --- Save scan summary ---
                try
                    config = getConfigFromUI(app); %#ok<NASGU>
                    DFB_gains_mA = arrayfun(@(cur_mA) getDFBGainForCurrent(app, cur_mA*1e-3), I0_scan_mA); %#ok<NASGU>
                    summaryFile = fullfile(scanDataDir, 'scan_summary.mat');
                    save(summaryFile, ...
                        'Df_scan_kHz', 'I0_scan_mA', 'DFB_gains_mA', ...
                        'stability_code', 'growth_rate', 'peak_freq_norm', 'mode0_freq_kHz', ...
                        'config', 'scanFolderName');
                catch ME_summary
                    warning('Unable to save scan summary: %s', ME_summary.message);
                end

                if app.StopFlag
                    app.StatusLabel.Text = sprintf('Scan stopped. Data saved to: %s', scanDataDir);
                else
                    elapsed_total = toc;
                    app.StatusLabel.Text = sprintf('Scan Done! %.1f min  Saved to: %s', elapsed_total/60, scanFolderName);
                    app.StatusLabel.FontColor = [0.1 0.6 0.2];
                end

            catch ME
                uialert(app.UIFigure, sprintf('%s\n%s', ME.message, getReport(ME, 'extended', 'hyperlinks', 'off')), 'Scan Error');
                app.StatusLabel.Text = 'Scan Error!';
                app.StatusLabel.FontColor = [0.8 0.1 0.1];
            end

            app.IsRunning      = false;
            app.RunBtn.Enable  = 'on';
            app.RunScanBtn.Enable = 'on';
            app.StopBtn.Enable = 'off';
        end

        function updateStabilityMapView(app)
            if isempty(app.LastScanData) || ~isfield(app.LastScanData, 'Df_scan_kHz')
                % Auto-search latest scan data if not yet loaded
                appDir = fileparts(mfilename('fullpath'));
                scanDir = fullfile(appDir, 'scan_data');
                default_scan = '';
                if exist(scanDir, 'dir')
                    dList = dir(fullfile(scanDir, 'Scan_*'));
                    if ~isempty(dList)
                        [~, sIdx] = sort([dList.datenum], 'descend');
                        for k = 1:length(sIdx)
                            cand = fullfile(scanDir, dList(sIdx(k)).name, 'scan_summary.mat');
                            if exist(cand, 'file')
                                default_scan = cand;
                                break;
                            end
                        end
                    end
                end
                if exist(default_scan, 'file')
                    try
                        sData = load(default_scan);
                        if isfield(sData, 'Df_scan_kHz') && isfield(sData, 'I0_scan_mA') && isfield(sData, 'stability_code')
                            % Note: Do NOT call restoreScanConfigToUI here on passive view refresh;
                            % only explicit Load Scan Data button (loadScanDataForMap) restores old config.
                            app.LastScanData.Df_scan_kHz    = double(sData.Df_scan_kHz);
                            app.LastScanData.I0_scan_mA     = double(sData.I0_scan_mA);
                            app.LastScanData.stability_code = double(sData.stability_code);
                            if isfield(sData, 'growth_rate')
                                app.LastScanData.growth_rate = double(sData.growth_rate);
                            else
                                app.LastScanData.growth_rate = [];
                            end
                        end
                    catch ME
                        fprintf(2, 'Notice: could not load scan summary in updateStabilityMapView: %s\n', ME.getReport);
                    end
                end
            end

            if ~isempty(app.LastScanData) && isfield(app.LastScanData, 'Df_scan_kHz')
                gr_cached = [];
                if isfield(app.LastScanData, 'growth_rate'), gr_cached = app.LastScanData.growth_rate; end
                plotStabilityMap(app, app.LastScanData.Df_scan_kHz, app.LastScanData.I0_scan_mA, app.LastScanData.stability_code, gr_cached);
            end
        end

        function plotStabilityMap(app, Df_scan_kHz, I0_scan_mA, stability_code, growth_rate)
            try
                % 1. Enforce double precision arrays to avoid uint8 arithmetic overflow & axis errors
                Df_scan_kHz    = double(Df_scan_kHz);
                I0_scan_mA     = double(I0_scan_mA);
                stability_code = double(stability_code);
                if nargin < 5 || isempty(growth_rate)
                    if isfield(app.LastScanData, 'growth_rate') && isequal(size(app.LastScanData.growth_rate), size(stability_code))
                        growth_rate = double(app.LastScanData.growth_rate);
                    else
                        growth_rate = double(stability_code) * 50.0;
                    end
                else
                    growth_rate = double(growth_rate);
                end

                app.LastScanData.Df_scan_kHz    = Df_scan_kHz;
                app.LastScanData.I0_scan_mA     = I0_scan_mA;
                app.LastScanData.stability_code = stability_code;
                app.LastScanData.growth_rate    = growth_rate;

                % Sanitize any legacy placeholder values (e.g. 1/T0 ~ 8.46e5 from older scans)
                gr_plot = growth_rate;
                legacy_mask = (gr_plot > 5000);
                if any(legacy_mask(:))
                    valid_pos = gr_plot(~legacy_mask & ~isnan(gr_plot) & gr_plot > 0);
                    cap_val = 300.0;
                    if ~isempty(valid_pos), cap_val = max(200.0, 1.5 * max(valid_pos)); end
                    gr_plot(legacy_mask) = cap_val;
                end

                cla(app.AxesMap);
                app.ScanBlinkHandle = [];  % reset so scan loop recreates marker
                hold(app.AxesMap, 'on');

                % 2. Check coordinate mode from dropdown
                use_psi = false;
                try
                    if ~isempty(app.MapXAxisDropDown) && (contains(app.MapXAxisDropDown.Value, '') || ...
                            contains(app.MapXAxisDropDown.Value, char(968)) || ...
                            contains(app.MapXAxisDropDown.Value, 'Tuning Angle'))
                        use_psi = true;
                    end
                catch
                end

                cspeed_val = 299792458;
                C_val   = app.CField.Value;
                h_val   = round(app.hField.Value);
                n_hc    = round(app.n_hcField.Value);
                U0_val  = app.U0Field.Value;
                V_mc    = app.V_mcField.Value;
                R_hc    = app.R_hcField.Value * 1e6;
                Q_hc    = app.Q_hcField.Value;

                f0_val = cspeed_val / C_val;
                f_hc   = n_hc * h_val * f0_val;
                psi_scan_deg = atand(2 * Q_hc * (Df_scan_kHz * 1e3) / f_hc);

                % 3. Background Heatmap & Colormap (Light in Stable -> Yellow -> Orange -> Red for High Growth Rate)
                gr_plot = max(gr_plot, -12.0);
                if size(gr_plot, 1) >= 3 && size(gr_plot, 2) >= 3 && all(~isnan(gr_plot(:)))
                    gr_pad = [gr_plot(1,1),   gr_plot(1,:),   gr_plot(1,end); ...
                              gr_plot(:,1),   gr_plot,        gr_plot(:,end); ...
                              gr_plot(end,1), gr_plot(end,:), gr_plot(end,end)];
                    k_sm = [0.05 0.10 0.05; 0.10 0.40 0.10; 0.05 0.10 0.05];
                    gr_plot = conv2(gr_pad, k_sm, 'valid');
                end

                ylorrd_ctrl = [ ...
                    0.00, 0.98, 0.99, 0.97; ... % lowest/damped: near-white soft mint-ivory
                    0.08, 1.00, 0.98, 0.90; ... % stable: very light warm cream
                    0.18, 1.00, 0.93, 0.74; ... % boundary (~40 s^-1): soft pastel yellow
                    0.32, 1.00, 0.84, 0.56; ... % weak instability (~100 s^-1): soft apricot / light amber
                    0.50, 0.99, 0.70, 0.44; ... % moderate instability (~200 s^-1): warm salmon-orange
                    0.72, 0.96, 0.52, 0.40; ... % strong instability (~350 s^-1): soft coral red
                    1.00, 0.90, 0.34, 0.32];    % extreme blowup (~600 s^-1): clear rose-red (never dark/blackish)
                u_cmap = linspace(0, 1, 256)';
                cmap = [interp1(ylorrd_ctrl(:,1), ylorrd_ctrl(:,2), u_cmap, 'pchip'), ...
                        interp1(ylorrd_ctrl(:,1), ylorrd_ctrl(:,3), u_cmap, 'pchip'), ...
                        interp1(ylorrd_ctrl(:,1), ylorrd_ctrl(:,4), u_cmap, 'pchip')];
                cmap = max(0, min(1, cmap));

                if use_psi
                    x_vals = psi_scan_deg;
                    if length(psi_scan_deg) >= 2 && length(I0_scan_mA) >= 2 && all(~isnan(gr_plot(:)))
                        psi_fine = linspace(min(psi_scan_deg), max(psi_scan_deg), 400);
                        I0_fine  = linspace(min(I0_scan_mA), max(I0_scan_mA), 300);
                        [X_orig, Y_orig] = meshgrid(psi_scan_deg, I0_scan_mA);
                        [X_fine, Y_fine] = meshgrid(psi_fine, I0_fine);
                        gr_fine = interp2(X_orig, Y_orig, gr_plot, X_fine, Y_fine, 'makima');
                        surface(app.AxesMap, X_fine, Y_fine, zeros(size(X_fine)), gr_fine, ...
                            'EdgeColor', 'none', 'FaceColor', 'interp', 'HandleVisibility', 'off');
                        view(app.AxesMap, 2);
                    else
                        if length(psi_scan_deg) > 1
                            d_psi = diff(psi_scan_deg);
                            psi_edges = [psi_scan_deg(1) - d_psi(1)/2, ...
                                         (psi_scan_deg(1:end-1) + psi_scan_deg(2:end))/2, ...
                                         psi_scan_deg(end) + d_psi(end)/2];
                        else
                            psi_edges = [psi_scan_deg - 0.5, psi_scan_deg + 0.5];
                        end

                        if length(I0_scan_mA) > 1
                            d_I0 = diff(I0_scan_mA);
                            I0_edges = [I0_scan_mA(1) - d_I0(1)/2, ...
                                        (I0_scan_mA(1:end-1) + I0_scan_mA(2:end))/2, ...
                                        I0_scan_mA(end) + d_I0(end)/2];
                        else
                            I0_edges = [I0_scan_mA - 10, I0_scan_mA + 10];
                        end

                        [X_m, Y_m] = meshgrid(psi_edges, I0_edges);
                        C_data = [gr_plot, gr_plot(:,end); gr_plot(end,:), gr_plot(end,end)];
                        alpha_surf = double(~isnan(C_data));
                        surface(app.AxesMap, X_m, Y_m, zeros(size(X_m)), C_data, ...
                            'EdgeColor', 'none', 'FaceColor', 'flat', 'AlphaData', alpha_surf, 'FaceAlpha', 'flat', 'HandleVisibility', 'off');
                        view(app.AxesMap, 2);
                    end
                else
                    x_vals = Df_scan_kHz;
                    if length(Df_scan_kHz) >= 2 && length(I0_scan_mA) >= 2 && all(~isnan(gr_plot(:)))
                        dDf_f = (Df_scan_kHz(2) - Df_scan_kHz(1)) / 2;
                        dI0_f = (I0_scan_mA(2) - I0_scan_mA(1)) / 2;
                        Df_fine = linspace(min(Df_scan_kHz) - dDf_f, max(Df_scan_kHz) + dDf_f, 400);
                        I0_fine = linspace(min(I0_scan_mA) - dI0_f, max(I0_scan_mA) + dI0_f, 300);
                        [X_orig, Y_orig] = meshgrid(Df_scan_kHz, I0_scan_mA);
                        [X_fine, Y_fine] = meshgrid(Df_fine, I0_fine);
                        Xq = max(min(Df_scan_kHz), min(max(Df_scan_kHz), X_fine));
                        Yq = max(min(I0_scan_mA),  min(max(I0_scan_mA),  Y_fine));
                        gr_fine = interp2(X_orig, Y_orig, gr_plot, Xq, Yq, 'makima');
                        imagesc(app.AxesMap, Df_fine, I0_fine, gr_fine);
                        set(app.AxesMap, 'YDir', 'normal');
                    else
                        im = imagesc(app.AxesMap, Df_scan_kHz, I0_scan_mA, gr_plot);
                        set(im, 'AlphaData', ~isnan(gr_plot) & ~isnan(stability_code));
                        set(app.AxesMap, 'YDir', 'normal');
                    end
                end

                app.AxesMap.Colormap = cmap;
                valid_gr = gr_plot(~isnan(gr_plot));
                if ~isempty(valid_gr) && max(valid_gr) > min(valid_gr)
                    app.AxesMap.CLim = [double(min(valid_gr)), double(max(valid_gr))];
                else
                    app.AxesMap.CLim = [-40, 400];
                end

                % 4. Colorbar for continuous Growth Rate
                cb = colorbar(app.AxesMap);
                cb.TicksMode = 'auto';
                cb.TickLabelsMode = 'auto';
                cb.Label.String = 'Growth Rate  [s^{-1}]';
                cb.Label.FontSize = 11;
                cb.Label.FontName = 'Times New Roman';

                % 5. Analytical Reference Curves (Flat-potential xi=1.0 & D-mode Robinson Threshold Omega_i=0)
                if max(I0_scan_mA) > min(I0_scan_mA)
                    I0_line = linspace(min(I0_scan_mA), max(I0_scan_mA), 60) * 1e-3;
                else
                    I0_line = linspace(max(50, I0_scan_mA(1) - 15), I0_scan_mA(1) + 15, 60) * 1e-3;
                end
                try
                    E0_map    = app.E0Field.Value;
                    ac_map    = app.alpha_cField.Value;
                    se0_map   = app.sigma_e0Field.Value;
                    tau_z_map = app.tau_zField.Value;
                    R_mc_map  = app.R_mcField.Value;
                    Q_mc_map  = app.Q_mcField.Value;
                    df_mc_fix = app.fre_shift_mcField.Value;
                    use_ff_map = app.FormFactorIterCheckbox.Value;

                    Df_opt_kHz = zeros(size(I0_line));
                    Df_105_kHz = zeros(size(I0_line));
                    psi_opt_deg = zeros(size(I0_line));
                    psi_105_deg = zeros(size(I0_line));
                    for idx_fp = 1:length(I0_line)
                        fp_100 = calc_bunch_form_factor_eq('detune_for_xi', I0_line(idx_fp), V_mc, U0_val, ...
                            E0_map, C_val, h_val, ac_map, se0_map, n_hc, R_hc, Q_hc, 1.0, use_ff_map);
                        fp_105 = calc_bunch_form_factor_eq('detune_for_xi', I0_line(idx_fp), V_mc, U0_val, ...
                            E0_map, C_val, h_val, ac_map, se0_map, n_hc, R_hc, Q_hc, 1.05, use_ff_map);
                        Df_opt_kHz(idx_fp)  = fp_100.fre_shift_hc * 1e-3;
                        Df_105_kHz(idx_fp)  = fp_105.fre_shift_hc * 1e-3;
                        psi_opt_deg(idx_fp) = fp_100.psi_hc_deg;
                        psi_105_deg(idx_fp) = fp_105.psi_hc_deg;
                    end

                    if use_ff_map
                        lbl_fp100 = 'Flat-pot. (\xi_{eff}=1.0, F_n)';
                        lbl_fp105 = '\xi_{eff}=1.05 (PTBL, F_n)';
                    else
                        lbl_fp100 = 'Flat-pot. (\xi=1.0)';
                        lbl_fp105 = '\xi=1.05 (PTBL)';
                    end

                    if use_psi
                        plot(app.AxesMap, psi_opt_deg, I0_line * 1e3, 'k--', 'LineWidth', 2.0, 'DisplayName', lbl_fp100);
                        plot(app.AxesMap, psi_105_deg, I0_line * 1e3, 'g:', 'LineWidth', 1.8, 'DisplayName', lbl_fp105);
                    else
                        plot(app.AxesMap, Df_opt_kHz, I0_line * 1e3, 'k--', 'LineWidth', 2.0, 'DisplayName', lbl_fp100);
                        plot(app.AxesMap, Df_105_kHz, I0_line * 1e3, 'g:', 'LineWidth', 1.8, 'DisplayName', lbl_fp105);
                    end

                    % Critical D-mode Robinson Threshold Curve (Omega_i = 0, T. He et al., PRAB 2024)
                    Df_dmode_kHz = zeros(size(I0_line));
                    for idx_dm = 1:length(I0_line)
                        Ib_dm = I0_line(idx_dm);
                        df_opt_dm = Df_opt_kHz(idx_dm) * 1e3;
                        if app.MCAutoDetuneCheckbox.Value
                            df_mc_dm = computeOptimalMCDetune(app, Ib_dm, V_mc, U0_val, R_mc_map, Q_mc_map, R_hc, Q_hc, df_opt_dm, n_hc);
                        else
                            df_mc_dm = df_mc_fix;
                        end
                        g_dfb_dm = getDFBGainForCurrent(app, Ib_dm);
                        [~, ~, thr_hz_dm] = check_dmode_robinson(Ib_dm, R_hc, Q_hc, n_hc, h_val, f0_val, ...
                            U0_val, V_mc, tau_z_map, max(df_opt_dm, 1e3), E0_map, ac_map, R_mc_map, Q_mc_map, df_mc_dm, g_dfb_dm, 0);
                        Df_dmode_kHz(idx_dm) = thr_hz_dm * 1e-3;
                    end
                    if any(Df_dmode_kHz > 0.1)
                        if use_psi
                            psi_dmode = atand(2 * Q_hc * (Df_dmode_kHz * 1e3) / f_hc);
                            plot(app.AxesMap, psi_dmode, I0_line * 1e3, 'r-.', 'LineWidth', 2.0, 'DisplayName', 'D-mode Thresh. (\Omega_i=0)');
                        else
                            plot(app.AxesMap, Df_dmode_kHz, I0_line * 1e3, 'r-.', 'LineWidth', 2.0, 'DisplayName', 'D-mode Thresh. (\Omega_i=0)');
                        end
                    end
                catch ME_opt
                    fprintf(2, 'Notice: could not plot analytical reference curves: %s\n', ME_opt.message);
                end

                % Overlay Fitted Physical Safety Boundary if scan data has boundary transitions
                try
                    crit_I0_a = []; crit_Df_a = [];
                    for iI_a = 1:size(stability_code, 1)
                        row_codes_a = stability_code(iI_a, :);
                        valid_cols_a = find(~isnan(row_codes_a));
                        if length(valid_cols_a) < 2, continue; end
                        row_s = (row_codes_a(valid_cols_a) == 0);
                        % Use rightmost unstable -> stable transition (outer LCBI safety boundary psi_crit)
                        idx_rel_a = find(diff(double(row_s)) == 1, 1, 'last');
                        if ~isempty(idx_rel_a)
                            idx1_a = valid_cols_a(idx_rel_a);
                            idx2_a = valid_cols_a(idx_rel_a + 1);
                            crit_I0_a(end+1) = I0_scan_mA(iI_a); %#ok<AGROW>
                            crit_Df_a(end+1) = 0.5 * (Df_scan_kHz(idx1_a) + Df_scan_kHz(idx2_a)); %#ok<AGROW>
                        end
                    end
                    if length(unique(crit_I0_a)) >= 2
                        p_a = polyfit(crit_I0_a, crit_Df_a, 1);
                        I0_dense = linspace(min(I0_scan_mA), max(I0_scan_mA), 100);
                        Df_dense = polyval(p_a, I0_dense);
                        sign_a = '+'; if p_a(2) < 0, sign_a = '-'; end
                        if use_psi
                            f0_a = 299792458 / C_val;
                            fhc_a = n_hc * h_val * f0_a;
                            psi_dense = atand(2 * Q_hc * (Df_dense * 1e3) / fhc_a);
                            plot(app.AxesMap, psi_dense, I0_dense, '-', 'Color', [0.05 0.30 0.85], 'LineWidth', 2.4, ...
                                'DisplayName', sprintf('Safety Boundary (\\psi_{crit})'));
                        else
                            plot(app.AxesMap, Df_dense, I0_dense, '-', 'Color', [0.05 0.30 0.85], 'LineWidth', 2.4, ...
                                'DisplayName', sprintf('Safety: \\Delta f > %.3f I_0 %s %.1f', p_a(1), sign_a, abs(p_a(2))));
                        end
                    end
                catch ME_sb
                    fprintf(2, 'Notice: could not overlay safety boundary: %s\n', ME_sb.message);
                end

                % 6. Proxy plots for clear Legend using valid character symbols
                plot(app.AxesMap, nan, nan, 'o', 'MarkerFaceColor', [0.18 0.72 0.22], 'MarkerEdgeColor', [0 0.28 0], 'MarkerSize', 7, 'DisplayName', 'Stable (o)');
                plot(app.AxesMap, nan, nan, 'x', 'Color', 'k', 'LineWidth', 2.2, 'MarkerSize', 8, 'DisplayName', 'PTBL (x)');
                plot(app.AxesMap, nan, nan, 'p', 'MarkerFaceColor', [0.10 0.55 0.98], 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'DisplayName', ['Coupled D-Q (' char(9733) ')']);
                plot(app.AxesMap, nan, nan, '^', 'MarkerFaceColor', [0.58 0.12 0.78], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'DisplayName', ['LCBI (' char(9650) ')']);
                plot(app.AxesMap, nan, nan, 'd', 'MarkerFaceColor', [0.45 0.45 0.45], 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'DisplayName', ['Other (' char(9670) ')']);

                legend(app.AxesMap, 'Location', 'northeast', 'FontSize', 10, 'FontName', 'Times New Roman', 'Box', 'on', 'Color', [1 1 1 0.85]);

                % 7. Overlay simulation points
                [X_grid, I0_grid] = meshgrid(x_vals, I0_scan_mA);
                m0 = (stability_code == 0);
                m1 = (stability_code == 1);
                m2 = (stability_code == 2);
                m3 = (stability_code == 3);
                m4 = (stability_code == 4);

                if any(m0(:)), plot(app.AxesMap, X_grid(m0), I0_grid(m0), 'o', 'MarkerFaceColor', [0.18 0.72 0.22], 'MarkerEdgeColor', [0 0.28 0], 'MarkerSize', 6, 'HandleVisibility', 'off'); end
                if any(m1(:))
                    plot(app.AxesMap, X_grid(m1), I0_grid(m1), 'x', 'Color', 'w', 'LineWidth', 3.6, 'MarkerSize', 9, 'HandleVisibility', 'off');
                    plot(app.AxesMap, X_grid(m1), I0_grid(m1), 'x', 'Color', 'k', 'LineWidth', 2.0, 'MarkerSize', 8, 'HandleVisibility', 'off');
                end
                if any(m2(:)), plot(app.AxesMap, X_grid(m2), I0_grid(m2), 'p', 'MarkerFaceColor', [0.10 0.55 0.98], 'MarkerEdgeColor', 'k', 'MarkerSize', 9, 'HandleVisibility', 'off'); end
                if any(m3(:)), plot(app.AxesMap, X_grid(m3), I0_grid(m3), '^', 'MarkerFaceColor', [0.58 0.12 0.78], 'MarkerEdgeColor', 'k', 'MarkerSize', 7, 'HandleVisibility', 'off'); end
                if any(m4(:)), plot(app.AxesMap, X_grid(m4), I0_grid(m4), 'd', 'MarkerFaceColor', [0.50 0.50 0.50], 'MarkerEdgeColor', 'k', 'MarkerSize', 6, 'HandleVisibility', 'off'); end

                % 8. Titles, Parameter Summary Subtitle, and Limits (all strictly double)
                pi_str = 'PI: Off';
                if app.PIEnableCheckbox.Value
                    pi_str = sprintf('PI(K_p=%.2g, K_i=%.0e)', app.KPField.Value, app.KIField.Value);
                end
                dfb_str = 'DFB: Off';
                if app.DFBEnableCheckbox.Value
                    if app.DFB_adaptiveCheckbox.Value
                        dfb_str = 'DFB: Adaptive';
                    elseif app.DFB_auto_scaleCheckbox.Value
                        dfb_str = sprintf('DFB: Auto(%.2g-%.2g)', app.DFB_G1Field.Value, app.DFB_G2Field.Value);
                    else
                        dfb_str = sprintf('DFB(g=%.2g)', app.DFB_gainField.Value);
                    end
                end
                hom_str = 'HOM: Off';
                if ~isempty(app.HOMsTable.Data)
                    if app.HOMAutoAlignCheckbox.Value
                        hom_str = 'HOM: Aligned';
                    else
                        hom_str = 'HOM: Nominal';
                    end
                end
                param_sub = sprintf('V_{mc}=%.2f MV, R_{hc}=%.2f M\\Omega, Q_{hc}=%.0f  |  %s  |  %s  |  %s', ...
                    V_mc * 1e-6, R_hc * 1e-6, Q_hc, pi_str, dfb_str, hom_str);

                if use_psi
                    title(app.AxesMap, {'2D Stability Map vs HC Tuning Angle \psi_{HC}', ['\fontsize{9}\rm' param_sub]});
                    xlabel(app.AxesMap, 'HHC Tuning Angle \psi_{HC}  [^\circ]');
                    psi_lo = min(65.0, floor(min(psi_scan_deg)) - 0.5);
                    psi_hi = min(90.0, max(85.0, ceil(max(psi_scan_deg)*2)/2 + 0.2));
                    xlim(app.AxesMap, [double(psi_lo), double(psi_hi)]);
                else
                    title(app.AxesMap, {'2D Stability Map vs \Delta f_{HC}', ['\fontsize{9}\rm' param_sub]});
                    xlabel(app.AxesMap, 'HHC Detuning \Delta f_{HC}  [kHz]');
                    dDf = 2.5;
                    if length(Df_scan_kHz) > 1, dDf = (Df_scan_kHz(2)-Df_scan_kHz(1))/2; end
                    xlim(app.AxesMap, [double(min(Df_scan_kHz)-dDf), double(max(Df_scan_kHz)+dDf)]);
                end

                ylabel(app.AxesMap, 'Beam Current I_0  [mA]');
                dI0 = 12.5;
                if length(I0_scan_mA) > 1, dI0 = (I0_scan_mA(2)-I0_scan_mA(1))/2; end
                ylim(app.AxesMap, [double(min(I0_scan_mA)-dI0), double(max(I0_scan_mA)+dI0)]);

                grid(app.AxesMap, 'on');
                hold(app.AxesMap, 'off');
                drawnow;

            catch ME
                fprintf(2, 'Error in plotStabilityMap: %s\n', ME.getReport);
                app.StatusLabel.Text = 'Map plot error';
                app.StatusLabel.FontColor = [0.8 0.1 0.1];
            end
        end
    end
end
