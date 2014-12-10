function oldSettings = fillPage(h, varargin)
IDX_PAPERTYPE = 1;
IDX_PAPERSIZE = 2;
IDX_PAPERPOSITION = 3;
IDX_PAPERPOSITIONMODE = 4;
if ishandle(h) 
    setprops = {'PaperType', 'PaperSize', 'PaperPosition', 'PaperPositionMode'};
    try 
        settings = get(h, setprops);
    catch 
        error('fillPage unable to get figure properties');
    end
    
    oldSettings = cell2struct(settings, setprops, 2); % to return to caller
    paramSettings = processArgs(varargin{:});
    if isstruct(paramSettings) 
        % process papersize/type
        if isfield(paramSettings, 'papertype')
            if ~isempty(paramSettings.papertype)
                settings{IDX_PAPERTYPE} = paramSettings.papertype;
                % surround with try/catch in case setting is invalid
                try
                   set(h, 'PaperType', paramSettings.papertype);
                catch
                    error('fillPage: invalid PaperType ''%s'' requested ', ...
                        paramSettings.papertype );
                end
                if isfield(paramSettings, 'papersize') && ...
                    ~isempty(paramSettings.papersize) 
                    settings{IDX_PAPERSIZE} = paramSettings.papersize;
                else
                    % since changing papertype might change papersize...
                    % re-read value
                    settings{IDX_PAPERSIZE} = get(h, 'PaperSize'); 
                end
            end
        end
        
        % process margins
        if isfield(paramSettings, 'margins') 
            margins = paramSettings.margins;
            papersize = settings{IDX_PAPERSIZE};
            settings{IDX_PAPERPOSITIONMODE} = 'manual';
            settings{IDX_PAPERPOSITION} = [ ...
                 margins.left ... % X = leftMargin
                 margins.bottom ... % Y = bottomMargin
                 papersize(1) - (margins.right + margins.left) ... % W =  paperwidth - rightMargin - leftMargin
                 papersize(2) - (margins.top + margins.bottom) ... % H = paperheight - topMargin - bottomMargin
                 ];
        end
    end
    set(h, setprops, settings);
else
    error('fillPage requires a handle to a figure');
end
end
function param_settings = processArgs(varargin) 
    margins.left   = .25;
    margins.right  = .25;
    margins.top    = .5;
    margins.bottom = .5;
    papertype = []; 
    papersize = []; 
    for i = 1 : 2 : length(varargin)-1 
        param_arg = varargin{i};
        
        % if caller specified margins to use
        if ischar(param_arg) && strcmpi(param_arg, 'margins')
           val_arg = varargin{i+1};
           if isnumeric(val_arg) && length(val_arg) == 4
               margins.left   = val_arg(1);
               margins.right  = val_arg(2);
               margins.top    = val_arg(3);
               margins.bottom = val_arg(4);
           else 
               warning('fillpage:InvalidMargin', ...
                   'fillpage ignoring invalid margin setting; using default');
           end
        end
        % if caller specified papersize to use
        if ischar(param_arg) && strcmpi(param_arg, 'papersize')
           val_arg = varargin{i+1};
           if ischar(val_arg) 
              papertype = val_arg;
            elseif isnumeric(val_arg) && length(val_arg) == 2
              papertype = '<custom>';
              papersize = val_arg;
           else
               warning('fillpage:InvalidPapersize', ...
                   'fillpage ignoring invalid papersize setting; using figure''s current papersize');
           end
        end
    end
    
    param_settings.margins = margins;
    param_settings.papertype = papertype;
    param_settings.papersize = papersize;
end
