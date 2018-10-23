%%
% install_guilayout();
fh = figure('Units', 'Pixels', ...
   'OuterPosition', [100 100 500 500], ...
   'Toolbar', 'none', 'Menu', 'none');
figSettings = get(fh);
fnames = fieldnames(figSettings);
vals = struct2cell(figSettings);
charIdx = cellfun('isclass', vals, 'char');
d = [fnames(charIdx), vals(charIdx)];

p = uiextras.TabPanel('Parent', fh, ...      % Tab Component
   'Padding', 5 );
axes('Parent', p);                           % 1st Tab
b1 = uiextras.HBox('Parent', p, ...          % 2nd Tab - Horiz Box
   'Spacing', 5);
% load table_data
uitable('Parent', b1, ...                    %   Left Box
   'Data', d, ...
   'RowName', '', ...
   'ColumnName', {'Property', 'Value'}, ...
   'ColumnWidth', {150,100});
b2 = uiextras.VBoxFlex('Parent', b1, ...     %   Right Box - Flex Vert Box
   'Spacing', 5);
uicontrol('Parent', b2, ...                  %     Top Box
   'Style', 'listbox');
uicontrol('Parent', b2, ...                  %     Middle Box
   'Style', 'listbox');
b3 = uiextras.HButtonBox('Parent', b2, ...   %     Bottom Box - Horz Button
   'ButtonSize', [100 25]);
uicontrol('Parent', b3, ...                  %       Left Button
   'String', 'Set');
uicontrol('Parent', b3, ...                  %       Right Button
   'String', 'Clear');

% Set Tab Names
set(p, 'TabNames', {'Main Plot', 'Settings'});

% Set sizes of boxes
set(b1, 'Sizes', [-1 205]);
set(b2, 'Sizes', [-1 -1 25]);