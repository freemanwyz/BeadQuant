function handles = beadGatherParam(handles)
handles.pp.timeRg0 = str2double(get(handles.time0,'String'));
handles.pp.timeRg1 = str2double(get(handles.time1,'String'));
handles.pp.timeStep = str2double(get(handles.timeGap,'String'));
handles.pp.thrDetect = str2double(get(handles.thr0,'String'));
handles.pp.thrMulti = str2double(get(handles.thrMulti,'String'));
handles.pp.usePara = get(handles.usePara,'Value');
handles.pp.useDff = get(handles.useDff,'Value');
handles.pp.radRg0 = str2double(get(handles.radRg0,'String'));
handles.pp.radRg1 = str2double(get(handles.radRg1,'String'));
handles.pp.bgch = get(handles.selbgch,'Value');
handles.pp.sigch = get(handles.selsigch,'Value');

% contents = get(handles.fitModel,'String'); 
% handles.pp.fitModel = contents{get(handles.fitModel,'Value')};
handles.pp.fitModel = get(handles.fitModel,'Value');
handles.pp.smoMethod = get(handles.smoMethod,'Value');
handles.pp.extractMethod = get(handles.extractMethod,'Value');