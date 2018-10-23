function beadSetParam(handles)
set(handles.img0Name,'String',handles.pp.fname0);
set(handles.img1Name,'String',handles.pp.fname1);
set(handles.time0,'String',num2str(handles.pp.timeRg0));
set(handles.time1,'String',num2str(handles.pp.timeRg1));
set(handles.timeGap,'String',num2str(handles.pp.timeStep));
set(handles.thr0,'String',num2str(handles.pp.thrDetect));
set(handles.thrMulti,'String',num2str(handles.pp.thrMulti));
set(handles.usePara,'Value',handles.pp.usePara);
set(handles.radRg0,'String',num2str(handles.pp.radRg0));
set(handles.radRg1,'String',num2str(handles.pp.radRg1));
set(handles.channelNum,'String',num2str(handles.pp.nChan));
set(handles.selbgch,'Value',handles.pp.bgch);  % FIXME when items in selbgch is less than pp.bgch
set(handles.selsigch,'Value',handles.pp.sigch);
set(handles.fitModel,'Value',handles.pp.fitModel);
set(handles.smoMethod,'Value',handles.pp.smoMethod);
set(handles.extractMethod,'Value',handles.pp.extractMethod);

if length(handles.pp.outPath)>40
    folder_name = ['..',handles.pp.outPath(end-20:end)];
else
    folder_name = handles.pp.outPath;
end
set(handles.outFolderName,'String',folder_name );

if handles.pp.doneStage>=1
    set(handles.done1,'String', 'Done');
end
if handles.pp.doneStage>=2
    set(handles.done2,'String', 'Done');
end
if handles.pp.doneStage>=3
    set(handles.done3,'String', 'Done');
end

set(handles.myinfo,'String','Session loaded, please continue');

% [X,map] = imread('./cfg/beads.png');
% if ~isempty(map)
%     Im = ind2rgb(X,map);
% end
% imshow(Im);

Im = imread('./cfg/beads.png');
imshow(Im);
[Nx,Ny,~] = size(Im);
text(round(Nx/5),round(Ny/2.5),'BEAD','FontSize',128,'Color',[1 0.5 0])
text(round(Nx/5)+15,round(Ny/1.7),'Version 0.7','FontSize',32,'Color',[1 0.5 0])

end