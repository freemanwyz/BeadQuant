function controlBeadClick(h,~)
controlCenters = getappdata(h,'test');
axesHandle  = get(h,'Parent');
point = get(axesHandle,'CurrentPoint');
row0 = round(point(1,2));
col0 = round(point(1,1));
switch get(ancestor(h,'figure'),'SelectionType')
    case 'normal' %left click, add circle
        controlCenters(end+1,:) = [row0 col0];
        setappdata(h,'test',controlCenters);
        fprintf(1,'Add X,Y = %.2f,%.2f\n',point(1,1),point(1,2));
        cd1 = get(h,'CData');
        cd1(row0,col0,3) = 1;
        set(h,'CData',cd1);
    case 'alt' %alternate click, remove circle
        nB = size(controlCenters,1);
        for ii=1:nB
            if controlCenters(ii,1)==row0 && controlCenters(ii,2)==col0
                controlCenters(ii,:) = [];
                setappdata(h,'test',controlCenters);
                fprintf(1,'Remove X,Y = %.2f,%.2f\n',point(1,1),point(1,2));
                cd1 = get(h,'CData');
                cd1(row0,col0,3) = 0;
                set(h,'CData',cd1);
                break
            end
        end
end

end