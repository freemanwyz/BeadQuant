function Click_CallBack(h,~,sigma)
userData = get(h,'userData'); %Store x,y in axis userData
userData1 = getappdata(h,'test');
switch get(ancestor(h,'figure'),'SelectionType')
    case 'normal' %left click 
        axesHandle  = get(h,'Parent');
        point = get(axesHandle,'CurrentPoint');
        userData(end+1,:) = [point(1,1) point(1,2)];
        userData1(end+1,:) = [point(1,1) point(1,2)];
        set(h,'userData',userData)
        setappdata(h,'test',userData1)
        fprintf(1,'X,Y = %.2f,%.2f\n',point(1,1),point(1,2));
        cd1 = rand(500,500)*sigma;
        set(h,'CData',cd1);
    case 'alt' %alternate click
        % Reset figure pointer
        %         set(ancestor(a,'figure'), 'Pointer','arrow');
        %Clear button down fcn to prevent errors later
        %         set(get(gca,'Children'),'ButtonDownFcn',[]);
        %Wipe out userData 
        %         set(h,'userData',[]);
        %         x = userData(:,1);
        %         y = userData(:,2);
        %         save('myMatFile', 'x', 'y'); %Save to MAT file ... replace name
end