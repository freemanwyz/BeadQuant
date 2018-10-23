function ImageClickCallback ( objectHandle , ~, sigma0, a )
axesHandle  = get(objectHandle,'Parent');
coordinates = get(axesHandle,'CurrentPoint'); 
coordinates = coordinates(1,1:2);
% message     = sprintf('x: %.1f , y: %.1f',coordinates (1) ,coordinates (2));
% helpdlg(message);
% cd0 = get(objectHandle,'CData');
cd1 = rand(500,500)*sigma0;
set(objectHandle,'CData',cd1);
set(a,'userData',cd1);
end