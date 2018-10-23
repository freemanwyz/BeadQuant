function Close_CallBack(h,~)
% userData = get(h,'userData'); %Store x,y in axis userData
cc=findall(h,'type','Image');
dd = get(cc,'userData');
getappdata(cc,'test')
assignin('base','dd',dd);
% userData
delete(h)
end