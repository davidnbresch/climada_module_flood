function ls_test2()
%function to test stuff... can be deleted afterwards

[X,Y,Z] = peaks(50);
figure('Position',[280 400 1200 450])
% Original surface with too many edges
subplot(1,2,1)
surf(X,Y,Z,'FaceColor','interp');
xlabel('X')
ylabel('Y')
zlabel('Z')
% Compare to:
subplot(1,2,2)
s = surf(X,Y,Z,'FaceColor','interp','EdgeColor','none');
xlabel('X')
ylabel('Y')
zlabel('Z')
%% Extract X,Y and Z data from surface plot
x=s.XData;
y=s.YData;
z=s.ZData;
% For R2014a and earlier:
% x=get(s,'XData');
% y=get(s,'YData');
% z=get(s,'ZData');
%% Create vectors out of surface's XData and YData
x=x(1,:);
y=y(:,1);
%% Divide the lengths by the number of lines needed
xnumlines = 10; % 10 lines
ynumlines = 10; % 10 partitions
xspacing = round(length(x)/xnumlines);
yspacing = round(length(y)/ynumlines);
%% Plot the mesh lines 
% Plotting lines in the X-Z plane
hold on
for i = 1:yspacing:length(y)
    Y1 = y(i)*ones(size(x)); % a constant vector
    Z1 = z(i,:);
    plot3(x,Y1,Z1,'-k');
end
% Plotting lines in the Y-Z plane
for i = 1:xspacing:length(x)
    X2 = x(i)*ones(size(y)); % a constant vector
    Z2 = z(:,i);
    plot3(X2,y,Z2,'-k');
end
hold off

disp('hier')


end