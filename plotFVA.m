function plotFVA(x,y1,y2)

% plot(x, y1, 'r', 'LineWidth', 2);
% hold on;
% plot(x, y2, 'b', 'LineWidth', 2);
% x2 = [x, fliplr(x)];

% patch([x fliplr(x)], [y1 fliplr(y2)], 'g')

% inBetween = [y1, fliplr(y2)];
% fill(x2, inBetween, 'g');

% y=[0,1;2,2;2,3;3,5];
% xlimit=y(:,1);
% ylimit=y(:,2);
% x=[0,5;1,5;2,5;3,5];
% boundryxlim=x(:,1);
% boundryylim=x(:,2);

ylimit = y1;
boundryylim = y2;
xlimit = x;
boundryxlim = x;

% figure
fill([xlimit; flipud(boundryxlim)], [ylimit; flipud(boundryylim)], 'g')
hold on
plot(xlimit,ylimit, 'b' )
plot(boundryxlim,boundryylim, 'r')
hold off


end