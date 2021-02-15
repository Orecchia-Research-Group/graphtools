function[] = TimeGraph(folderPath, runNumber);
[A] = TimeTestFolder(folderPath, runNumber);
dimA = size(A);
lenA = dimA(1);
B = zeros(lenA/2,4);
for i=1:lenA/2
    B(i, 1) = A(i, 1);
    B(i, 2) = A(i, 3);
    B(i, 3) = A(i+lenA/2, 1);
    B(i, 4) = A(i+lenA/2, 3);
end
B = sort(B, 'ascend');
for i=1:lenA/2
    B(i, 2) = B(i,2)/B(i,1);
    B(i, 3) = B(i,3)/B(i,1);
    B(i, 4) = B(i,4)/B(i,1);
end
B = B(:,2:end)
graph = bar(B);
hold on;
title('Comparative Runtimes for Dynamic Trees');
ylabel('Relative Time to No Matching');
xlabel('Algorithm Used');
set(gca, 'box', 'off');
set(gca,'XTickLabel',{'Karate','Amazon', 'MinedDBLP', 'Youtube'});
for k1 = 1:size(B,2)
    ctr(k1,:) = bsxfun(@plus, graph(1).XData, graph(k1).XOffset');    % Note: ‘XOffset’ Is An Undocumented Feature, This Selects The ‘bar’ Centres
    ydt(k1,:) = graph(k1).YData;                                     % Individual Bar Heights
    text(ctr(k1,:), ydt(k1,:), sprintfc('%.3f', ydt(k1,:)), 'HorizontalAlignment','center', 'VerticalAlignment','bottom', 'FontSize',8, 'Color','black')
end
y=1;
line([0, 5], [y,y], 'Color', 'red', 'LineStyle', '--');
legend('master no matching', 'dynamic trees no matching','dynamic trees matching');
hold off;

