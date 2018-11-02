function fig_map(indx_cluster)
coor = [
    26.403873	-97.727356;
    26.494133	-97.630707;
    26.989964	-97.642059;
    27.127579	-97.864624;
    27.385967	-98.852478;
    27.470905	-98.976196;
    27.954189	-97.402191;
    29.381058	-100.377167;
    30.881897	-101.927979;
    30.926517	-102.861115;
    31.196365	-100.514923;
    31.188366	-102.230103;
    31.539753	-95.628265;
    31.611843	-98.647949;
    31.843525	-100.985474;
    31.99437	-100.118378;
    32.060436	-99.84375;
    32.135536	-100.212769;
    32.111942	-101.391663;
    32.067677	-102.654785;
    32.311188	-98.266968;
    32.279575	-100.243256;
    32.235004	-101.466003;
    32.358738	-100.097534;
    32.285675	-101.836121;
    32.421467	-100.704803;
    32.444733	-100.210541;
    32.446957	-100.555634;
    32.454082	-100.405182;
    32.750698	-99.385345;
    32.710697	-100.746704;
    32.720074	-100.920593;
    32.908356	-101.151581;
    32.99371	-100.549316;
    32.978638	-101.831238;
    33.186367	-98.187744;
    33.199898	-98.384216;
    33.372696	-98.739777;
    33.522728	-98.592712;
    33.712151	-97.371185;
    33.662434	-101.385803;
    34.284252	-99.375488;
    35.349194	-101.390228;
    35.434502	-101.1716;
    35.717995	-100.65036;
    ];
sid = [
    820;
    900;
    1249;
    1342;
    1575;
    1622;
    2061;
    4816;
    8289;
    8476;
    8705;
    8845;
    8979;
    9047;
    9411;
    9572;
    9641;
    9773;
    9842;
    9937;
    10069;
    10148;
    10191;
    10334;
    10365;
    10480;
    10494;
    10526;
    10527;
    10986;
    11008;
    11038;
    11462;
    11614;
    11726;
    11873;
    11913;
    12236;
    12445;
    12708;
    12851;
    14342;
    20481;
    21055;
    22235;
    ];
figure
axesm('mercator','Grid', 'on', 'MapLonLimit',[-107 -93], 'MapLatLimit',[25 37]);
alabamahi = shaperead('usastatehi', 'UseGeoCoords', true,...
    'Selector',{@(name) strcmpi(name,'Texas'), 'Name'});
geoshow(alabamahi, 'FaceColor', [1, 1, 1]);
tightmap;
setm(gca,'mlinelocation',5);
setm(gca,'plinelocation',5);
plabel('PLabelLocation',5);
mlabel('MLabelLocation',5);
setm(gca, 'MlabelParallel', 'south');
framem('FlineWidth',0.5,'FEdgeColor','k');
if nargin == 0
    plotm(coor(:, 1), coor(:, 2), 'or', 'MarkerSize', 3, 'MarkerFaceColor', 'r');
else
    K = unique(indx_cluster);
    listcolor = {'r','b','g','k','y',[.5 .6 .7],[.8 .2 .6]};
    for k = 1: length(K)
        if k<= 7
            tmp_color = listcolor{k};
        else
            tmp_color = rand(1, 3);
        end
        plotm(...
            coor(indx_cluster==K(k), 1), coor(indx_cluster==K(k), 2), 'o',...
            'MarkerSize', 3,...
            'MarkerFaceColor', tmp_color,...
            'MarkerEdgeColor', tmp_color);
    end
end
end

function [fig] = old_GS_style()
% Plot the Gulf Stream style graph, cv is color value
% Returns a MATLAB struct
id = figure();
fig = struct('fig', nan, 'c1', nan, 'h1', nan, 'c2', nan, 'h2', nan,...
    'h_img', nan, 'h_cbar', nan);
fig.id = id;
[Z,R] = arcgridread('TX.asc');
[row col] = size(Z);
x11 = R(3, 1); y11 = R(3, 2);
dx = R(2, 1); dy = R(1, 2);
x = 1: col; y = 1: row;
x = x11 + x.*dx; y = y11 + y.*dy;
% [fig.c1,fig.h1] = contour(x, y, Z, [0 0]);
[fig.c1,fig.h1] = contourf(x, y, Z, 0:500:5000, 'FaceColor', [.8, .8, .8]);
hold on;
set(fig.h1, 'LineWidth', 2, 'Color', 'k');
[fig.c2,fig.h2] = contour(x, y, Z, [-100 -1000 -2000 -3000], 'k');
grid on;
% fig.h_img = imagesc(lon_range, lat_range, cv);
% set(gca,'YDir', 'normal', 'FontSize', 16);
% set(gcf,'Color', 'white');
% set(fig.h_img,'AlphaData',0.9.*~isnan(cv));
% fig.h_cbar = colorbar('FontSize', 16);
% colormap jet;
hold off;
end
