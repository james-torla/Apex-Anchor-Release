%Space_Elevator_Graphing
% By James Torla

%% Full Cycle
clear 
clc
load('/Users/james/Documents/MATLAB/Mars Interceptions/MATLAB Script and Functions/Cycle_Heatmap_v3.mat')
tab = cell2table(intdata,'VariableNames',{'phi','Time_of_Launch','Time_of_Intercept','Time_of_Flight','Delta_V','Radiuses','Planetary_Radius','Mars_Radius'});
tab.Time_of_Flight = round(tab.Time_of_Flight,0);
for tt = 1:890
Time_of_Flight(tt,:) = round(intdata{tt,4},0);
Delta_V(tt,:) = intdata{tt,:};
launch(tt,:) = [convertCharsToStrings(num2str(intdata{tt,2}(1,3))),convertCharsToStrings(num2str(intdata{tt,2}(1,2))),convertCharsToStrings(num2str(intdata{tt,2}(1,1)))];
end
Delta_V(82,:) = 0;
Delta_V(83,:) = 0;
Delta_V(84,:) = 0;
Delta_V(85,:) = 0;
Delta_V(87,:) = 0;
Delta_V(88,:) = 0;
Delta_V(89,:) = 0;
Delta_V(90,:) = 0;
Delta_V(91,:) = 0;
Delta_V(92,:) = 0;
Delta_V(93,:) = 0;
Delta_V(94,:) = 0;
Delta_V(95,:) = 0;
Delta_V(96,:) = 0;
Delta_V(97,:) = 0;
Delta_V(98,:) = 0;
Delta_V(99,:) = 0;
%tab.Delta_V = Delta_V;
tab.Time_of_Launch = join(launch(1:890,:),"-");
Time_of_Launch = join(launch(1:890,:),"-");
figure;
j = heatmap(tab,'Time_of_Launch','Time_of_Flight','colorvariable','Delta_V','Colormap',jet)
j.MissingDataColor = [0.8 0.8 0.8];
j.MissingDataLabel = 'No Data';
j.XDisplayData = unique(Time_of_Launch,'rows','stable');
j.YDisplayData = flipud(unique(num2str(Time_of_Flight),'rows'));
j.YDisplayLabels = fliplr(["75","","","","","","","","","","","","","","100","","","","","","","","","","","","","","","125","","","","","","","","","","","","","","","150","","","","","","","","","","","","","","","175","","","","","","","","","","","","","","","200","","","","","","","","","","","","","","","225","","","","","","","","","","","","","","","250","","","","","","","","","","","","","","","275","","","","","","","","","","","","","","","300","","","","","","","","","","","","","","325","","","","","","","","","","","","","","350","","","","","","","","","","","","","","365"]);
j.XDisplayLabels = ["","","","","","","Jan 2035","","","","","","","","","","","","","","","","","","Feb 2035","","","","","","","","","","","","","","","","","","","","","","","","","","Mar 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Apr 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","May 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jun 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jul 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Aug 2035","","","","","","","","","","","","","","","","","Sep 2035","","","","","","","","","Oct 2036","","","","","","","","","","","","Nov 2036","","","","","","","","","","","","","","","Jan 2037","","","","","","","","","","","","","","","Feb 2037"",","","","",""];
j.Title = 'Delta V for Mars Intercept TOFs and Launch times, Jan 2035 to March 2037';
j.GridVisible = 'off'

%% July
clear 
clc
load('/Users/james/Documents/MATLAB/Mars Interceptions/MATLAB Script and Functions/Month_Heatmap_v3.mat')
tab = cell2table(intdata,'VariableNames',{'phi','Time_of_Launch','Time_of_Intercept','Time_of_Flight','Delta_V','Radiuses','Planetary_Radius','Mars_Radius'});
tab.Time_of_Flight = round(tab.Time_of_Flight,0);
for tt = 1:277
Time_of_Flight(tt,:) = round(intdata{tt,4},0);
Delta_V(tt,:) = intdata{tt,5};
Time_of_Launch(tt,:) = datestr(intdata{tt,2});
% launch(tt,:) = [convertCharsToStrings(num2str(intdata{tt,2}(1,3))),convertCharsToStrings(num2str(intdata{tt,2}(1,2))),convertCharsToStrings(num2str(intdata{tt,2}(1,1)))];
end
Delta_V(82,:) = 0;
Delta_V(83,:) = 0;
Delta_V(84,:) = 0;
Delta_V(85,:) = 0;
Delta_V(87,:) = 0;
Delta_V(88,:) = 0;
Delta_V(89,:) = 0;
Delta_V(90,:) = 0;
Delta_V(91,:) = 0;
Delta_V(92,:) = 0;
Delta_V(93,:) = 0;
Delta_V(94,:) = 0;
Delta_V(95,:) = 0;
Delta_V(96,:) = 0;
Delta_V(97,:) = 0;
Delta_V(98,:) = 0;
Delta_V(99,:) = 0;
% tab.Delta_V = Delta_V;
tab.Time_of_Launch = Time_of_Launch;
% tab.Time_of_Launch = join(launch(1:277,:),"-");
% Time_of_Launch = join(launch(1:277,:),"-");
figure;
j = heatmap(tab,'Time_of_Launch','Time_of_Flight','colorvariable','Delta_V','Colormap',jet)
j.MissingDataColor = [0.8 0.8 0.8];
j.MissingDataLabel = 'No Data';
j.XDisplayData = unique(Time_of_Launch,'rows','stable');
j.YDisplayData = flipud(unique(num2str(Time_of_Flight),'rows'));
j.YDisplayLabels = fliplr(["","75","","","","","","100","","","","","","125","","","","","","150","","","","","","175","","","","","","200","","","","","","225","","","","","","250","","","","","","","275","","","","","","300","","","","","","325","","","","","","350","","","","365"]);
j.XDisplayLabels = datestr(datetime(2035,7,1,'InputFormat','dd/mm/yyyy'):datetime(2035,7,31,'InputFormat','dd/mm/yyyy'))
j.Title = 'Delta V for Mars Intercept TOFs and Launch times, Jul 2035';
j.GridVisible = 'off'



%% Lambert
clear
clc
load('/Users/james/Documents/MATLAB/Mars Interceptions/MATLAB Script and Functions/Lambert_v2_2.mat')
tab = cell2table(intdata,'VariableNames',{'phi','Time_of_Launch','Time_of_Intercept','Time_of_Flight','Delta_V'});
tab.Time_of_Flight = round(tab.Time_of_Flight,0);
Time_of_Flight = tab.Time_of_Flight;
for tt = 1:803 
Time_of_Launch(tt,:) = datestr(intdata{tt,2});
if tt > 1
    if Time_of_Launch(tt,:) == Time_of_Launch(tt-1,:)
        Delta_V(tt,:) = 0;
    end
end
Delta_V(tt,:) = intdata{tt,5};
end
tab.Delta_V = Delta_V;
test2 = 4;
tab.Time_of_Launch = Time_of_Launch;
% Time_of_Launch = join(launch(1:803,:),"-");
figure;
mycolor(1,:) = [.8,.8,.8];
mycolor(2:65,:) = jet;
j = heatmap(tab,'Time_of_Launch','Time_of_Flight','colorvariable','Delta_V','Colormap',mycolor)
j.ColorMethod = 'sum';
test2 = 1;
j.MissingDataColor = [0.8 0.8 0.8];
test2 = 3;
j.MissingDataLabel = 'No Data';
test2 = 2;
j.XDisplayData = unique(Time_of_Launch,'rows','stable');
j.XDisplayLabels = ["","","","","","","","","","","","","","","Jan 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Feb 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Mar 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Apr 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","May 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jun 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jul 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Aug 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Sep 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Oct 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Nov 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Dec 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jan 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Feb 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Mar 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Apr 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","May 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jun 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jul 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Aug 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Sep 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Oct 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Nov 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Dec 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jan 2037","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Feb 2037","","",""];
j.YDisplayLabels = (["","60","","","","","","","","","","","","","","75","","","","","","","","","","","","","","","100","","","","","","","","","","","","","","","","125","","","","","","","","","","","","","","","","150","","","","","","","","","","","","","","","","","175","","","","","","","","","","","","","","","","200","","","","","","","","","","","","","","","","225","","","","","","","","","","","","","","","","250","","","","","","","","","","","","","","","","275","","","","","","","","","","","","","","","300","","","","","","","","","","","","","","325","","","","","","","","","","","","","","350","","","","","","","","","","","","","","365"]);
j.YDisplayData = flipud(unique(num2str(Time_of_Flight),'rows'));
j.Title = 'Delta V for Mars Intercept TOFs and Launch times, Jan 2035 to March 2037, Lamberts Method';
j.GridVisible = 'off'


%% Total

clear
clc
load('/Users/james/Documents/MATLAB/Mars Interceptions/MATLAB Script and Functions/Cycle_Heatmap_v3.mat')
table = cell2table(intdata,'VariableNames',{'phi','Time_of_Launch','Time_of_Intercept','Time_of_Flight','Delta_V','Radiuses','Planetary_Radius','Mars_Radius'});
table.Time_of_Flight = round(table.Time_of_Flight,0);
for tt = 1:890
Time_of_Flight(tt,:) = round(intdata{tt,4},0);
Delta_V(tt,:) = intdata{tt,:};
launch(tt,:) = [convertCharsToStrings(num2str(intdata{tt,2}(1,3))),convertCharsToStrings(num2str(intdata{tt,2}(1,2))),convertCharsToStrings(num2str(intdata{tt,2}(1,1)))];
end
table.Time_of_Launch = datestr(join(launch(1:890,:),"-"));
Time_of_Launch_1 = datestr(join(launch(1:890,:),"-"));
Delta_V_1 = table.Delta_V;
min = min(Delta_V);
[row,col] = find(Delta_V_1==min);
phi_V_1 = table.phi;
Time_of_Intercept_V_1 = table.Time_of_Intercept;
Time_of_Flight_1 = table.Time_of_Flight;

% for y = 1:890
%     t{y,1} = table.phi(y,:);
%     t{y,2} = table.Time_of_Launch(y,:);
%     t{y,3} = table.Time_of_Intercept(y,:);
%     t{y,4} = table.Time_of_Flight(y,:);
%     t{y,5} = table.Delta_V(y,:);
% end
% table = cell2table(t,'VariableNames',{'phi','Time_of_Launch','Time_of_Intercept','Time_of_Flight','Delta_V'});

load('/Users/james/Documents/MATLAB/Mars Interceptions/MATLAB Script and Functions/Lambert_v2_2.mat')
tab = cell2table(intdata,'VariableNames',{'phi','Time_of_Launch','Time_of_Intercept','Time_of_Flight','Delta_V'});
tab.Time_of_Flight = round(tab.Time_of_Flight,0);
Time_of_Flight = tab.Time_of_Flight;
Delta_V = tab.Delta_V;
for tt = 1:803 
Time_of_Launch(tt,:) = datestr(intdata{tt,2});
if tt > 1
    if Time_of_Launch(tt,:) == Time_of_Launch(tt-1,:)
        Delta_V(tt-1,:) = 0;
    end
end
end
tab.Time_of_Launch = Time_of_Launch;
tab.Delta_V = Delta_V;
Delta_V_2 = tab.Delta_V;
phi_V_2 = tab.phi;

Time_of_Intercept_V_2 = tab.Time_of_Intercept;
Time_of_Flight_2 = tab.Time_of_Flight;
% for y = 1:803
%     b(y,:) = tab.phi(y,:);
%     b(y,:) = tab.Time_of_Launch(y,:);
%     b(y,:) = tab.Time_of_Intercept(y,:);
%     b(y,:) = tab.Time_of_Flight(y,:);
%     b(y,:) = tab.Delta_V(y,:);
% end

p = [(1:803)';phi_V_1];
p2 = [Time_of_Launch;Time_of_Launch_1];
p3 = [Time_of_Flight_2;Time_of_Flight_1];
p4 = [Delta_V_2;Delta_V_1];

for n = 1:length(p2)
    if n <= 803
        h{n,1} = phi_V_2;
    else
        h{n,1} = p(n,:);
    end
    h{n,2} = p2(n,:);
    h{n,3} = p3(n,:);
    h{n,4} = p4(n,:);
end


tab = cell2table(h,'VariableNames',{'phi','Time_of_Launch','Time_of_Flight','Delta_V'});

Time_of_Launch = tab.Time_of_Launch;
Time_of_Flight = tab.Time_of_Flight;

% Time_of_Launch = join(launch(1:803,:),"-");
% total = [table,tab];
figure;
mycolor(1,:) = [.8,.8,.8];
mycolor(2:65,:) = jet;
j = heatmap(tab,'Time_of_Launch','Time_of_Flight','colorvariable','Delta_V','Colormap',mycolor)
j.ColorMethod = 'mean';
j.GridVisible = 'off'
j.MissingDataColor = [0.8 0.8 0.8];
j.MissingDataLabel = 'No Data';
j.XDisplayData = unique(Time_of_Launch,'rows','stable');
j.XDisplayLabels = ["","","","","","","","","","","","","","","","","Jan 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Feb 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Mar 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Apr 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","May 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jun 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jul 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Aug 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Sep 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Oct 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Nov 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Dec 2035","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jan 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Feb 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Mar 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Apr 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","May 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jun 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jul 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Aug 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Sep 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Oct 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Nov 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Dec 2036","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Jan 2037","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","","Feb 2037","","","","","",""];
j.YDisplayLabels = (["","","60","","","","","","","","","","","","","","","","","75","","","","","","","","","","","","","","","","","","100","","","","","","","","","","","","","","","","","","125","","","","","","","","","","","","","","","","","","","150","","","","","","","","","","","","","","","","","","","","175","","","","","","","","","","","","","","","","","","","200","","","","","","","","","","","","","","","","","","","","225","","","","","","","","","","","","","","","","","","","","250","","","","","","","","","","","","","","","","","","","","275","","","","","","","","","","","","","","","","","","","300","","","","","","","","","","","","","","","","","","325","","","","","","","","","","","","","","350","","","","","","","","","","","","","","","","","","","","","365","",""]);
j.YDisplayData = flipud(unique(num2str(Time_of_Flight),'rows'));
j.Title = 'Delta V for Mars Intercept TOFs and Launch times, Jan 2035 to March 2037, Lambets and Free Release';


%%

figure;
for tt = 1:1693
temp(tt,1) = intdata{tt,5};
end

plot(1:803,temp);
%%
figure;
scatter(tab.Time_of_Launch,tab.Delta_V);

