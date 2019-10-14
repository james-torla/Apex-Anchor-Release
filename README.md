# Apex-Anchor-Release
Matlab Files developed to determine Low Delta V and Low Time of Flight interplanetary launch options from a Space elevator Apex Anchor.
There are 4 Matlab files within this code: 
Both the TOF3D_v3_cycle_heatmap.m and the TOF3D_v3_month_heatmap.m contain the Free Release Algorithm for time periods of the full Earth/Mars synodic period and a month long period, respectivly. The Free Release Algorithm calculates launch trajectories with no Delta v required other than transition to the ecliptic plane and orbital injection upon arrival at Mars.
The TOF3D_Lambert_v2.m contains the Lambert's Method Algorithm which utilizes a modified form of Richard Battin's Lambert’s Method to find the absolute Delta V required as a function of TOF and release time. The algorithm searches for the minimum Delta v over all release times for a given day and over all times-of-flight up to 365 days. 
The Space_Elevator_Graphing.m file calls the saved results from all three previous files and produces Matlab table heatmaps.
