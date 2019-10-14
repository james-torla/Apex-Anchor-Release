%%% 3-D TOF to Mars for Space Elevator Releases
% James Torla & Nicholas Iannacone
clear
clc

%% Read me
%{
Version changes:
    %{
    1)Utilized MATLAB's datetime function. Fixes month overflow where intercept 
      months would be larger than 12 for long-way intercepts, and months having
      more days than possible.
    2)Increased resolution of true anomaly for finding intercepts.
    3)'break' logic at end of script was improved to allow for longer TOF
      intercepts < 365 days.
    %}

Purpose: To estimate release dates, TOF, and delta V required for
    insertion to a low mars orbit (LMO). Increased complexity from 
    preliminary TOF analysis with planetary positions in 3-D.

Notes: 
    1) Two intercept types can be found; those that require an inclination
       phase change at release (PC), and those that don't (NPC). PC intercepts
       assume orbits are coplanar due to a phase change at release. NPC
       trajectories are rare and may still require a phase change
       depending on mission requirements.
    2) Release angles to reach Mars should have some prograde velocity
       component, i.e. angles from [0,90] and [270,360].
    3) The script is computationally expensive, so lower-end system users
       may want to download the relevant data files for analysis.
    4) For quaternion rotations, the Aerospace Toolbox must be installed to
       utilize the 'quatmultiply' function.
    5) The total frequency of intercepts depends on user constraints for
       intercept. There are more possible intercepts if we increase the 
       allowable distance between Mars and the spacecraft.

Uses: Use this script and relevant data files to analyze trajectories and
    answer relevant research questions including:
        1) What is the minimum TOF for possible and actual intercepts (NPC
           and PC). Compare to NASA's Interplanetary Mission Handbook. Do
           minimum TOF release dates correspond to minimum energy years? Most
           NPC trajectories have a similar TOF to a Hohmann Transfer.
           Minimum TOF for an actual PC intercept I've found is ~75 days.
%}

%% Main script
tic
global mu
mu=1.327E11; % Standard gravitational parameter of the Sun
intdata=[]; % Intercept data [rel angle, dates, rel/int position, orbits, mars positions, TOF, delta V]
inttype='PC'; % PC = Phase change, NPC = No phase change

% Space elevator release variables
phi=[0:359]; % Release angle wrt Earth's velocty, deg.
w_earth=7.2921e-5; % Angular velocity of Earth, rad/s
releaseHeight=106378; % Elevator release height - 106378, 156378 km
Vrmag=w_earth*releaseHeight; % Velocity magnitude at release

% Release date variables
years=2035:2037;
months=1:12;
days=1:31;

n=1; % 
y=1; % Iteration counters
q=1; %
d=1; %

nu=0:.1:360; % True anomaly for plotting
for t=1:length(years)
for u=1:length(months)
for x=1:length(days)
depdate=datetime(years(t),months(u),days(x));
[year,month,day]=ymd(depdate);
hour=0; minute=0; second=0;

% Earth at departure
planet_id=3;
[earthcoe, r, v, jd] = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second);
h_earth=earthcoe(1);
e_earth=earthcoe(2);
RA_earth=earthcoe(3);
i_earth=earthcoe(4);
w_earth=earthcoe(5);
a_earth=earthcoe(7);

% Mars at departure
planet_id=4;
[marscoe, rmars, vmars, jdmars] = planet_elements_and_sv(planet_id, year, month, day, hour, minute, second);
h_mars=marscoe(1);
e_mars=marscoe(2);
RA_mars=marscoe(3);
i_mars=marscoe(4);
w_mars=marscoe(5);
a_mars=marscoe(7);
mSOI=0.576E6;
mu_mars=42820;
rpmars=a_mars*(1-e_mars);
ramars=a_mars*(1+e_mars);

for k=1:length(phi)
% Heliocentric rotations
Rz=[cosd(-phi(k)) -sind(-phi(k)) 0; % R3 rotation matrix
    sind(-phi(k)) cosd(-phi(k)) 0;
    0 0 1];
uv=v./norm(v); % Earth velocity unit vector
uvprime=Rz*uv'; % Unit vector in direction of release, relative to Earth velocty
uvperp=null(uvprime.'); uvperp=uvperp(:,1); % Perpendicular axis for inclination rotation

% Quaternion rotations of 23.5 deg about perpendicular axis
Q1=[0, uvprime'];
Q2=[cosd(-23.5/2), uvperp(1)*sind(-23.5/2), uvperp(2)*sind(-23.5/2), uvperp(3)*sind(-23.5/2)];
Q2star=[cosd(-23.5/2), -uvperp(1)*sind(-23.5/2), -uvperp(2)*sind(-23.5/2), -uvperp(3)*sind(-23.5/2)];

if month<6 || month==12
    Q3=quatmultiply(Q2,Q1); Q3=quatmultiply(Q3,Q2star);
else
    Q3=quatmultiply(Q2star,Q1); Q3=quatmultiply(Q3,Q2);
end
ur=Q3(2:end); % Release velocity unit vector in heliocentric
if ur(3)>=0
    ur(3)=-ur(3);
end
Vs=v+ur*Vrmag; % Total release velocity in heliocentric

% Trajectory orbital elements in heliocentric
coe=coe_from_sv(r,Vs,mu);
h=coe(1); % Angular momentum
e=coe(2); % Eccentricity
RA=coe(3)*180/pi; % Right ascension of the ascending node
incl=coe(4)*180/pi; % Inclination
w=coe(5)*180/pi; % Argument of periapse
a=coe(7); % Semi-major axis
f=coe(6)*180/pi; % True anomaly at release
if rpmars>a*(1+e)
    continue
end

if string(inttype)=='PC'
    dv1=2*norm(Vs)*sind(incl/2); % Delta V required for phase change at release
end

% Orbit propagation
% Transfer orbit
rt=a*(1-e^2)./(1+e.*cosd(nu));
[rt3, vt3]=PQWtoECI(rt,h,e,nu,RA,w,incl);

rrelease=a*(1-e^2)/(1+e*cosd(f));
rrel3d=PQWtoECI(rrelease,h,e,f,RA,w,incl);
rrel3d2=PQWtoECI(a*(1-e^2)/(1+e*cosd(f+10)),h,e,f+10,RA,w,incl);

% Earth's orbit
re=a_earth*(1-e_earth^2)./(1+e_earth.*cosd(nu));
[re3,ve3]=PQWtoECI(re,h_earth,e_earth,nu,RA_earth,w_earth,i_earth);
    
% Mars' orbit
rm=a_mars*(1-e_mars^2)./(1+e_mars.*cosd(nu));
[rm3,vm3]=PQWtoECI(rm,h_mars,e_mars,nu,RA_mars,w_mars,i_mars);

rrelease=a*(1-e^2)/(1+e*cosd(f));
rrel3d=PQWtoECI(rrelease,h,e,f,RA,w,incl);

rrel3d2=PQWtoECI(a*(1-e^2)/(1+e*cosd(f+10)),h,e,f+10,RA,w,incl);

% Find intercepts
for i=1:length(nu)
    foundint=0;
    if (norm(rt3(i,:))<0.8*rpmars || norm(rt3(i,:))>1.2*ramars)
        continue
    end
    
    rdist=rt3(i,:)-rm3;
    if string(inttype)=='NPC'
        rdistnorms=sqrt(rdist(:,1).^2+rdist(:,2).^2+rdist(:,3).^2);
        intindex=find(rdistnorms<=3*mSOI); % Allow for larger distance
    elseif string(inttype)=='PC'
        rdistnorms=sqrt(rdist(:,1).^2+rdist(:,2).^2);
        intindex=find(rdistnorms<=1.25*mSOI); % In 2-D, if an intercept occurs it will be at rSOI of Mars
    end 
        if isempty(intindex)==0        
            rintercept(n,:)=rt3(i,:);
            vintercept(n,:)=vt3(i,:);
            foundint=1;
            if dot(rt3(i,:),vt3(i,:))>=0 % Short way
                TOF(n)=lambert2(rrel3d,rintercept(n,:),a,mu)/84000;
            else % Long way orbit propagation
                ta=real(360-acosd(((a*(1-e^2)/norm(rintercept(n,:)))-1)/e));
                Erelease=2*atan(sqrt((1-e)/(1+e))*tand(f/2));
                if Erelease<0
                    Erelease=Erelease+2*pi;
                end
                Eintercept=2*atan(sqrt((1-e)/(1+e))*tand(ta/2));
                if Eintercept<0
                    Eintercept=Eintercept+2*pi;
                end
                Mrelease=Erelease+e*sin(Erelease);
                if Mrelease<0
                    Mrelease=Mrelease+2*pi;
                end
                Mintercept=Eintercept+e*sin(Eintercept);
                if Mintercept<0
                    Mintercept=Mintercept+2*pi;
                end
                if Mintercept<Mrelease
                     Mintercept=Mintercept+2*pi;
                end
                TOF(n)=(Mintercept-Mrelease)/sqrt(mu/a^3)/84000;
            end
            n=n+1;
            clear intindex
        end
    if foundint==1 % Find position of Mars at intercept
        intdate=depdate+caldays(floor(TOF(n-1)));
        [yearint,monthint,dayint]=ymd(intdate);
        
        [marscoe, rmarsint, vmarsint, jdmarsint] = planet_elements_and_sv(planet_id, yearint, monthint, dayint, hour, minute, second);
        if string(inttype)=='NPC'
            g(q)=norm(rmarsint-rintercept(n-1,:));
            vinf=norm(vmarsint-vintercept(n-1,:)); % This is almost certainly wrong
            amarssoi = -42830/vinf^2;
            dv2 = norm(sqrt(2*42830/3389.5-42830/amarssoi) - 5.03);
            deltav=abs(dv1)+abs(dv2);
        elseif string(inttype)=='PC'
            g(q)=norm(rmarsint(1:2)-rintercept(n-1,1:2));
            vinf=norm(vmarsint(1:2)-vintercept(n-1,1:2));     % Talk about vinf (how big should it be?)           
            amarssoi = -42830/vinf^2;
            dv2 = norm(sqrt(2*42830/3389.5-42830/amarssoi) - 5.03);
            deltav=abs(dv1)+abs(dv2);
        end
        
        if g(q)<=3*mSOI % Check if Mars is close to the spacecraft at intercept
            if TOF(n-1)>365
                break
            end
            intdata{y,1}=phi(k); intdata{y,2}=datetime([year month day]);
            intdata{y,3}=datetime([yearint monthint dayint]);
            launch(y,:) = [convertCharsToStrings(num2str(day)),convertCharsToStrings(num2str(month)),convertCharsToStrings(num2str(year))];
            Time_of_Launch(y,:) = join(launch(y,:),"/");
            intdata{y,6}=[rrel3d; rintercept(n-1,:); rrel3d2];
            intdata{y,7}=[rt3, re3, rm3];
            intdata{y,8}=[rmars; rmarsint];
            intdata{y,4}=[TOF(n-1)];
            intdata{y,5}=[deltav];
            y=y+1;
            q=q+1;
            break
        end
        q=q+1;
    end
end
end
if d >= 803
    break
end
d = d+1;
end
if d >= 803
    break
end
end
if d >= 803
    break
end
end
save Cycle_Heatmap_v3.mat intdata % Shows an example intercept trajectory if one exists
toc                          
%% Graphing

% figure;
% Delta_V = intdata{1,7}';
% Time_of_Flight = round(intdata{1,6}',0);
% tab = table(Time_of_Launch,Time_of_Flight,Delta_V);
% torb = [Time_of_Launch,Time_of_Flight,Delta_V];
% 
% j = heatmap(tab,'Time_of_Launch','Time_of_Flight','colorvariable','Delta_V')
% j.MissingDataColor = [0.8 0.8 0.8];
% j.MissingDataLabel = 'No Data';
% j.XDisplayData = unique(Time_of_Launch,'rows','stable');
% j.YDisplayData = flipud(unique(num2str(Time_of_Flight),'rows'));
% j.Title = 'Delta V for Mars Intercept TOFs and Launch times, July 2035';
% j.GridVisible = 'off'
% 
% figure
% function [h] = plotintercept(intdatain1,intdatain2,intdatain3,intdatain4,intdatain5,intdatain6,intdatain7)
%     phi=intdatain1;
%     rt3=intdatain4(:,1:3);
%     re3=intdatain4(:,4:6);
%     rm3=intdatain4(:,7:9);
%     rrel3d=intdatain3(1,:);
%     rintercept=intdatain3(2,:);
%     rrel3d2=intdatain3(3,:);
%     rmars=intdatain5(1,:);
%     rmarsint=intdatain5(2,:);
%     
%     
%     h=figure;
%     plot3(re3(:,1),re3(:,2),re3(:,3),'b')
%     hold on
%     plot3(rm3(:,1),rm3(:,2),rm3(:,3),'r')
%     plot3(rt3(:,1),rt3(:,2),rt3(:,3),'k')
%     plot3(rrel3d(1),rrel3d(2),rrel3d(3),'m*')
%     plot3(rrel3d2(1),rrel3d2(2),rrel3d2(3),'y*')
%     plot3(rintercept(1),rintercept(2),rintercept(3),'g*')
%     plot3(rmars(1),rmars(2),rmars(3),'c*')
%     plot3(rmarsint(1),rmarsint(2),rmarsint(3),'g*')
%     hold off
%     legend('Earth','Mars','Trajectory','Release','Transfer direction','Intercept','Mars at release')
%     axis equal; grid on
% end
% 
