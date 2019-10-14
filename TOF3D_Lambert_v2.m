%%% 3-D TOF to Mars for Space Elevator Releases
%  James Torla
clear all
clc
tic
%% Read me
%{
Version changes:
    %{
    1)Used Lamber_universal to calculate a lambert's interception of mars
    for every phi for every day.
    2)Increased resolution of true anomaly for finding intercepts.
    3)'break' logic at end of script was improved to allow for longer TOF
      intercepts < 500 days.
    %}

Purpose: To estimate release dates, TOF, and delta V required for
    insertion to a low mars orbit (LMO) by utilizing Richard Battin's
    method for Lambert's Problem

Notes: 
    1) Release angles to reach Mars must have some prograde velocity
       component, i.e. angles from [0,90] and [270,360].
    2) The script is computationally expensive, so lower-end system users
       may want to download the relevant data files for analysis.
    3) For quaternion rotations, the Aerospace Toolbox must be installed to
       utilize the 'quatmultiply' function.

%}
loop = 0;

%% Main script
tic
global mu
mu=1.327E11; % Standard gravitational parameter of the Sun
intdata=[]; % Intercept data [release angle, dates, TOF, delta V]

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
    loop = loop + 1;
    Vr_initial(k,:) = Vs;
end

c = abs(rmars-r);
s = (c+r+rmars)/2;
amin = s/2;
T_O_F_min = sqrt(2)/3*sqrt(s.^3/mu).*(1-((s-c)./s).^(3/2));
T_O_F_min = norm(sqrt(abs(T_O_F_min.^2)))/(60*60*24);

T_O_F_max = 365;            % this is chosen because the calculated tof max is multiple years

TOF = T_O_F_min:.1:T_O_F_max;
for z = 1:length(TOF)
    
    %setting time after TOF
    depdate_2=depdate + TOF(z);
    [year_2,month_2,day_2]=ymd(depdate_2);
    hour=0; minute=0; second=0;
    
    % Mars sfter TOF
    planet_id=4;
    [marscoe, rmars_2, vmars_2, jdmars] = planet_elements_and_sv(planet_id, year_2, month_2, day_2, hour, minute, second);
    
    %lambert with r1,r2,TOF
    if 0 <= incl <= 90
        dm = 'pro';
    elseif 90 < incl <= 180
        dm = 'retro';
    end
    
    t2 = TOF(z);
    [v1(z,:), v2(z,:)] = LAMBERTBATTIN(r, rmars_2, dm, TOF(z));

    for k=1:length(phi)
        DV_1(z,k) = norm(v1(z,:)-Vr_initial(k,:));
        vinf_mars(z,:) = norm(v2(z,:) - vmars_2);
        a_hypo_mars(z,:) = -42830/vinf_mars(z,:).^2;
        DV_2(z,:) = norm(sqrt(2*42830/3389.5-42830/a_hypo_mars(z,:)) - 5.03); 
    end
    
    DV_T(z,:) = DV_1(z,:) + DV_2(z,:);

end

[DV_min(n),idx]=min(min(DV_T(1:z,:)));
[row,col] = find(DV_T==DV_min(n));

depdate_2=depdate + TOF(min(row));
[yearint,monthint,dayint]=ymd(depdate_2);

intdata{n,1}=phi(col); 
intdata{n,2}=datetime([year month day]);
intdata{n,3}=datetime([yearint monthint dayint]);
intdata{n,4}=TOF(min(row));
intdata{n,5}=DV_min(n);
if n >= 803
    break
end
n = n+1;
end
if n >= 803
    break
end
end
if n >= 803
    break
end
end
save Lambert_v2_3.mat intdata
toc

% figure;
% for tt = 1:803
% temp(tt,1) = intdata{tt,5};
% end
% plot(1:803,temp);