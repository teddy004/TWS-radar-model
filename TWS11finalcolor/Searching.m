function [pos_x,pos_y,pos_z]=Searching(stR,dt,stT,theta_Azi,theta_Eli)
%sTR=stR(1);
Radar_Pos=stR.Radar_Pos;
min_power=stR.min_power;%sensitivity of the radar receiver
Pt=stR.Pt;
Gt=stR.Gt;
Gr=stR.Gr;

lamda=stR.lamda;
Rmax=stR.Rmax;%maximum range with given power with by taking RCS=1
Rmin=stR.Rmin;
%target_pos_x=stT.target_pos_x;

%target_pos_y=stT.target_pos_y;

%target_pos_z=stT.target_pos_z;
%target_velocity=stT.target_velocity;
RCS=stT.RCS;

target_numbers=stT.target_numbers;
pos_x=[];
pos_y=[];
pos_z=[];

%Generating targets
for i=1:target_numbers
% p{i}={{target_pos_x(i),target_pos_y(i),target_pos_z(i)},target_velocity,RCS(i)};
end
%update target parameters
%Calculating range of the targets
xt=[];
yt=[];
zt=[];
power=[];
Azimuth=[];
Elevation=[];
Range=[];
%for t=0:dt:T
    for i=1:target_numbers
        %k=p{i};
        %x=k{1,1}{1}+dt*target_velocity;
        stT.target_pos_x(i)=stT.target_pos_x(i)+dt*stT.target_velocity;
        xt=[xt stT.target_pos_x(i)];
        %y=k{1,1}{2}+dt*target_velocity;
        stT.target_pos_y(i)=stT.target_pos_y(i)+dt*stT.target_velocity;
        yt=[yt stT.target_pos_y(i)];
        %z=k{1,1}{3};
        stT.target_pos_z(i)=stT.target_pos_z(i);
        zt=[zt stT.target_pos_z(i)];
        %v=k{1,2}+randn;
        v=stT.target_velocity;
        %rcs=k{1,3};
        rcs=RCS(i);
        R=sqrt((Radar_Pos(1)-stT.target_pos_x(i))^2+(Radar_Pos(2)-stT.target_pos_y(i))^2+(Radar_Pos(3)-stT.target_pos_z(i))^2);
        Range=[Range R];
        PR=(Pt*Gt*Gr*lamda^2*rcs)./((4*pi)^3*Range(i)^4);
        power=[power PR];
        Az=atan((stT.target_pos_x(i)-Radar_Pos(1))/(stT.target_pos_y(i)-Radar_Pos(2)))*57.3;
        
        El=atan((stT.target_pos_z(i)-Radar_Pos(3))/sqrt((stT.target_pos_y(i)-Radar_Pos(2)).^2+(stT.target_pos_x(i)-Radar_Pos(1))^.2))*57.3;
        Elevation=[Elevation El];
        if((stT.target_pos_x(i)>=0)&&(stT.target_pos_y(i)>=0))
            Az=90-Az;
            Azimuth=[Azimuth Az];
        elseif(((stT.target_pos_x(i)>=0)&&(stT.target_pos_y(i)<0)))
                Az=90+Az;
                Azimuth=[Azimuth Az];
        elseif((stT.target_pos_x(i)<0)&&(stT.target_pos_y(i)<=0))
            Az=270-Az;
            Azimuth=[Azimuth Az];
        elseif(((stT.target_pos_x(i)<0)&&(stT.target_pos_y(i)>=0)))
                Az=270+Az;
            Azimuth=[Azimuth Az];
        end           
        
    end
    
%Detection for every target
for i=1:target_numbers
    if(power(i)>min_power)   
        if((-2<theta_Azi-Azimuth(i)<2) && (Rmin< Range(i) <Rmax) && (-2<theta_Eli-Elevation(i)<2) )%consider minimum range
        %pause(0.05)
        pos_x=[pos_x xt(i)];
        pos_y=[pos_y yt(i)];
        pos_z=[pos_z zt(i)];
        end
    end
end
%end
