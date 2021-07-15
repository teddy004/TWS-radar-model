function varargout = PPI_Screen(varargin)
% PPI_SCREEN MATLAB code for PPI_Screen.fig
%      PPI_SCREEN, by itself, creates a new PPI_SCREEN or raises the existing
%      singleton*.
%
%      H = PPI_SCREEN returns the handle to a new PPI_SCREEN or the handle to
%      the existing singleton*.
%
%      PPI_SCREEN('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PPI_SCREEN.M with the given input arguments.
%
%      PPI_SCREEN('Property','Value',...) creates a new PPI_SCREEN or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PPI_Screen_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PPI_Screen_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PPI_Screen

% Last Modified by GUIDE v2.5 01-Aug-2019 02:33:41

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PPI_Screen_OpeningFcn, ...
                   'gui_OutputFcn',  @PPI_Screen_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before PPI_Screen is made visible.
function PPI_Screen_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PPI_Screen (see VARARGIN)

% Choose default command line output for PPI_Screen
axes(handles.axes1);
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes PPI_Screen wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = PPI_Screen_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
set(gca,'Color','k');
hold on
teta=0:0.01:2*pi;
r=360;%maximum range in kilometer
m=2;
for i=90:90:r
    m=m+1;
    x=i*cos(teta);
    y=i*sin(teta);
    plot(x,y,'r');
    hold on
end

xlim([-r r])
ylim([-r r])


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%[x,y]=scanning();
%[x,y,z]=Searching(0,0);
% prev=0;
chk = false;
global xp_n
global pminus_n
axes(handles.axes1)
teta=45:0.01:2*pi/4;
r=360;%maximum range in kilometer
for i=90:90:r
    x=i*cos(teta);
    y=i*sin(teta);
    p=plot(handles.axes1,x,y,'k');
    x=get(p);
    hold(handles.axes1,'on')
end
hold(handles.axes1,'on')
grid(handles.axes1,'on')
L1=[];
 for x=0:0.01:254.558
     y=x; 
     L1=[L1;x y];
   
 end
 L2=[];
for x=0:0.01:3.6
    y2=0.9656*x;
    L2=[L2;x y2];
end
L3=[];
for x=0:0.01:3.7
    y2=0.9324*x;
    L3=[L3; x y2];
end
L4=[];
for x=0:0.01:3.7
    y2=0.9003*x;
    L4=[L4;x y2];
end
L5=[];
for x=0:0.01:3.78
    y2=0.8692*x;
    L5=[L5;x y2];
end
L6=[];
for x=0:0.01:3.85
    y2=0.8390*x;
    L6=[L6;x y2];
end
%------------------------------Radar Init-----------------------------
%Simulation parameter
T=10;%duration of simulation
dt=0.1;%update time of simulation



%Radar parameters
Radar_Pos=[0 0 0];
min_power=1.59*10^-19;% 1.59*10^-19 sensitivity of the radar receiver
%Pt=10*10^3;
Pt=str2num(get(handles.Radar_power,'String'))*10^3;
Gt=39;
Gr=39;
lamda=0.1;
RCs=1;%radar cross section considered to calculate maximum range
Rmax=((Pt*Gt*Gr*lamda^2*RCs)/((4*pi)^3*min_power))^(1/4);
%Rmax=240;%maximum range with given power with by taking RCS=1
Ptmin=10*10^3;%minimu power correspondes to minimum range
Rmin=300;%((Ptmin*Gt*Gr*lamda^2*RCs)/((4*pi)^3*min_power))^(1/4);
stR=struct('Radar_Pos',Radar_Pos,'min_power',min_power,'Pt',Pt,'Gt',Gt,'Gr',Gr,'lamda',lamda,'Rmax',Rmax,'Rmin',Rmin);
%Target parameters
target_numbers=2000;
target_pos_x=-Rmin+(Rmax-Rmin)*randn(target_numbers,1);
target_pos_y=-Rmin+(Rmax-Rmin)*randn(target_numbers,1);
target_pos_z=-Rmin+(Rmax-Rmin)*randn(target_numbers,1);
target_velocity=100;
RCS=5*rand(target_numbers,1)+10;%target radar cross section
stT=struct('RCS',RCS,'target_velocity',target_velocity,'target_pos_x',target_pos_x,'target_pos_y',target_pos_y,'target_pos_z',target_pos_z,'target_numbers',target_numbers);
%--------------------------------------------------------------------

P1=plot(handles.axes1,L1(:,1),L1(:,2),'LineWidth',2.5,'Color','b');
hold(handles.axes1,'on')
%P2=plot(L2(:,1),L2(:,2),'LineWidth',2.5,'Color','g');

%P3=plot(L3(:,1),L3(:,2),'LineWidth',2.5,'Color','g');

%P4=plot(L4(:,1),L4(:,2),'LineWidth',2.5,'Color','g');

%P5=plot(L5(:,1),L5(:,2),'LineWidth',2.5,'Color','g');

%P6=plot(L6(:,1),L6(:,2),'LineWidth',2.5,'Color','g');%this make line thicker
xlim(handles.axes1,[-r r])
ylim(handles.axes1,[-r r])
P = [0,0,1];  % Rotation vector
Q = [0,0,-1];
Rp = [0 0]; % Point about which to rotate.

plot(handles.axes1,Rp(1),Rp(2),'.r')

tetaA=45;%starting antenna beam in azimuth
tetaE=30;%antenna elevation angle
test=[];
% kk=[];
% prevX=zeros(6,1);
% prevP=zeros(6,6);
% crentX=zeros(6,1);
% crentP=zeros(6,6);
xpn_c=xp_n;
pminn_c=pminus_n;
id_c=0;
xpn_p=xp_n;
pminn_p=pminus_n;
id_p=0;
sct= struct('xpnext',xpn_c,'pminnext',pminn_c,'id_c',id_c);
spt= struct('xpnext_p',xpn_p,'pminnext_p',pminn_p,'id_p',id_p);
isit=struct([]);
A=([struct('xpnext',xpn_c,'pminnext',pminn_c,'id_c',id_c)]);
% A=([struct('xpnext',0,'pminnext',0,'id_c',0)]);
getidi=0;
% counter=0;
fdggdi=0;
% T = table(id,position_X,Velocity_X,position_Y,Velocity_Y,position_Z,Velocity_Z,Azimuth,Elevation,Range);
while(1) 
     
     rotate(P1,P,5) % rotate 5 degrees at a time    
    
%     rotate(P1,P,0.072)
%     rotate(P3,P,2)
%     rotate(P4,P,2)
%     rotate(P5,P,2)
%     rotate(P6,P,2)
    
    tetaA=tetaA+5;
    pause(0.05)%0.069
    
   % plot(x1,y1,'ro');
   
     if(tetaA>=135)
         while(tetaA>45)
         rotate(P1,P,-5)
         tetaA=tetaA-5;
         pause(0.05)
         end
       % P = [0,0,1];
        %tetaA=0;
      
    end
    Radar_Pos=stR.Radar_Pos;
    min_power=stR.min_power;%sensitivity of the radar receiver
    Pt=stR.Pt;
    Gt=stR.Gt;
    Gr=stR.Gr;
    lamda=stR.lamda;
    Rmax=stR.Rmax;%maximum range with given power with by taking RCS=1
    Rmin=stR.Rmin;
    RCS=stT.RCS;
    target_numbers=stT.target_numbers;
    pos_x=[];
    pos_y=[];
    pos_z=[];
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
        R=sqrt((Radar_Pos(1)-stT.target_pos_x(i))^2+(Radar_Pos(2)-stT.target_pos_y(i))^2+(Radar_Pos(3)-stT.target_pos_z(i))^2);% to make it in kilo meter
        Range=[Range R];
        PR=(Pt*Gt*Gr*lamda^2*rcs)./((4*pi)^3*Range(i)^4);
        power=[power PR];
        Az=atan((stT.target_pos_x(i)-Radar_Pos(1))/(stT.target_pos_y(i)-Radar_Pos(2)))*57.3;
        
        El=atan((stT.target_pos_z(i)-Radar_Pos(3))/sqrt((stT.target_pos_y(i)-Radar_Pos(2)).^2+(stT.target_pos_x(i)-Radar_Pos(1))^.2))*57.3;
        Elevation=[Elevation El];
        if((stT.target_pos_x(i)>=0)&&(stT.target_pos_y(i)>=0))
            Azi=90-Az;
            Azimuth=[Azimuth Azi];
        elseif((stT.target_pos_x(i)<=0)&&(stT.target_pos_y(i)>=0))
            Azi=90+abs(Az);
            Azimuth=[Azimuth Azi];
        elseif((stT.target_pos_x(i)<0)&&(stT.target_pos_y(i)<=0))
            Azi=270-abs(Az);
            Azimuth=[Azimuth Azi];
        elseif((((stT.target_pos_x(i)>=0)&&(stT.target_pos_y(i)<=0))))
            Azi=270+abs(Az);
            Azimuth=[Azimuth Azi];
        end           
        
    end
    %Detection for every target
    pos_x=[];
    pos_y=[];
    for i=1:target_numbers
        if(power(i)>min_power)
            if((-2<tetaA-Azimuth(i))&&(tetaA-Azimuth(i)<2) && (Rmin< Range(i))&& (Range(i)<Rmax) && (-2<tetaE-Elevation(i))&&(tetaE-Elevation(i)<2) )%consider minimum range
                           
                pos_x=[pos_x xt(i)];
                %pos_x= xt(i);
                pos_y=[pos_y yt(i)];
                set(handles.positionin_X,'String',xt(i));
                set(handles.positionin_Y,'String',yt(i));
                set(handles.positionin_Z,'String',zt(i));
                set(handles.Target_Range,'String',Range(i));
                set(handles.Target_Azimuth,'String',Azimuth(i));
                set(handles.Target_Elevation,'String',Elevation(i));
                            
                
                %pos_y= yt(i);
                %pos_z=[pos_z zt(i)]; sins_x=sign(pos_x);
                sins_y=sign(pos_y);
                xd=pos_x./1000;
                yd=pos_y./1000; 
    
                plot(handles.axes1,xd,yd,'y*');   
                
               
                
                if(Range(i)<=360000)
                
                    o=i;
%                     test=[];
                   
                    T=1;
                    %nvar = 0.5;
                    trajectory = [xt(i);yt(i);zt(i)] ;
                    x0 = [xt(1),300,yt(1),300,zt(1),300]';
                    
%                     x0 = zeros(6,1);
                    P0 = diag([.5 .5 .5 0 0 0]);
                    phi = [1 T 0 0 0 0; 0 1 0 0 0 0 ; 0 0 1 T 0 0 ; 0 0 0 1 0 0  ; 0 0 0 0 1 T ; 0 0 0 0 0 1];
                    R = diag([.5^2 .5^2 .5^2]);
                    Q = diag([.035^2 .035^2 .035^2 0 0 0]);
                    Azi=Azimuth(i);
                    El=Elevation(i);
%                     if(prev~=i

                      if(ismember(i,test)==false)
%                         if(chk == false)
                            [filtered, residuals , covariances, kalmgains,xhat,range1,Azimuth1,Elevation1] = kalfilt1(trajectory, x0, P0, phi, R, Q,Azi,El);
%                             prev=i;
%                               prev=[prev i];
%                                 keyy=i
%                                 prev=xp_n
%                                 pp=pminus_n
% sct= struct('xpnext',xpn_c,'pminnext',pminn_c,'id_c',id_c);
% spt= struct('xpnext_p',xpn_p,'pminnext_p',pminn_p,'id_p',id_p);
                                sct.xpn_c=xp_n;
                                sct.pminn_c=pminus_n;
                                sct.id_c=i;
                             
                                A=sct;
%                               prev=[prev pminus_n]
%                                 maybe=cat(2,keyy,xp_n)
                               
%                            
%                             chk=true;
                      else
% %                            if (spt.id_c==i)
%                              if(exist('spt.id_c')==true)
%                            if (i==getfield(A,'id_c'))
                             disp('hello')
                               for ii = 1:length(A)
                                   disp('hello2')
                                   disp(length(A))
                                if (A(ii).id_c==i) 
                                    disp('hello3')
                                   xp_n=A(ii).xpn_c
                                   pminus_n=A(ii).pminn_c
                                end   
                                
                               end
%           
                                
                                [filtered, residuals , covariances, kalmgains,xhat,range1,Azimuth1,Elevation1] = kalfilt1(trajectory, xp_n, pminus_n, phi, R, Q,Azi,El)
                                sct.xpn_c=xp_n;
                                sct.pminn_c=pminus_n;
                                sct.id_c=i;
%                                 fdggdi==getfield(spt,'id_c')
%                             end    
%                             end
%                             keyy=i
%                             crentX=xp_n
%                             crentP=pminus_n
          
                            
                            
                         
%                     end
                      end
              
%                     global data


%                    prevX=crentX;
%                    prevP=crentP;
                   
                   
                  
                   test=[test i]
                   id=int32(o);
                  
                   hope=cat(2,id,filtered');
                   dataa=cat(2,hope,Azi,El,range1);

                   data=get(handles.uitable1,'data');
                   data(end+1,:)=dataa;
              

                   set(handles.uitable1,'data',data );
                  
                   
                   pushbutton2_Callback(hObject, eventdata, handles,filtered,trajectory)
%                    qq={'id','xp_n','pminus_n'}
%                    ss=[i];
%                    pp=[uint8(xp_n)];
%                    ll=[uint8(pminus_n)];
%                    T = table(ss,pp,ll);
                  
%                    spt=sct;
%                    isit=[spt spt]
                   
                   A=[A sct];
                
%                    [ outStruct ] = table2structofarrays( uitable1 )
%                   S = table2struct(data)
%                   S = table2struct(data,'ToScalar',true)
                end
                
             end
        end
%         counter=counter +1
%       spt=sct;
%       isit=[spt spt]
%       T = struct2table(spt)
%       spt=cat(2,spt,sct)
     
    end
      
    G = cell(1,length(A));
for ii = 1:length(A)
   G{ii} = A(ii);
end
G = [G{:}];
spt=G;
% getidi=getfield(spt,'id_c')

end




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles,filtered,trajectory)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes6);
plot(handles.axes6,filtered(1),filtered(3),'r*')

% plot3(filtered(1,loop),filtered(3,loop),filtered(5,loop),'ro');
% text(filtered(1,loop)+2,filtered(3,loop)+2,filtered(5,loop)+2,num2str(loop))
% hold on
% plot3(xhatminus_next(1,loop),xhatminus_next(3,loop),xhatminus_next(5,loop),'b.');
% hold on
plot(handles.axes6,trajectory(1),trajectory(2),'k.');
% text(z(1,loop)+3,z(2,loop)+3,z(3,loop)+3,num2str(loop));
hold(handles.axes6,'on')
grid(handles.axes6,'on')
xlabel(handles.axes6,'x-position');
ylabel(handles.axes6,'y-position');
zlabel(handles.axes6,'z-position');
title(handles.axes6,'Estimated Target 3-D visualization');
xlim(handles.axes1,[-360 360])
ylim(handles.axes1,[-360 360])


% --- Executes during object deletion, before destroying properties.
function axes1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
grid(gca,'Color','g');
set(handles,'Color','g');

function positionin_X_Callback(hObject, eventdata, handles)
% hObject    handle to positionin_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of positionin_X as text
%        str2double(get(hObject,'String')) returns contents of positionin_X as a double


% --- Executes during object creation, after setting all properties.
function positionin_X_CreateFcn(hObject, eventdata, handles)
% hObject    handle to positionin_X (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function positionin_Y_Callback(hObject, eventdata, handles)
% hObject    handle to positionin_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of positionin_Y as text
%        str2double(get(hObject,'String')) returns contents of positionin_Y as a double


% --- Executes during object creation, after setting all properties.
function positionin_Y_CreateFcn(hObject, eventdata, handles)
% hObject    handle to positionin_Y (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function positionin_Z_Callback(hObject, eventdata, handles)
% hObject    handle to positionin_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of positionin_Z as text
%        str2double(get(hObject,'String')) returns contents of positionin_Z as a double


% --- Executes during object creation, after setting all properties.
function positionin_Z_CreateFcn(hObject, eventdata, handles)
% hObject    handle to positionin_Z (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Target_Range_Callback(hObject, eventdata, handles)
% hObject    handle to Target_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Target_Range as text
%        str2double(get(hObject,'String')) returns contents of Target_Range as a double


% --- Executes during object creation, after setting all properties.
function Target_Range_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Target_Range (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Target_Azimuth_Callback(hObject, eventdata, handles)
% hObject    handle to Target_Azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Target_Azimuth as text
%        str2double(get(hObject,'String')) returns contents of Target_Azimuth as a double


% --- Executes during object creation, after setting all properties.
function Target_Azimuth_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Target_Azimuth (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function Target_Elevation_Callback(hObject, eventdata, handles)
% hObject    handle to Target_Elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Target_Elevation as text
%        str2double(get(hObject,'String')) returns contents of Target_Elevation as a double


% --- Executes during object creation, after setting all properties.
function Target_Elevation_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Target_Elevation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
z=get(hObject,'value');
set(handles.Radar_power,'String',z);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function Radar_power_Callback(hObject, eventdata, handles)
% hObject    handle to Radar_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Radar_power as text
%        str2double(get(hObject,'String')) returns contents of Radar_power as a double


% --- Executes during object creation, after setting all properties.
function Radar_power_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Radar_power (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function uitable1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)












% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see MATLAB.UI.CONTROL.TABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
% a={}
% b=[]



function edit52_Callback(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit52 as text
%        str2double(get(hObject,'String')) returns contents of edit52 as a double


% --- Executes during object creation, after setting all properties.
function edit52_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit52 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit53_Callback(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit53 as text
%        str2double(get(hObject,'String')) returns contents of edit53 as a double


% --- Executes during object creation, after setting all properties.
function edit53_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit53 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit54_Callback(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit54 as text
%        str2double(get(hObject,'String')) returns contents of edit54 as a double


% --- Executes during object creation, after setting all properties.
function edit54_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit54 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object deletion, before destroying properties.
function uitable1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
