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

% Last Modified by GUIDE v2.5 23-Aug-2017 16:19:02

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
teta=45:0.01:2*pi/4;
r=360;%maximum range in kilometer
for i=90:90:r
    x=i*cos(teta);
    y=i*sin(teta);
    p=plot(x,y,'c');
    x=get(p);
    hold on
end
hold on
grid on
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
Pt=str2num(get(handles.edit8,'String'))*10^3;
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
target_numbers=200;
target_pos_x=-Rmin+(Rmax-Rmin)*randn(target_numbers,1);
target_pos_y=-Rmin+(Rmax-Rmin)*randn(target_numbers,1);
target_pos_z=-Rmin+(Rmax-Rmin)*randn(target_numbers,1);
target_velocity=100;
RCS=5*rand(target_numbers,1)+10;%target radar cross section
stT=struct('RCS',RCS,'target_velocity',target_velocity,'target_pos_x',target_pos_x,'target_pos_y',target_pos_y,'target_pos_z',target_pos_z,'target_numbers',target_numbers);
%--------------------------------------------------------------------

P1=plot(L1(:,1),L1(:,2),'LineWidth',2.5,'Color','b');
hold on
%P2=plot(L2(:,1),L2(:,2),'LineWidth',2.5,'Color','g');

%P3=plot(L3(:,1),L3(:,2),'LineWidth',2.5,'Color','g');

%P4=plot(L4(:,1),L4(:,2),'LineWidth',2.5,'Color','g');

%P5=plot(L5(:,1),L5(:,2),'LineWidth',2.5,'Color','g');

%P6=plot(L6(:,1),L6(:,2),'LineWidth',2.5,'Color','g');%this make line thicker
xlim([-r r])
ylim([-r r])
P = [0,0,1];  % Rotation vector
Q = [0,0,-1];
Rp = [0 0]; % Point about which to rotate.

plot(Rp(1),Rp(2),'.r')

tetaA=45;%starting antenna beam in azimuth
tetaE=30;%antenna elevation angle

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
                %pause(0.05)
%                 tetaA
%                 Azimuth(i)
%                 
%                 tetaA-Azimuth(i)
                
               
%                 Rmin
%                 Rmax
%                 
                Range(i)
                
%                 tetaE
%                 Elevation(i)
%                 
%                 tetaE-Elevation(i)
%                 
                pos_x=[pos_x xt(i)];
                %pos_x= xt(i);
                pos_y=[pos_y yt(i)];
                set(handles.edit1,'String',xt(i));
                set(handles.edit2,'String',yt(i));
                set(handles.edit3,'String',zt(i));
                set(handles.edit4,'String',Range(i));
                set(handles.edit5,'String',Azimuth(i));
                set(handles.edit6,'String',Elevation(i));
                %pos_y= yt(i);
                %pos_z=[pos_z zt(i)]; sins_x=sign(pos_x);
                sins_y=sign(pos_y);
                xd=pos_x./1000;
                yd=pos_y./1000; 
    
                plot(xd,yd,'.y');   
                
                
                if(Range(i)<=360000)
                    
                    
                    figure(1)
                    plot(i,Range(i),'b--o');
                    hold on
                    % end
                    grid on
                    
                    
                    T=1;
                    %nvar = 0.5;
                    trajectory = [xt(i);yt(i);zt(i)] ;
                    x0 = [xt(1),300,yt(1),300,zt(1),300]';
                    
%                     x0 = zeros(6,1);
                    P0 = diag([.5 .5 .5 0 0 0]);
                    phi = [1 T 0 0 0 0; 0 1 0 0 0 0 ; 0 0 1 T 0 0 ; 0 0 0 1 0 0  ; 0 0 0 0 1 T ; 0 0 0 0 0 1];
                    R = diag([.5^2 .5^2 .5^2]);
                    Q = diag([.035^2 .035^2 .035^2 0 0 0]);
                    [filtered, residuals, xhat] = kalfilt1(trajectory, x0, P0, phi, R, Q );
                    %set(handles.oputput,'string',xhat);
                    %xhat;
                    %trajectory;
                    %for i=1:target_numbers
%                     figure(1)
%                     plot(Range(i),'b--o');
%                     hold on
%                     % end
%                     grid on
%                     axis([-0.5 2.5 0 2000]) ;
%                     hold off
%                     figure(2)
%                     plot(filtered,'k*');
                    %                    grid on
                    %                     axis([0 2000 -0.5 2.5]) ;
                    % figure(3)
                    %                     plot(residuals(1,:))
                    %                     grid on
%                     axis([-0.1 0.1 0 2000 ]);
                
                
                end
                
             end
        end
    end
        
end




% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes during object deletion, before destroying properties.
function axes1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
grid(gca,'Color','g');
set(handles,'Color','g');

function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
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
set(handles.edit8,'String',z);

% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
