function varargout = lbox2(varargin)
% LBOX2 Application M-file for lbox2.fig
%   LBOX2, by itself, creates a new LBOX2 or raises the existing
%   singleton*.
%
%   H = LBOX2 returns the handle to a new LBOX2 or the handle to
%   the existing singleton*.
%
%   LBOX2('CALLBACK',hObject,eventData,handles,...) calls the local
%   function named CALLBACK in LBOX2.M with the given input arguments.
%
%   LBOX2('Property','Value',...) creates a new LBOX2 or raises the
%   existing singleton*.  Starting from the left, property value pairs are
%   applied to the GUI before lbox2_OpeningFunction gets called.  An
%   unrecognized property name or invalid value makes property application
%   stop.  All inputs are passed to lbox2_OpeningFcn via varargin.
%
%   *See GUI Options - GUI allows only one instance to run (singleton).
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help lbox2

% Last Modified by GUIDE v2.5 02-Oct-2013 08:47:06

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',          mfilename, ...
                   'gui_Singleton',     gui_Singleton, ...
                   'gui_OpeningFcn',    @lbox2_OpeningFcn, ...
                   'gui_OutputFcn',     @lbox2_OutputFcn, ...
                   'gui_LayoutFcn',     [], ...
                   'gui_Callback',      []);
if nargin & isstr(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    varargout{1:nargout} = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before lbox2 is made visible.
function lbox2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to lbox2 (see VARARGIN)

% Choose default command line output for lbox2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

if nargin == 3,
    initial_dir = pwd;
elseif nargin > 4
    if strcmpi(varargin{1},'dir')
        if exist(varargin{2},'dir')
            initial_dir = varargin{2};
        else
            errordlg('Input argument must be a valid directory','Input Argument Error!')
            return
        end
    else
        errordlg('Unrecognized input argument','Input Argument Error!');
        return;
    end
end
% Populate the listbox
load_listbox(initial_dir,handles)
% Return figure handle as first output argument
    
% UIWAIT makes lbox2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = lbox2_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% ------------------------------------------------------------
% Callback for list box - open .fig with guide, otherwise use open
% ------------------------------------------------------------
function varargout = listbox1_Callback(h, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1

% ------------------------------------------------------------
% Read the current directory and sort the names
% ------------------------------------------------------------
%% Load listbox function
function load_listbox(dir_path, handles)
% function adapated from that found in MATLAB\User Guide\Creating Graphical
% User Interfaces\Creating GUIs with GUIDE\Examples of GUIDE GUIs\List Box
% Directory Reader
cd (dir_path) % change directory to selected path
dir_struct = dir(dir_path); 
% get a list of the files/folders in the specified folder
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
% sort the file/folder names
handles.file_names = sorted_names;
% save the file/folder names as handles
guidata(handles.figure1,handles)
% save the handles structure
set(handles.listbox1,'String',handles.file_names,'Value',1)
% display the file/folder names
% setting the "Value" (index to item) to 1 ensures that there are never
% more files/folders than there are file/folder names

listEntries = get(handles.listbox1,'String');

% get a list of all the files/folders listed

match = strfind(listEntries,'files_meanFFT.xls');

% is there a list of selected files?

if sum(cellfun('length',match))~=0

    mac = 1; %using mac?

    if mac

        T=readtable('files_meanFFT.txt');

        T=table2cell(T);

    else

        [N,T]=xlsread('files_meanFFT');

    end

    for t=1:length(T)

        match2 = strfind(listEntries,T{t});

        selected(t) = find(cellfun('length',match2));

    end

    set(handles.listbox1,'Value',selected)

    delete files_meanFFT.xls

end



%% Listbox create function
% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO) eventdata  reserved - to be
% defined in a future version of MATLAB handles    empty - handles not
% created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'),...
        get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in savebutton.
function savebutton_Callback(hObject, eventdata, handles)
% hObject    handle to savebutton (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB

% handles    structure with handles and user data (see GUIDATA)

 

indexSelected = get(handles.listbox1,'Value');

% find the indexes for all the items selected from the list

length = size(indexSelected,2);

% determine how many items were selected

listEntries = get(handles.listbox1,'String');

% get a list of all the files/folders listed

for i=1:length

fileNames(i) = {listEntries{indexSelected(i)}};

end

% create an array of all the filenames selected

 

mac = 1; %if using mac

 

if mac

    writetable(cell2table(fileNames'),'files_meanFFT.txt')

else

    xlswrite('files_meanFFT',fileNames') %save a list of files used

end

delete(handles.figure1) %close the list box figure

