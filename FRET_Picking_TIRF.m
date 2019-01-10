function Picked=FRET_Picking(M,N,CUT)

% This script will automatically pick out data sets that contain FRET transitions. 
% 
% 
% M - Size in frames of initial smoothing window - 5 is good - less for
%   shorter frame rate if acting up. More if very noisy
% N - Size in frames of the second window used to measure fret transition
%           20 is good - less for faster transitions. More if noisy 
% CUT - Cutoff value used to determine if transition has occured - 200 is
%           good - decreasing will increase sensitivity, increasing will
%           decrease sensitivity




% 
% I - First channel selected - 1 for Cy3, 2 for Cy3.5
% J - Second channel selected - 3 for Cy5, 4 for Cy5.5

I=input('Select first channel: ');
J=input('Select second channel: ');

M = M;
N = N;
CUT = CUT;

% % File Selection
% Looks for .mat files that start with "Run" and end with an I.mat where I
% is an integer

% % % Dir = dir;
% % % DirName = {Dir.name};
% % % m=1;
% % % for i=1:size(DirName,2)
% % %     if (size(DirName{i},2)>4)
% % %         if (strcmp(DirName{i}(1:3),'Run')==1 && (size(str2num(DirName{i}(1,size(DirName{i},2)-4)),1)==1))
% % %             C{m}=DirName{i};
% % %             m=m+1;
% % %         end
% % %     end
% % % end
% % % clear Dir DirName

% END File Selection
% If automatic file selection is not working, or you only want to select a
% subset of the data sets to analyze, comment out the above File Selection
% code and hardcode names of data sets into C below
%
%C={'Run_001.mat','Run_002.mat','Run_003.mat'};
C={'Run_001.mat'};

% This is the portion of the script that does all the heavy lifting
% Input file read and channels of interest are grabbed from ttotal in
%       FRET_GrabChannels
% FRET_Pick calls the FRET_Measure script. FRET_Measure is what is
%   computing the actual measure. FRET_Pick then looks at this measure to
%   determine which traces meet the specified criteria and should be picked.
% The remainder of the code is doing the actual picking of the traces from
%   the initial ttotal 

Folder = strcat('PickedTraces_Ch',num2str(I),'_Ch',num2str(J),'_',num2str(M),'_',num2str(N),'_',num2str(CUT));
mkdir(Folder)
for i=1:size(C,2)
    load(C{i});
    Y=FRET_GrabChannels_TIRF(ttotal,I,J);
    fprintf(char(strcat('Input File:',{' '}, C{i}(1:size(C{i},2)-4),'\r')));
    Picked{:,i}=FRET_Pick(Y,M,N,CUT);
    ttotal_picked=ttotal(:,1);
    ttotal_unfilt_picked=ttotal(:,1);
    m=1;
    ttotal_xy_picked = [];
    if (size(Picked{:,i}>0))
        for j=Picked{:,i}'
            ttotal_unfilt_picked(:,4*m-2:4*m+1)=ttotal(:,4*j-2:4*j+1);
            ttotal_picked(:,4*m-2:4*m+1)=ttotal(:,4*j-2:4*j+1);
            ttotal_xy_picked(m,:)=ttotal_xy(j,:);
            m=m+1;
        end
    end
    seq=Picked{:,i}';
    keep Picked C i FrameRate holeXY seq ttotal_unfilt_picked ttotal_picked ttotal_xy_picked CUT M N I J Folder;
    ttotal=ttotal_picked;
    ttotal_unfilt=ttotal_unfilt_picked;
    ttotal_xy=ttotal_xy_picked;
    clear ttotal_picked ttotal_xy_picked ttotal_unfilt_picked;
    clear ii;
    save(strcat(Folder,'/',C{i}(1:end-4),'-p','.mat'));
end
