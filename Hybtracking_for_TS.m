% HybTrack: A hybrid single-particle tracking software using manual and automatic detection of dim signals
%
%     Copyright (C) 2017  Byung Hun Lee, bhlee1117@snu.ac.kr
%     Copyright (C) 2017  Hye Yoon Park, hyeyoon.park@snu.ac.kr
%
%     Website : http://neurobiophysics.snu.ac.kr/
%
% HybTrack is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     HybTrackis distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with HybTrack.  If not, see <http://www.gnu.org/licenses/>.
%% Load Image and Parameter Setting 
[FileName, PathName] = uigetfile('*.tif', 'Choose the mRNA image(.tif) file');
PathFileName=[PathName,FileName];
TXN_number=4;
%scan_row and col should be odd
Scan_row=15;
Scan_col=15;
window_size=9;
Threshold=99; % Validation
Th2=95; % Detection
SaveAVI=1;
% swch_input.tp=2; % 1 = Linear motion, 2 = Manual
swch_input.fit=2; % 1 = Gaussian fitting, 2 = Centroid fitting
kk=1; deno=PathName(1,1);
while deno~='\' && deno~='/'
    kk=kk+1;
    deno=PathName(1,kk);
end
outputpath=[uigetdir(PathName,'Choose the output path'),deno];
%% Default
Colorset=[1 0 1;0.3 0.2 0.1;0 1 0;0 0 1;1 0 0; 0 1 1; 0.2 0 0.8;0.4 0 0.6; 0.6 0 0.4; 0.8 0 0.2;0 0.2 0.8;0 0.4 0.6;0 0.6 0.4;0 0.8 0.2]; 
maxima_avg=3;
man_scan=15;
imginfo=imfinfo([PathFileName]);
frmmax=numel(imginfo);
frm_i=1; frm_f=frmmax;
%% Tracking Start
%if SaveAVI==1
%immax = averagestack([PathFileName],frmmax);
%end
%% Std Max Projection
if SaveAVI==1
imdevmax = devmaxstack([PathFileName],frmmax);
end
%% Particle selection from maximum projected image
figure(1);
%LBH = imagesc(imread([PathFileName],'tif',1));
CHY = imagesc(imdevmax)
axis equal
for i=1:TXN_number
pts{i} = ginput(1);
end
close;
%%
Name=split(PathFileName,deno);
filename=split(FileName,'.');
 if SaveAVI==1
aviobj = VideoWriter([outputpath,char(filename(1,1)),'.avi']);
 fig=figure(2);
 hold all;
colormap('gray')
open(aviobj);
 end
 tr=zeros(1,6);
 for ii=1:TXN_number

 for i=frm_i:frm_f
     im = imread([PathFileName],'tif',i);
     if i==1
         [fitted, temp1(ii,:,i) , swch,error]= Hybtrack_localization_TS_erf(im,pts{ii},Scan_col,Scan_row,window_size,Threshold,zeros(1,3),man_scan,swch_input,Th2);
     
     else
         i
     [fitted, temp1(ii,:,i) , swch,error]= Hybtrack_localization_TS_erf(im,pts{ii},Scan_col,Scan_row,window_size,Threshold,tr(1:i-1,3*ii-2:3*ii),man_scan,swch_input,Th2);
     end
     
     if error.swch==1
        error.message;
     end
     
     if swch==1
        
        break;
     end
      %for k=1:RNA_number
      tr(i,3*ii-2:3*ii)=[ fitted(1,2) fitted(1,1) fitted(1,3)]; 
    if isnan(temp1(ii,2,i))
    else
        pts{ii}(1,1)=temp1(ii,2,i); pts{ii}(1,2)=temp1(ii,1,i);
        prm_pts(1,1)=floor(fitted(1,1)); prm_pts(1,2)=floor(fitted(1,2));
    end
%end
    if SaveAVI==1   
        outputimage=im;
        subplot(4,2,[1 3]);
        colormap('gray')
        h=imagesc(outputimage);
        axis equal
        hold all;

        plot(tr(i,3*ii-2),tr(i,3*ii-1),'o','Color',[abs(ii/TXN_number),abs(1-ii/TXN_number),0],'markersize',15)
        if isnan(temp1(ii,2,i))
            line([pts{ii}(1,1)-floor(Scan_row/2) pts{ii}(1,1)-floor(Scan_row/2)],[pts{ii}(1,2)-floor(Scan_row/2) pts{ii}(1,2)+floor(Scan_row/2)],'color','y')
            hold all
            line([pts{ii}(1,1)-floor(Scan_row/2) pts{ii}(1,1)+floor(Scan_row/2)],[pts{ii}(1,2)-floor(Scan_row/2) pts{ii}(1,2)-floor(Scan_row/2)],'color','y')
            line([pts{ii}(1,1)+floor(Scan_row/2) pts{ii}(1,1)+floor(Scan_row/2)],[pts{ii}(1,2)-floor(Scan_row/2) pts{ii}(1,2)+floor(Scan_row/2)],'color','y')
         line([pts{ii}(1,1)+floor(Scan_row/2) pts{ii}(1,1)-floor(Scan_row/2)],[pts{ii}(1,2)+floor(Scan_row/2) pts{ii}(1,2)+floor(Scan_row/2)],'color','y')
        end
%title(char(Name(size(Name,1),1)))
        subplot(4,2,[2 4]);
        colormap('gray')
        xlim([0,size(im,2)])
        ylim([0,size(im,1)])
        imagesc(imdevmax);
        axis equal
        hold all;

        plot(tr(1:i,3*ii-2),tr(1:i,3*ii-1),'Color',[abs(ii/TXN_number),abs(1-ii/TXN_number),0]);

        axis off
        writeVideo(aviobj,getframe(fig));
        if mod(i,100)==0
            close(fig)
            fig=figure(2);
            hold all;
            colormap('gray')
            open(aviobj);
        end
    end
    hold all
    subplot(4,2,[5:8]);
    if ii==1
        plot([1:1:i],tr(:,3*ii),'Color',[abs(ii/TXN_number),abs(1-ii/TXN_number),0])
        xlim([1 frm_f])
        ylim([0 max(tr(:,3*ii))])
        hold all
 
    else
        plot([1:1:i],tr(1:i,3*ii),'Color',[abs(ii/TXN_number),abs(1-ii/TXN_number),0])
        hold all
        %plot([1:1:i],tr(1:i,6),'g')
        xlim([1 frm_f])
        ylim([0 max([tr(:,3);tr(:,3*ii)])])
    end
 end

 end
if SaveAVI==1 
    close(fig);
    close(aviobj);
end

%% Plot kymograph
% Hybtrack_kymooverlap_script(PathFileName,tr,'x');
% hold all;
% figure(4)
% imagesc(imaverage);
% colormap('gray');
% axis equal
% hold all
% position=tr;
% for i=1:size(position,2)/3
%     plot(position(:,3*i-2),position(:,3*i-1),'Color',Colorset(i,:));
%     hold all
% end
%% Save
Name=split(FileName,'.');
tmpname=strcat(outputpath,char(Name(1,1)),'.mat'); 
data.tr=tr;
data.temp1=temp1;
save(char(tmpname),'data');
type(char(tmpname));