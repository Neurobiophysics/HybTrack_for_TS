% HybTrack: A hybrid single-particle tracking software using manual and automatic detection of dim signals
%
%     Copyright (C) 2017  Byung Hun Lee, bhlee1117@snu.ac.kr
%     Copyright (C) 2017  Hye Yoon Park, hyeyoon.park@snu.ac.kr
%
%     Website : http://neurobiophysics.snu.ac.kr/

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
function [fitted ,temp ,swch,error] = Hybtrack_localization_TXN(im,pts,scan_col,scan_row,ptl_sz,thres,tr,man_scan,swch_input,Th2)
%scan_row는 세로로 얼마나 scan할것인가
%scan_col는 가로로 얼마나 scan할것인가
%%
N=round(sqrt(ptl_sz));
fit_x=2*N+3;
fit_y=2*N+3;
tp=0; wx=0; wy=0;
%%%%%%%%%%%%%%%%%%%%
if tr(1,1)==0 %First frame error prevention
    tr_t=pts;
    tr_t(:,3)=10000;
    tr=reshape(tr_t',1,3*size(pts,1));
end
swch=0;
error.swch=0;
%%
cc=floor(pts(1,1)); rr=floor(pts(1,2));
fitted(1,1)=pts(1,2);
fitted(1,2)=pts(1,1);

scan_size = scan_row*scan_col;
max=0; max2=0; min2=65535;
rr_min = rr-(scan_row-1)/2;
rr_max = rr+(scan_row-1)/2;
cc_min = cc-(scan_col-1)/2;
cc_max = cc+(scan_col-1)/2;

Scanwindow=im(rr_min:rr_max,cc_min:cc_max); %tmp is croped square image

% 1d array sorted by pixel intensities from 2d croped image and inverse error function fitting
Scanwindow_sort = double(sort(reshape(Scanwindow,[],1))); 
    Scanwindow_sort_moved = Scanwindow_sort - Scanwindow_sort((scan_size+1)/2);
    X = linspace(-1,1,scan_size+2);
    X(1)=[];
    X(scan_size+1)=[];
    X = transpose(X);
        ft = fittype('a*erfinv(x)');
        [fitobj,gof,numobs] = fit(X(1:(scan_size+1)/2),Scanwindow_sort_moved(1:(scan_size+1)/2),ft,'StartPoint', 1000);
        error_rms_left = rms(Scanwindow_sort_moved(1:(scan_size-1)/2)-fitobj(X(1:(scan_size-1)/2)));
        error_rms_right = rms(Scanwindow_sort_moved((scan_size+3)/2:scan_size)-fitobj(X((scan_size+3)/2:scan_size)));
        error_ratio = error_rms_right/error_rms_left; %ratio between left and right side of sorted intensity

for j=1:(scan_col*scan_row)  %포인트에서 행 열만큼 스캔
        
    [rr2,cc2]=ind2sub(scan_row,j);
    tmp=im(floor(rr-rr2+1+(scan_row-1)/2)+(N-1)/2-2:floor(rr-rr2+1+(scan_row-1)/2)+(N-1)/2,floor(cc-cc2+1+(scan_col-1)/2)+(N-1)/2-2:floor(cc-cc2+1+(scan_col-1)/2)+(N-1)/2); %tmp is croped square image
m=mean(mean(tmp));

    if m>=max
        max=m;
        temp(1,1)=floor(rr-rr2+1+scan_row/2); temp(1,2)=floor(cc-cc2+1+scan_col/2);
        ttmp(1,1)=floor(rr-rr2+1+scan_row/2); ttmp(1,2)=floor(cc-cc2+1+scan_col/2);
    end
end
    %%
    if isnan(tr(size(tr,1),1))
        
            temp(1,1:2)=NaN;
    else
    %particle detection with error ratio threshold
     if error_ratio < 5
         error.swch=1;
         numbering={'TXN sites'};
                error.message=[char(numbering(1,1)),' cannot be detected.'];
        swch =menu ([error.message,' Do you want to stop at this frame?'],'Stop','Manual detection','No TXN');
        if swch==1
             temp(1,:)=NaN;
        end
        if swch==3
                    temp(1,:)=NaN;    
        end

       if swch==2       
            figure(99)
            hold all;
            colormap('gray')
            subplot(3,1,1:2)
            imagesc(im);     
            axis equal
            hold all;
            plot(tr(size(tr,1),1),tr(size(tr,1),2),'ro')
        temp2 = ginput(1);
        close;
        c=floor(temp2(1,1)); r=floor(temp2(1,2));
        max=0;
        for a=1:(man_scan*man_scan)
            [r2,c2]=ind2sub(man_scan,a);
            m=0;
            for b=1:ptl_sz
                
                [r3,c3]=ind2sub(N,b);
                tmp2(r3,c3)=im(floor(r-r2+1+(man_scan-1)/2)-r3+1+(N-1)/2,floor(c-c2+1+(man_scan-1)/2)-c3+1+(N-1)/2);
                m=mean(mean(tmp2));
            end
            if m>=max
                max=m;
                temp(1,1)=floor(r-r2+1+(man_scan-1)/2);
                temp(1,2)=floor(c-c2+1+(man_scan-1)/2);
                
            end
        end
       end
        end
        
    end
    %if abs(max2/max)<=Th2/100
     if error_ratio > 15
        if isnan(tr(size(tr,1),1))
         numbering={'TXN sites'};
                error.message=[char(numbering(1,1)),' detected.'];
        swch =menu ([error.message,' Do you want to stop at this frame?'],'Stop','Manual detection','No TXN');
        if swch==1
             temp(1,:)=NaN;
        end
        if swch==3
                    temp(1,:)=NaN;    
        end

       if swch==2       
            figure(99)
            hold all;
            colormap('gray')
            subplot(3,1,1:2)
            imagesc(im);    
            axis equal
            hold all;
            plot(tr(size(tr,1),1),tr(size(tr,1),2),'ro')
        temp2 = ginput(1);
        close;
        c=floor(temp2(1,1)); r=floor(temp2(1,2));
        max=0;
        for a=1:(man_scan*man_scan)
            [r2,c2]=ind2sub(man_scan,a);
            m=0;
            for b=1:ptl_sz
                
                [r3,c3]=ind2sub(N,b);
                tmp2(r3,c3)=im(floor(r-r2+1+(man_scan-1)/2)-r3+1+(N-1)/2,floor(c-c2+1+(man_scan-1)/2)-c3+1+(N-1)/2);
                m=mean(mean(tmp2));
            end
            if m>=max
                max=m;
                temp(1,1)=floor(r-r2+1+(man_scan-1)/2);
                temp(1,2)=floor(c-c2+1+(man_scan-1)/2);
                
            end
        end
       end
        end
    end
    close(figure);

    %%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 2D GAUSSIAN FITTING %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    temp_final = zeros(1,2);
    if isnan(temp(1,1)) && swch_input.fit==1
        fitted(1,1:3)=NaN;
    else
    if isnan(temp(1,1))
        temp_final = ttmp;
    else
        temp_final = temp;
    end
    for c=1:fit_x*fit_y
        
        [r4,c4]=ind2sub(fit_y,c);
        %crop(r4,c4)=im(ttmp(1,1)+r4-round(fit_y/2),ttmp(1,2)+c4-round(fit_x/2));
        
        crop(r4,c4)=im(temp_final(1,1)+r4-round(fit_y/2),temp_final(1,2)+c4-round(fit_x/2));
    end
    crop=double(crop);
    
    switch swch_input.fit
        case 1  % 2d Gaussian fitting
            
    [xc , yc, Amp ,wx ,wy] = gaussian_fit_2d(crop);
    if tr(1,3*1)==10000
        tr(size(tr,1),3*1)=Amp*pi*wx*wy;
    end
    background=median(reshape(crop,fit_x*fit_y,1));
       if yc<N-1 || xc<N-1  || xc>fit_x-N+1 || yc>fit_y-N+1 || isnan(xc) || isnan(yc) %fitting error -> centroid
       
        cent_x=0; cent_y=0;             %
        
        for h=1:size(crop)       % Centroid
            for g=1:size(crop)
                cent_x=cent_x + crop(h,g)*h;
                cent_y=cent_y + crop(h,g)*g; 
            end
        end
        xc=cent_x/sum(sum(crop));
        yc=cent_y/sum(sum(crop));  %Centroid fitting end
         if yc<N-1 || xc<N-1 || xc>fit_x-N+1 || yc>fit_y-N+1 || isnan(xc) || isnan(yc) % Centroid error -> Local maxima // rarely happen
        yc=round(fit_y/2);
        xc=round(fit_x/2);
         end  
           fitted(1,3)=sum(sum(im(round(fitted(1,1))-ceil(N/2):round(fitted(1,1))+ceil(N/2),round(fitted(1,2))-ceil(N/2):round(fitted(1,2))+ceil(N/2))))-background*(2*ceil(N/2)+1)^2;
       else 
        fitted(1,3)=Amp*2*pi*wx*wy;
       end
       
    fitted(1,1)=ttmp(1,1)-round(fit_y/2)+yc;
    fitted(1,2)=ttmp(1,2)-round(fit_x/2)+xc;
    
    fitted(1,1)=temp_final(1,1)-round(fit_y/2)+yc;
    fitted(1,2)=temp_final(1,2)-round(fit_x/2)+xc;
    
        case 2 % Centroid Fitting
           
          cent_x=0; cent_y=0;            
        background=min(reshape(crop,fit_x*fit_y,1));
        for h=1:size(crop)       
            for g=1:size(crop)
                cent_x=cent_x + crop(h,g)*h;
                cent_y=cent_y + crop(h,g)*g; 
            end
        end
        xc=cent_x/sum(sum(crop));
        yc=cent_y/sum(sum(crop));
         if yc<N-1 || xc<N-1 || xc>fit_x-N+1 || yc>fit_y-N+1 || isnan(xc) || isnan(yc) 
        yc=round(fit_y/2);
        xc=round(fit_x/2);
         end
         
    fitted(1,1)=temp_final(1,1)-round(fit_y/2)+yc;
    fitted(1,2)=temp_final(1,2)-round(fit_x/2)+xc;

    fitted(1,3)=sum(sum(im(round(fitted(1,1))-ceil(N/2):round(fitted(1,1))+ceil(N/2),round(fitted(1,2))-ceil(N/2):round(fitted(1,2))+ceil(N/2))))-background*(2*ceil(N/2)+1)^2;
    if isnan(temp(1,1))
        fitted(1,1:2)=[NaN NaN];
    end
    end
    end

end

function gcbf=menucallback2(btn, evd, index)                                 %#ok
set(gcbf,swch, index);
end
