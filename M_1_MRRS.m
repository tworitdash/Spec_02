
video_file=[imgDir_video1,video_1,'.avi']; 

writerObj = VideoWriter(video_file); 

open(writerObj); 

for j=1:1 
    
    
    N_Doppler=512; 

     

    start_time=1+N_Doppler*(j-1); 

    x=Data_out(start_time:start_time+N_Doppler-1,:); 

    RD=fftshift(fft(x, N_Doppler),1); 

    frequency=[-500:1000/(N_Doppler+1):500]; % how this has to be changed for diff PRF? 
    hfig=figure;
    imagesc(frequency,range,db(abs(RD'))) 
    colorbarset(gca,'ydir','norm') 


    set(gca,'clim',[10,70]) % If you do not see the range-Doppler plane similar to slide 10,  

    % comment (or edit) the codelineset(gca,'clim',[10,70])
    xlabel('Doppler frequency, ms');
    ylabel('Range, m');
    title(['{',title_str,' 1ms, burst ',num2str(j),'}']);
    
    
    frame = getframe(hfig); 
    
    writeVideo(writerObj,frame); 
    
    close all 
    
end

close(writerObj);