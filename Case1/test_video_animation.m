%-------------------------------------------------------------------------%
%------------                    14th April 2022            --------------%
%------------                   Autopilot design            --------------%
%------------            Straight line path following    -----------------%
%-------------------------------------------------------------------------%

close all;clear all; clc

data = readtable('pva.csv');
inputVid=VideoReader('VJ.mp4');
mergedobj = VideoWriter('compositevid','Motion JPEG AVI');
mergedobj.FrameRate = inputVid.FrameRate;  %match same framerate 
mergedobj.Quality=100;
open(mergedobj); 
%start the stitch
hfig = figure;
t = 0:(1/60):(162/60);
%while loop until there are no more frames
i = 1;
while hasFrame(inputVid)
    %read in frame
    singleFrame = readFrame(inputVid);   
    % display frame
    subplot(3,2,[1 3 5]),imagesc(singleFrame), axis off, axis equal;
    %my gen of dummy data or whatever you want to do
    subplot(3,2,2);
    plot(t(i),data.Position(i),'o','MarkerFaceColor','red');
    hold on;
    plot(t(1:i),data.Position(1:i),'-k');
    axis([0 t(end) 0 max(data.Position)+0.2])
    subplot(3,2,4);
    plot(t(i),data.Velocity(i),'o','MarkerFaceColor','red');
    hold on;
    plot(t(1:i),data.Velocity(1:i),'-k');
    axis([0 t(end) -3 3])
    subplot(3,2,6);
    plot(t(i),data.Accleration(i),'o','MarkerFaceColor','red');
    hold on;
    plot(t(1:i),data.Accleration(1:i),'-k');
    axis([0 t(end) -18 18])
    pause(0.1)
    if i ~= length(t)
        clf
        i= i+1;
    else
        close(mergedobj)
    end   
end