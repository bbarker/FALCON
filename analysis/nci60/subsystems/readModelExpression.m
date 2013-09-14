%remove header of model_expression and this function reads it into a cell
%array
function C = readModelExpression ()
%INPUT
%OUTPUT
%C is a 1x7 containing the columns of model_expression
%Narayanan Sadagopan September 2013

%read file
fileID = fopen('model_expression.csv');
modelExpression = textscan(fileID,'%s %s %s %s %f %f %f', 11879);
modelExpression2 = textscan(fileID,'%s %s %s %s %s %f %f %f', 1480);
modelExpression3 = textscan(fileID,'%s %s %s %s %f %f %f', 14828);
modelExpression4 = textscan(fileID,'%s %s %s %s %s %f %f %f', 1482);
modelExpression5 = textscan(fileID,'%s %s %s %s %f %f %f', 16321);
modelExpression6 = textscan(fileID,'%s %s %s %s %s %f %f %f', 1493);
modelExpression7 = textscan(fileID,'%s %s %s %s %f %f %f', 14826);
modelExpression8 = textscan(fileID,'%s %s %s %s %s %f %f %f', 1483);
modelExpression9 = textscan(fileID,'%s %s %s %s %f %f %f', 19293);
modelExpression10 = textscan(fileID,'%s %s %s %s %s %f %f %f', 1489);
modelExpression11 = textscan(fileID,'%s %s %s %s %f %f %f');
fclose(fileID);
C = modelExpression;

%combine all the segments
count = 1;
for x = 11880:11879+1480
    C{1,1}(x) = modelExpression2{1,1}(count);
    C{1,5}(x) = modelExpression2{1,6}(count);
    C{1,6}(x) = modelExpression2{1,7}(count);
    C{1,7}(x) = modelExpression2{1,8}(count);
    count = count + 1;
end
count = 1;
for x = 11880+1480:11879+1480+14828
    C{1,1}(x) = modelExpression3{1,1}(count);
    C{1,5}(x) = modelExpression3{1,5}(count);
    C{1,6}(x) = modelExpression3{1,6}(count);
    C{1,7}(x) = modelExpression3{1,7}(count);
    count = count + 1;
end
count = 1;
for x = 11880+1480+14828:11879+1480+14828+1482
    C{1,1}(x) = modelExpression4{1,1}(count);
    C{1,5}(x) = modelExpression4{1,6}(count);
    C{1,6}(x) = modelExpression4{1,7}(count);
    C{1,7}(x) = modelExpression4{1,8}(count);
    count = count + 1;
end
count = 1;
for x = 11880+1480+14828+1482:11879+1480+14828+1482+16321
    C{1,1}(x) = modelExpression5{1,1}(count);
    C{1,5}(x) = modelExpression5{1,5}(count);
    C{1,6}(x) = modelExpression5{1,6}(count);
    C{1,7}(x) = modelExpression5{1,7}(count);
    count = count + 1;
end
count = 1;
for x = 11880+1480+14828+1482+16321:11879+1480+14828+1482+16321+1493
    C{1,1}(x) = modelExpression6{1,1}(count);
    C{1,5}(x) = modelExpression6{1,6}(count);
    C{1,6}(x) = modelExpression6{1,7}(count);
    C{1,7}(x) = modelExpression6{1,8}(count);
    count = count + 1;
end
count = 1;
for x = 11880+1480+14828+1482+16321+1493:11879+1480+14828+1482+16321+1493+14826
    C{1,1}(x) = modelExpression7{1,1}(count);
    C{1,5}(x) = modelExpression7{1,5}(count);
    C{1,6}(x) = modelExpression7{1,6}(count);
    C{1,7}(x) = modelExpression7{1,7}(count);
    count = count + 1;
end
count = 1;
for x = 11880+1480+14828+1482+16321+1493+14826:11879+1480+14828+1482+16321+1493+14826+1483
    C{1,1}(x) = modelExpression8{1,1}(count);
    C{1,5}(x) = modelExpression8{1,6}(count);
    C{1,6}(x) = modelExpression8{1,7}(count);
    C{1,7}(x) = modelExpression8{1,8}(count);
    count = count + 1;
end
count = 1;
for x = 11880+1480+14828+1482+16321+1493+14826+1483:11879+1480+14828+1482+16321+1493+14826+1483+19293
    C{1,1}(x) = modelExpression9{1,1}(count);
    C{1,5}(x) = modelExpression9{1,5}(count);
    C{1,6}(x) = modelExpression9{1,6}(count);
    C{1,7}(x) = modelExpression9{1,7}(count);
    count = count + 1;
end
count = 1;
for x = 11880+1480+14828+1482+16321+1493+14826+1483+19293:11879+1480+14828+1482+16321+1493+14826+1483+19293+1489
    C{1,1}(x) = modelExpression10{1,1}(count);
    C{1,5}(x) = modelExpression10{1,6}(count);
    C{1,6}(x) = modelExpression10{1,7}(count);
    C{1,7}(x) = modelExpression10{1,8}(count);
    count = count + 1;
end
count = 1;
for x = 11880+1480+14828+1482+16321+1493+14826+1483+19293+1489:11879+1480+14828+1482+16321+1493+14826+1483+19293+1489+2966
    C{1,1}(x) = modelExpression11{1,1}(count);
    C{1,5}(x) = modelExpression11{1,5}(count);
    C{1,6}(x) = modelExpression11{1,6}(count);
    C{1,7}(x) = modelExpression11{1,7}(count);
    count = count + 1;
end