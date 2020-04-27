clc
clear all
close all
% assigning x as a symbol
syms x;
%function to maximize
y= 2*sin(x)-0.1*x^2;
%number of loop iterations
NoOfIterations=3;
%starting initial guess
InitialGuess= 2.5;
X=InitialGuess;

%Newton method
fprintf('Newton Method\n')
%iterator on number of iterations
for i=1:NoOfIterations
    %vpa (variable percision floating-point arithmetic
    %subs (symbolic substitution)
    %diff (calculates differences between adjacent elements)
    %x(i+1)=x(i)-first diff/second diff
    %so we subs the first diff of y over the 2nd diff of y
    X= X- vpa(subs(diff(y)/diff(y,2),x,X));
    fprintf('i= %d\tx=%.5f\ty=%.5f\n',i,X,vpa(subs(y,x,X)));
end  

%Golden Method
fprintf('\n Golden Ratio Method\n')
%specified interval in the problem
%need to detect single peak
%we need to compare values of the function at 3 points
StartOfInterval=0;
EndOfInterval=4;
%calculated in lecture
GoldenRatio=0.618;
%no of iterations
NoOfIterations=8;
xl=StartOfInterval;
xu=EndOfInterval;

%coordinates of the two interior points
x1=xl+GoldenRatio*(xu-xl);
x2=xu-GoldenRatio*(xu-xl);
%table printing initialization
fprintf('i\txl\t\tx2\t\tx1\t\txu\t\ty2\t\ty1\t\t\t\t\t   Error Bound\n');
for i=1:NoOfIterations
    %same as newton in terms of functions described above
    y1= vpa(subs(y,x,x1));
    y2= vpa(subs(y,x,x2));
    %use 4 decimal places in all your calculations
    fprintf('%d\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t\t   %.4f\n',i,xl,x2,x1,xu,y2,y1,(xu-xl),(x2-xl));
    %if fx1>fx2 then xl and x2 can be eliminated because they don't contain
    %maximum
    if y1>y2
        xl=x2;
        x2=x1;
        x1=xl+GoldenRatio*(xu-xl);
    %if not then we eliminate x1,xu in that case
    else
        xu=x1;
        x1=x2;
        x2=xu-GoldenRatio*(xu - xl);
    end
end    