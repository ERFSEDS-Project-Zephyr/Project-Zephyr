%{

Name: Luke Farmer
Email: farmerl2@my.erau.edu
Date: 1/29/2023
Zephyr
Program Description: Fin flutter calculator
%}

clc
clear all

%INPUTS
%All variables have apogee representation listed after, also a test value
%of kleos's fins

thickness = input('Input fin thickness (in) '); %t thickness = 0.125 ;Zephyr = 

rootChord = input('Input root chord (in) '); %cr  7.618;

tipChord = input('Input tip chord (in) '); %ct  2.738;

semiSpan = input('Input semi-span (in) '); %b  5.248;

shearModulus = input('Input shear modulus of fin material (psi) '); %G  580151;

altitude = input('Input the altitude at the rocket''s max velocity (ft) '); %h  850;

%Calculations

%Fin area
S = 0.5*(rootChord + tipChord) * semiSpan; 

%Aspect Ratio
AR = semiSpan^2 / S;

%Taper Ratio (lambda)
taperRatio = tipChord / rootChord;

%Temperature (F)
temp = 59 - (0.00356 * altitude);

%Pressure (lbs/ft^2)
pressure = (2116 / 144) * ((temp + 459.7) / 518.6)^5.265;

%Speed of Sound
sound = (1.4 * 1716.59 * (temp + 460))^(1/2);

%Flutter Speed (ft/s)
middle = 1.337 * AR^3 * pressure * (taperRatio + 1);
bottom = 2 * (AR + 2) * (thickness / rootChord)^3;
bot2 = middle / bottom; %wacky order of operations needs to divide middle and bottom terms before shear mod is divided
flutterSpeed =sound * sqrt(shearModulus / bot2 );

fprintf('%f', flutterSpeed);
%Should equal 842.44 for kleos's fins





