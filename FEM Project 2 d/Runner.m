clear all
close all

%x = [0 0; 2 0; 1 1];      %%%%%% COORDINATES OF ELEMENT %%%%%%%

x = [0 0; 2 0; 1 2];

GPE = 1;                  %%%%%% GAUSS POINTS PER ELEMENT %%%%%%%

element_type = "D2TR3N";

Stiffness(x,GPE, element_type)

