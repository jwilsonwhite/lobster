%load habitat

load('vector_reef_polygon.txt');
Hab=vector_reef_polygon;
%Nation=textread('vector_nation.txt','%s');
Nation = importdata('nation.csv');
Nation = Nation.data;

