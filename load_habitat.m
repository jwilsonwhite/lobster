%load habitat

load('vector_reef_polygon_no_grounds.txt');
Hab=vector_reef_polygon_no_grounds;
%Nation = importdata('nation_no_grounds.csv');
Nation = importdata('nation_ecoregions_no_grounds.csv');
Nation = Nation.data;

