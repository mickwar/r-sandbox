options ls=120;
options FORMDLIM="#";
title 'MANOVA for Sowing-Date and Variety';

data snap;
    infile 'c:\Users\mickfwop\Desktop\SNAPBEAN.DAT';
    input sowdate variety replicate x1 x2 x3 x4;
run;

proc glm data = snap;
    class sowdate variety;
    model x1 x2 x3 x4 = sowdate | variety / nouni;
    contrast 'contrast on variety' variety 1 -2 1;
    contrast 'linear on sowdate' sowdate -3 -1 1 3;
    contrast 'quadratic on sowdate' sowdate 1 -1 -1 1;
    contrast 'cubic on sowdate' sowdate -1 3 -3 1;
    manova h = sowdate | variety / printe printh;
    means variety / bon alpha = 0.0125;
run;

proc glm data = snap;
    class sowdate variety;
    model x1 x2 = sowdate | variety / nouni;
    manova h = sowdate | variety / printe printh;
run;

title  'Repeated measures for Mandible';
data mandible;
    infile 'c:\Users\mickfwop\Desktop\MANDIBLE.DAT';
    input subject group a1t1 a1t2 a1t3 a2t1 a2t2 a2t3 a3t1 a3t2 a3t3 ;
run;

proc glm data = mandible;
    class group;
    model a1t1 a1t2 a1t3 a2t1 a2t2 a2t3 a3t1 a3t2 a3t3 = group / nouni;
    repeated activator 3 contrast(1), time 3 polynomial / printm printe summary;
    manova h = group / printe printh;
run;
