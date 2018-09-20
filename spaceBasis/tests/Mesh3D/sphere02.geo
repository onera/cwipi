DefineConstant[ R1 = {1, Min 0.5, Max 2.0, Step 1, Name "Parameters/R1"} ];
DefineConstant[ R2 = {2, Min 2.0, Max 3.0, Step 1, Name "Parameters/R2"} ];
DefineConstant[  h = {1, Min 1.0, Max 5.0, Step 1, Name "Parameters/h"} ];

Printf("Paramètres du maillage");

Printf("R1= %g",R1);
Printf("R2= %g",R2);
Printf("h = %g",h);


// Cercle intérieur R1

Point(1)={  0,  0,  0, h};
Point(2)={ R1,  0,  0, h};
Point(3)={  0, R1,  0, h};
Point(4)={-R1,  0,  0, h};
Point(5)={  0,-R1,  0, h};
Point(6)={  0,  0, R1, h};
Point(7)={  0,  0,-R1, h};

Circle (1) ={2 ,1 ,3};
Circle (2) ={3 ,1 ,4};
Circle (3) ={4 ,1 ,5};
Circle (4) ={5 ,1 ,2};
Circle (5) ={6 ,1 ,3};
Circle (6) ={3 ,1 ,7};
Circle (7) ={7 ,1 ,5};
Circle (8) ={5 ,1 ,6};

Line Loop(1) = {7, 4, 1, 6};
Ruled Surface(1) = {1};

Line Loop(2) = {4, 1, -5, -8};
Ruled Surface(2) = {2};

Line Loop(3) = {8, 5, 2, 3};
Ruled Surface(3) = {3};

Line Loop(4) = {3, -7, -6, 2};
Ruled Surface(4) = {4};

Surface Loop(5)     = {1, 2, 3, 4};
Physical Surface(1) = {1, 2, 3, 4};

// Cercle extérieur R2

Point(11)={  0,  0,  0, h};
Point(12)={ R2,  0,  0, h};
Point(13)={  0, R2,  0, h};
Point(14)={-R2,  0,  0, h};
Point(15)={  0,-R2,  0, h};
Point(16)={  0,  0, R2, h};
Point(17)={  0,  0,-R2, h};

Circle (11) ={12 ,11 ,13};
Circle (12) ={13 ,11 ,14};
Circle (13) ={14 ,11 ,15};
Circle (14) ={15 ,11 ,12};
Circle (15) ={16 ,11 ,13};
Circle (16) ={13 ,11 ,17};
Circle (17) ={17 ,11 ,15};
Circle (18) ={15 ,11 ,16};

Line Loop(11) = {17, 14, 11, 16};
Ruled Surface(11) = {11};

Line Loop(12) = {14, 11, -15, -18};
Ruled Surface(12) = {12};

Line Loop(13) = {18, 15, 12, 13};
Ruled Surface(13) = {13};

Line Loop(14) = {13, -17, -16, 12};
Ruled Surface(14) = {14};

Surface Loop(15)    = {11, 12, 13, 14};
Physical Surface(2) = {11, 12, 13, 14};



Volume(1) = {5, 15};
Physical Volume(1) = {1};
