// Localisation cwipi

DefineConstant[ R1 = {1.0, Min 0.5, Max 2.0, Step 1.0, Name "Parameters/R1"} ];
DefineConstant[  h = {1.0, Min 0.5, Max 5.0, Step 0.5, Name "Parameters/h"} ];


Printf("Paramètres du maillage");

Printf("R1= %g",R1);
Printf("h = %g",h);


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
Line Loop(2) = {4, 1, -5, -8};
Line Loop(3) = {8, 5, 2, 3};
Line Loop(4) = {3, -7, -6, 2};

Ruled Surface(1) = {1};
Ruled Surface(2) = {2};
Ruled Surface(3) = {3};
Ruled Surface(4) = {4};

Surface Loop(17)    = {1,2,3,4};

Coherence;


Physical Surface("CondLim") = {1,2,3,4};


/*
Transfinite Line {6, 7, 8, 5, 1, 2, 3, 4} = 10 Using Progression 1;
Recombine Surface {1,2,3,4};
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Recombine Surface {4, 3, 2, 1};
*/

Coherence;
