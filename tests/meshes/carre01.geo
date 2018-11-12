/******************************
    Maillage carre
******************************/

DefineConstant[ alpha = {13, Min 1, Max 500, Step 1, Name "Parameters/alpha"} ];

DefineConstant[ Lx = {1, Min 1, Max 50, Step 1, Name "Parameters/Lx"} ];
DefineConstant[ Ly = {1, Min 1, Max 50, Step 1, Name "Parameters/Ly"} ];

h=1;

Printf("Paramètres du maillage cartésien");
Printf("Lx= %g",h);
Printf("Ly= %g",h);
Printf("alpha= %g",alpha);

Point(1) = { 0,  0, 0, h};
Point(2) = {Lx,  0, 0, h};
Point(3) = {Lx, Ly, 0, h};
Point(4) = { 0, Ly, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Physical Line(1) = {1};
Physical Line(2) = {2};
Physical Line(3) = {3};
Physical Line(4) = {4};

Transfinite Line {1,3} = alpha*Lx+1 Using Progression 1;
Transfinite Line {2,4} = alpha*Ly+1 Using Progression 1;

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Transfinite Surface {1};
Recombine Surface {1}; //quadCoherence;
Physical Surface(1) = {1};
