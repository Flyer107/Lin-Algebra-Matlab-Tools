%%1
%Memory allocation (Speed)
A = zeros(5); 
%Loop to form matrix with A(i,j) = (i^2 + j^3 + 1)*(i + j^2 + 1)
for i = 1:1:5 
    for j = 1:1:5   
        A(i,j) = (i^2 + j^3 + 1)*(i + j^2 + 1);   
    end
end
%%2
syms a b c 
A = [ 1 2 4 ; 0 1 5 ; -2 -4 -3 ];
B = [ -2; 2; 9 ];
% Ax = b , therefore x = (A^-1)b || X = A\B 
X = A\B;
%Solves the Augmented matrix [A B]

%%3
%LU Decomposition
A = zeros(8); 
for i = 1:1:8
    for j = 1:1:8
        A(i,j) = min(A(i,i), A(j,j));
        A(i,i) = factorial(i+1)/(factorial(2)*factorial(i-1));   
    end 
end
%Intentionally looped twice for proper display
for i = 1:1:8    
    for j = 1:1:8        
        A(i,j) = min(A(i,i), A(j,j));        
        A(i,i) = factorial(i+1)/(factorial(2)*factorial(i-1));    
    end
end
[L, U] = lu(A);
%To display L and U matrices 
L
U 
%To calculate with given b
%[L, U] = lu(A2);
%disp('L = ');disp(L) 
%disp('U = ');disp(U) 
%d = L\b; 
%x = U\d; 
%disp('x = '); disp(x)

%%4
%Plotting unsigned Area
%Remember: (1/2) b/c triangle = half parallelogram visually 
%Area(s) = (1/2)*abs(det(a)) 
%Area( T*S ) = (1/2)*abs(det(a))*area(s) 
clear all; 
syms t
A = [ 7*cos((3/5)*t + (1/5)*sin(t)); 7*sin((3/5)*t + (1/5)*sin(t)) ];
B = [ 5*cos((t + 1)*(2*t)); 5*sin((t + 1)*(2*t)) ];
% AB = Concatenation of A & B 
AB = [7*cos((3/5)*t + (1/5)*sin(t)), 5*cos((t + 1)*(2*t));  7*sin((3/5)*t + (1/4)*sin(t)), 5*sin((t + 1)*(2*t)) ];
% ezplot(FUN,[A,B]) plots FUN(X) over A < X < B. 
ezplot(((1/2)*(abs(det(AB)))), [0, 2]); 
% fplot( ((1/2)*(abs(det(AB)))), [0 2] ); Works as well

%%5
%Graphing Nullspace and Columnspace
clear all; 
figure
%Creates another space for another graph 
%C = [ 1 3 11; 4 2 14; -2 -2 -10]; 
C = [ 1 3 11; 4 2 14; -2 -2 -10]; 

%Nullspace/ Kernel 
Z = null(C); 
%Start-> End 
% x   y   z 
%Z(1), Z(2), Z(3) retrieve the element in the matrix by the passed in 
%parameter 
%From zero to vector's "point" in space. 
plot3([0 Z(1)], [0 Z(2)], [0 Z(3)])
%Columnspace figure 
%Creates another space for another graph
A = sym([ 1 3 11; 4 2 14; -2 -2 -10]); 
B = colspace(A);
%z = -2.5(x+y) 
[x, y] = meshgrid(-10:1:10);
z = -2.5*(x+y); 
surf(x,y,z);
view([-50 45 20]);

%%6
%Compute eigenvalues and eigenvectors, given a matrix
H = [ 0 -1 0 -4 0 9; -1 0 -4 0 9 0; 0 -1 0 4 0 -9; -1 0 4 0 -9 0; 0 1 0 -4 0 -9; 1 0 -4 0 -9 0 ];
%... - Is row continuation 
[eig_vector_Matrix, eig_value_Matrix] = eig(H); 
%Separate Vectors 
eig_vectors = { eig_vector_Matrix(:,1) , eig_vector_Matrix(:,2) ...
    eig_vector_Matrix(:,3) ,  eig_vector_Matrix(:,4)...
    eig_vector_Matrix(:,5) ,  eig_vector_Matrix(:,6) };
%Separate Values 
eig_values = { eig_value_Matrix(1,1), eig_value_Matrix(2,2)... 
    eig_value_Matrix(3,3), eig_value_Matrix(4,4)... 
    eig_value_Matrix(5,5), eig_value_Matrix(6,6) } ;
% A = (P)(D)(P^-1) (P = e-vectrs, D = e-values diagonally)
