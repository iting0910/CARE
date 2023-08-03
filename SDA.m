function [it,X] = SDA(A,G,H)
format long e;
AA = A;
GG = G;
HH = H;

it=1; %iteration
max_it=1.0e+20;

n = size(G,1);
W = eye(n)+G*H;
V1 = inv(W)*A;
V2 = G*inv(W');
A_hat = A*V1;
G_hat = G+A*V2*A';
H_hat = H+(V1')*H*A;


A = A_hat;
G = G_hat;
H = H_hat;

while it<=max_it && norm(A)>=1.0e-14
    W = eye(n)+G*H;
    V1 = inv(W)*A;
    V2 = G*inv(W');
    A_hat = A*V1;
    G_hat = G+A*V2*A';
    H_hat = H+(V1')*H*A;

    G_hat = (G_hat+G_hat')/2;
    H_hat = (H_hat+H_hat')/2;
    
    A = A_hat;
    G = G_hat;
    H = H_hat;
    it=it+1;
end
X = H;

%rsdl = norm(-X*GG*X+AA'*X+X*AA+HH)

%rsdl = norm(AA'*H*inv(eye(n)+GG*H)*AA+HH-H,"fro");



