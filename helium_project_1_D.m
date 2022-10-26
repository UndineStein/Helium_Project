clc
clear all
close all

# Atomic number
z = 2;
n_max = 19;
# Maximum principal quantum number
y = zeros(n_max,1);

for n_max = 1:n_max

# Size of Hamiltonian
size_H = 2*n_max - 1;

# Initialization
H = zeros(size_H);
V2 = zeros(size_H);

# Function returning the energy of separated electrons
function e = ener_sep ( z, n1, n2 )
  # z = atomic number
  # n1,n2 = principal quantum numbers of electrons
  # e = energy in atomic units (hartrees)
  e = -z^2/2. * ( 1./n1^2 + 1./n2^2 );
endfunction

#ener_sep(2,1,1)


##### Functions dedicated to the calculation of <1/r12> #####

# Various quantities
function a = A_n ( z, n )
  a = -sqrt(z/n) * factorial(n-1);
endfunction

function b = B_nk ( z, n, k )
  b = 1. / ( factorial(k) * factorial(k+1) * factorial(n-k-1) );
endfunction

function c = C_kk ( z, np, kp, n, k )
  b_n  = 2.*z/n;
  b_np = 2.*z/np;
  c = (-1.)^(k+kp) * (2./(b_n+b_np))^(k+kp+2) * B_nk(z,np,kp) * B_nk(z,n,k) * b_n^(k+1) * b_np^(kp+1);
endfunction

#C_kk(2,1,0,1,0)

# Function returning the Coulombic energy
# <n1p s, n2p s| 1/r12 |n1 s, n2 s>
function e = ener_r12 ( z, n1p, n2p, n1, n2 )
  # z = atomic number
  # n1,n2,n1p,n2p = principal quantum numbers of electrons
  # e = energy in atomic units (hartrees)

  a  = A_n(z,n1)*A_n(z,n1p)*A_n(z,n2)*A_n(z,n2p);
  b1 = 2.*z*(1./n1+1./n1p);
  b2 = 2.*z*(1./n2+1./n2p);

  e = 0.;
  for k1p = 0:(n1p-1)
    for k1 = 0:(n1-1)
      c1 = C_kk ( z, n1p, k1p, n1, k1 );
      for k2p = 0:(n2p-1)
        for k2 = 0:(n2-1)
          c2 = C_kk ( z, n2p, k2p, n2, k2 );
          e = e + c1*c2*2./b2*factorial(k1+k1p+1)*factorial(k2+k2p+2) ...
                * ( 1. - (b1/(b1+b2))^(k1+k1p+2) );
          for l2 = 0:(k2+k2p)
            e = e + c1*c2*2./b1*factorial(k2+k2p+1) ...
                  * (b1/(b1+b2))^(k1+k1p+3) * (b2/(b1+b2))^l2 ...
                  * factorial(k1+k1p+2+l2) / factorial(l2) ...
                  * ( 1. - (k2+k2p+2.)/(l2+1.) );
          end
        end
      end
    end
  end
  e = e * a;
endfunction

#ener_r12(2,1,1,1,1)


#####  Builds the Hamiltonian  #####


# 1st line and 1st column corresponds
# to n1=n2=1 or n1p=n2p

H(1,1) = ener_sep ( z, 1., 1. ) + ener_r12 ( z, 1., 1., 1., 1. );

# 1st line: <1s1s|H|1s n2 s> in col 2,4, ...
for n2 = 2:n_max
  j = 2*n2-2;
  H(1,j) = ener_r12 ( z, 1., 1., 1., n2 );
end
# 1st line: <1s1s|H|n1 s 1s> in col 2,4, ...
for n1 = 2:n_max
  j = 2*n1-1;
  H(1,j) = ener_r12 ( z, 1., 1., n1, 1. );
end

# 1st col.: <1s n2p s|H|1s1s> in line 2,4, ...
for n2p = 2:n_max
  i = 2*n2p-2;
  H(i,1) = ener_r12 ( z, 1., n2p, 1., 1. );
end
# 1st col.: <n1p s 1s|H|1s1s> in line 2,4, ...
for n1p = 2:n_max
  i = 2*n1p-1;
  H(i,1) = ener_r12 ( z, n1p, 1., 1., 1. );
end
# Rest of the Hamiltonian
#
# Scan lines
for i = 2:size_H
  # < 1s n2p s |
  if ( mod(i,2) == 0 )
    n1p = 1;
    n2p = (i+2)/2;
  # < n1p s 1s |
  else
    n1p = (i+1)/2;
    n2p = 1;
  endif
  # Scan columns
  for j = 2:size_H
    # | 1s n2 s >
    if ( mod(j,2) == 0 )
      n1 = 1;
      n2 = (j+2)/2;
    # | n1 s 1s >
    else
      n1 = (j+1)/2;
      n2 = 1;
    endif
    H(i,j) = ener_r12 ( z, n1p, n2p, n1, n2 );
    if ( i == j )
    H(i,j) = H(i,j) + ener_sep ( z, n1, n2 );
    endif
  end
end
# up to the first order

##### Eigenvalues of the Hamiltonian
a = eig(H);
[V,D] = eig(H);
y_h(n_max) = D(1,1);

##### Adding corrections to the ground state energy #######

ener_pert = ener_sep(z, 1, 1) + ener_r12 ( z, 1., 1., 1., 1. );

for i = 2:n_max
  ener_pert = ener_pert - (ener_r12(z, i, 1, 1, 1))**2/(ener_sep ( z, i, 1 ) - ener_sep ( z, 1, 1 ));
endfor

for i = 2:n_max
  ener_pert = ener_pert -(ener_r12(z, 1, i, 1, 1))**2/(ener_sep ( z, 1, i ) - ener_sep ( z, 1, 1 ));
endfor

##################

ener_pert_s = ener_sep(z, 1, 2) + (ener_r12 ( z, 1, 2, 1, 2 ) + ener_r12(z,1,2,2,1));
ener_pert_a = ener_sep(z, 1, 2) + (ener_r12 ( z, 1, 2, 1, 2 ) - ener_r12(z,1,2,2,1));

####Symmetrical part
for i = 3:n_max
     ener_pert_s = ener_pert_s - 0.5*(ener_r12(z, i, 1, 1, 2) + ener_r12(z, i, 1, 2, 1))**2/(ener_sep ( z, i, 1 ) - ener_sep ( z, 1, 2 ));
endfor

for i = 3:n_max
      ener_pert_s = ener_pert_s - 0.5*(ener_r12(z, 1, i, 1, 2) + ener_r12(z, 1, i, 2, 1))**2/(ener_sep ( z, 1, i ) - ener_sep ( z, 1, 2 ));
endfor

####Antisymmetrical part
for i = 3:n_max
     ener_pert_a = ener_pert_a - 0.5*(ener_r12(z, i, 1, 1, 2) - ener_r12(z, i, 1, 2, 1))**2/(ener_sep ( z, i, 1 ) - ener_sep ( z, 1, 2 ));
endfor

for i = 3:n_max
      ener_pert_a = ener_pert_a - 0.5*(ener_r12(z, 1, i, 1, 2)-ener_r12(z, 1, i, 2, 1))**2/(ener_sep ( z, 1, i ) - ener_sep ( z, 1, 2 ));
endfor

#### Add ground state energy impact
ener_pert_s = ener_pert_s - 0.5*(ener_r12(z,1.,1.,1.,2.) + ener_r12(z,1.,1.,2.,1.))**2/(ener_sep ( z, 1, 1 ) - ener_sep ( z, 1, 2 ));

##### Initialization of the values we're going to graph
y(n_max) = ener_pert;
y_a(n_max) = ener_pert_a;
y_s(n_max) = ener_pert_s;
endfor

##### Interpolation ######
x = (1:n_max)';
xnew = linspace(max(x),min(x),100);
vq = interp1(x,y_h,xnew,'pchip');

######### Plot ######
plot(x,y_h,'o', xnew, vq)
xlabel('nmax')
ylabel('Energy in atomic unites')
title('Commparison of the ground state energy without taking into account the corrections depending on nmax')

ener_pert
ener_pert_a
ener_pert_s

#### Eigenvalues of the first excited state matrix
E1 = [H(2,2) H(2,3); H(3,2) H(3,3)]
b = eig(E1)
[V1,D1] = eig(E1)


##### Eigenvalues of the Hamiltonian
#{
a = eig(H)
[V,D] = eig(H)
#}
