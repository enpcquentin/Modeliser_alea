// Notations
// A : matrice de transition de S (taille 2x2)
// B : proba conditionnelle de Y sachant S (taille 2x4)
// filt : proba conditionnelle de S sachant Y (taille 4x2)
// base : nombre entre 1 et 4 indexant la base (ACGT)

// Lit le fichier ADN (NE PAS MODIFIER)
function [adn]=lecture_adn(toto)
  A = mgetl(toto);
  A=[A;'tccggtgatc cgacaggtta cg'];
  B = ascii(A);
  B(find(B==32))=[]; 
  // replace ACGT by 1234
  code= "acgt" ;  
  acode = ascii(code);
  val = [1,2,3,4]; 
  for i=1:length(code);
  B(find(B==acode(i))) = val(i);
  end;
  adn=B;
endfunction;

// Programme principal (NE PAS MODIFIER)
// Retourne la probabilite d'etre en chaque valeur possible de l'etat adjoint
// connaissant la sequence d'ADN, ainsi que
// l'evolution des matrices A et B au cours des iterations
function [region,Aevol,Bevol]=result(adn,nbre_etat,nbre_iteration,A,B,p0)
   
  nbre_base=4;
  Ai=A; Bi=B; p0i=p0;
  Aevol=zeros(nbre_etat,nbre_iteration);
  Bevol=zeros(nbre_base*nbre_etat,nbre_iteration);
  n_adn=length(adn);
  l1=zeros(nbre_etat,n_adn);
  
  for i=1:nbre_iteration,
    disp(i);
    [Ai,Bi,p0i,region]=Estimation(Ai,Bi,p0i,adn,nbre_base,nbre_etat);
    
    Aevol(:,i)=diag(Ai);
    
    Bevol(:,i)=matrix(Bi',nbre_etat*nbre_base,1);
    
  end;
endfunction;

// Estimation des parametres (NE PAS MODIFIER)
// A chaque iteration, le programme principal
// fait appel Ã  cette fonction, qui se charge
// d'appeler les fonctions de prevision, de filtrage
// et de lissage
function [Aestim,Bestim,p0estim,l1]=Estimation(A,B,p0,adn,nbre_base,nbre_etat)
  Aestim=zeros(A);
  Bestim=zeros(B);
  [prediction,filtering]=PrevFilt(A,B,p0,adn,nbre_base,nbre_etat);
  [smooth1,smooth2]=Smooth(A,prediction,filtering,adn,nbre_base,nbre_etat);
  for x=1:nbre_base,
	Bestim(:,x)=sum(smooth1(:,(adn==x)),'c');
  end;
  Bestim=Bestim./(sum(Bestim,'c')*ones(1,nbre_base));
  p0estim=sum(smooth1,'c')/sum(smooth1);
  z=sum(smooth2,'c');
  z=matrix(z,nbre_etat,nbre_etat);
  Aestim=z./(sum(z,'c')*ones(1,nbre_etat));
  l1=smooth1;
endfunction;

// NE PAS MODIFIER (Fonction resolvant les Ã©quations de prevision et de filtrage en parallele)
function [prediction,filtering]=PrevFilt(A,B,p0,adn,nbre_base,nbre_etat);
  n=length(adn);
  filtering=zeros(nbre_etat,n);
  prediction=zeros(nbre_etat,n)
  z=Filtering(B,p0,adn(1));
  filtering(:,1)=z;
  for i=2:n,
	z=Prediction(A,z);
	prediction(:,i)=z;
	z=Filtering(B,z,adn(i));
	filtering(:,i)=z;
  end;
endfunction;

// NE PAS MODIFIER (Fonction resolvant les deux Ã©quations de lissage en parallÃ¨le)
function [smooth1,smooth2]=Smooth(A,prediction,filtering,adn,nbre_base,nbre_etat)
  n=length(adn);
  smooth1=zeros(nbre_etat,n);
  smooth2=zeros(nbre_etat^2,n);
  smooth1(:,n)=filtering(:,n);
  for i=n-1:-1:1,
	y=smooth1(:,i+1);
	filt=filtering(:,i);
	prev=prediction(:,i+1);
	z=Smooth1(y,A,prev,filt);
	smooth1(:,i)=z;
	smooth2(:,i+1)=Smooth2(y,A,prev,filt);
  end;
endfunction;

// Notations
// filt, prev ont mÃªme nature
// Ce sont les probas que Sn soit dans l'Ã©tat i
// Ã  l'issue du filtrage, de la prevision

// Ã‰quations de prediction
// Renvoie la proba que Sn soit dans l'etat i (Lemme 3)
function [y]=Prediction(A,filt)
  // A FAIRE
  y = A' * filt;
endfunction;

// equations de filtrage
// Renvoie la proba que Sn soit dans l'etat i (Lemme 4)
function [y]=Filtering(B,prev,base)
  // A FAIRE
  y = ( B(:,base) .* prev ) ./(sum(B(:,base) .* prev ));
endfunction;

// Ã‰quations de lissage
function [y1]=Smooth1(y0,A,prev,filt)
  // A FAIRE
  // Deuxieme equation de lissage du lemme 5
  inv = ones(2,1) ./ prev;
  y1 = (A *  (inv .* y0))  .* filt 

  //y1 = zeros(2,1);
  //for j=1:2
  //  y1(j,1) = ( A(j,:) * (inv .* y0)) .* filt(j) ;
  //end;
endfunction;

function [y2]=Smooth2(y1,A,prev,filt)
  // A FAIRE
  // Premiere equation de lissage du lemme 5
  y2 = zeros(4,1);
  for i=1:2
    for j=1:2
      y2(i+(j-1)*2,1) = A(i,j) * filt(i) * y1(j) / (prev(j)) ;
    end;
  end;
endfunction;

adn = lecture_adn('seq_lambda.txt');
nbre_etat = 2;
nbre_iteration = 1000;
matA0 = [0.28,0.72;0.19,0.81];
matB0 = [0.21,0.36,0.37,0.06;0.27,0.27,0.26,0.20];
vecp0 = [0.5;0.5] ;

[region,Aevol,Bevol]=result(adn,nbre_etat,nbre_iteration,matA0,matB0,vecp0);


figure(0);
plot2d(1:length(adn),region(1,:));
xtitle('Estimation des zones homogènes dans la séquence d ADN');


figure(1);
for i=1:8
  s=i;
  if i==8
    s=10;
  end;
  plot2d(1:nbre_iteration,Bevol(i,:),style=s);
end;
legends(['b(1,1)','b(2,1)','b(1,2)','b(2,2)','b(1,3)','b(2,3)','b(1,4)','b(2,4)'],[1,2,3,4,5,6,7,10],'ur') ;
xtitle('Evolution des coefficients de la matrice b au cours de l algorithme EM');

figure(2);
for i=1:2
  plot2d(1:nbre_iteration,Aevol(i,:),style=i)
end;
legends(['a(1,1)','a(2,1)'],[1,2],'ur') ;
xtitle('Evolution des coefficients diagonaux de la matrice a au cours de l algorithme EM');

