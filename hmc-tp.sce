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
// Retourne la probabilitÃ© d'Ãªtre en chaque valeur possible de l'Ã©tat adjoint
// connaissant la sÃ©quence d'ADN, ainsi que
// l'Ã©volution des matrices A et B au cours des itÃ©rations
function [region,Aevol,Bevol]=result(adn,nbre_etat,nbre_iteration,A,B,p0)
   
  nbre_base=4;
  Ai=A; Bi=B; p0i=p0;
  Aevol=zeros(nbre_etat,nbre_iteration);
  Bevol=zeros(nbre_base*nbre_etat,nbre_iteration);
  n_adn=length(adn);
  l1=zeros(nbre_etat,n_adn);
  
  for i=1:nbre_iteration,
    
    [Ai,Bi,p0i,region]=Estimation(Ai,Bi,p0i,adn,nbre_base,nbre_etat);
    
    Aevol(:,i)=diag(Ai);
    
    Bevol(:,i)=matrix(Bi',nbre_etat*nbre_base,1);
    
  end;
endfunction;

// Estimation des paramÃ¨tres (NE PAS MODIFIER)
// A chaque itÃ©ration, le programme principal
// fait appel Ã  cette fonction, qui se charge
// d'appeler les fonctions de prÃ©vision, de filtrage
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

// NE PAS MODIFIER (Fonction rÃ©solvant les Ã©quations de prÃ©vision et de filtrage en parallÃ¨le)
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

// NE PAS MODIFIER (Fonction rÃ©solvant les deux Ã©quations de lissage en parallÃ¨le)
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
// Ã  l'issue du filtrage, de la prÃ©vision

// Ã‰quations de prÃ©diction
// Renvoie la proba que Sn soit dans l'Ã©tat i (Lemme 3)
function [y]=Prediction(A,filt)
  // A FAIRE
  y = A' * filt;
endfunction;

// Ã‰quations de filtrage
// Renvoie la proba que Sn soit dans l'Ã©tat i (Lemme 4)
function [y]=Filtering(B,prev,base)
  // A FAIRE
  y = ( B(:,base) .* prev ) ./(sum(B(:,base) .* prev ));
endfunction;

// Ã‰quations de lissage
function [y1]=Smooth1(y0,A,prev,filt)
  // A FAIRE
  // DeuxiÃ¨me Ã©quation de lissage du lemme 5
  inv = ones(2,1) ./ prev;
  y1 = zeros(2,1);
  for j=1:2
    y1(j,1) = ( A(j,:) * (inv .* y0)) .* filt(j) ;
  end;
endfunction;

function [y2]=Smooth2(y1,A,prev,filt)
  // A FAIRE
  // PremiÃ¨re Ã©quation de lissage du lemme 5
  y2 = zeros(4,1);
  for i=1:2
    for j=1:2
      y2(i+(j-1)*2,1) = A(i,j) * filt(i) * y1(j) / (prev(j)) ;
    end;
  end;
endfunction;

adn = lecture_adn('seq_lambda2.txt');
nbre_etat = 2;
nbre_iteration = 1000;
matA0 = [0.45 ,0.55;0.51 ,0.49];
matB0 = [0.2, 0.2, 0.3, 0.3;0.2, 0.3, 0.2, 0.3];
vecp0 = [0.5 ; 0.5] ;
[region,Aevol,Bevol]=result(adn,nbre_etat,nbre_iteration,matA0,matB0,vecp0);


figure(0);
plot(1:length(adn),region(1,:));


figure(1);
for i=1:8
  plot(1:nbre_iteration,Bevol(i,:))
end;


figure(2);

for i=1:2
  plot(1:nbre_iteration,Aevol(i,:))
end;


//crabe-duchemin-oreistein
//crabe-duchemin_oreistein