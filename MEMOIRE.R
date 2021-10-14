rm(list = ls())

#____________________________ Simulation d'un mvmt Brownien ______________________________________________________________
simuX = function(n, p, J = 50, lambda = (1:J) ^ (-2), loi = 'norm',...) {
  if (loi == 'norm')
    ksi = rnorm(n * J)
  else
    ksi = runif(n * J, -sqrt(3), sqrt(3))
  ksi = matrix(ksi, n, J) * matrix(rep(sqrt(lambda), n), n, J, byrow = TRUE)
  x = matrix(rep(seq(0, 1, length.out = p), J), J, p, byrow = TRUE)
  Jm = matrix(rep(seq(0.5, J - 0.5, by = 1), p), J, p)
  base = sqrt(2) * sin(pi * Jm * x)
  pertinit <- matrix(rep(rnorm(n), p), n, p, byrow = FALSE)
  ksi %*% base + pertinit
}
MvB <- function(n, p) {
  X = rnorm(n * p) / sqrt(p)
  X = matrix(X, ncol = p)
  MvB = matrix(0, ncol = p, nrow = n)
  for (i in 1:n)
  {
    MvB[i, ] = cumsum(X[i, ])
  }
  MvB
}

#______________________ Implementation des bases histogrammes et Fourier de L^2([0,1])______________________________________
base_histogramme <- function(n, t) {
  L <- matrix(0, n, length(t))
  for (i in 1:n) {
    L[i, ] = sqrt(n) * (t >= (i - 1) / n & t < (i) / n)
  }
  return(L)
}

fourier_pair <- function(t, k) {
  return(sqrt(2) * cos(k * pi * t))
}
fourier_impair <- function(t, k) {
  return(sqrt(2) * sin((k - 1) * pi * t))
}
base_fourier <- function(n, t) {
  L <- matrix(1, n, length(t))
  for (i in 2:n) {
    if (i %% 2 == 0) {
      L[i, ] = fourier_pair(t, i)
    }
    else{
      L[i, ] = fourier_impair(t, i)
    }
  }
  return(L)
}

#____________________ Fonction donnant la matrice Ghat caractérisant l'opérateur de covariance____________________________
M_ACPF <-  function(X, L, d, p) {
  X.cent <- scale(X, scale = FALSE) #On centre les données
  K_hat <- crossprod(X.cent) / (d * p)
  G_hat <- (L %*% K_hat %*% t(L)) / p
  return(G_hat) 
}

#__________________________________ Application sur données réelles et simulées______________________________________________________________
library("fda")

###Données réelles Canada Weather
Donneesreelles <- CanadianWeather # Chargement des données
DonneesreellesTemp <- t(Donneesreelles[["monthlyTemp"]]) # Chargement des données températures
n <- nrow(DonneesreellesTemp)
p <- ncol(DonneesreellesTemp)
t <- seq(0, 1, by = 1 / (p - 1))
L <- base_fourier(n, t) #On utilise la base Fourier ici


X_ACPF <- M_ACPF(DonneesreellesTemp, L, n, p) # Calcul de GHat
eigenACPF <- eigen(X_ACPF)$vector # Recuperation des fonctions propres de GHat
A <- crossprod(eigenACPF, L) # Reconstruction des fonctions propres de l'estimateur de Gamma

#Affichage des premières composantes principales
plot(
  A[1, ],
  type = "b",
  col = "blue",
  ylim = c(-5, 1),
  ylab = "Première composante principale",
  xlab = "Temps"
)
plot(
  A[2, ],
  type = "b",
  col = "blue",
  ylim = c(-3, 4),
  ylab = "Deuxième composante principale",
  xlab = "Temps"
)
plot(
  A[3, ],
  type = "b",
  col = "blue
",
  ylim = c(-4, 4),
  ylab = "Troisème composante principale",
  xlab = "Temps"
)
plot(
  A[4, ],
  type = "b",
  col = "blue",
  ylim = c(-4, 4),
  ylab = "Quatrième composante principale",
  xlab = "Temps"
)

# Caractéristiques données
###########################
eigenACPF_valeur <-
  eigen(X_ACPF)$values # Récuperation des valeurs propres de l'ACPF

V <- numeric(length(eigenACPF_valeur))
for (i in 1:length(eigenACPF_valeur)) {
  eigen = 0
  valeur = (eigen + eigenACPF_valeur[i]) / sum(eigenACPF_valeur)
  eigen = valeur
  V[i] = valeur
}
matrice_V <- t(matrix(V,  nrow = 1))
barplot(V[1:5], xlab = c("Vp_1       Vp_2       Vp_3        Vp_4       Vp_5"))

Scores <- DonneesreellesTemp %*% t(A) #Fonctions scores associées à toutes les composantes principales

#Affichage des scores des premières composantes principales
plot(
  Scores[, 1],
  type = "p",
  col = "blue",
  xlab = "Données",
  ylab = "Score de la deuxième composante principale",
  main = "Scores des stations météorologiques"
)
plot(
  Scores[, 2],
  type = "p",
  col = "blue",
  xlab = "Données",
  ylab = "Score de la deuxième composante principale",
  main = "Les scores des stations météorologiques"
)
plot(
  Scores[, 1],
  Scores[, 2],
  type = "p",
  col = "blue",
  xlab = "Score de la première composante principale",
  ylab = "Score de la deuxième composante principale",
  main = "Les scores des stations météorologiques"
)




###################################################
###################################################
##################################################





#Données réelles Pinch
Donneesreelles <- t(pinchraw) #Chargement des données
n <- nrow(Donneesreelles)
p <- ncol(Donneesreelles)
t <- seq(0, 1, by = 1 / (p - 1))
L <- base_fourier(n, t) #On utilise la base Fourier

X_ACPF <- M_ACPF(Donneesreelles, L, n, p) # Calcul de GHat
eigenACPF <- eigen(X_ACPF)$vector  # Recuperation des fonctions propres de GHat
A <- crossprod(eigenACPF, L) # Reconsctruction des fonctions propres de l'estimateur de Gamma

matplot(t(Donneesreelles),
        type = "l",
        xlab = "Temps en secondes" ,
        ylab = "Force exercée exprimée en Newton")

# Affichages des composantes principales
plot(
  A[1, ],
  type = "b",
  col = "blue",
  ylab = "Première composante principale",
  xlab = "Temps"
)
plot(
  A[2, ],
  type = "b",
  col = "blue",
  ylab = "Deuxième composante principale",
  xlab = "Temps"
)
plot(
  A[3, ],
  type = "b",
  col = "blue",
  ylab = "Troisème composante principale",
  xlab = "Temps"
)
plot(
  A[4, ],
  type = "b",
  col = "blue",
  ylab = "Quatrième composante principale",
  xlab = "Temps"
)

# Caractéristiques données
###########################
eigenACPF_valeur <- eigen(X_ACPF)$values # Récuperation des valeurs propres de l'ACPF
V <- numeric(length(eigenACPF_valeur))
for (i in 1:length(eigenACPF_valeur)) {
  eigen = 0
  valeur = (eigen + eigenACPF_valeur[i]) / sum(eigenACPF_valeur)
  eigen = valeur
  V[i] = valeur
}
matrice_V <- t(matrix(V,  nrow = 1))
barplot(V[1:5], xlab = c("Vp_1       Vp_2       Vp_3        Vp_4       Vp_5"))

Scores <- Donneesreelles %*% t(A) #Fonctions scores associées à toutes les composantes principales

#Score des premières composantes principales
plot(
  Donneesreelles[, 1],
  Scores[, 1],
  type = "p",
  col = "blue",
  xlab = "Données",
  ylab = "Score de la première composante principale"
)
plot(
  Donneesreelles[, 2],
  Scores[, 2],
  type = "p",
  col = "blue",
  xlab = "Données",
  ylab = "Score de la deuxième composante principale"
)
plot(
  A[1, ],
  A[2, ],
  type = "p",
  col = "blue",
  xlab = "Score de la première composante principale",
  ylab = "Score de la deuxième composante principale"
)
###################################################
###################################################
###################################################
##################################################


#Données simulées
n <- 1000
p <- 500
X <- simuX(n, p) #Chargement des données
Y <- MvB(n, p) #Chargement des données
t <- seq(0, 1, by = 1 / (p - 1))
L <- base_fourier(n, t) #On utilise la base fourier 

X_ACPF <- M_ACPF(X, L, n, p) # Calcul de GHat
X_ACPFB <- M_ACPF(Y, L, n, p)
eigenACPF <- eigen(X_ACPF)$vector  # Recuperation des fonctions propres de GHat
A <- crossprod(eigenACPF, L)# Reconstruction des fonctions propres de l'estimateur de Gamma
B <- crossprod(eigen(X_ACPFB)$vector, L)

for (i in 1:4) {
  matplot(
    t,
    A[i, ],
    type = "l",
    col = "blue",
    ylim = c(-3, 4),
    ylab = "Données simulées",
    xlab = "Temps"
  )
  par(new = T)
  matplot(
    t,
    B[i, ],
    type = "l",
    col = "green",
    ylim = c(-3, 4),
    ylab = "Données simulées",
    xlab = "Temps"
  )
}

###################################################

#Calculer importance axes
##################################################
eigenACPF_valeur <-
  eigen(X_ACPF)$values # Récuperation des valeurs propres de l'ACPF
V <- numeric(length(eigenACPF_valeur))
for (i in 1:length(eigenACPF_valeur)) {
  eigen = 0
  valeur = (eigen + eigenACPF_valeur[i]) / sum(eigenACPF_valeur)
  eigen = valeur
  V[i] = valeur
}
matrice_V <- t(matrix(V,  nrow = 1))
barplot(V[1:5], xlab = c("Vp_1       Vp_2       Vp_3        Vp_4       Vp_5"))
Scores <-
  X %*% t(A) #Fonctions scores associées à toutes les composantes principales
###################################################



#Convergence de ı = base=sqrt(2)*sin(pi*Jm*x)
###################################################
p = 500
J = 50
Jm = matrix(rep(seq(0.5, J - 0.5, by = 1), p), J, p)
x = matrix(rep(seq(0, 1, length.out = p), J), J, p, byrow = TRUE)
base = sqrt(2) * sin(pi * Jm * x)

n_1 = c(10, 50, 100, 500)
n2 = c("green", "red", "orange", "blue")
t <- seq(0, 1, by = 1 / (p - 1))
plot.new()
#Modifier
matplot(
  base[2, ],
  type = "l",
  ylim = c(-5, 5),
  ylab = "Données simulées",
  xlab = "Temps"
)
for (i in 1:4) {
  L <- base_fourier(n_1[i], t)
  X_ACPF <-
    M_ACPF(t(subset(t(X), select = seq(1, n_1[i], by = 1))), L, n, p) # Calcul de GHat
  eigenACPF <-
    eigen(X_ACPF)$vector  # Recuperation des fonctions propres de GHat
  A <-
    crossprod(eigenACPF, L)# Reconsctruction des fonctions propres de l'estimateur de Gamma
  par(new = T)
  matplot(
    A[2, ],
    type = "l",
    col = n2[i],
    ylab = "Données simulées",
    ylim = c(-5, 5),
    xlab = "Temps"
  )
  par(new = T)
}

############################################################################

eigenacpf_valeur <- eigen(x_acpf)$values # Récuperation des valeurs propres de l'acpf

v <- numeric(length(eigenacpf_valeur))

for (i in 1:length(eigenacpf_valeur)) {
  
  eigen <- 0
  
  valeur <- (eigen + eigenacpf_valeur[i]) / sum(eigenacpf_valeur)
  
  eigen <- valeur
  
  v[i] <- valeur
  
}

matrice_v <- t(matrix(v, nrow = 1))

barplot(v[1:5], xlab = c("vp_1       vp_2       vp_3        vp_4       vp_5"))

Scores <- DonneesreellesTemp %*% t(a) # Fonctions scores associées à toutes les composantes principales

eps<- matrix(0,p,n)

for(i in 1:p){
  for(j in 1:n){
    eps[i,j] = crossprod(DonneesreellesTemp[i,],a[j,])
  }
}
for(i in 1:p){
  for(j in 1:n){
    eps[i,j] = eps[i,j]/matrice_v[j,1]
  }
}
espace<-numeric(p)
for(i in 1:p){
  s=0
  for(j in 1:n){
    s = s + sum(matrice_v[j,1]*eps[i,j]*a[j,])
  }
  espace[i]=s
}
matplot(espace,type="l")