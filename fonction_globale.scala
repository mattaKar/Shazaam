object fonction_globale extends App{

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                   Class Complex                          /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // Les complexes ne sont pas définis en Scala, or la FTT qui est l'outil privilégié pour calculer les 
  // fréquences d'un signal donné nécessite des complexes dans ses calculs, nous allons donc implementer une classe
  // Complex qui nous sera très utile par la suite.
  
  
//on importe certaines fonctions de base du module scala.math pour réaliser la classe Complex
  import scala.math.{ Pi, cos, sin, cosh, sinh, abs }
  
case class Complex(reel:Double,imag:Double) {
  //définition des opérations de base (nécessaires pour le calcul de la fft)
  def + (x:Complex):Complex = Complex(reel+x.reel,imag+x.imag)
  def - (x:Complex):Complex = Complex(reel-x.reel,imag-x.imag)
  def * (x:Double):Complex = Complex(reel*x,imag*x)
  def * (x:Complex):Complex = Complex(reel*x.reel - imag*x.imag,reel*x.imag + imag* x.reel)
  def / (x:Double):Complex = Complex(reel/x,imag/x)
  
  //définition de l'affichage de la classe Complex
  override def toString():String = {
    val x="%1.4f" format reel
    val y="%1.4f" format abs(imag)
    (x,y) match {
      case(_,"0,0000") => x
      case("0,0000",_) => y
      case(_,_)        => if (imag>0) {
    	  				        	  x + "+" + "i" + y
                  				}else{
                  					x + "-" + "i" + y
                  				}    
    }
  }
  
//définition des exponentiels complexes -- exp(a+ib)=exp(a)*exp(ib)=(cosh(a)+sinh(a))*(cos(b)+isin(b)) --
def expo() : Complex = {
    val r = (cosh(reel) + sinh(reel))
    Complex(cos(imag), sin(imag)) * r

  }
def module():Double={
  math.sqrt(math.pow(reel,2)+math.pow(imag,2))
}
  }
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 Fonction butterfly                       /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

 // la fonction butterfly est la partie principale de la FFT, elle scinde le calcul en différentes étapes en 
 // supposant que le nombre de données à traiter soit une puissance de 2, elle réalise ensuite des calculs par
 // paire et obtient la FFT complexe

//définition de omega , fonction intermédiaire qui correspond au calcul de w(k,n) dans la formule G(k)+-w(k,n)F(k)
   def omega(o:Double,i:Int):Complex={
       val t=new Complex(0,-2*Pi*o/i)
       t.expo() 
    } 

//fonction principale
   
  def butterfly(A:Array[Complex]):Array[Complex] = {
    var B=A                                   // On construit un nouveau tableau à partir de celui donné en entrée
    var b:Int=1                               // b est le compteur d'étapes
    var k:Int=0                               // k est l'indice qui parcourt le tableau
    var c=0                                   // c est l'indice qui parcourt les sous paquets de calculs 
    var w:Double=0.0                          // l'indice de l'exponentiel complexe
    val n:Int=A.length 
    val u= math.log(n)/math.log(2)            // nombre d'étapes nécessaires au calcul
    var Bbis= Array.fill(n)(Complex(0,0))     // on construit un tableau de 0 que l'on remplira ensuite
    while (b <= u) {                          // on boucle sur les différentes étapes 
      k=0                         
      while (k < n) {                         //on boucle tant que l'on n'a pas atteint le bout du tableau
        c=0                       
        w=0
       while (c < math.pow(2,b-1).toInt){     //on boucle sur le compteur de sous tableau selon une première méthode de calcul
          Bbis(k)= B(k) + omega(w,n)*B(k+math.pow(2,b-1).toInt)
          c+=1
          k+=1
          w= w + math.pow(2,u-b)

        }
        c=0
        w=0
      while (c < math.pow(2,b-1).toInt){      //on boucle sur le compteur de sous tableau selon une seconde méthode de calcul          
          Bbis(k)= Complex(0,0) - (omega(w,n))*B(k) + B(k-math.pow(2,b-1).toInt)
          c+=1
         k+=1
          w= w + math.pow(2,u-b)

        }
        }
      b+=1
      for(i<-0 to n-1){
      B(i)=Bbis(i)
      }                                       // on écrase l'ancien tableau par le nouveau avant de changer d'étape 
      
    }
    return B
  }
  
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 Fonction makeComplex                     /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // définition de makeComplex, une fonction qui prend un tableau de décimaux en entrée et le transforme 
  // en un tableau de type Complex ie: x(i):reel -> x(i):(reel,0)
  
  def makeComplex(A:Array[Double]):Array[Complex] ={
     val n:Int=A.length
     var B=Array.fill(n)(Complex(0,0))
     for (i<-0 to n-1) {
       B(i)=Complex(A(i),0)
     }
     B
   }
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 Fonction inver                     /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // La fonction inver est une sorte de tri fusion : à chaque appel récursif, 
  // elle met les indices pairs dans  un sous tableau et les indices impairs dans un autre.
  // Et lorsqu'elle ne peut plus séparer les tableaux en 2 elle met tous les tableaux obtenus 
  // les uns à la suite des autres.
  
  def inver(x : Array[Complex]) : Array[Complex] = {  // A utiliser seulement avec des tableaux de taille pow(2,p)
 
     var n = x.length
     if (n==1) {
       x
     }
     else {
       var t1:Array[Complex]= new Array[Complex](n/2)
       var t2:Array[Complex]= new Array[Complex](n/2)
       for (i<-0 to n-1) {
         if (i%2==0) {
           t1(i/2)=x(i)
         }
         else{
           t2(i/2)=x(i)
         }
       }
       inver(t1)++inver(t2)
     }
  }
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 appel de la fonction                     /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
  // La fonction "FFT" utilise les différentes fonctions que nous avons définis avant pour enfin réaliser la ...
  // *roulement de tambour...* transformée de Fourier rapide d'un tableau donné en entrée et ce sous la forme
  // d'un tableau de décimaux (Double).
  
 def FFT(X:Array[Double]):Array[Double]={
   var V=X
   var N=makeComplex(V)
   N=inver(N)
   var G=butterfly(N):Array[Complex]
   var Q= Array.fill(G.length)(0.0) 
   for(i<-0 to G.length-1){ 
     var s=G(i).module()
     Q(i)=s
   }
   Q
   }
  
  
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 Conversion Stéréo à Mono                 /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  
 // Tous les fichiers de la base de données ne sont pas fournit en Mono, on est capable d'étudier seulement les 
 // fichiers en Mono, la fonction StoM permet donc de convertir les fichiers Stereo en Mono!
 // NB : on part du principe qu'un fichier à au maximum 2 canaux.
 
   def StoM (A : Array[Array[Int]]) : Array[Array[Int]] = {
    if (A(0)(1)==1) { // On regarde si on est en Mono ou pas,
      A               // comme ça on applique la fonction à chaque fois sans se poser de question.
    }
    else {
      A(0)(1)=1
      val n = A(1).length
      var C = new Array[Int](n)
      for (i<-0 to n-1) {
        C(i) = A(1)(i)+A(2)(i)  //On fait la somme des canaux
      }
      Array(A(0),C)
    }
  }
  

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                    harmonisation                         /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 
 // tous les fichiers musicaux fournis ne sont pas enregistrés de la même manière, ie selon la même fréquence 
 // d'échantillonage, la fonction harmonisation a pour but de ramener tous les fichiers sous un même format où
 // fe = 11025Hz
 // NB : on part du principe qu'un fichier est soit en 11025Hz soit en 22050Hz soit en 44100Hz
   
 def harmonisation(fichier:Array[Array[Int]]):Array[Array[Int]]={
   if (fichier(0)(0)==11025){     // Si jamais le fichier est déjà dans le bon échantillonage on le conserve
     fichier
   }else{
     var A=Array(Array(fichier(0)(0),fichier(0)(1),fichier(0)(2)))
    A(0)(0)=11025
     if (fichier(0)(0)==22050){     // on le modifie dans les autres cas
     A(0)(2)=A(0)(2)/2            // et on remplace aussi l'information dans le tableau caractéristique 
     var B=Array.fill(fichier(1).length/2)(0)
       for (i<-0 to B.length-1){
         B(i)=(fichier(1)(2*i) + fichier(1)(2*i+1))/2
       }
     Array(A(0),B)
   }
     else {                       // le cas 44Khz
     A(0)(2)=A(0)(2)/4
      var B=Array.fill(fichier(1).length/4)(0)
       for (i<-0 to B.length-1){
         B(i)=(fichier(1)(4*i) + fichier(1)(4*i+1) + fichier(1)(4*i+2) + fichier(1)(4*i+3))/4  
       }
     Array(A(0),B)
   }
 }
 }



///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                 création d'empreintes                    /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
import com.tncy.top.files.WavWrapper
import com.tncy.top.files.Utils;
 // La fonction empreintes permet de prendre une musique issue de notre base de données dont on fournit le 
 // chemin d'accès, on indique aussi la taille des échantillons que l'on souhaite prélever et la fonction fournit
 // un fichier similaire à celui entré où les canaux sont convertis en mono puis separés en autant d'échantillons
 // possibles.

// NB : on perd donc un peu d'information puisqu'on ne traite pas les valeurs en excès qui se trouve à la fin (celles 
// qui ne suiffisent pas pour créer un nouvel échantillon... sachant qu'un échantillon de taille 1024 ne correspond qu'à 0,1s
// on considère que cette perte d'information est négligeable. 

 def empreintes(file:String,hashsize:Int)={
   var wrappedWav : WavWrapper = new WavWrapper(file)          //on saisit le fichier desiré
   var wav2D1 : Array[Array[Int]] = wrappedWav.getWav()
   var wav2D=StoM(wav2D1)                          //conversion en Mono si le fichier ne l'est pas
   wav2D=harmonisation(wav2D) 
   var y=(math.pow(2,hashsize))                                //taille des échantillons desirée
   var W0=Array(Array(wav2D(0)(0).toDouble,wav2D(0)(1).toDouble,wav2D(0)(2).toDouble))   //initialisation du tableau d'échantillons desiré  
//var W0 : Array[Array[Double]] = Array(wav2D(0)) 
   var t=(wav2D(1).length/y).toInt                             //nombre d'échantillons en fonction de la taille de la chanson écoutée
// ici on perd des informations puisqu'on troncature le morceau en convertissant en entiers (on peut avoir length=t*y+quelquechose) 
   for (x<-0 to t-1){                                          //on remplie W0 d'emplacements pour fournir les différents échantillons  
     W0 = W0 ++ Array(Array.fill(y.toInt)(0.0))			// (ex:[wav2D(0),FFT(empreinte1),FFT(empreinte2),....])
   }
   for(i<-1 to t){                                             //on remplie réellement les échantillons
     for (j<-0 to y.toInt-1){
     W0(i)(j)=wav2D(1)((i-1)*y.toInt +j).toDouble
     }
   }
   for (i<-1 to W0.length-1){                                  //on applique la FFT sur chaque échantillon
     W0(i)=FFT(W0(i))
   }
   W0
 }
 
 

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////             Déterminer les points d'intérêts            /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// cette fonction est la première pour permettre de définir les marqueurs : 
// elle parcourt un tableau d'amplitudes en fonction de la fréquence (fourni par la FFT) à un temps donné et ressort
// les fréquences où les amplitudes sont remarquables. 
// On scinde la recherche des points en 6 bandes logarithmiques de fréquences pour pallier les conséquences de notre 
// audition logarithmique. 
// (Pour l'instant) On ne garde que la fréquence où l'amplitude est maximale et on ne garde la fréquence de la bande
// que si l'amplitude corresponde est supérieure à la moyenne des amplitudes retenues.

  import scala.math.pow 
  val Te=0.1 //seconde   
  val fe=11025 //Hz 
  val N=1024
  val factf=fe/N
  
  def points_interet(A:Array[Double]): Array[Double]={ // on entre un tableau d'amplitude (FFT) 
// les indices permettent d'avoir les fréquences (f=fe*k/N=k*factf)  
    var F : Array[Double] = Array() // tableau des fréquences à renvoyer 
    var B = Array(Array(A(0),0), Array(A(10),10*factf), Array(A(20),20*factf), Array(A(40),40*factf), Array(A(80),80*factf), Array(A(160),160*factf))
// tableau de 6 couples représentant (amplitudemax,freq) pour chaque bande de fréquence... 
// on initialise à la première valeur de chaque bande pour faire après la comparaison et chercher le maximum
    var AmplMoyTot : Double =0  // permettra de définir le seuil pour savoir si la bande est à regarder
    
// B(0) Bande faible (0 à 10)
    for (k<-1 to 9) {
      if (A(k)>B(0)(0))
        B(0)=Array(A(k),k*factf)
    }
    AmplMoyTot+=B(0)(0)         // on ajoute la valeur de l'amplitude retenue 

// B(1) Bande basse (10 à 20)  B(2) Bande moyenne basse (20 à 40)  B(3) Bande moyenne (40 à 80)  B(4) Bande moyenne haute(80 à 160)
    for (p <-0 to 3){
      for (k <- 10*math.pow(2,p).toInt+1 to 20*math.pow(2,p).toInt-1){
        if (A(k)>B(p+1)(0))
          B(p+1)=Array(A(k),k*factf)
      }
      AmplMoyTot+=B(p+1)(0)     // on ajoute la valeur de l'amplitude retenue
    }
    
// B(5) Bande élevée (160 à 511)
    for (k<-81 to 511) {
      if (A(k)>B(5)(0))
        B(5)=Array(A(k),k*factf)
    }
    AmplMoyTot+=B(5)(0)         // on ajoute la valeur de l'amplitude retenue

    AmplMoyTot=AmplMoyTot/6     // calcul de la moyenne des amplitudes retenues 
    for (p <- 0 to 5){
      if (B(p)(0)>AmplMoyTot){
        F=F++Array(B(p)(1))     // si on retient la bande, on préserve uniquemment la valeur de la fréquence 
      }
    }
    F
  }

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                       Marqueurs                          /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// cette fonction définie les marqueurs temporels dans un échantillon de musique 
// ce sont avec ces derniers qu'on va voir s'il y a une corrélation entre échantillon à tester et base de données 
// un marqueur temporel = (tps de l'ancrage, freq ancrage, freq d'un point, temps entre le point et l'ancrage)
// on défini que la taille de la zone cible, ie de la fenêtre est de 0,4 seconde

  
def marqueurs(spectro:Array[Array[Double]]) : Array[Array[Double]]={ 
// spectro est le spectrogramme de l'échantillon (issue des FFT successives) 

	var A = new Array[Array[Double]](0)
	for (i <- 1 to spectro.length-1){
		A=A++Array(points_interet(spectro(i)))
	}
// A est alors le tableau contenant les points d'interêt par temps 

	var M = new Array[Array[Double]](0)
	var x=0
	while (A(x).length==0){   // on se debarasse des premieres valeur de A qui sont parfois nulles
	  x+=1                    // (cas ou la musique ne commence pas au tout debut de l'enregistrement)
	}
	
	var tancr : Double = x*Te
	var fancr : Double = A(x)(0)
	var deltat : Double = 0
	
	for (t <- x to A.length-5){   // on parcourt le tableau des premières valeurs existantes jusqu'à la dernière "fenêtre"
	  if (A(t).length!=0){
		for(i<- 1 to 4){						// pour chaque temps de la zone cible 
		  deltat=i
			for (f <- 0 to A(t+i).length-1)	{	// pour chaque fréquence de ce temps on crée un marqueur 
				M= M++Array(Array(tancr,fancr,A(t+i)(f),deltat/10))	
			}
			
		}
		fancr = A(t)(0)           // puis on se déplace à la prochaine valeur que l'on défini comme fréquence d'ancrage
	  }
		tancr = t*Te 							// et on redéfinit la valeur temporelle correspondante 
	}
	// cette seconde boucle permet de réduire la taille de la fenêtre pour faire des marqueurs sur les dernières valeurs
	var decompteur=4
	for (t <- A.length-5 to A.length-1){
	  if (A(t).length!=0){
		for(i<- 1 to decompteur){						    // pour chaque temps de la zone cible 
		  deltat=i
			for (f <- 0 to A(t+i).length-1)	{			// pour chaque fréquence de ce temps on crée un marqueur 
				M= M++Array(Array(tancr,fancr,A(t+i)(f),deltat/10))	
			}
			
		}
		decompteur -=1
		fancr = A(t)(0)           // puis on se déplace à la prochaine valeur que l'on défini comme fréquence d'ancrage
	  }
		tancr = t*Te 							// et on redéfinit la valeur temporelle correspondante  
	}
	  
	 M
}


///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                    Base de données                        ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

var NomsChansons = new Array[String](10)
var M = new Array[Array[Array[Double]]](0) //Tableaux des marqueurs
val hashsize = 10   //échantillons de 1024 valeurs (=2^(10))

// La fonction BasedeDonnees va stocker tous les marqueurs des chansons connues dans un tableau de marqueurs M 
// qui est une variable globale

def BasedeDonnees (path : String) = {
  NomsChansons = Utils.listFiles(path)
  println("building the database")
  for (i<- 1 to NomsChansons.length-1) {
    println(((i.toDouble/NomsChansons.length.toDouble)*100).toInt + "%")
    var Ei = empreintes(path+"/"+NomsChansons(i),hashsize) // Donne les FFT sur la chanson i
    var Marqueursi = marqueurs(Ei) // Tableau où l'on va stocker les marqueurs de la chanson i 
    val n = Ei.length
    M = M ++ Array(Marqueursi)
  }
  println("100%")
  println("database is built")
  M
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                    match sample                           ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// La fonction match_sample procède en deux étapes, tout d'abord, elle répertorie dans un tableau l'ensemble des 
// marqueurs qui appartiennent communément à la musique et à l'echantillon testé, puis une fois ce tableau 
// construit, elle compte le nombre de fois qu'on retrouve les mêmes différence de temps(retard) entre les marqueurs
// échantillon et musique. Enfin la fonction sélectionne le maximum et donne ainsi le plus grande concordance
// de marqueurs cohérant en temps.


def match_sample(marq : Array[Array[Double]],i: Int, M:Array[Array[Array[Double]]]) : Int = {
  var Mi = M(i)   // on restocke les tableaux
  var tps = new Array[Array[Double]](0)
  val n = marq.length
  val m = Mi.length
  var c=0
  for (k<- 0 to n-1) {    //Les boucles for regardent les marqueurs identiques et remplissent le tableau tps avec les temps qui correspondent à ses marqueurs 
    for (j<- 0 to m-1) {
      if (Mi(j)(1)==marq(k)(1) && Mi(j)(2)==marq(k)(2) && Mi(j)(3)==marq(k)(3)){
        tps = tps ++ Array(Array(Mi(j)(0),marq(k)(0)))
      }
    }
  }
  if (tps.length == 0){  //S'il n'y a aucun marqueurs identiques : ce n'est pas la bonne musique 
    0
  }
  else {
    var U1=Array.fill(7000)(0)  // nous avons choisit cette taille de tableau pour indexé l'ensemble des différences en temps, c'est grossier mais ça marche
    for(t<-0 to tps.length-1){
    var O=((tps(t)(0)-tps(t)(1))*10).toInt   // comme on est sur des dixièmes de seconde, on multiplie par 10 la différence ce qui nous donne un indice du tableau à indenté
    U1(O+3500)= U1(O+3500) + 1               // ce qui correspond à l'indentation
    }
    var w=U1(0)                              // les prochaines lignes correspondent à la recherche du maximum
    for (p<-1 to U1.length-1){
      if (U1(p)>w){
        w=U1(p)
      }
    }

    w   // on renvoi donc le maximum !
    
  }  
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                    Reconnaissance                         ////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// La fonction reconnaissance, à un échantillon de musique donné (par un chemin) va regarder
// si celui-ci correspond à une musique de la base de données (en appelant match_sample) et renvoie
// une phrase en fonction du résultat.

import scala.io.Source

var NomsChansonsbase = new Array[String](NomsChansons.length)

def reconnaissance(path : String,basepath:String) = {
  var E = empreintes(path,hashsize)
  var marq = marqueurs(E)
  var L1=Array.fill(NomsChansons.length)(0)
  var nbmarqueur=Array.fill(NomsChansons.length)(0)
  var c1=0
  for(line <- Source.fromFile("/Users/matta/Desktop/ecriture/nbmarqueur.txt").getLines()){
   nbmarqueur(c1)=line.toInt
     c1+=1
   }
   var tabl:Array[Array[Array[Double]]]=Array(Array.fill(nbmarqueur(0))(Array(0.0,0.0,0.0,0.0)))
   for (i<-1 to NomsChansons.length-2){
   tabl= tabl ++ Array(Array.fill(nbmarqueur(i))(Array(0.0,0.0,0.0,0.0)))
   }
   for (i<-0 to NomsChansons.length-2){ 
   var cbis=0
   var marqactu=0
   for(line <- Source.fromFile(basepath+"/chanson"+(i+1).toString+".txt").getLines()){
    if (cbis!=3){
  tabl(i)(marqactu)(cbis)=line.toDouble
  cbis+=1
  }
  else{
    tabl(i)(marqactu)(cbis)=line.toDouble
    cbis=0
    marqactu+=1
  }
}
   }
  for(i<-0 to NomsChansons.length-2) { 
    L1(i) = match_sample(marq,i,tabl)
  }
  var maxi=L1(0)
  var indiceclef=0
  for(i<-0 to L1.length-1){
    if (maxi<L1(i)){  // on cherche la chanson qui a la meilleure correspondance 
      maxi=L1(i)
      indiceclef=i
    }
    
  }
  if (maxi<=60){     // on s'est fixé 60 comme valeur limite, on a en effet remarqué que les chansons qui n'étaient pas la bonne avait rarement plus de 30 valeurs tandis que les autres en avait rarement moins de 80
      if(maxi>=30){
        println("La chanson que vous écoutez est peut etre " + NomsChansons(indiceclef+1) + " ?")
    
      }else{
        println("Désolé, nous n'avons pas reconnu la chanson...")
      }
  }
  else {
    println("La chanson que vous écoutez est " + NomsChansons(indiceclef+1) + " !")
  }
}

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////           ecriture de la base de donnée en machine       /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

import java.io._
 def ecriturebase(base:String)={
var BD=BasedeDonnees(base)
for(i<-0 to BD.length-1){
  val writer = new PrintWriter(new File("/Users/matta/Desktop/ecriture/chanson"+(i+1).toString+".txt" ))
  for(j<-0 to BD(i).length-1){
    for(k<-0 to BD(i)(j).length-1){
      writer.println(BD(i)(j)(k))
    }
  }
  writer.close()
}
val writer = new PrintWriter(new File("/Users/matta/Desktop/ecriture/nbmarqueur.txt" ))
for(i<-0 to BD.length-1){
  writer.println(BD(i).length)
}
writer.close
}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////                         tests                            /////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//ecriturebase("/Users/matta/Desktop/test")
NomsChansons = Utils.listFiles("/Users/matta/Desktop/test")
println("TEST CONNU MONO 11KHz")
println("testing laurie ")
reconnaissance("/Users/matta/Desktop/musiques_top/WAV_Echantillons_connus_Mono_11025Hz/Je_serai_Mono-11025Hz.wav","/Users/matta/Desktop/ecriture")
println("testing little talks")
reconnaissance("/Users/matta/Desktop/musiques_top/WAV_Echantillons_connus_Mono_11025Hz/Little_Talks_Mono-11025Hz.wav","/Users/matta/Desktop/ecriture")
 println("testing raiders march")
 reconnaissance("/Users/matta/Desktop/musiques_top/WAV_Echantillons_connus_Mono_11025Hz/The_Raiders_March_Mono-11025Hz.wav","/Users/matta/Desktop/ecriture")

 println("TEST CONNU STEREO 44KHz")
 println("testing laurie ")
reconnaissance("/Users/matta/Desktop/musiques_top/WAV_Echantillons_connus_Stereo_44100Hz/Je_serai.wav","/Users/matta/Desktop/ecriture")
println("testing little talks")
reconnaissance("/Users/matta/Desktop/musiques_top/WAV_Echantillons_connus_Stereo_44100Hz/Little_Talks.wav","/Users/matta/Desktop/ecriture")
 println("testing raiders march")
 reconnaissance("/Users/matta/Desktop/musiques_top/WAV_Echantillons_connus_Stereo_44100Hz/The_Raiders_March.wav","/Users/matta/Desktop/ecriture")

 println("TEST INCONNU STEREO 44KHz")
 println("testing melissa ")
reconnaissance("/Users/matta/Desktop/musiques_top/WAV_Echantillons_inconnus_Stereo_44100Hz/Ma_Melissa.wav","/Users/matta/Desktop/ecriture")
println("testing ponies")
reconnaissance("/Users/matta/Desktop/musiques_top/WAV_Echantillons_inconnus_Stereo_44100Hz/Ponies.wav","/Users/matta/Desktop/ecriture")
 println("testing pirates")
 reconnaissance("/Users/matta/Desktop/musiques_top/WAV_Echantillons_inconnus_Stereo_44100Hz/The_Beast_of_Pirates_bay.wav","/Users/matta/Desktop/ecriture")
}