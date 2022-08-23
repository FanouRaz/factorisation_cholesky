#include<iostream>
#include<vector>
#include<fstream>
#include<iomanip>
#include<cmath>

using namespace std;

/*
    Factorisation de Cholesky d'une matrice A symétrique définie positive.
    ie décomposer A sous forme A = B*(B)t
    où B est une matrice triangulaire supérieure
*/

vector<vector<float>> cholesky(vector<vector<float>> A);
void afficher(vector<vector<float> > matrice); 
void lireData(int &n,vector<vector<float>> &A, vector <float> &b);
vector<float> remplirLigne(string ligne, int n);
void afficherVect(vector<float>tab);
vector<vector<float>> initialisation(int n);
void resolution(vector<vector<float>> A,vector<float> b);
void afficheEquation(vector<vector<float>> A,vector<float> B, int dim);
void transpositionTriangulaire(vector<vector<float>> &M,int n);

int main(){
    int n = 0;
    vector<vector<float>> A({});
    vector<float> b({});
    lireData(n,A,b);
    cout<<"Ce programme permet de résoudre un système d'equation "
    <<"sous forme d'équation matricielle Ax = B par la méthode de Cholesky, où A est une matrice symétrique définie positive . L'utilisateur inserera les matrices A et B dans le fichier data.txt et le programme"
    <<" affichera les solutions"<<endl;
    cout<<"_________________________________________________________________________________________________________________"<<endl<<endl;
    cout<<"Le système entré par l'utilisateur est:"<<endl;
    afficheEquation(A,b,n);
    cout << "La matrice triangulaire B obtenue par la factorisation de Cholesky de A est: "<<endl;
    afficher(cholesky(A));
    resolution(A,b);
    return 0;
}

void transpositionTriangulaire(vector<vector<float>> &M,int n){
    for(int i(0);i<n;i++){
        for(int j(0);j<n;j++){
            if(j<i){
                M[i][j] = 0;
            }
            else{
                M[i][j] = M[j][i];
            }
        }
    }
}

void lireData(int &n,vector<vector<float>> &A, vector <float> &b){
    fstream donnee;    
    string line;
    int taille(0);
    int tailleA=0;
    int tailleB=0;
    //ouverture du fichier en mode lecture
    donnee.open("data.txt", ios::in);                                
    if(!donnee){
        cout<<"Erreur lors du chargement des données. Veuillez placer le fichier data.txt dans le répertoire courant."<<endl;
        exit(-1);                                                                                       
    }
    
    while(getline(donnee,line)){
        if(taille==0)  n=stoi("0"+line);
        else if(tailleA < n){
            A.push_back(remplirLigne(line , n));
            tailleA++;
        }
        else{
            b.push_back(stof(line));
            tailleB++;
        }
        taille++;
    }

    donnee.close();
    if(tailleA!= n){
        cout<<"La dimension de la matrice A dans le fichier 'data.txt' n'est pas valide"<<endl;
        exit(-1);
    }
    if(tailleA != tailleB){
        cout<<"La dimension du vecteur B dans le fichier 'data.txt' n'est pas valide"<<endl;
        exit(-1);
        }
    }

//Construire les matrices A et b à partir des données recupérées dans le fichier data.txt 
vector<float> remplirLigne(string ligne, int n){
    vector <float> row;
    int i=0;
    int colonne=0;
    ligne+=" 0";
    int len =ligne.length();
    string valeurString="";
    while(i<len){
        if((ligne.at(i)>='0'&& ligne.at(i)<='9')||(ligne.at(i)=='.')||(ligne.at(i)=='+')||(ligne.at(i)=='-')){
            valeurString+=ligne.at(i);
        }
        else{
            if(valeurString.length()){
                row.push_back(stof(valeurString));
                colonne++;
            }
            valeurString="";
        }
            i++;
    }
    if(colonne != n){
        cout<<"La dimension de la matrice dans le fichier 'data.txt' n'est pas valide"<<endl;
        exit(-1);
    }
    return row;
}

//Determination de la matrice B obtenu par la factorisation de Cholesky de A 
vector<vector<float>> cholesky(vector<vector<float>> A){
    int n (A.size());
    vector<vector<float>> B(initialisation(n));
    float s(0);

    for(int i(0);i<n;i++){
        for(int j(0);j<=i;j++){
            if(i == j){
               for(int k(0);k<i;k++){
                 s += B[i][k]*B[i][k];
               } 
               B[i][j] = sqrt(A[i][i] - s);
               s=0;
            }
            else{
                for(int k(0);k<j;k++){
                    s += B[i][k]*B[j][k];
                } 
               B[i][j] = (A[i][j] - s) / B[j][j];
               s=0;
            }
        }
    }

    return B;
}

void afficheEquation(vector<vector<float>> A,vector<float> B, int n){
	
	for(int i=0; i< n ; i++){
		for(int j=0 ; j < n ; j++){
		
			cout << setw(10)<< A[i][j] << " x" << j << "\t+";
		}
		cout << setw(10)<<"= \t" << B[i] <<setw(2) << endl << endl;
	}
		
}

//Résolution de l'équation A.x = b , on remplace A par B(B)t puis on pose  y= (B)t x
void resolution(vector<vector<float>> A,vector<float> b){
    int n = A.size();
    vector<vector<float>>B(cholesky(A));
    float y[n] = {0};
	float x[n] = {0} ;
	float s;
    /*
        On a B.y = b    et (B)t . x = y
        B étant triangulaire inférieure, (B)t est triangulaire supérieure
    */

   //on calcule d'abord les yi en itérant de haut en bas
	for(int i(0);i<n;i++){
		s = 0;
		for(int j(0);j<i;j++){
			s += B[i][j] * y[j];
		}		
 		y[i]= ((b[i]-s)/B[i][i]);
	}

    //on transpose ensuite B et on calcule les xi à l'aide de l'equation (B)t.x = y, On itère cette fois vers le bas
	transpositionTriangulaire(B,n);

	for(int i(n-1);i>=0;i--){
		s=0;
		for(int j(i);j<n;j++){
			s += B[i][j] * x[j];
		}		
 		x[i]= ((y[i]-s)/B[i][i]);
 	}
    cout<<endl<<"Ainsi, comme A = B.(B)t, Ax = B.(B)t.x. Si on pose y=(B)t.x on a le système:"
        <<"B.y = b et (B)t.x = y. B et (B)t étant triangulaire, on peut résoudre facilement cette équation."<<endl<<"les solutions de cette équation sont:"<<endl;

    for(int i(0);i<n;i++){
        cout<<"x"<<i<<" = "<<x[i]<<endl;
    }
}

//Initialisation matrice
vector<vector<float>> initialisation(int n){
    vector<float> ligne;
    vector<vector<float>> M;
    for(int j(0);j<n;j++){
            ligne.push_back(0);
    }
    
    for(int i(0);i<n;i++){
        M.push_back(ligne);
    }
    return M;
}

//Afficher une matrice n x n
void afficher(vector<vector<float> > matrice){
    int size = matrice.size();
    cout<<setw(10);
    for(int i=0;i<size;i++){
        for(float ligne: matrice[i]){
            cout << ligne << setw(10);
        }
        cout<< endl;
    }
}

//Fonction pour afficher un vecteur colonne
void afficherVect(vector<float> tab){
    for(float elt : tab){
        cout<<elt<<endl;
    }
    cout<<endl;
}
