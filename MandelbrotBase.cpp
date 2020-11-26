# include <iostream>
# include <cstdlib>
# include <string>
# include <chrono>
# include <cmath>
# include <vector>
# include <fstream>
#include <mpi.h>

/** Une structure complexe est d�finie pour la bonne raison que la classe
 * complex propos�e par g++ est tr�s lente ! Le calcul est bien plus rapide
 * avec la petite structure donn�e ci--dessous
 **/
struct Complex
{
    Complex() : real(0.), imag(0.)
    {}
    Complex(double r, double i) : real(r), imag(i)
    {}
    Complex operator + (const Complex& z)
    {
        return Complex(real + z.real, imag + z.imag);
    }
    Complex operator * (const Complex& z)
    {
        return Complex(real * z.real - imag * z.imag, real * z.imag + imag * z.real);
    }
    double sqNorm() { return real * real + imag * imag; }
    double real, imag;
};

std::ostream& operator << (std::ostream& out, const Complex& c)
{
    out << "(" << c.real << "," << c.imag << ")" << std::endl;
    return out;
}

/** Pour un c complexe donn�, calcul le nombre d'it�rations de mandelbrot
 * n�cessaires pour d�tecter une �ventuelle divergence. Si la suite
 * converge, la fonction retourne la valeur maxIter
 **/
int iterMandelbrot(int maxIter, const Complex& c)
{
    Complex z{ 0.,0. };
    // On v�rifie dans un premier temps si le complexe
    // n'appartient pas � une zone de convergence connue :
    // Appartenance aux disques  C0{(0,0),1/4} et C1{(-1,0),1/4}
    if (c.real * c.real + c.imag * c.imag < 0.0625)
        return maxIter;
    if ((c.real + 1) * (c.real + 1) + c.imag * c.imag < 0.0625)
        return maxIter;
    // Appartenance � la cardio�de {(1/4,0),1/2(1-cos(theta))}    
    if ((c.real > -0.75) && (c.real < 0.5)) {
        Complex ct{ c.real - 0.25,c.imag };
        double ctnrm2 = sqrt(ct.sqNorm());
        if (ctnrm2 < 0.5 * (1 - ct.real / ctnrm2)) return maxIter;
    }
    int niter = 0;
    while ((z.sqNorm() < 4.) && (niter < maxIter))
    {
        z = z * z + c;
        ++niter;
    }
    return niter;
}

/**
 * On parcourt chaque pixel de l'espace image et on fait correspondre par
 * translation et homoth�tie une valeur complexe c qui servira pour
 * it�rer sur la suite de Mandelbrot. Le nombre d'it�rations renvoy�
 * servira pour construire l'image finale.

 Sortie : un vecteur de taille W*H avec pour chaque case un nombre d'�tape de convergence de 0 � maxIter
 MODIFICATION DE LA FONCTION :
 j'ai supprim� le param�tre W �tant donn� que maintenant, cette fonction ne prendra plus que des lignes de taille W en argument.
 **/
void
computeMandelbrotSetRow(int W, int H, int maxIter, int num_ligne, int* pixels)
{
    // Calcul le facteur d'�chelle pour rester dans le disque de rayon 2
    // centr� en (0,0)
    double scaleX = 3. / (W - 1);
    double scaleY = 2.25 / (H - 1.);
    //
    // On parcourt les pixels de l'espace image :
    for (int j = 0; j < W; ++j) {
        Complex c{ -2. + j * scaleX,-1.125 + num_ligne * scaleY };
        pixels[j] = iterMandelbrot(maxIter, c);
    }
}

std::vector<int>
computeMandelbrotSet(int W, int H, int maxIter)
{
    std::chrono::time_point<std::chrono::system_clock> start, end;
    start = std::chrono::system_clock::now();

  



    // On parcourt les pixels de l'espace image :

    int nbproc;
    MPI_Comm_size(MPI_COMM_WORLD, &nbproc);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    

    std::vector<int> pixels(W * H);

    if (H % nbproc != 0) { std::cerr << "H n'est pas divisible par nbproc" << std::endl; MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE); }

    int imin_loc = H / nbproc * rank;
    int imax_loc = H / nbproc * (rank + 1);
    int H_loc = H / nbproc;

    std::vector<int> pixels_loc(W * H_loc);

    for (int i = imin_loc; i < imax_loc; ++i) {
        computeMandelbrotSetRow(W, H, maxIter, i, pixels_loc.data() + W * (H_loc - i - 1));
    }

 
   
        
    MPI_Gather(pixels_loc.data(), W*H_loc, MPI_INT, pixels.data(), W*H_loc, MPI_INT, 0, MPI_COMM_WORLD);

   
    // else
    //   MPI_Send(pixels_loc, ...) 

    end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end - start;
    std::cout << "Temps calcul ensemble mandelbrot : " << elapsed_seconds.count()
        << std::endl;
    return pixels;
}

/** Construit et sauvegarde l'image finale **/
void savePicture(const std::string& filename, int W, int H, const std::vector<int>& nbIters, int maxIter)
{
    double scaleCol = 1. / maxIter;//16777216
    std::ofstream ofs(filename.c_str(), std::ios::out | std::ios::binary);
    ofs << "P6\n"
        << W << " " << H << "\n255\n";
    for (int i = 0; i < W * H; ++i) {
        double iter = scaleCol * nbIters[i];
        unsigned char r = (unsigned char)(256 - (unsigned(iter * 256.) & 0xFF));
        unsigned char b = (unsigned char)(256 - (unsigned(iter * 65536) & 0xFF));
        unsigned char g = (unsigned char)(256 - (unsigned(iter * 16777216) & 0xFF));
        ofs << r << g << b;
    }
    ofs.close();
}

int main(int argc, char* argv[])
{
    const int W = 800;
    const int H = 600;
    // Normalement, pour un bon rendu, il faudrait le nombre d'it�rations
    // ci--dessous :
    //const int maxIter = 16777216;
    const int maxIter = 8 * 65536;

    MPI_Init(&argc, &argv);

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    auto iters = computeMandelbrotSet(W, H, maxIter);
    if (rank == 0) {
        savePicture("mandelbrot.tga", W, H, iters, maxIter);
    }

    MPI_Finalize();

    return EXIT_SUCCESS;
}
