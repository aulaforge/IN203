// Produit matrice-vecteur
# include <cassert>
# include <vector>
# include <iostream>
# include <mpi.h>

// ---------------------------------------------------------------------
class Matrix : public std::vector<double>
{
public:
    Matrix(int dim);
    Matrix(int nrows, int ncols);
    Matrix(const Matrix& A) = delete;
    Matrix(Matrix&& A) = default;
    ~Matrix() = default;

    Matrix& operator = (const Matrix& A) = delete;
    Matrix& operator = (Matrix&& A) = default;

    double& operator () (int i, int j) {
        return m_arr_coefs[i + j * m_nrows];
    }
    double  operator () (int i, int j) const {
        return m_arr_coefs[i + j * m_nrows];
    }

    std::vector<double> operator * (const std::vector<double>& u) const;

    std::ostream& print(std::ostream& out) const
    {
        const Matrix& A = *this;
        out << "[\n";
        for (int i = 0; i < m_nrows; ++i) {
            out << " [ ";
            for (int j = 0; j < m_ncols; ++j) {
                out << A(i, j) << " ";
            }
            out << " ]\n";
        }
        out << "]";
        return out;
    }
private:
    int m_nrows, m_ncols;
    std::vector<double> m_arr_coefs;
};
// ---------------------------------------------------------------------
inline std::ostream&
operator << (std::ostream& out, const Matrix& A)
{
    return A.print(out);
}
// ---------------------------------------------------------------------
inline std::ostream&
operator << (std::ostream& out, const std::vector<double>& u)
{
    out << "[ ";
    for (const auto& x : u)
        out << x << " ";
    out << " ]";
    return out;
}
// ---------------------------------------------------------------------
std::vector<double>
Matrix::operator * (const std::vector<double>& u) const
{
    const Matrix& A = *this;
    assert(u.size() == unsigned(m_ncols));
    std::vector<double> v(m_nrows, 0.);
    for (int i = 0; i < m_nrows; ++i) {
        for (int j = 0; j < m_ncols; ++j) {
            v[i] += A(i, j) * u[j];
        }
    }
    return v;
}

// =====================================================================
Matrix::Matrix(int dim) : m_nrows(dim), m_ncols(dim),
m_arr_coefs(dim* dim)
{
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            (*this)(i, j) = (i + j) % dim;
        }
    }
}
// ---------------------------------------------------------------------
Matrix::Matrix(int nrows, int ncols) : m_nrows(nrows), m_ncols(ncols),
m_arr_coefs(nrows* ncols)
{
    int dim = (nrows > ncols ? nrows : ncols);
    for (int i = 0; i < nrows; ++i) {
        for (int j = 0; j < ncols; ++j) {
            (*this)(i, j) = (i + j) % dim;
        }
    }
}
// =====================================================================
int main(int nargs, char* argv[])
{
    MPI_Init(&nargs, &argv);


    int rang, nproc;

    MPI_Comm_rank(MPI_COMM_WORLD, &rang);

    MPI_Comm_size(MPI_COMM_WORLD, &nproc);





    const int N = 8;

    if (N % nproc != 0) { std::cerr << "N n'est pas divisible pas le nombre de processus" << std::endl; MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);}
    else {
        const int Taille = N / nproc;
        Matrix A(N);
        //std::cout << "A : " << A << std::endl;
        
        Matrix Mpart (N, Taille);
        int imin = Taille * rang;
        int imax = Taille * (rang + 1);
        for (int i = imin; i < imax; i++) {
            for (int j = 0; j < N; j++) {
                Mpart(j, i) = A(j, i);
            }
        }


        std::vector<double> u(N);
        for (int i = 0; i < N; ++i) u[i] = i + 1;
        //std::cout << " u : " << u << std::endl;
        
        std::vector<double> vpart(N);
        for (int i = 0; i < N; i++) {
            for (int j = Taille*rang; j < Taille*(rang+1); j++) {
                vpart[i] += Mpart(i, j) * u[j];
            }
        }
        
        double result[N];

        MPI_Allreduce(vpart.data(), result, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        
        
        
        std::vector<double> v = A * u;
        std::cout << "Je suis le processus " << rang << ", A.u = " << v << " et j'ai obtenu [";
        
        for (int i = 0; i < N - 1; i++) {

            std::cout << result[i] << " ; ";
        }

        std::cout << result[N - 1] << "]" << std::endl;


       
    }
    MPI_Finalize();

    return EXIT_SUCCESS;
}


