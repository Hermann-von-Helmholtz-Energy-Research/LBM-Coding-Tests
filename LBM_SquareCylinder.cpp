#include <iostream>
#include <vector>

using namespace std;

// Funções que serão definidas posteriormente
void collition(int nx, int ny, vector<vector<double>>& u, vector<vector<double>>& v,
               vector<int>& cx, vector<int>& cy, double omega,
               vector<vector<vector<double>>>& f, vector<vector<double>>& rho);

void stream(int nx, int ny, vector<vector<vector<double>>>& f);

void boundary(int nx, int ny, vector<vector<vector<double>>>& f, double uo,
              vector<vector<double>>& rho);

void obstc(int nx, int ny, vector<vector<vector<double>>>& f, double uo,
           vector<vector<double>>& rho);

void ruv(int nx, int ny, vector<vector<vector<double>>>& f,
         vector<vector<double>>& rho, vector<vector<double>>& u, vector<vector<double>>& v);

// Função para salvar os resultados em um arquivo de texto
void result(int nx, int ny, vector<double>& x, vector<double>& y,
            vector<vector<double>>& u, vector<vector<double>>& v, double uo,
            vector<vector<double>>& rho, vector<double>& count, vector<double>& utim) {
    // Nome do arquivo onde os resultados serão salvos
    string filename = "resultados.txt";

    // Abre o arquivo para escrita
    ofstream outfile(filename);

    if (!outfile.is_open()) {
        cerr << "Erro ao abrir o arquivo " << filename << " para escrita." << endl;
        return;
    }

    // Escreve os resultados no arquivo
    outfile << "Resultados da simulação:" << endl;
    outfile << "Contagem de iterações:" << endl;
    for (int i = 1; i <= 10; ++i) {  // Escrevendo as primeiras 10 iterações
        outfile << "Iteração " << i << ": " << count[i] << ", Rho médio: " << utim[i] << endl;
    }

    // Exemplo de escrita dos dados de u no arquivo
    outfile << endl << "Exemplo de dados de u:" << endl;
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            outfile << "u[" << i << "][" << j << "] = " << u[i][j] << endl;
        }
    }

    // Fechar o arquivo
    outfile.close();

    cout << "Resultados salvos em " << filename << endl;
}

int main() {
    // Definindo variáveis e parâmetros
    int nx = 501;
    int ny = 81;
    double uo = 0.1;
    double c2 = 1.0 / 3.0;
    double dx = 1.0;
    double dy = 1.0;
    double xl = static_cast<double>(nx - 1) / (ny - 1);
    double yl = 1.0;
    double alpha = 0.01;
    double ReH = uo * (ny - 1) / alpha;
    double ReD = uo * 10.0 / alpha;
    double omega = 1.0 / (3.0 * alpha + 0.5);

    // Inicializando vetores
    vector<vector<vector<double>>> f(nx, vector<vector<double>>(ny, vector<double>(9, 0.0)));
    vector<vector<vector<double>>> feq(nx, vector<vector<double>>(ny, vector<double>(9, 0.0)));
    vector<vector<double>> u(nx, vector<double>(ny, 0.0)); // Inicializa u com zeros
    vector<vector<double>> v(nx, vector<double>(ny, 0.0));
    vector<vector<double>> rho(nx, vector<double>(ny, 2.0));
    vector<double> x(nx, 0.0);
    vector<double> y(ny, 0.0);
    vector<double> count(8001, 0.0); // Inicializa com 8001 elementos
    vector<double> utim(8001, 0.0);  // Inicializa com 8001 elementos
    vector<double> Tm(nx, 0.0);
    vector<double> Tvm(nx, 0.0);

    // Pesos
    vector<double> w = {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0,
                        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0};

    // Vetores unitários ao longo das direções de streaming da rede de Boltzmann
    vector<int> cx = {1, 0, -1, 0, 1, -1, -1, 1, 0};
    vector<int> cy = {0, 1, 0, -1, 1, 1, -1, -1, 0};

    // Loop principal
    for (int kk = 1; kk <= 8000; ++kk) {
        // Collitions
        collition(nx, ny, u, v, cx, cy, omega, f, rho);

        // Streaming
        stream(nx, ny, f);

        // Boundary condition
        boundary(nx, ny, f, uo, rho);

        // Obsticale
        obstc(nx, ny, f, uo, rho);

        // Calculate rho, u, v
        ruv(nx, ny, f, rho, u, v);

        // Atualização de count e utim
        count[kk] = kk;
        utim[kk] = rho[(nx - 1) / 2][(ny - 1) / 2]; // Ajustar índices conforme necessário
    }

    // Chamada da função result para salvar os dados em um arquivo
    result(nx, ny, x, y, u, v, uo, rho, count, utim);
    
    return 0;
}

