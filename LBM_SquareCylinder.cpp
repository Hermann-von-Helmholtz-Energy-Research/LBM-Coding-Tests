#include <iostream>
#include <vector>

using namespace std;

// Função boundary
void boundary(int nx, int ny, vector<vector<vector<double>>>& f, double uo, vector<vector<double>>& rho) {
    // Right hand boundary
    for (int j = 0; j < ny; ++j) {
        f[nx-1][j][2] = f[nx-2][j][2];
        f[nx-1][j][6] = f[nx-2][j][6];
        f[nx-1][j][5] = f[nx-2][j][5];
    }

    // Bottom and top boundary, bounce back
    for (int i = 0; i < nx; ++i) {
        f[i][0][1] = f[i][0][3];
        f[i][0][4] = f[i][0][6];
        f[i][0][5] = f[i][0][7];
        f[i][ny-1][3] = f[i][ny-1][1];
        f[i][ny-1][6] = f[i][ny-1][4];
        f[i][ny-1][7] = f[i][ny-1][5];
    }

    // Left boundary, velocity is given= uo
    for (int j = 1; j < ny-1; ++j) {
        f[0][j][0] = f[0][j][2] + 2.0 * rho[0][j] * uo / 3.0;
        f[0][j][4] = f[0][j][6] - 0.5 * (f[0][j][1] - f[0][j][3]) + rho[0][j] * uo / 6.0;
        f[0][j][7] = f[0][j][5] + 0.5 * (f[0][j][1] - f[0][j][3]) + rho[0][j] * uo / 6.0;
    }
}

// Função collition
void collition(int nx, int ny, vector<vector<double>>& u, vector<vector<double>>& v,
               vector<int>& cx, vector<int>& cy, double omega,
               vector<vector<vector<double>>>& f, vector<vector<double>>& rho) {
    double t1, t2;
    vector<vector<vector<double>>> feq(nx, vector<vector<double>>(ny, vector<double>(9, 0.0)));

    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            t1 = u[i][j] * u[i][j] + v[i][j] * v[i][j];
            for (int k = 0; k < 9; ++k) {
                t2 = u[i][j] * cx[k] + v[i][j] * cy[k];
                feq[i][j][k] = rho[i][j] * w[k] * (1.0 + 3.0 * t2 + 4.5 * t2 * t2 - 1.5 * t1);
                f[i][j][k] = (1.0 - omega) * f[i][j][k] + omega * feq[i][j][k];
            }
        }
    }
}
// Função obstc
void obstc(int nx, int ny, vector<vector<vector<double>>>& f, double uo, vector<vector<double>>& rho) {
    int nxb = (nx - 1) / 5;
    int nxe = nxb + 10;
    int nyb = ((ny - 1) - 10) / 2;
    int nye = nyb + 10;

    // Replace at the entrance, Back Fase Flow
    for (int i = nxb; i <= nxe; ++i) {
        f[i][nyb][3] = f[i][nyb][1];
        f[i][nyb][7] = f[i][nyb][5];
        f[i][nyb][6] = f[i][nyb][8];
        f[i][nye][1] = f[i][nye][3];
        f[i][nye][5] = f[i][nye][7];
        f[i][nye][8] = f[i][nye][6];
    }

    // Bottom and top boundary, bounce back
    for (int j = nyb; j <= nye; ++j) {
        f[nxb][j][2] = f[nxb][j][4];
        f[nxb][j][7] = f[nxb][j][5];
        f[nxb][j][6] = f[nxb][j][8];
        f[nxe][j][4] = f[nxe][j][2];
        f[nxe][j][5] = f[nxe][j][7];
        f[nxe][j][8] = f[nxe][j][6];
    }

    // Set u and v to 0 inside the obstacle
    for (int i = nxb; i <= nxe; ++i) {
        for (int j = nyb; j <= nye; ++j) {
            u[i][j] = 0.0;
            v[i][j] = 0.0;
        }
    }
}

// Função ruv
void ruv(int nx, int ny, vector<vector<vector<double>>>& f,
         vector<vector<double>>& rho, vector<vector<double>>& u, vector<vector<double>>& v) {
    // Calcular a densidade rho
    for (int i = 0; i < nx; ++i) {
        for (int j = 0; j < ny; ++j) {
            rho[i][j] = f[i][j][0] + f[i][j][1] + f[i][j][2] + f[i][j][3] + f[i][j][4]
                      + f[i][j][5] + f[i][j][6] + f[i][j][7] + f[i][j][8];
        }
    }

// Função para salvar os resultados em um arquivo de texto
void result(int nx, int ny, vector<double>& x, vector<double>& y,
            vector<vector<double>>& u, vector<vector<double>>& v, double uo,
            vector<vector<double>>& rho, vector<double>& count, vector<double>& utim) {
    // Nome do arquivo onde os resultados serão salvos
    string filename = "results.txt";
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

