#include <iostream>
#include <vector>

using namespace std;

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
    vector<vector<double>> u(nx, vector<double>(ny, uo));
    vector<vector<double>> v(nx, vector<double>(ny, 0.0));
    vector<vector<double>> rho(nx, vector<double>(ny, 2.0));
    vector<double> x(nx, 0.0);
    vector<double> y(ny, 0.0);
    vector<double> count(1001, 0.0);
    vector<double> utim(1001, 0.0);
    vector<double> Tm(nx, 0.0);
    vector<double> Tvm(nx, 0.0);

    // Pesos
    vector<double> w = {1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 1.0 / 9.0, 
                        1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 1.0 / 36.0, 4.0 / 9.0};

    // Vetores unitários ao longo das direções de streaming da rede de Boltzmann
    vector<int> cx = {1, 0, -1, 0, 1, -1, -1, 1, 0};
    vector<int> cy = {0, 1, 0, -1, 1, 1, -1, -1, 0};

    return 0;
}
