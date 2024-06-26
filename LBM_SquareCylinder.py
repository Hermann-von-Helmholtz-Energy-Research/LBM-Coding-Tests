import numpy as np
import matplotlib.pyplot as plt

# Inicialização e parâmetros
nx, ny = 501, 81 # (número de pontos na direção x e y)
uo = 0.1 # velocidade inicial

# Inicialização das matrizes
f = np.zeros((nx, ny, 9))
feq = np.zeros((nx, ny, 9))

# Vetores para armazenar dados temporais
utim = np.zeros(1001)
count = np.zeros(1001)

# Inicialização das matrizes bidimensionais de velocidade nas direções de x e y
u = uo * np.ones((nx, ny))
v = np.zeros((nx, ny))

# Matriz bidimensional de densidade
rho = 2 * np.ones((nx, ny)) # Armazena a densidade do fluido em cada ponto

x = np.zeros(nx) # coordenadas x
y = np.zeros(ny) # coordenadas y

# Inicialização das Matrizes temporais para cálculos intermediários
Tm = np.zeros(nx) 
Tvm = np.zeros(nx)

# Pesos
w = np.zeros(9)
w = np.array([1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36, 4/9])

# Vetores unitários nas direções da grade
cx = np.array([1, 0, -1, 0, 1, -1, -1, 1, 0])
cy = np.array([0, 1, 0, -1, 1, 1, -1, -1, 0])

c2 = 1.0 / 3.0  # constante velocidade ao quadrado
dx, dy = 1.0, 1.0  # incremento de espaços

xl = (nx - 1) / (ny - 1)  # comprimento da grade na direção x
yl = 1.0  # comprimento da grade na direção y

x = np.arange(nx)
y = np.arange(ny)
"""
np.arange() em Python
np.arange() é uma função da biblioteca NumPy que cria um array contendo uma sequência de números. Ela é similar à função range() do Python padrão, 
mas retorna um array NumPy em vez de uma lista padrão do Python.

Por que x e y usam np.arange(nx) e np.arange(ny)?
No código Matlab original, x e y são inicializados com vetores que representam as coordenadas ao longo das direções x e y da grade. 
Em Matlab, 0:1:nx-1 cria um vetor de 0 até nx-1 com um passo de 1.

Para replicar isso em Python com NumPy, usamos np.arange(nx) e np.arange(ny), que cria vetores de 0 até nx-1 e de 0 até ny-1, respectivamente. Essa função é conveniente porque:

Sintaxe Simples: np.arange(start, stop, step) permite controlar o início, o fim e o passo da sequência.
Compatibilidade com NumPy: Retorna um array NumPy, que é eficiente para operações matriciais e vetoriais.
Portanto, x = np.arange(nx) cria um vetor de números de 0 a 500 (para nx=501), e y = np.arange(ny) cria um vetor de números de 0 a 80 (para ny=81), 
representando as coordenadas ao longo das direções x e y da grade, respectivamente.
"""

alpha = 0.01  # coeficiente de viscosidade cinemática
ReH = uo * (ny - 1) / alpha  # número de Reynolds baseado na altura do canal
ReD = uo * 10.0 / alpha  # número de Reynolds baseado no tamanho do obstáculo.
omega = 1.0 / (3.0 * alpha + 0.5) # fator de relaxação 

count[0] = 0

# setting velocity
for j in range(1, ny-1):
    u[0, j] = uo

# Função boundary
def boundary(nx, ny, f, uo, rho):
    # Right hand boundary
    for j in range(ny):
        f[nx-1, j, 2] = f[nx-2, j, 2]
        f[nx-1, j, 5] = f[nx-2, j, 5]
        f[nx-1, j, 6] = f[nx-2, j, 6]
    
    # Bottom and top boundary, bounce back
    for i in range(nx):
        f[i, 0, 4] = f[i, 1, 4]
        f[i, 0, 7] = f[i, 1, 7]
        f[i, 0, 8] = f[i, 1, 8]
        
        f[i, ny-1, 2] = f[i, ny-2, 2]
        f[i, ny-1, 5] = f[i, ny-2, 5]
        f[i, ny-1, 6] = f[i, ny-2, 6]
        
        u[i, 0] = 0.0
        v[i, 0] = 0.0
        u[i, ny-1] = 0.0
        v[i, ny-1] = 0.0
    
    # Left boundary, velocity is given = uo
    for j in range(1, ny-1):
        f[0, j, 1] = f[0, j, 3] + 2.0 * rho[0, j] * uo / 3.0
        f[0, j, 5] = f[0, j, 7] - 0.5 * (f[0, j, 2] - f[0, j, 4]) + rho[0, j] * uo / 6.0
        f[0, j, 8] = f[0, j, 6] + 0.5 * (f[0, j, 2] - f[0, j, 4]) + rho[0, j] * uo / 6.0
        u[0, j] = uo
        v[0, j] = 0.0
    
    return f

# Função collition
def collition(nx, ny, u, v, cx, cy, omega, f, rho, w):
    feq = np.zeros((nx, ny, 9))  # Inicializa matriz feq com zeros
    
    for j in range(ny):
        for i in range(nx):
            t1 = u[i, j]**2 + v[i, j]**2
            for k in range(9):
                t2 = u[i, j] * cx[k] + v[i, j] * cy[k]
                feq[i, j, k] = rho[i, j] * w[k] * (1.0 + 3.0 * t2 + 4.5 * t2**2 - 1.5 * t1)
                f[i, j, k] = (1.0 - omega) * f[i, j, k] + omega * feq[i, j, k]
    
    return f

# Objeto
def obstc(nx, ny, f, uo, rho):
    # Length of obstacle = nx/5, and has sides of 10 units
    nxb = int((nx-1) / 5)
    nxe = nxb + 10
    nyb = int(((ny-1) - 10) / 2)
    nyb = 35
    nye = nyb + 10
    
    for i in range(nxb, nxe+1):
        f[i, nyb, 2] = f[i, nyb, 4]
        f[i, nyb, 5] = f[i, nyb, 7]
        f[i, nyb, 6] = f[i, nyb, 8]
        f[i, nye, 4] = f[i, nye, 2]
        f[i, nye, 7] = f[i, nye, 5]
        f[i, nye, 8] = f[i, nye, 6]
    
    # Bottom and top boundary, bounce back
    for j in range(nyb, nye+1):
        f[nxb, j, 3] = f[nxb, j, 1]
        f[nxb, j, 7] = f[nxb, j, 5]
        f[nxb, j, 6] = f[nxb, j, 8]
        f[nxe, j, 1] = f[nxe, j, 3]
        f[nxe, j, 5] = f[nxe, j, 7]
        f[nxe, j, 8] = f[nxe, j, 6]
    
    for i in range(nxb, nxe+1):
        for j in range(nyb, nye+1):
            u[i, j] = 0.0
            v[i, j] = 0.0
    
    return f

# result
def result(nx, ny, x, y, u, v, uo, rho, count, utim):
    Tm1 = np.zeros(ny)
    Tm2 = np.zeros(ny)
    Tm3 = np.zeros(ny)
    Tm4 = np.zeros(ny)
    
    for j in range(ny):
        Tm1[j] = u[51, j] / uo
        Tm2[j] = u[101, j] / uo
        Tm3[j] = u[261, j] / uo
        Tm4[j] = u[301, j] / uo
    
    umx = np.zeros(nx)
    vmx = np.zeros(nx)
    
    for i in range(nx):
        umx[i] = u[i, int((ny-1)/2)] / uo
        vmx[i] = v[i, int((ny-1)/2)] / uo
    
    plt.figure()
    plt.plot(x / (nx-1), umx, label='Ux', linewidth=1.5)
    plt.plot(x / (nx-1), vmx, label='Uy', linewidth=1.5)
    plt.xlabel('X')
    plt.ylabel('Velocity')
    plt.legend()
    
    plt.figure()
    plt.plot(Tm1, y, label='Tm1', linewidth=1.5)
    plt.plot(Tm2, y, label='Tm2', linewidth=1.5)
    plt.plot(Tm3, y, label='Tm3', linewidth=1.5)
    plt.plot(Tm4, y, label='Tm4', linewidth=1.5)
    plt.xlabel('U')
    plt.ylabel('Y')
    plt.legend()
    
    plt.figure()
    plt.plot(count, utim)
    plt.xlabel('Iterations')
    plt.ylabel('Ut')
    
    plt.show()

# função ruv
def ruv(nx, ny, f):
    rho = np.sum(f, axis=2)  # soma ao longo da terceira dimensão (índice 2)
    
    for i in range(nx):
        rho[i, ny-1] = f[i, ny-1, 8] + f[i, ny-1, 0] + f[i, ny-1, 2] + 2.0 * (f[i, ny-1, 1] + f[i, ny-1, 5] + f[i, ny-1, 4])
    
    # calculate velocity components
    u = ( np.sum(f[:, :, [0, 4, 7]], axis=2) - np.sum(f[:, :, [2, 5, 6]], axis=2) ) / rho
    v = ( np.sum(f[:, :, [1, 4, 5]], axis=2) - np.sum(f[:, :, [3, 6, 8]], axis=2) ) / rho
    
    return rho, u, v

# função stream
def stream(f):
    # Direções de streaming (shifts)
    shifts = [(+1, 0), (0, +1), (-1, 0), (0, -1), (+1, +1), (-1, +1), (-1, -1), (+1, -1)]
    
    # Aplicando os deslocamentos para cada direção de streaming
    for k, shift in enumerate(shifts):
        f[:, :, k] = np.roll(f[:, :, k], shift, axis=(0, 1))
    
    return f

# def stream_function_calculation(nx, ny, x, y, u):
    sx = np.zeros((nx, ny))
    sy = np.zeros((nx, ny))
    str = np.zeros((nx, ny))
    
    for j in range(ny):
        sx[:, j] = x[:]
    
    for i in range(nx):
        sy[i, :] = y[:]
    
    for i in range(nx):
        for j in range(1, ny):
            str[i, j] = str[i, j-1] + 0.5 * (u[i, j] + u[i, j-1])
    
    plt.figure()
    plt.contour(sx, sy, str)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Stream Function Contour')
    
    plt.figure()
    plt.contour(sx, sy, u, linewidths=1.0)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.title('Velocity Contour')
    
    plt.show()

# Main Loop
for kk in range(1, 8001):  # Python range é exclusivo, então vai até 8000
    # Collitions
    f = collition(nx, ny, u, v, cx, cy, omega, f, rho, w)
    
    # Streaming
    f = stream(f)
    
    # Boundary condition
    f = boundary(nx, ny, f, uo, rho)
    
    # Obstaculo
    f = obstc(nx, ny, f, uo, rho)
    
    # Calculate rho, u, v
    rho, u, v = ruv(nx, ny, f)
    
    count[kk-1] = kk  # kk-1 porque Python usa índices baseados em zero
    utim[kk-1] = rho[(nx-1)//2, (ny-1)//2]  # Usamos // para divisão inteira em Python
    # Ao final do loop, count e utim contêm os resultados desejados


# Plotting data
result(nx, ny, x, y, u, v, uo, rho, count, utim)

# Calculando a função de corrente
sx = np.zeros((nx, ny))
sy = np.zeros((nx, ny))
str = np.zeros((nx, ny))

for j in range(ny):
    sx[:, j] = x[:]
for i in range(nx):
    sy[i, :] = y[:]

for i in range(nx):
    for j in range(1, ny):
        str[i, j] = str[i, j-1] + 0.5 * (u[i, j] + u[i, j-1])

# Plotando a função de corrente
plt.figure()
plt.contour(sx, sy, str)

# Plotando o contorno da velocidade u
plt.figure()
plt.contour(sx, sy, u, linewidths=1.0)

# Exibindo os gráficos
plt.show()
