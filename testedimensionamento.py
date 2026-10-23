import math

FS = 1.5          # fator de seguranca
esp = 0.3           # mm por camada
n_max = 3   # limite da tabela de propriedades
n_min = 2

# geometria (mm)
alavanca_tdp, roda_tdp, L_tdp = 52.2, 40.0, 450.0
alavanca_tdn, roda_tdn, L_tdn, angulo_tdn = 115.5, 30.0, 136.5, 7

# diametros externos comerciais (mm)
diametros = [12.7, 15.88, 19.05, 22.2, 25.4, 28.57, 31.75]

# propriedades do laminado (MPa), na ordem: 1, 2 e 3 camadas
rho = 1150e-9   # kg/mm3
E_lam   = [13400, 20520, 36900]    # modulo de elasticidade
G_lam   = [1200, 3100, 3800]       # modulo de cisalhamento
tra_lam = [254.5, 285.7, 361.2]    # resistencia a tracao
com_lam = [203.6, 228.56, 288.96]  # resistencia a compressao
cis_lam = [74.4, 96.2, 112.1]      # resistencia ao cisalhamento

# casos de carga: vertical, horizontal, lateral (N)
casos_tdp = [
    (168.6, 31.6, 0.0),     # 3 pontos
    (181.3, 33.9, 0.0),     # 2 pontos / 1 ponto
    (114.2, 0.0, 67.9),     # 2 pontos, carga lateral
    (114.2, 0.0, 44.8),     # 2 pontos, carga lateral
    (132.8, 0.0, 105.5),    # taxi: decolagem abortada
    (76.4, 0.0, 59.6),      # taxi: zigue-zague
]
casos_tdn = [
    (25.3, 4.7, 0.0),       # 3 pontos
    (9.6, 0.0, 6.3),        # taxi: decolagem abortada
    (8.0, 0.0, 3.5),        # taxi: zigue-zague
]

nomes_tdp = [
    "3 pontos",
    "2 pontos / 1 ponto",
    "2 pontos, lateral 67,9 N",
    "2 pontos, lateral 44,8 N",
    "taxi: decolagem abortada",
    "taxi: zigue-zague",
]
nomes_tdn = [
    "3 pontos",
    "taxi: decolagem abortada",
    "taxi: zigue-zague",
]

# furo na bequilha (mm) e curva K do grafico (placa com furo circular, K x 2r/D)
d_furo = 3.0     # furo do parafuso do mancal (mm)
n_furos = 2      # parafuso atravessa as duas paredes do tubo
x_grafico = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]
K_grafico = [3.0, 2.72, 2.50, 2.36, 2.27, 2.18, 2.11, 2.05]


def fator_K(x):
    for i in range(len(x_grafico) - 1):
        if x <= x_grafico[i + 1]:
            x1 = x_grafico[i]
            x2 = x_grafico[i + 1]
            k1 = K_grafico[i]
            k2 = K_grafico[i + 1]
            return k1 + (k2 - k1) * (x - x1) / (x2 - x1)
    return K_grafico[-1]

def verifica_tdp(casos, alavanca, R, d, n, comp):
    Sc = com_lam[n - 1]
    St = tra_lam[n - 1]
    Ss = cis_lam[n - 1]
    D = d + 2 * n * esp
    I = math.pi / 64 * (D**4 - d**4)
    A = math.pi / 4 * (D**2 - d**2)
    c = D / 2 #ponto mais externo até o centro

    th_max = 0.0
    caso_c = ""

    for i, (V, H, L) in enumerate(casos):
        My = (V * alavanca) + (L * R)               # vertical + momento da lateral
        Mz = H * alavanca                           # arrasto
        M = FS * math.sqrt(My**2 + Mz**2)
        N = FS * L
        Vres = FS * math.sqrt(V**2 + H**2)
        tau =  Vres / A
        sig = -(M * c / I) - N/A
        sig_1 = sig/2 + math.sqrt((sig/2)**2 + tau**2)
        sig_2 = sig/2 - math.sqrt((sig/2)**2 + tau**2)
        tau12 = (sig_1 - sig_2)/2

        TH_tracao = (max(sig_1, 0) / St)**2 + (tau12 / Ss)**2
        TH_compressao = (max(-sig_2, 0) / Sc)**2 + (tau12 / Ss)**2

        TH = max(TH_tracao, TH_compressao)

        if TH > th_max:
            th_max = TH
            caso_c = nomes_tdp[i]

    ms_th = 1 / math.sqrt(th_max) - 1
    massa = A * comp * rho * 1000
    return D, ms_th, massa, caso_c


def verifica_tdn(casos, s, R, d, n, comp):
    Sc = com_lam[n - 1]
    Ss = cis_lam[n - 1]
    D = d + 2 * n * esp
    I = math.pi / 64 * (D**4 - d**4)
    A = math.pi / 4 * (D**2 - d**2)
    c = D / 2
    t = n * esp
    r = (D - t) / 2

    tdn_ang = math.radians(angulo_tdn)
    sen = math.sin(tdn_ang)
    cos = math.cos(tdn_ang)

    W = 2 * math.pi * r             # largura "da placa" = circunferencia media do tubo
    x_furo = d_furo / W        
    K = fator_K(x_furo)
    x = n_furos * x_furo            # parte da largura perdida pelos furos

    th_max = 0.0
    caso_c = ""

    for i, (V, H, Lat) in enumerate(casos):
        for Hs in (H, -H):                                # arrasto para um lado e para o outro
            N = FS * (V * cos - Hs * sen)           # forca axial (ao longo do tubo)
            F_in = FS * (V * sen + Hs * cos)        # forca transversal no plano
            M_in = FS * (Hs * (R + (s * cos)) + (V * s * sen))   # momento de flexao no plano
            M_lat = FS * Lat * (s + (R * cos))        # flexao causada pela lateral
            M = math.sqrt(M_in**2 + M_lat**2)
            F = math.sqrt(F_in**2 + (FS * Lat)**2)         # forca transversal total
            sig = (M * c / I) + ((N) / A)
            tau = F / A
            sig_semfuro = sig / (1 - x)             
            sig_conc = K * sig_semfuro              # tensao concentrada: sigma_max = K * sigma
            sig_x = -sig_conc
            sig_2 = sig_x/2 - math.sqrt((sig_x/2)**2 + tau**2)
            sig_1 = sig_x/2 + math.sqrt((sig_x/2)**2 + tau**2)
            tau12 = (sig_1 - sig_2)/2
            TH = (max(-sig_2, 0) / Sc)**2 + (tau12 / Ss)**2
            if TH > th_max:
                th_max = TH
                caso_c = nomes_tdn[i]
                if Hs < 0:
                    caso_c = caso_c + " (arrasto para tras)"

    ms_th = 1 / math.sqrt(th_max) - 1
    massa = A * comp * rho * 1000
    return D, ms_th, massa, caso_c


def dimensiona(nome, tipo, casos, alavanca, R, comp):
    print(f"\n{nome} | FS = {FS}")

    achou = False
    melhor_massa = 1000
    melhor_d = 0
    melhor_n = 0
    melhor_D = 0
    melhor_th = 0
    melhor_caso = ""

    for d in diametros:
        for n in range(n_min, n_max + 1):
            if tipo == "tdp":
                D, ms_th, massa, caso = verifica_tdp(casos, alavanca, R, d, n, comp)
            else:
                D, ms_th, massa, caso = verifica_tdn(casos, alavanca, R, d, n, comp)

            if ms_th >= 0.1:                      
                if massa < melhor_massa:           # criterio: MENOR MASSA
                    achou = True
                    melhor_massa = massa
                    melhor_d = d
                    melhor_n = n
                    melhor_D = D
                    melhor_th = ms_th
                    melhor_caso = caso
                

    if not achou:
        print("nao tanka")
        return

    print(f"ACHEI! (menor massa, {melhor_massa:.1f} g): molde {melhor_d} mm | {melhor_n} camadas | D final = {melhor_D:.2f} mm")
    print(f"MS GERAL = {melhor_th * 100:.2f}%  (tsai-hill)")
    print(f"CASO CRITICO: {melhor_caso}")

dimensiona("trem de pouso", "tdp", casos_tdp, alavanca_tdp, roda_tdp, L_tdp)
dimensiona("bequilha", "tdn", casos_tdn, alavanca_tdn, roda_tdn, L_tdn)